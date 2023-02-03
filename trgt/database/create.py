import sys
import pysam
import duckdb
import truvari
import pandas as pd

import trgt

def pull_alleles(data, sample_id=0):
    """
    Builds rows for Alelle and SampleAlleleProperties tables
    """
    allele_rows = []
    sap_rows = []
    
    allele_id = 0 # PRIMARY KEY
    for idx, row in data.iterrows():
        # For deduping alleles per-locus
        cur_locus_alleles = {row['REF']: allele_id}
        
        allele_rows.append([allele_id, row["ID"], len(row['REF']), None])
        
        for idx, alt in enumerate(row["ALT"]):
            if alt not in cur_locus_alleles:
                allele_id += 1
                cur_locus_alleles[alt] = allele_id
                allele_rows.append([allele_id, row["ID"], len(alt), trgt.dna_encode(alt)])
            sap_rows.append([sample_id, cur_locus_alleles[alt], row["SD"][idx]])
        allele_id += 1
                
    allele_rows = pd.DataFrame(allele_rows, columns=["AlleleID", "LocusID", "allele_length", "sequence"])
    sap_rows = pd.DataFrame(sap_rows, columns=["SampleID", "AlleleID", "spanning_reads"])

    return allele_rows, sap_rows


def create_main(args):
    """
    Create a new duckdb and fill it with a VCF
    """
    # TODO: argparse
    in_vcf = args[0]
    out_db = args[1]

    data = truvari.vcf_to_df(in_vcf, with_info=True, with_fmt=True, no_prefix=True, with_seqs=True)
    data['ID'] = range(len(data))

    locus_df = data[["ID", "chrom", "start", "end"]]
    con = duckdb.connect(database=out_db ,read_only=False)
    with open("schema.sql") as fh:
        con.execute(fh.read())

    con.execute("INSERT INTO Locus SELECT * FROM locus_df")

    sample_name = pysam.VariantFile(in_vcf).header.samples[0]
    con.execute(f"INSERT INTO Sample (SampleID, name) VALUES (0, '{sample_name}')")

    allele_df, sap_df = pull_alleles(data)
    con.execute("INSERT INTO Allele SELECT * FROM allele_df")
    con.execute("INSERT INTO SampleAlleleProperties SELECT * FROM sap_df")

    con.close()
