import os
import sys
import argparse

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
    parser = argparse.ArgumentParser(prog="trgt db create", description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("vcf", metavar="VCF", type=str,
                        help="Input TRGT VCF")
    parser.add_argument("dbname", metavar="DB", type=str,
                        help="Output TRGT DB")
    args = parser.parse_args(args)

    if not os.path.exists(args.vcf):
        raise RuntimeError(f"input {args.vcf} does not exist")
    if os.path.exists(args.dbname):
        raise RuntimeError(f"output {args.dbname} already exists")

    data = truvari.vcf_to_df(args.vcf, with_info=True, with_fmt=True, no_prefix=True, with_seqs=True)
    data['ID'] = range(len(data))

    locus_df = data[["ID", "chrom", "start", "end"]]

    con = duckdb.connect(database=args.dbname, read_only=False)
    here = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(here, "schema.sql"), 'r') as fh:
        con.execute(fh.read())

    con.execute("INSERT INTO Locus SELECT * FROM locus_df")

    sample_name = pysam.VariantFile(args.vcf).header.samples[0]
    con.execute(f"INSERT INTO Sample (SampleID, name) VALUES (0, '{sample_name}')")

    allele_df, sap_df = pull_alleles(data)
    con.execute("INSERT INTO Allele SELECT * FROM allele_df")
    con.execute("INSERT INTO SampleAlleleProperties SELECT * FROM sap_df")

    con.close()
