import os
import sys
import argparse
import joblib

import pysam
import duckdb
import logging
import truvari
import numpy as np
import pandas as pd

import trgt

def pull_alleles(data):
    """
    Turn alleles and their sample allele properties into tables
    """
    alleles = pd.DataFrame(data["ALT"].to_list(), columns=["ALT1", "ALT2"], index=data.index)
    sap = pd.DataFrame(data["SD"].to_list(), columns=["SD1", "SD2"], index=data.index)

    data = pd.concat([data, alleles, sap], axis=1)

    data['a1'] = np.where(data['REF'] == data['ALT1'], None, data['ALT1'])
    data['a2'] = np.where(data['REF'] == data['ALT2'], None, data['ALT2'])
    part1 = data[['ID', 'a1']].dropna().copy()
    part1.columns = ["LocusID", "sequence"]
    part2 = data[['ID', 'a2']].dropna().copy()
    part2.columns = ["LocusID", "sequence"]
    all_alleles = pd.concat([part1, part2]).sort_values(["LocusID", "sequence"]).drop_duplicates()
    all_alleles['allele_number'] = all_alleles.groupby('LocusID')['LocusID'].rank(method='first').astype(int)
    all_alleles['allele_length'] = all_alleles['sequence'].str.len()

    aidx = all_alleles.set_index(["LocusID", "sequence"])

    part1 = data[["ID", "a1", "SD1"]]
    part1.columns = ["LocusID", "sequence", "spanning_reads"]
    part2 = data[["ID", "a2", "SD2"]]
    part2.columns = ["LocusID", "sequence", "spanning_reads"]
    all_sap = pd.concat([part1, part2]).set_index(["LocusID", "sequence"]).join(aidx, how='left')
    all_sap["allele_number"] = all_sap['allele_number'].fillna(0)
    all_sap = all_sap.reset_index().drop(columns=['sequence'])
    return all_alleles.reset_index(drop=True), all_sap.reset_index(drop=True)

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

    truvari.setup_logging()
    if not os.path.exists(args.vcf):
        raise RuntimeError(f"input {args.vcf} does not exist")
    if os.path.exists(args.dbname):
        raise RuntimeError(f"output {args.dbname} already exists")

    logging.info("Loading VCF")
    data = truvari.vcf_to_df(args.vcf, with_info=True, with_fmt=True, no_prefix=True, with_seqs=True)
    logging.info("Parsed %d Loci", len(data))
    data['ID'] = range(len(data))

    con = duckdb.connect(database=args.dbname, read_only=False)
    here = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(here, "schema.sql"), 'r') as fh:
        con.execute(fh.read())

    locus_df = data[["ID", "chrom", "start", "end"]]
    con.execute("INSERT INTO Locus SELECT * FROM locus_df")

    sample_name = pysam.VariantFile(args.vcf).header.samples[0]
    con.execute(f"INSERT INTO Sample (SampleID, name) VALUES (0, '{sample_name}')")

    logging.info("Wrangling Alleles")
    allele_df, sap_df = pull_alleles(data)
    sap_df["SampleID"] = 0
    #allele_df = allele_df[['LocusID', 'allele_number', 'allele_length', 'sequence']]
    #sap_df = sap_df[["SampleID", "LocusID", "allele_number", "spanning_reads"]]
    
    logging.info("Encoding Alleles")
    allele_df['sequence'] = allele_df['sequence'].apply(trgt.dna_encode)
    logging.info("Parsed %d Alleles", len(allele_df))

    con.execute("INSERT INTO Allele SELECT LocusID, allele_number, allele_length, sequence FROM allele_df")
    con.execute("INSERT INTO SampleAlleleProperties SELECT SampleID, LocusID, allele_number, spanning_reads FROM sap_df")

    con.close()
    logging.info("Finished DuckDB")

    jl_out = {'l':locus_df, 's':sample_name, 'a':allele_df, 'sap':sap_df}
    joblib.dump(jl_out, args.dbname + '.jl')
    locus_df.to_parquet(args.dbname + '.l.pq')
    allele_df.to_parquet(args.dbname + '.a.pq')
    sap_df.to_parquet(args.dbname + '.s.pq')

