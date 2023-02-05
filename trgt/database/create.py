import os
import sys
import argparse
import joblib

import pysam
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

    # Split VCF line alleles and dedup
    data['a1'] = np.where(data['REF'] == data['ALT1'], None, data['ALT1'])
    data['a2'] = np.where(data['REF'] == data['ALT2'], None, data['ALT2'])
    part1 = data[['LocusID', 'a1']].dropna().rename(columns={'a1':'sequence'})
    part2 = data[['LocusID', 'a2']].dropna().rename(columns={'a2':'sequence'})
    all_alleles = pd.concat([part1, part2], ignore_index=True).sort_values(["LocusID", "sequence"]).drop_duplicates()
    all_alleles['allele_number'] = all_alleles.groupby('LocusID')['LocusID'].rank(method='first')
    all_alleles['allele_length'] = all_alleles['sequence'].str.len()
    # Adding reference allele
    partr = data[["LocusID"]].copy()
    partr["allele_number"] = 0
    partr["allele_length"] = data["end"] - data["start"]
    all_alleles = pd.concat([all_alleles, partr], ignore_index=True)
    all_alleles['allele_number'] = all_alleles['allele_number'].astype(int)

    # allele number lookup
    aidx = all_alleles.set_index(["LocusID", "sequence"])

    # Split VCF line sample allele properties and tie to their allele number
    part1 = data[["LocusID", "a1", "SD1"]].rename(columns={'a1':"sequence", "SD1":"spanning_reads"})
    part2 = data[["LocusID", "a2", "SD2"]].rename(columns={'a2':"sequence", "SD2":"spanning_reads"})
    all_sap = pd.concat([part1, part2]).set_index(["LocusID", "sequence"]).join(aidx['allele_number'], how='left')
    # all_sap["allele_number"] = all_sap['allele_number'].fillna(0).astype(int)
    all_sap = all_sap.reset_index().drop(columns=['sequence'])

    return all_alleles.reset_index(drop=True), all_sap

def create_main(args):
    """
    Create a new database and fill it with a VCF
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

    os.mkdir(args.dbname)

    logging.info("Loading VCF")
    data = truvari.vcf_to_df(args.vcf, with_info=True, with_fmt=True, no_prefix=True, with_seqs=True)
    logging.info("Parsed %d Loci", len(data))
    data['LocusID'] = range(len(data))

    pq_fns = trgt.get_tdb_files(args.dbname)
    data[["LocusID", "chrom", "start", "end"]].reset_index(drop=True).to_parquet(pq_fns['locus'])

    logging.info("Wrangling Alleles")
    allele_df, sap_df = pull_alleles(data)

    sample_name = pysam.VariantFile(args.vcf).header.samples[0]
    s_fn = os.path.join(args.dbname, f'sample.{sample_name}.pq')
    sap_df[["LocusID", "allele_number", "spanning_reads"]].to_parquet(s_fn)

    logging.info("Encoding Alleles")
    allele_df['sequence'] = allele_df[~allele_df['sequence'].isna()]['sequence'].apply(trgt.dna_encode)
    #np.where(allele_df['sequence'].isna(), None, allele_df['sequence'].apply(trgt.dna_encode)
    a_fn = os.path.join(args.dbname, f'allele.pq')
    allele_df[['LocusID', 'allele_number', 'allele_length', 'sequence']].to_parquet(pq_fns['allele'])
    logging.info("Parsed %d Alleles", len(allele_df))

    logging.info("Finished")
