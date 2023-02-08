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
    Turn alleles into a table
    """
    alleles = pd.DataFrame(data["ALT"].to_list(), columns=["ALT1", "ALT2"], index=data.index)
    gt = pd.DataFrame(data["GT"].to_list(), columns=["GT1", "GT2"], index=data.index)
    alleles = pd.concat([data[["LocusID", "REF"]], alleles, gt], axis=1)
    alleles = alleles.melt(id_vars=["LocusID", "REF", "ALT1", "ALT2"], value_vars=["GT1", "GT2"], value_name="allele_number")

    def seq_chooser(x):
        """
        What sequence does the genotype point to
        """
        if x["allele_number"] == 0:
            return x["REF"]
        if x["allele_number"] == 1:
            return x["ALT1"]
        if x["allele_number"] == 2:
            return x["ALT2"]

    alleles['sequence'] = alleles.apply(seq_chooser, axis=1)
    alleles = alleles.drop(columns=["REF", "ALT1", "ALT2", "variable"])
    alleles.insert(2, "allele_length", alleles['sequence'].str.len())
    return alleles

def pull_saps(data):
    """
    Turn sample allele properties into a table
    """
    gt = pd.DataFrame(data["GT"].to_list(), columns=["GT1", "GT2"], index=data.index)
    span = pd.DataFrame(data["SD"].to_list(), columns=["SD1", "SD2"], index=data.index)
    alci = pd.DataFrame(data["ALCI"].to_list(), columns=["ALCI1", "ALCI2"], index=data.index)
    sap = pd.concat([data[["LocusID"]], span, alci, gt], axis=1)

    renamer = {"SD1": "spanning_reads", "SD2": "spanning_reads",
               "ALCI1":"ALCI", "ALCI2":"ALCI", "GT1":"allele_number", "GT2":"allele_number"}
    sap1 = sap[["LocusID", "GT1", "SD1", "ALCI1"]].rename(columns=renamer)
    sap2 = sap[["LocusID", "GT2", "SD2", "ALCI2"]].rename(columns=renamer)
    sap = pd.concat([sap1, sap2], axis=0)
    sap[["ALCI_lower", "ALCI_upper"]] = sap["ALCI"].str.split('-', expand=True).astype(int)
    sap = sap.drop(columns=["ALCI"])
    return sap.reset_index(drop=True)

def vcf_to_tdb(vcf_fn, dbname):
    """
    Turn a vcf into a TRGT database
    """
    if not os.path.exists(vcf_fn):
        raise RuntimeError(f"input {vcf_fn} does not exist")
    if os.path.exists(dbname):
        raise RuntimeError(f"output {dbname} already exists")

    os.mkdir(dbname)

    logging.info("Loading VCF")
    data = truvari.vcf_to_df(vcf_fn, with_info=True, with_fmt=True, no_prefix=True, with_seqs=True)
    logging.info("Parsed %d Loci", len(data))
    data['LocusID'] = range(len(data))

    pq_fns = trgt.get_tdb_files(dbname)
    data[["LocusID", "chrom", "start", "end"]].reset_index(drop=True).to_parquet(pq_fns['locus'])

    logging.info("Wrangling Alleles")
    allele_df = pull_alleles(data)
    logging.info("Encoding Alleles")
    allele_df['sequence'] = allele_df[~allele_df['sequence'].isna()]['sequence'].apply(trgt.dna_encode)
    a_fn = os.path.join(dbname, f'allele.pq')
    allele_df.to_parquet(pq_fns['allele'])
    logging.info("Parsed %d Alleles", len(allele_df))

    logging.info("Pulling Sample Information")
    sap_df = pull_saps(data)
    sample_name = pysam.VariantFile(vcf_fn).header.samples[0]
    s_fn = os.path.join(dbname, f'sample.{sample_name}.pq')
    sap_df.to_parquet(s_fn)
    

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
    vcf_to_tdb(args.vcf, args.dbname)
    logging.info("Finished")
