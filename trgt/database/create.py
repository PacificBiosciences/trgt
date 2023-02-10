"""
Create a tdb
"""
import os
import sys
import logging
import argparse

import pysam
import truvari

import trgt

def get_samples(file):
    """
    Gets the sample name from vcf or tdb inputs
    """
    if file.endswith((".vcf", ".vcf.gz")):
        return list(pysam.VariantFile(file).header.samples)
    return trgt.get_tdb_samplenames(file)

def check_args(args):
    """
    Preflight checks on arguments. Returns True if there is a problem
    """
    check_fail = False

    if os.path.exists(args.output):
        logging.error(f"output {args.output} already exists")
        check_fail = True
    if not args.output.endswith(".tdb"):
        logging.error(f"output {args.output} must end with `.tdb`")
        check_fail = True
    seen_samples = {}
    for i in args.inputs:
        if not os.path.exists(i):
            logging.error(f"input {i} does not exist")
            check_fail = True
        if not i.endswith((".vcf", ".vcf.gz", ".tdb")):
            logging.error("unrecognized file extension on {i}")
            logging.error("expected .vcf .vcf.gz or .tdb")
            check_fail = True
        else: # can only check sample of valid file names
            for s in get_samples(i):
                if s in seen_samples:
                    logging.error(f"input {i} has redundant sample with {seen_samples[i]}")
                seen_samples[s] = i
    return check_fail

def create_main(args):
    """
    Create a new trgt.db from multiple input calls
    """
    parser = argparse.ArgumentParser(prog="trgt db create", description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", metavar="OUT", required=True,
                        help="Output trgt.db directory")
    parser.add_argument("inputs", metavar="IN", nargs="+",
                        help="VCF or trgt.db files")
    args = parser.parse_args(args)

    truvari.setup_logging()
    if check_args(args):
        logging.error("cannot create database. exiting")
        sys.exit(1)

    m_data = None
    for i in args.inputs:
        logging.info("Loading %s", i)
        n_data = trgt.load_tdb(i) if i.endswith(".tdb") else trgt.vcf_to_tdb(i)
        m_data = n_data if m_data is None else trgt.tdb_consolidate(m_data, n_data)

    logging.info("Writing parquet files")
    trgt.dump_tdb(m_data, args.output)
    logging.info("Finished")
