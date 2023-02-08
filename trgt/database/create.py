import os
import sys
import argparse

import pysam
import truvari

import trgt

def get_sample(file):
    """
    Gets the sample name from vcf or tdb inputs
    """
    if file.endswith((".vcf", ".vcf.gz")):
        return pysam.VariantFile(vcf_fn).header.samples[0]
    return trgt.get_tdb_samplename(file)       

def check_args(args):
    """
    Preflight checks on arguments. Returns True if there is a problem
    """
    check_fail = False
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
            samp = get_sample(i)
            if samp in seen_samples:
                logging.error(f"input {i} has redundant sample with {seen_samples[i]}")
            seen_samples[samp] = i
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
    if os.path.exists(args.output):
        logging.error(f"output {dbname} already exists")
        sys.exit(1)
    
    if check_args(args):
        logging.error("cannot create database. exiting")
        sys.exit(1)

    logging.info("Loading %d files", len(args.inputs))
    m_data = None
    for i in inputs:
        n_data = trgt.load_tdb(i) if inputs[0].endswith(".tdb") else trgt.vcf_to_tdb(i)
        m_data = n_data if m_data is None else trgt.tdb_combine(m_data, n_data)

    pq_fns = get_tdb_files(dbname)
    logging.info("Writing parquet files")
    m_data['locus'].to_parquet(pq_fns['locus'])
    m_data['allele'].to_parquet(pq_fns['allele'])
    for sample, value in m_data['sample'].items():
        value.to_parquet(os.path.join(dbname, f"sample.{sample}.pq"))
    logging.info("Finished")
