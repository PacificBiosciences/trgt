"""
De-identify a TRGTdb
"""
import os
import sys
import shutil
import logging
import argparse
import pandas as pd
import truvari
import trgt

def check_args(args):
    """
    Preflight checks on arguments. Returns True if there is a problem
    """
    check_fail = False
    if not os.path.exists(args.input):
        logging.error(f"input {args.input} does not exists")
        check_fail = True
    if not args.input.endswith(".tdb"):
        logging.error(f"input {args.input} must end with `.tdb`")
        check_fail = True
    if os.path.exists(args.output):
        logging.error(f"output {args.output} exists")
        check_fail = True
    if not args.output.endswith(".tdb"):
        logging.error(f"output {args.output} must end with `.tdb`")
        check_fail = True
    return check_fail

def deid_main(args):
    """
    Remove genotypes from a TRGTdb
    """
    parser = argparse.ArgumentParser(prog="trgt db append", description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
    # these could be positional arguments, but a little bit of user friction
    # will help prevent unintentional overwriting
    parser.add_argument("-i", "--input", metavar="DB", type=str,
                        help="Input DB")
    parser.add_argument("-o", "--output", metavar="FR", type=str,
                        help="Output DB")
    parser.add_argument("-s", "--remove-sequences", action="store_true",
                        help="Remove sequences from allele table")
    args = parser.parse_args(args)
    truvari.setup_logging()

    if check_args(args):
        logging.error("argument error. exiting")
        sys.exit(1)
    
    in_file_names = trgt.get_tdb_filenames(args.input)
    out_file_names = trgt.get_tdb_filenames(args.output)
    os.mkdir(args.output)
    shutil.copy(in_file_names['locus'], out_file_names['locus'])

    alleles = pd.read_parquet(in_file_names['allele'])
    if args.remove_sequences:
        alleles["sequence"] = ""
    ref = alleles[alleles["allele_number"] == 0].copy()
    alt = alleles[alleles["allele_number"] != 0].copy().sort_values(["LocusID", "allele_length"])
    alt["allele_number"] = alt.groupby(["LocusID"]).cumcount() + 1
    out = pd.concat([ref, alt]).sort_values(["LocusID", "allele_number"])
    out.to_parquet(out_file_names['allele'], index=False, compression='gzip')
    logging.info("Finished")
