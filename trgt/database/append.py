"""
Add new VCF to existing database
"""
import os
import argparse
import pandas as pd
import truvari
import logging
import trgt
import numpy as np
import math

def append_main(args):
    """
    Create a new database and fill it with a VCF
    """
    parser = argparse.ArgumentParser(prog="trgt db create", description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dbname", metavar="DB", type=str,
                        help="Existing TRGT DB")
    parser.add_argument("vcf", metavar="VCF", type=str,
                        help="Input TRGT VCF")
    args = parser.parse_args(args)
    truvari.setup_logging()

    if not os.path.exists(args.dbname):
        raise RuntimeError(f"output {args.dbname} does not exists")
    if not os.path.exists(args.vcf):
        raise RuntimeError(f"input {args.vcf} does not exist")

    tmp_name = truvari.make_temp_filename(suffix=".tdb")
    trgt.vcf_to_tdb(args.vcf, tmp_name)
    new_db = trgt.tdb_to_pd(tmp_name)
    exist_db = trgt.tdb_to_pd(args.dbname)

    tdb_combine()
    logging.info("Finished")
