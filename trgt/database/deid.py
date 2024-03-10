"""
De-identify a TRGTdb
"""
import os
import sys
import shutil
import logging
import argparse
import truvari
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
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
    parser.add_argument("-S", "--shuffle-samples", action="store_true",
                        help="Shuffle sample tables together and split (experimental)")
    args = parser.parse_args(args)
    truvari.setup_logging()

    if check_args(args):
        logging.error("argument error. exiting")
        sys.exit(1)
    
    in_file_names = trgt.get_tdb_filenames(args.input)
    out_file_names = trgt.get_tdb_filenames(args.output)
    os.mkdir(args.output)
    logging.info("Copying locus")
    shutil.copy(in_file_names['locus'], out_file_names['locus'])

    alleles = pd.read_parquet(in_file_names['allele'])
    if args.remove_sequences:
        alleles["sequence"] = ""
    logging.info("Changing allele_numbers")
    ref = alleles[alleles["allele_number"] == 0].copy()
    alt = alleles[alleles["allele_number"] != 0].copy().sort_values(["LocusID", "allele_length"])
    alt["allele_number"] = alt.groupby(["LocusID"]).cumcount() + 1
    out = pd.concat([ref, alt]).sort_values(["LocusID", "allele_number"])
    out.to_parquet(out_file_names['allele'], index=False, compression='gzip')

    s_schema = pa.schema([('LocusID', pa.uint32()),
                          ('allele_number', pa.uint16()),
                          ('spanning_reads', pa.uint16()),
                          ('length_range_lower', pa.uint16()),
                          ('length_range_upper', pa.uint16()),
                          ('average_methylation', pa.uint16())
                        ])

    if args.shuffle_samples:
        logging.info("Shuffling samples (experimental)")
        parts = []
        n_samps = 0
        for i in in_file_names['sample'].values():
            parts.append(pd.read_parquet(i))
            n_samps += 1
        samps = pd.concat(parts).sample(frac=1).sort_values(["LocusID"])
        for idx in range(n_samps):
            o_fn = os.path.join(args.output, f"sample.{idx}.pq")
            value = samps.iloc[idx:len(samps):n_samps]
            # Awful code duplication. Need to separate out dump_tdb
            m_table = pa.Table.from_pandas(value)
            n_table = []
            n_names = []
            for col, dtype in zip(s_schema.names, s_schema.types):
                n_table.append(pa.compute.cast(m_table[col], dtype))
                n_names.append(col)
            n_table = pa.Table.from_arrays(n_table, names=n_names)
            writer = pq.ParquetWriter(o_fn, s_schema, compression='gzip')
            writer.write_table(n_table)
            writer.close()

    logging.info("Finished")
