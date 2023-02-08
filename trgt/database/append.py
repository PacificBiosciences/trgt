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

    # Join Locus tables
    logging.info("Consolidating Loci")
    el = exist_db["locus"].set_index(["chrom", "start", "end"])
    nl = new_db["locus"].set_index(["chrom", "start", "end"])
    union = el.join(nl, rsuffix='_new', how='outer').sort_values("LocusID")

    # Need new LocusIDs for anything that exists in the new but not in the original
    new_ids = np.array(range(len(el), len(el) + int(union["LocusID"].isna().sum())))
    union["LocusID"] = np.hstack([union[~union["LocusID"].isna()]["LocusID"].values, np.array(new_ids)])
    union["LocusID_new"] = union["LocusID_new"].fillna(-1).astype(int)
    union["LocusID"] = union["LocusID"].astype(int)
    # Safe to write this out now
    new_locus = union.reset_index()[["LocusID", "chrom", "start", "end"]]
    logging.info("Total of %d Loci", len(new_locus))

    # Will need to update the other tables' LocusID
    new_locusids = dict(zip(union["LocusID_new"], union["LocusID"]))
    new_db["allele"]["LocusID"] = new_db["allele"]["LocusID"].map(new_locusids)
    # Assume there's only one sample
    sample_name = list(new_db["sample"].keys())[0]
    new_db["sample"][sample_name]["LocusID"] = new_db["sample"][sample_name]["LocusID"].map(new_locusids)

    # Join allele tables
    logging.info("Consolidating Alleles")
    ea = exist_db["allele"].set_index(["LocusID", "sequence"])
    na = new_db["allele"].set_index(["LocusID", "sequence"])
    new_allele = (ea.join(na, rsuffix='_new', how='outer')
                    .reset_index()
                    .sort_values(["LocusID", "allele_number", "allele_number_new"])
                    .drop_duplicates())

    # Identical Locus/sequence will use existing allele_numbers
    # New alleles in the locus will need next available allele_numbers
    def allele_num_consolidate(x):
        l_val = x.iloc[0]["allele_number"]
        l_val = int(x.iloc[0]["allele_number_new"]) if math.isnan(l_val) else int(l_val)
        return np.arange(l_val, l_val + len(x), dtype='int')
    new_allele["n_an"] = np.hstack(new_allele.groupby(["LocusID"]).apply(allele_num_consolidate)).astype(int)
    new_allele["allele_length"] = new_allele["allele_length"].fillna(new_allele["allele_length_new"]).astype(int)
    # hold this for samples
    union = new_allele.dropna(subset="allele_number_new").copy()
    new_allele["allele_number"] = new_allele["n_an"]
    # Write new_allele[["LocusID", "allele_number", "allele_length", "sequence"]]
    logging.info("Total of %d Alleles", len(new_allele))

    # Pass the new allele numbers to the new_db's sample
    logging.info("Consolidating Sample Information")
    union["allele_number_new"] = union["allele_number_new"].astype(int)
    lookup = (union[["LocusID", "allele_number_new", "n_an"]]
                .rename(columns={"allele_number_new":"allele_number"})
                .set_index(["LocusID", "allele_number"])["n_an"])
    nsi = new_db["sample"][sample_name].set_index(["LocusID", "allele_number"])
    new_sample = nsi.join(lookup, how='left').reset_index()
    new_sample["allele_number"] = new_sample["n_an"].fillna(new_sample["allele_number"]).astype(int)
    # Write new_sample[["LocusID", "allele_number", "spanning_reads", "ALCI_lower", "ALCI_upper"]]
    logging.info("Finished")
