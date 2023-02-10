"""
Utilities for interacting with a trgt database
"""
import os
import glob
import math
import logging

import pysam
import truvari
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

import trgt


def get_tdb_samplenames(file):
    """
    Parses the sample name from a tdb sample parquet files
    """
    ret = []
    for i in glob.glob(os.path.join(file, "sample.*.pq")):
        ret.append(os.path.basename(i)[len('sample.'):-len('.pq')])
    return ret

def get_tdb_files(dbname):
    """
    Return names of parquet table files in a tdb
    """
    l_fn = os.path.join(dbname, 'locus.pq')
    a_fn = os.path.join(dbname, 'allele.pq')
    s_files = glob.glob(os.path.join(dbname, 'sample.*.pq'))
    s_names = get_tdb_samplenames(dbname)
    s_dict = dict(zip(s_names, s_files))

    return {'locus': l_fn, 'allele': a_fn, 'sample': s_dict}


def load_tdb(dbname):
    """
    Loads files into dataframes
    return dict with keys of table name and values of the table(s)
    Note that 'sample' will have a value sample_name:DataFrame
    """
    names = get_tdb_files(dbname)
    ret = {}
    ret['locus'] = pd.read_parquet(names['locus'])
    ret['allele'] = pd.read_parquet(names['allele'])
    ret['sample'] = {}
    for samp, fn in names['sample'].items():
        ret['sample'][samp] = pd.read_parquet(fn)
    return ret

def set_tdb_types(d):
    """
    Sets tdb datatypes of table columns in place
    """
    l_types = {"LocusID": np.uint32,
               "chrom": str,
               "start": np.uint32,
               "end": np.uint32
               }
    a_types = {"LocusID": np.uint32,
               "allele_number": np.uint16,
               "allele_length": np.uint16,
               # "sequence": np.bytes_ this truncates for some reason
               }
    s_types = {"LocusID": np.uint32,
               "allele_number": np.uint16,
               "spanning_reads": np.uint16,
               "ALCI_lower": np.uint16,
               "ALCI_upper": np.uint16}
    d['locus'] = d['locus'].astype(l_types)

    d['allele'] = d['allele'].astype(a_types)

    for samp, val in d['sample'].items():
        d['sample'][samp] = val.astype(s_types)

def dump_tdb(data, output):
    """
    Write tdb data to output folder. output folder must already exist.

    WARNING: will overwrite existing data
    """
    if not os.path.exists(output):
        os.mkdir(output)
    pq_fns = get_tdb_files(output)
    set_tdb_types(data)
    # We have options here
    # See https://stackoverflow.com/questions/35789412/spark-sql-difference-between-gzip-vs-snappy-vs-lzo-compression-formats
    # Possibly https://arrow.apache.org/docs/python/generated/pyarrow.parquet.ParquetWriter.html
    # And help(df.to_parquet)
    data['locus'].to_parquet(pq_fns['locus'], index=False, compression='gzip')
    data['allele'].to_parquet(pq_fns['allele'], index=False, compression='gzip')
    for sample, value in data['sample'].items():
        value.to_parquet(os.path.join(output, f"sample.{sample}.pq"), index=False, compression='gzip')

def _dump_tdb_experimental(data, output):
    """
    Write tdb data to output folder. output folder must already exist.

    WARNING: will overwrite existing data
    """
    if not os.path.exists(output):
        os.mkdir(output)
    pq_fns = get_tdb_files(output)
    set_tdb_types(data)
    # Possibly https://arrow.apache.org/docs/python/generated/pyarrow.parquet.ParquetWriter.html
    # And help(df.to_parquet)
    data['locus'].to_parquet(pq_fns['locus'], index=False, compression='gzip')

    allele = pa.Table.from_pandas(data['allele'])
    a_schema = pa.schema([('LocusID', pa.uint32()),
                          ('allele_number', pa.uint16()),
                          ('allele_length', pa.uint16()),
                          # somehow the dtype gets messed up by.. consol_allele..?
                          # it was (a tiny bit) smaller, though
                          ('sequence', pa.Int8Array())
                         ])
    writer = pq.ParquetWriter(pq_fns['allele'], a_schema, compression='gzip')
    writer.write_table(allele)
    writer.close()

    #data['allele'].to_parquet(pq_fns['allele'], index=False, compression='gzip')
    for sample, value in data['sample'].items():
        value.to_parquet(os.path.join(output, f"sample.{sample}.pq"), index=False, compression='gzip')

def pull_alleles(data, encode=False):
    """
    Turn alleles into a table
    """
    alleles = (pd.DataFrame(data["alleles"].to_list(), index=data.index)
               .reset_index()
               .melt(id_vars='hash', value_name="sequence")
               .drop(columns="variable")
               .dropna()
               .set_index('hash'))
    alleles["LocusID"] = data["LocusID"]
    alleles["allele_number"] = alleles.groupby(["LocusID"]).cumcount()
    alleles = (alleles.sort_values(["LocusID", "allele_number"])
                    .reset_index(drop=True)
                    .drop_duplicates(subset=["LocusID", "sequence"]))
    alleles["allele_length"] = alleles["sequence"].str.len()
    alleles.loc[alleles["allele_number"] == 0, "sequence"] = ""
    alleles = (alleles.sort_values(["LocusID", "allele_number"])
                    [["LocusID", "allele_number", "allele_length", "sequence"]]
                    .reset_index(drop=True))
    if encode:
        #alleles['sequence'] = alleles['sequence'].where(~alleles['sequence'].isna(), [])
        alleles['sequence'] = alleles[~alleles['sequence'].isna()]['sequence'].apply(trgt.dna_encode)

    return alleles

def pull_saps(data, sample):
    """
    Turn sample allele properties into a table
    """
    gt = pd.DataFrame(data[f"{sample}_GT"].to_list(), columns=["GT1", "GT2"], index=data.index)
    span = pd.DataFrame(data[f"{sample}_SD"].to_list(), columns=["SD1", "SD2"], index=data.index)
    alci = pd.DataFrame(data[f"{sample}_ALCI"].to_list(), columns=["ALCI1", "ALCI2"], index=data.index)
    sap = pd.concat([data[["LocusID"]], span, alci, gt], axis=1)

    renamer = {"SD1": "spanning_reads", "SD2": "spanning_reads",
               "ALCI1":"ALCI", "ALCI2":"ALCI", "GT1":"allele_number", "GT2":"allele_number"}
    sap1 = sap[["LocusID", "GT1", "SD1", "ALCI1"]].rename(columns=renamer)
    sap2 = sap[["LocusID", "GT2", "SD2", "ALCI2"]].rename(columns=renamer)
    sap = pd.concat([sap1, sap2], axis=0)
    sap[["ALCI_lower", "ALCI_upper"]] = sap["ALCI"].str.split('-', expand=True).astype(int)
    sap = sap.drop(columns=["ALCI"])
    return sap.reset_index(drop=True)

def vcf_to_tdb(vcf_fn):
    """
    Turn a vcf into an in-memory TRGT database
    """
    if not os.path.exists(vcf_fn):
        raise RuntimeError(f"input {vcf_fn} does not exist")

    ret = {}
    data = truvari.vcf_to_df(vcf_fn, with_info=True, with_fmt=True, alleles=True)
    logging.info("locus count:\t%d", len(data))
    data["LocusID"] = range(len(data))

    ret["locus"] = data[["LocusID", "chrom", "start", "end"]].reset_index(drop=True).copy()  # pylint: disable=unsubscriptable-object # pylint/issues/3139

    logging.info("Wrangling alleles")
    allele_df = pull_alleles(data, encode=False)
    ret["allele"] = allele_df
    logging.info("allele count:\t%d", len(allele_df))

    logging.info("Pulling samples")
    ret["sample"] = {}
    gt_count = 0
    for sample in pysam.VariantFile(vcf_fn).header.samples:
        ret['sample'][sample] = pull_saps(data, sample)
        gt_count += len(ret['sample'][sample])
    logging.info("genotype count:\t%d", gt_count)
    return ret

def locus_consolidator(exist_db, new_db):
    """
    Consolidate Locus tables
    """
    el = exist_db["locus"].set_index(["chrom", "start", "end"])
    nl = new_db["locus"].set_index(["chrom", "start", "end"])
    union = el.join(nl, rsuffix='_new', how='outer').sort_values("LocusID")

    # Need new LocusIDs for anything that exists in new but not in original
    new_ids = np.array(range(len(el), len(el) + int(union["LocusID"].isna().sum())))
    union["LocusID"] = (np.hstack([union[~union["LocusID"].isna()]["LocusID"].values,
                        np.array(new_ids)]))
    union["LocusID_new"] = (union["LocusID_new"]
                            .fillna(-1)
                            .astype(int))
    union["LocusID"] = union["LocusID"].astype(int)
    ret = union.reset_index()[["LocusID", "chrom", "start", "end"]].copy()
    return ret, union

def allele_consolidator(exist_db, new_db, consol_locus):
    """
    Consolidate allele tables
    """
    new_locusids = dict(zip(consol_locus["LocusID_new"], consol_locus["LocusID"]))
    new_db["allele"]["LocusID"] = new_db["allele"]["LocusID"].map(new_locusids)

    for table in new_db["sample"].values():
        table["LocusID"] = table["LocusID"].map(new_locusids)

    ea = exist_db["allele"].set_index(["LocusID", "sequence"])
    na = new_db["allele"].set_index(["LocusID", "sequence"])
    new_allele = (ea.merge(na, how='outer', left_index=True, right_index=True)
                        .reset_index()
                        .sort_values(["LocusID", "allele_number_x", "allele_number_y"]))

    new_allele["allele_number"] = (new_allele.groupby(["LocusID"]).cumcount())
    new_allele["allele_length"] = (new_allele["allele_length_x"]
                                    .fillna(new_allele["allele_length_y"])
                                    .astype(int))
    allele_lookup = (new_allele[["LocusID", "allele_number_y", "allele_number"]]
                        .dropna()
                        .set_index(["LocusID", "allele_number_y"]))
    ret = new_allele[["LocusID", "allele_number", "allele_length", "sequence"]].copy()

    return ret, allele_lookup

def sample_consolidator(exist_db, new_db, allele_lookup):
    """
    Consolidate sample tables
    """
    ret = exist_db['sample']
    gt_count = len(exist_db['sample'])

    # Update new_db's allele numbers
    for sample, table in new_db["sample"].items():
        new_sample = table.set_index(["LocusID", "allele_number"])
        new_sample["allele_number"] = allele_lookup
        new_sample = new_sample.reset_index(level=0).reset_index(drop=True)
        ret[sample] = new_sample[["LocusID", "allele_number", "spanning_reads", "ALCI_lower", "ALCI_upper"]].copy()
        gt_count += len(ret[sample])
    return ret, gt_count

def tdb_consolidate(exist_db, new_db):
    """
    Combine two tdbs. returns new in-memory tdb
    """
    ret = {}
    logging.info("Consolidating loci")
    ret["locus"], consol_locus = locus_consolidator(exist_db, new_db)
    logging.info("locus count:\t%d", len(ret['locus']))

    logging.info("Consolidating alleles")
    ret["allele"], allele_lookup = allele_consolidator(exist_db, new_db, consol_locus)
    logging.info("allele count:\t%d", len(ret["allele"]))

    logging.info("Consolidating samples")
    ret['sample'], gt_count = sample_consolidator(exist_db, new_db, allele_lookup)
    logging.info("genotype count:\t%d", gt_count)

    return ret
