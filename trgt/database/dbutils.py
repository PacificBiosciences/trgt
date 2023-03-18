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


def load_tdb(dbname, samples=None, lfilters=None, afilters=None, sfilters=None):
    """
    Loads tdb into DataFrames
    returns dict of {'locus': DataFrame, 'allele': DataFrame, 'sample': {'sname': DataFrame}}
    
    If samples is provided, only a subset of sample tables are loaded.

    The (l)ocus, (a)llele, and (s)ample filters are passed to pyarrow.parquet.read_table
    filters during loading. Filters allow pulling subsets of data and have structure of
      List[Tuple] or List[List[Tuple]] or None (default)

    From their documentation:
      Each tuple has format: (key, op, value) and compares the key with the value. The
    supported op are: = or ==, !=, <, >, <=, >=, in and not in. If the op is in or not in,
    the value must be a collection such as a list, a set or a tuple.

    If a subset of loci are loaded via lfilters, then a filter of ('LocusID', 'in', loaded_locusids)
    is added to afilters and sfilters
    """
    def add_filter(filts, n_filt):
        """
        Add a new filter
        """
        if filts is None:
            return n_filt
        if isinstance(filts[0], tuple):
            return n_filt + filts
        filts.insert(0, n_filt)
        return filts
    names = get_tdb_files(dbname)
    ret = {}
    ret['locus'] = pq.read_table(names['locus'], filters=lfilters).to_pandas()
    if lfilters:
        loci = [("LocusID", "in", ret['locus']['LocusID'].values)]
        afilters = add_filter(afilters, loci)
        sfilters = add_filter(sfilters, loci)

    ret['allele'] = pq.read_table(names['allele'], filters=afilters).to_pandas()
    ret['allele']['sequence'] = ret['allele']['sequence'].apply(bytes)
    ret['sample'] = {}
    samp_to_fetch = samples if samples is not None else names["sample"].keys()
    for samp in samp_to_fetch:
        ret['sample'][samp] = pq.read_table(names['sample'][samp], filters=sfilters).to_pandas()
    return ret

def set_tdb_types(d):
    """
    Sets tdb datatypes of table columns in place
    """
    l_types = {"LocusID": np.uint32,
               "chrom": str,
               "start": np.uint32,
               "end": np.uint32}
    a_types = {"LocusID": np.uint32,
               "allele_number": np.uint16,
               "allele_length": np.uint16}
    s_types = {"LocusID": np.uint32,
               "allele_number": np.uint16,
               "spanning_reads": np.uint16,
               "length_range_lower": np.uint16,
               "length_range_upper": np.uint16}
    d['locus'] = d['locus'].astype(l_types)

    d['allele'] = d['allele'].astype(a_types)

    for samp, val in d['sample'].items():
        d['sample'][samp] = val.astype(s_types)

def dump_tdb(data, output):
    """
    Write tdb data to output folder

    WARNING: will overwrite existing data
    """
    if not os.path.exists(output):
        os.mkdir(output)
    pq_fns = get_tdb_files(output)
    set_tdb_types(data)
    data['locus'].to_parquet(pq_fns['locus'], index=False, compression='gzip')

    allele = pa.Table.from_pandas(data['allele'][["LocusID", "allele_number", "allele_length"]])
    seq = pa.Int8Array.from_pandas(data['allele']['sequence'], type=pa.list_(pa.uint8()))
    allele = allele.set_column(3, 'sequence', seq)
    a_schema = pa.schema([('LocusID', pa.uint32()),
                          ('allele_number', pa.uint16()),
                          ('allele_length', pa.uint16()),
                          ('sequence', pa.list_(pa.uint8()))
                         ])
    writer = pq.ParquetWriter(pq_fns['allele'], a_schema, compression='gzip')
    writer.write_table(allele)
    writer.close()

    for sample, value in data['sample'].items():
        value.to_parquet(os.path.join(output, f"sample.{sample}.pq"), index=False, compression='gzip')

def pull_alleles(data):
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
    alleles['sequence'] = alleles[~alleles['sequence'].isna()]['sequence'].apply(trgt.dna_encode)
    alleles = (alleles.sort_values(["LocusID", "allele_number"])
                    [["LocusID", "allele_number", "allele_length", "sequence"]]
                    .reset_index(drop=True))
    return alleles

def pull_saps(data, sample):
    """
    Turn sample allele properties into a table
    """
    # Remove sites with uninformative genotype information
    data = data[~data[f'{sample}_GT'].isin([(None,)])]
    gt = pd.DataFrame(data[f"{sample}_GT"].to_list(), columns=["GT1", "GT2"], index=data.index)
    span = pd.DataFrame(data[f"{sample}_SD"].to_list(), columns=["SD1", "SD2"], index=data.index)
    alci = pd.DataFrame(data[f"{sample}_ALLR"].to_list(), columns=["LR1", "LR2"], index=data.index)
    sap = pd.concat([data[["LocusID"]], span, alci, gt], axis=1)

    renamer = {"SD1": "spanning_reads", "SD2": "spanning_reads",
               "LR1":"LR", "LR2":"LR", "GT1":"allele_number", "GT2":"allele_number"}
    sap = pd.concat([sap[["LocusID", "GT1", "SD1", "LR1"]].rename(columns=renamer),
                     sap[["LocusID", "GT2", "SD2", "LR2"]].rename(columns=renamer)],
                     axis=0)
    sap[["length_range_lower", "length_range_upper"]] = sap["LR"].str.split('-', expand=True).astype(int)
    return sap.drop(columns=["LR"]).reset_index(drop=True)

def vcf_to_tdb(vcf_fn):
    """
    Turn a vcf into an in-memory TRGT database
    """
    if not os.path.exists(vcf_fn):
        raise RuntimeError(f"input {vcf_fn} does not exist")

    ret = {}
    data = truvari.vcf_to_df(vcf_fn, with_info=True, with_format=True, alleles=True)
    logging.info("locus count:\t%d", len(data))
    data["LocusID"] = range(len(data))

    ret["locus"] = data[["LocusID", "chrom", "start", "end"]].reset_index(drop=True).copy()  # pylint: disable=unsubscriptable-object # pylint/issues/3139

    logging.info("Wrangling alleles")
    allele_df = pull_alleles(data)
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
    ret = (union.reset_index()[["LocusID", "chrom", "start", "end"]]
                .sort_values(["chrom", "start", "end"])).copy()
    return ret, union

def allele_consolidator(exist_db, new_db, consol_locus):
    """
    Consolidate allele tables
    """
    new_locusids = dict(zip(consol_locus["LocusID_new"], consol_locus["LocusID"]))
    new_db["allele"]["LocusID"] = new_db["allele"]["LocusID"].map(new_locusids)

    for table in new_db["sample"].values():
        table["LocusID"] = table["LocusID"].map(new_locusids)

    ea = exist_db["allele"].set_index(["LocusID", "sequence", "allele_length"])
    na = new_db["allele"].set_index(["LocusID", "sequence", "allele_length"])
    new_allele = (ea.merge(na, how='outer', left_index=True, right_index=True)
                        .reset_index()
                        .sort_values(["LocusID", "allele_number_x", "allele_number_y"]))

    new_allele["allele_number"] = (new_allele.groupby(["LocusID"]).cumcount())

    allele_lookup = (new_allele[["LocusID", "allele_number_y", "allele_number"]]
                        .dropna()
                        .drop_duplicates()
                        .rename(columns={"allele_number": "n_an", "allele_number_y":"allele_number"})
                        .set_index(["LocusID", "allele_number"]))
    ret = new_allele[["LocusID", "allele_number", "allele_length", "sequence"]].copy()

    return ret, allele_lookup

def sample_consolidator(exist_db, new_db, allele_lookup):
    """
    Consolidate sample tables
    """
    ret = exist_db['sample']
    gt_count = len(exist_db['sample'])

    for sample, table in new_db["sample"].items():
        ret[sample] = (table.set_index(["LocusID", "allele_number"])
                        .join(allele_lookup)
                        .reset_index()
                        .drop(columns="allele_number")
                        .rename(columns={"n_an":"allele_number"})
                    )[["LocusID", "allele_number", "spanning_reads", "length_range_lower", "length_range_upper"]].copy()
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
