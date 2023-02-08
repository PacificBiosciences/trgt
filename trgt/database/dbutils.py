"""
Utilities for interacting with a trgt database
"""
import os
import glob
import logging

import pysam
import truvari
import pandas as pd

import trgt


def get_tdb_samplename(file):
    """
    Parses the sample name from a tdb sample parquet files
    """
    return os.path.basename(file)[len('sample.'):-len('.pq')]

def get_tdb_files(dbname):
    """
    Return names of parquet table files in a tdb
    """
    l_fn = os.path.join(dbname, 'locus.pq')
    a_fn = os.path.join(dbname, 'allele.pq')
    s_files = glob.glob(os.path.join(dbname, 'sample.*.pq'))
    s_names = [get_tdb_samplename(_) for _ in s_files]
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

def vcf_to_tdb(vcf_fn):
    """
    Turn a vcf into an in-memory TRGT database
    """
    if not os.path.exists(vcf_fn):
        raise RuntimeError(f"input {vcf_fn} does not exist")
    if os.path.exists(dbname):
        raise RuntimeError(f"output {dbname} already exists")

    ret = {}

    logging.info("Loading VCF %s", vcf_fn)
    data = truvari.vcf_to_df(vcf_fn, with_info=True, with_fmt=True, no_prefix=True, with_seqs=True)
    logging.info("Parsed %d loci", len(data))
    data['LocusID'] = range(len(data))

    ret['locus'] = data[["LocusID", "chrom", "start", "end"]].reset_index(drop=True).copy()

    logging.info("Wrangling alleles")
    allele_df = pull_alleles(data)
    logging.info("Encoding alleles")
    allele_df['sequence'] = allele_df[~allele_df['sequence'].isna()]['sequence'].apply(trgt.dna_encode)
    allele_df.to_parquet(pq_fns['allele'])
    logging.info("Parsed %d alleles", len(allele_df))

    logging.info("Pulling sample information")
    sap_df = pull_saps(data)
    sample_name = pysam.VariantFile(vcf_fn).header.samples[0]
    s_fn = os.path.join(dbname, f'sample.{sample_name}.pq')
    sap_df.to_parquet(s_fn)

def tdb_combine(exist_db, new_db):
    """
    Combine two tdbs
    """
    # Join Locus tables
    logging.info("Consolidating loci")
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
    logging.info("Total of %d loci", len(new_locus))

    # Will need to update the other tables' LocusID
    new_locusids = dict(zip(union["LocusID_new"], union["LocusID"]))
    new_db["allele"]["LocusID"] = new_db["allele"]["LocusID"].map(new_locusids)
    # Assume there's only one sample
    sample_name = list(new_db["sample"].keys())[0]
    new_db["sample"][sample_name]["LocusID"] = new_db["sample"][sample_name]["LocusID"].map(new_locusids)

    # Join allele tables
    logging.info("Consolidating alleles")
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
    logging.info("Total of %d alleles", len(new_allele))

    # Pass the new allele numbers to the new_db's sample
    logging.info("Consolidating sample information")
    union["allele_number_new"] = union["allele_number_new"].astype(int)
    lookup = (union[["LocusID", "allele_number_new", "n_an"]]
                .rename(columns={"allele_number_new":"allele_number"})
                .set_index(["LocusID", "allele_number"])["n_an"])
    nsi = new_db["sample"][sample_name].set_index(["LocusID", "allele_number"])
    new_sample = nsi.join(lookup, how='left').reset_index()
    new_sample["allele_number"] = new_sample["n_an"].fillna(new_sample["allele_number"]).astype(int)
    # Write new_sample[["LocusID", "allele_number", "spanning_reads", "ALCI_lower", "ALCI_upper"]]
