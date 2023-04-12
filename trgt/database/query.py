"""
Basic queries on a tdb
"""
import os
import sys
import argparse

import joblib
import numpy as np
import pandas as pd
import trgt

def tdb_opener(foo, *args, **kwargs):
    """
    Decorator for turning a tdb file name into a loaded tdb
    Allows queries to be used from command line or programmatically
    """
    def wrapper(data):
        if isinstance(data, str):
            data = trgt.load_tdb(data)
        return foo(data)
    wrapper.__doc__ = foo.__doc__
    return wrapper

@tdb_opener
def allele_count(data):
    """
    Allele counts and frequency
    """
    all_alleles = pd.concat([_[["LocusID", "allele_number"]] for _ in data['sample'].values()])
    lcnts = all_alleles["LocusID"].value_counts()
    acnts = (all_alleles.groupby(["LocusID"])["allele_number"]
                .value_counts()
                .rename("AC")
                .reset_index(level=1))
    acnts['AF'] = acnts['AC'] / lcnts
    return (data['locus'].set_index("LocusID")
                .join(acnts)
                .fillna(0)
                .astype({'allele_number':np.uint16, 'AC':np.uint16}))

@tdb_opener
def allele_count_length(data):
    """
    Allele counts and frequency by allele length
    """
    all_alleles = pd.concat([_[["LocusID", "allele_number"]] for _ in data['sample'].values()])
    lcnts = all_alleles["LocusID"].value_counts()
    all_alleles = all_alleles.set_index(["LocusID", "allele_number"])
    all_alleles["allele_length"] = data['allele'].set_index(['LocusID', "allele_number"])['allele_length']
    acnts = (all_alleles.reset_index()
                .groupby(["LocusID"])["allele_length"]
                .value_counts()
                .rename("AC")
                .reset_index(level=1))
    acnts['AF'] = acnts['AC'] / lcnts
    ret = (data['locus'].set_index("LocusID")
                .join(acnts)
                .fillna(0)
                .astype({'allele_length':np.uint16, 'AC':np.uint16}))
    # need a way to determine if an allele is the reference allele or not
    is_ref = (data['allele'].sort_values(["LocusID", "allele_length", "allele_number"])
              .drop_duplicates(["LocusID", "allele_length"])
              .set_index(["LocusID", "allele_length"])['allele_number'] == 0)
    ret = ret.reset_index().set_index(["LocusID", "allele_length"])
    ret['is_ref'] = is_ref
    return ret.reset_index(level=1)[["chrom", "start", "end", "is_ref", "allele_length", "AC", "AF"]]

def variant_length(allele_table):
    """
    Calculate variant length as allele's length minus locus' reference allele length
    """
    alleles = allele_table.set_index(["LocusID"])
    reflen = alleles[alleles["allele_number"] == 0]
    return (alleles["allele_length"].astype(int) - reflen["allele_length"].astype(int)).values

def allele_seqs(dbname):
    """
    Allele sequence, length, and difference from reference
    """
    tdb_fns = trgt.get_tdb_files(dbname)
    alleles = pd.read_parquet(tdb_fns["allele"])
    #alleles["sequence"] = alleles.apply(trgt.dna_decode_df, axis=1)
    alleles["ref_diff"] = variant_length(alleles)
    return alleles.reset_index()[["LocusID", "allele_number", "ref_diff", "sequence"]].dropna()

@tdb_opener
def monref(data):
    """
    Monozygotic reference sites per-sample and overall
    """
    out_table = []
    for samp,table in data["sample"].items():
        table["is_ref"] = table["allele_number"] == 0
        out_table.append([samp,
                          table["LocusID"].nunique(),
                          table.groupby(["LocusID"])["is_ref"].all().sum()
                         ])

    all_sap = pd.concat(data["sample"].values())
    all_sap["is_ref"] = all_sap["allele_number"] == 0
    out_table.append(['all',
                      len(data['locus']),
                      all_sap.groupby(["LocusID"])["is_ref"].all().sum()
                     ])

    out_table = pd.DataFrame(out_table, columns=["sample", "loci", "mon_ref"])
    out_table['pct'] = out_table['mon_ref'] / out_table['loci']
    return out_table

@tdb_opener
def gtmerge(data):
    """
    Collect per-locus genotypes
    """
    loci = data['locus'].set_index('LocusID')
    snames = {}
    gt_parts = []
    for idx, (samp, table) in enumerate(data['sample'].items()):
        table['gt_id'] = (table.groupby(["LocusID"]).cumcount())
        gts = (table.pivot(index="LocusID", columns="gt_id", values="allele_number")
                    .apply((lambda x: f"{x[0]:.0f}/{x[1]:.0f}"), axis=1))
        snames[idx] = samp
        gt_parts.append(gts)
    out = loci.join(pd.concat(gt_parts, axis=1, names=snames)).fillna('./.')
    return out.rename(columns=snames).sort_values(["chrom", "start", "end"])

def metadata(dbname):
    """
    Get table properties e.g. row counts and memory/disk sizes (mb)
    """
    def sizes(table, fname, df):
        denom = 1.0e6 #mega
        dsize = round(os.path.getsize(fname) / denom, 2)
        msize = round(df.memory_usage().sum() / denom, 2)
        shape = df.shape[0]
        return [table, dsize, msize, shape]

    fnames = trgt.get_tdb_files(dbname)
    data = trgt.load_tdb(dbname)
    header = ['table', 'disk', 'mem', 'rows']
    rows = [sizes("locus", fnames['locus'], data['locus']),
            sizes("allele", fnames['allele'], data['allele'])]
    for samp in fnames['sample']:
        rows.append(sizes(samp, fnames['sample'][samp], data['sample'][samp]))
    return pd.DataFrame(rows, columns=header)

@tdb_opener
def methyl(data):
    """
    Allele length, methylation, and CpG stats (PMID:3656447)
    """
    def cpg_stats(seq):
        #seq = df[0]
        obs = seq.count("CG")
        c_count = seq.count("C")
        g_count = seq.count("G")
        exp = (c_count * g_count) / len(seq) if len(seq) else 0
        density = obs * 2 / len(seq) if len(seq) else 0
        gc_pct = (c_count + g_count) / len(seq) if len(seq) else 0
        return obs, exp, density, gc_pct
    allele = data['allele'].set_index(["LocusID", "allele_number"])
    new_cols = ["CpG_obs", "CpG_exp", "CpG_density", "GC_pct"]
    allele["CpG_obs"], allele["CpG_exp"], allele["CpG_density"], allele["GC_pct"] = zip(*allele['sequence'].apply(cpg_stats))

    parts = []
    for samp in data['sample']:
        parts.append((data['sample'][samp]
                        .set_index(["LocusID", "allele_number"])
                        .join(allele)
                        .reset_index()
                        [["LocusID", "allele_number", "allele_length", "average_methylation"] + new_cols]
                    ))

    return pd.concat(parts)

QS = {"allele_cnts": allele_count,
      "allele_cnts_bylen": allele_count_length,
      "allele_seqs": allele_seqs,
      "monref": monref,
      "gtmerge": gtmerge,
      "metadata": metadata,
      "methyl": methyl,
}

USAGE = "TRGT queries:\n" + "\n".join([f"    {k:9}: {t.__doc__.strip()}" for k,t in QS.items()])

def query_main(args):
    """
    Main entrypoint
    """
    parser = argparse.ArgumentParser(prog="trgt db", description=USAGE,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("query", metavar="Q", choices=QS.keys(), type=str,
                        help="query to run")
    parser.add_argument("dbname", metavar="TDB", type=str,
                        help="TRGT db name")
    parser.add_argument("--index", action="store_true",
                        help="Write index with tsv/csv outputs")
    parser.add_argument("-o", "--output", type=str, default='/dev/stdout',
                        help="Output destination (stdout)")
    parser.add_argument("-O", "--output-type", default='t', choices=['t', 'c', 'p', 'j'],
                        help="Output type of [t]sv, [c]sv, [p]arquet, [j]oblib (%(default)s)")
    args = parser.parse_args(args)

    if args.output_type == 'p' and os.path.exists(args.output):
        sys.stderr.write("Cannot write parquet to existing file\n")
        sys.exit(1)

    result = QS[args.query](args.dbname)

    if args.output_type == 't':
        result.to_csv(args.output, sep='\t', index=args.index)
    elif args.output_type == 'c':
        result.to_csv(args.output, index=args.index)
    elif args.output_type == 'p':
        result.to_parquet(args.output)
    elif args.output_type == 'j':
        joblib.dump(result, args.output)
