"""
Basic queries on a tdb
"""
import os
import sys
import argparse
import joblib
import pandas as pd
import trgt

def allele_count(dbname):
    """
    Locus - allele number - allele count
    """
    data = trgt.load_tdb(dbname, decode=False)

    # For a single sample, get how many times an allele is found
    ac = data['allele'][["LocusID", "allele_number"]].copy()
    ac['allele_count'] = 0
    ac = allele_count.set_index(["LocusID", "allele_number"])

    for samp_data in data['sample'].values():
        num_samps = samp_data.reset_index().groupby(["LocusID", "allele_number"]).size()
        ac["allele_count"] += num_samps

    ac = allele_count.reset_index()
    view = data['locus'].join(ac, on='LocusID', rsuffix="_")
    view['allele_count'] = view['allele_count'].fillna(0)
    return view[["chrom", "start", "end", "allele_number", "allele_count"]]

def allele_seqs(dbname):
    """
    Locus - allele_number - sequence
    """
    tdb_fns = trgt.get_tdb_files(dbname)
    alleles = pd.read_parquet(tdb_fns["allele"])
    alleles['sequence'] = alleles.apply(trgt.dna_decode_df, axis=1)
    return alleles[["LocusID", "allele_number", "sequence"]].dropna()

def monz_ref(dbname):
    """
    Monozygotic reference sites per-sample and overall
    """
    data = trgt.load_tdb(dbname, decode=False)

    out_table = []
    for samp,table in data["sample"].items():
        table["is_ref"] = table["allele_number"] == 0
        out_table.append([samp,
                          table["LocusID"].nunique(),
                          table.groupby(["LocusID"])["is_ref"].all().sum()
                         ])

    # loci that are monozygotic across all samples
    all_sap = pd.concat(data["sample"].values())
    all_sap["is_ref"] = all_sap["allele_number"] == 0
    out_table.append(['all',
                      len(data['locus']),
                      all_sap.groupby(["LocusID"])["is_ref"].all().sum()
                     ])

    out_table = pd.DataFrame(out_table, columns=["sample", "loci", "mon_ref"])
    out_table['pct'] = out_table['mon_ref'] / out_table['loci']
    return out_table

def gtmerge(dbname):
    """
    Collect per-locus genotypes
    """
    data = trgt.load_tdb(dbname, decode=False)
    loci = data['locus'].set_index('LocusID')
    snames = {}
    gt_parts = []
    for idx, (samp, table) in enumerate(data['sample'].items()):
        table['gt_id'] = (table.groupby(["LocusID"]).cumcount())
        gts = (table.pivot(index="LocusID", columns="gt_id", values="allele_number")
                    .apply((lambda x: f"{x[0]:.0f}/{x[1]:.0f}"), axis=1))
        snames[idx] = samp
        gt_parts.append(gts)
    out = loci.join(pd.concat(gt_parts, axis=1, names=snames).fillna('./.'))
    return out.rename(columns=snames).sort_values(["chrom", "start", "end"])

def metadata(dbname):
    """
    Get table properties e.g. row counts and memory/disk sizes (mb)
    """
    def sizes(table, fname, df):
        dsize = round(os.path.getsize(fname) / 1.0e6, 1)
        msize = round(df.memory_usage().sum() / 1.0e6, 1)
        shape = df.shape[0]
        return [table, dsize, msize, shape]

    fnames = trgt.get_tdb_files(dbname)
    data = trgt.load_tdb(dbname, decode=False)
    header = ['table', 'disk', 'mem', 'rows']
    rows = [sizes("locus", fnames['locus'], data['locus']),
            sizes("allele", fnames['allele'], data['allele'])]
    for samp in fnames['sample']:
        rows.append(sizes(samp, fnames['sample'][samp], data['sample'][samp]))
    return pd.DataFrame(rows, columns=header)

QS = {"ac": allele_count,
      "as": allele_seqs,
      "monref": monz_ref,
      "gtmerge": gtmerge,
      "metadata": metadata,
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
        result.to_csv(args.output, sep='\t', index=False)
    elif args.output_type == 'c':
        result.to_csv(args.output, index=False)
    elif args.output_type == 'p':
        result.to_parquet(args.output)
    elif args.output_type == 'j':
        joblib.dump(result, args.output)
