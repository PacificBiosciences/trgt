"""
Standard, basic queries
"""
import os
import argparse
import pandas as pd
import trgt

def allele_count(dbname):
    """
    Locus - allele number - allele count
    """
    data = trgt.load_tdb(dbname)

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
    view[["chrom", "start", "end", "allele_number", "allele_count"]].to_csv('/dev/stdout', sep='\t', index=False)

def allele_seqs(dbname):
    """
    Locus - allele_number - sequence
    """
    tdb_fns = trgt.get_tdb_files(dbname)
    alleles = pd.read_parquet(tdb_fns["allele"])
    def deseq(x):
        return trgt.dna_decode(x['sequence'], x['allele_length'])
    alleles['sequence'] = alleles[~alleles["sequence"].isna()].apply(deseq, axis=0)
    # Need fasta fetching for reference alleles
    alleles[["LocusID", "allele_number", "sequence"]].dropna().to_csv("/dev/stdout", sep='\t', index=False)

def monz_ref(dbname):
    """
    Monozygotic reference sites per-sample and overall
    """
    data = trgt.load_tdb(dbname)

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
    out_table.to_csv('/dev/stdout', sep='\t', index=False)

def gtmerge(dbname):
    """
    Collect per-locus genotypes
    """
    data = trgt.load_tdb(dbname)
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
    out.rename(columns=snames).sort_values(["chrom", "start", "end"]).to_csv("/dev/stdout", sep='\t', index=False)

# Can also pass parameters as **kwargs?
# Though that gets difficult to auto format if we want to expose them
# unless we override so that -h shows one-liner where --help shows full docs
# and --help in conjunction with a Q only shows the single query's full docs
def metadata(dbname):
    """
    Get table properties e.g. row counts and memory/disk sizes (mb)
    """
    fnames = trgt.get_tdb_files(dbname)
    data = trgt.load_tdb(dbname)
    print(f"table\tdisk\tmem\trows")
    dsize = round(os.path.getsize(fnames['locus']) / 1.0e6, 1)
    msize = round(data['locus'].memory_usage().sum() / 1.0e6, 1)
    shape = data['locus'].shape
    print(f"locus\t{dsize}\t{msize}\t{shape[0]}")

    dsize = round(os.path.getsize(fnames['allele']) / 1.0e6, 1)
    msize = round(data['allele'].memory_usage().sum() / 1.0e6, 1)
    shape = data['allele'].shape
    print(f"allele\t{dsize}\t{msize}\t{shape[0]}")
    for samp in fnames['sample']:
        dsize = round(os.path.getsize(fnames['sample'][samp]) / 1.0e6, 1)
        msize = round(data['sample'][samp].memory_usage().sum() / 1.0e6, 1)
        shape = data['sample'][samp].shape
        print(f"{samp}\t{dsize}\t{msize}\t{shape[0]}")
       
QS = {
    "ac": allele_count,
    "as": allele_seqs,
    "monref": monz_ref,
    "gtmerge": gtmerge,
    "metadata": metadata,
    # copy numbers...?
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
    args = parser.parse_args(args)
    QS[args.query](args.dbname)
