"""
Standard, basic queries
"""
import argparse
import pandas as pd
import trgt

def allele_count(dbname):
    """
    Locus - allele number - allele count
    """
    data = trgt.load_tdb(dbname)

    # For a single sample, get how many times an allele is found
    allele_count = data['allele'][["LocusID", "allele_number"]].copy()
    allele_count['allele_count'] = 0
    allele_count = allele_count.set_index(["LocusID", "allele_number"])

    for samp_name, samp_data in data['sample'].items():
        num_samps = data['sample'][samp_name].reset_index().groupby(["LocusID", "allele_number"]).size()
        allele_count["allele_count"] += num_samps

    allele_count = allele_count.reset_index()
    view = data['locus'].join(allele_count, on='LocusID', rsuffix="_")
    view['allele_count'] = view['allele_count'].fillna(0)
    view[["chrom", "start", "end", "allele_number", "allele_count"]].to_csv('/dev/stdout', sep='\t', index=False)

def allele_seqs(dbname):
    """
    Locus - allele_number - sequence
    """
    tdb_fns = trgt.get_tdb_files(dbname)
    alleles = pd.read_parquet(tdb_fns["allele"])
    deseq = lambda x: trgt.dna_decode(x['sequence'], x['allele_length'])
    alleles['sequence'] = alleles[~alleles["sequence"].isna()].apply(deseq, axis=0)
    # Need fasta fetching for reference alleles
    alleles[["LocusID", "allele_number", "sequence"]].dropna().to_csv("/dev/stdout", sep='\t', index=False)

def monz_ref(dbname):
    """
    monozygotic reference per-sample and overall
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

QS = {
    "ac": ("Locus - allele number - allele count", allele_count),
    "as": ("Locus - allele number - allele sequence", allele_seqs),
    "monref": ("How many loci are monozygotic reference", monz_ref),
    # genotypes,
    # copy numbers...?
}

USAGE = "TRGT queries:\n" + "\n".join([f"    {k:9} {t[0]}" for k,t in QS.items()])

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
    QS[args.query][1](args.dbname)
