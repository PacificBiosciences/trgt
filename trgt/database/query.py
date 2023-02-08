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
    data = trgt.tdb_to_pd(dbname)

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
    view[["chrom", "start", "end", "allele_number", "allele_count"]].to_csv('/dev/stdout', index=False)

def allele_seqs(dbname):
    """
    Locus - allele_number - sequence
    """
    tdb_fns = trgt.get_tdb_files(dbname)
    alleles = pd.read_parquet(tdb_fns["allele"])
    deseq = lambda x: trgt.dna_decode(x['sequence'], x['allele_length'])
    alleles['sequence'] = alleles[~alleles["sequence"].isna()].apply(deseq, axis=0)
    # Need fasta fetching for reference alleles
    alleles[["LocusID", "allele_number", "sequence"]].dropna().to_csv("/dev/stdout", index=False)

QS = {
    'ac': ("Locus - allele number - allele count", allele_count),
    'as': ("Locus - allele sequence", allele_seqs),
    # genotypes,
    # copy numbers
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

"""
these are operations on all the data
we'll have them for everything we 'load' i.e. columns from VCF that are pulled during 'create'

But we won't make them for anything 'derived' e.g. if we make the MotifCount table described elsewhere
allele_length ---
--region / --Region (I got this logic somewhere)

Would need the sample.*pq
spanning_reads --- 
--samples (inline csv or a file)
"""
