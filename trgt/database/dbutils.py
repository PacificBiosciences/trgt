"""
Utilities for interacting with a trgt database
"""
import os
import glob
import pandas as pd


def get_tdb_files(dbname):
    """
    Return names of parquet table files in a tdb
    """
    l_fn = os.path.join(dbname, 'locus.pq')
    a_fn = os.path.join(dbname, 'allele.pq')

    snamestrip = lambda x: os.path.basename(x)[len('sample.'):-len('.pq')]
    s_dict = dict([(snamestrip(x), x) for x in glob.glob(os.path.join(dbname, 'sample.*.pq'))])

    return {'locus': l_fn, 'allele': a_fn, 'sample': s_dict}


def tdb_to_pd(dbname):
    """
    Loads  files into dataframes
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
