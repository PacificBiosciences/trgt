#!/bin/env python

import sys
import gzip
import itertools
import numpy as np
from numpy.typing import ArrayLike
from collections import Counter, defaultdict

# original method
def get_kmers(seq1, seq2, kmer_len = 5, min_freq = 5):
    def collect(seq):
        kmers = []
        index = 0
        while index + kmer_len <= len(seq):
            kmer = seq[index: index + kmer_len]
            kmers.append(kmer)
            index += 1
        return kmers

    kmers1 = collect(seq1)
    kmers2 = collect(seq2)

    frequent1 = set(k for k, c in Counter(kmers1).items() if c >= min_freq)
    frequent2 = set(k for k, c in Counter(kmers2).items() if c >= min_freq)
    frequent = frequent1.union(frequent2)

    kmers1 = set(kmers1).intersection(frequent)
    kmers2 = set(kmers2).intersection(frequent)

    return kmers1, kmers2

def compare(seq1, seq2, kmer_len=5, min_freq=5):
    kmers1, kmers2 = get_kmers(seq1, seq2, kmer_len, min_freq)

    if len(kmers1) + len(kmers2) == 0:
        return None
    intersection = kmers1.intersection(kmers2)
    union = kmers1.union(kmers2)
    return len(intersection) / len(union)


def get_dists(alleles):
    dists = []
    for index1 in range(len(alleles)):
        for index2 in range(index1 + 1, len(alleles)):
            seq1, count1 = alleles[index1]
            seq2, count2 = alleles[index2]

            dist = compare(seq1, seq2)
            if dist == None:
                continue
            dists.extend([dist] * (count1 * count2))
    return dists

def test_original():
    fn = "test_comp.txt"
    with open(fn, "r") as file:
        for locus, alleles in itertools.groupby(file, key=lambda rec: rec.split()[0]):
            locus = locus#.decode("utf8")
            alleles = [a.strip().split()[1:] for a in alleles]
            alleles = [a for a in alleles if len(a) == 2]
            alleles = [(seq, int(count)) for count, seq in alleles]
            #alleles = [(seq, count) for seq, count in alleles if count > 1]
            # Don't forget this minimum allele length check
            alleles = [(seq, count) for seq, count in alleles if len(seq) >= 10]
            if len(alleles) == 0:
                continue

            if len(alleles) == 1:
                print(locus, 0.0, 1.0, sep="\t")
                continue

            median_len = np.median([len(seq) for seq, count in alleles])
            dists = get_dists(alleles)
            mean_dist = np.mean(dists) if len(dists) > 0 else "NA"
            print(locus, median_len, mean_dist, sep="\t")

def old_method(alleles, counts):
    median_len = np.median([len(seq) for seq in alleles])
    dists = get_dists([_ for _ in zip(alleles, counts)])
    mean_dist = np.mean(dists) if len(dists) > 0 else "NA"
    return median_len, mean_dist

def make_kmer_sets(seq, kmer_len=5, min_freq=5):
    """
    Makes the sets of all kmers and set of kmers over min_freq
    """
    kmers = Counter([seq[i:i + kmer_len] for i in range(len(seq) - kmer_len + 1)])
    kmers_freq = {k for k, c in kmers.items() if c >= min_freq}
    return set(kmers.keys()), kmers_freq

def jaccard_compare_seqs(seq1, seq2, kmer_len=5, min_freq=5):
    """
    return the jaccard distance of two sequences
    """
    kmers1, kmers1_freq = make_kmer_sets(seq1, kmer_len, min_freq)
    kmers2, kmers2_freq = make_kmer_sets(seq1, kmer_len, min_freq)

    return jaccard_compare_kmers(kmers1, kmers1_freq, kmers2, kmers2_freq)

def jaccard_compare_kmers(kmers1, kmers1_freq, kmers2, kmers2_freq):
    """
    return the jacard distance of two kmer featurization arrays
    """
    # Will be True for freq1.union(freq2)
    frequent = kmers1_freq | kmers2_freq

    # Intersection of kmers
    k1_idx = kmers1 & frequent
    k2_idx = kmers2 & frequent

    intersection = len(k1_idx & k2_idx)
    union = len(k1_idx | k2_idx)

    return intersection / union if union else None

def alleles_jaccard_dist(alleles, counts, kmer_len=5, min_freq=5):
    """
    Given a list of alleles and their observed counts,
    return their mean jaccard dist
    """
    allele_cnt = len(alleles)
    # Only kfeaturize once - only compare high frequency alleles
    all_kfeats = [make_kmer_sets(_, kmer_len, min_freq) for _ in alleles]
    dist_total = 0
    tot_pairs = 0
    for idx1 in range(allele_cnt - 1):
        for idx2 in range(idx1 + 1, allele_cnt):
            dist = jaccard_compare_kmers(*all_kfeats[idx1], *all_kfeats[idx2])
            if dist is None:
                continue
            pair_cnt = counts[idx1] * counts[idx2]
            dist_total += dist * pair_cnt
            tot_pairs += pair_cnt

    return dist_total / tot_pairs if tot_pairs else None

def test():
    # Assumes a single locus
    fn = "test_comp.txt"
    fn = "delme.txt"
    alleles = []
    counts = []
    with open(fn, "r") as fh:
        for line in fh:
            data = line.strip().split('\t')
            counts.append(int(data[1]))
            alleles.append(data[2])
    print('orig', old_method(alleles, counts))
    print('new ', alleles_jaccard_dist(alleles, counts))

def test2():
    fn = "/Users/english/code/trgt/test_files/databases/hprc_105.tdb"#all_vcf.tdb/"
    import trgt
    data = trgt.load_tdb(fn, lfilters=[("LocusID", "<", 2000)])
    a_cnts = trgt.allele_count(data).reset_index().set_index(["LocusID", "allele_number"])
    allele = data['allele'].set_index(["LocusID", "allele_number"])
    a_cnts['sequence'] = allele['sequence']
    view = a_cnts.reset_index().groupby(['LocusID'])[["sequence", "AC"]]
    import time

    then = time.time()
    dists_o = view.apply(lambda x: old_method(x["sequence"].values, x["AC"].values))
    print('old', time.time() - then)

    then = time.time()
    dists_n = (view.apply(lambda x: alleles_jaccard_dist(x["sequence"].values, x["AC"].values)))
    print('new', time.time() - then)

    #print(a_cnts[['sequence', 'AC']])
    import joblib
    joblib.dump([dists_o, dists_n], 'test_jd2.jl')

if __name__ == '__main__':
    test2()
