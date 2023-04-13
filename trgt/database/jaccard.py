"""
Methods to calculate jaccard indexes
"""
import sys
import itertools
from collections import Counter

def make_kmer_sets(seq, kmer_len=5, min_freq=5):
    """
    Makes the sets of all kmers and set of kmers over min_freq
    """
    kmers = Counter([seq[i:i + kmer_len] for i in range(len(seq) - kmer_len + 1)])
    kmers_freq = {k for k, c in kmers.items() if c >= min_freq}
    return set(kmers.keys()), kmers_freq

def jaccard_compare_kmers(kmers1, kmers1_freq, kmers2, kmers2_freq):
    """
    return the jacard similarity of two kmer sets
    """
    frequent = kmers1_freq | kmers2_freq

    k1_set = kmers1 & frequent
    k2_set = kmers2 & frequent

    intersection = len(k1_set & k2_set)
    union = len(k1_set | k2_set)

    return intersection / union if union else None

def jaccard_compare_seqs(seq1, seq2, kmer_len=5, min_freq=5):
    """
    return the jaccard similarity of two sequences
    """
    return jaccard_compare_kmers(*make_kmer_sets(seq1, kmer_len, min_freq),
                                 *make_kmer_sets(seq2, kmer_len, min_freq))

def alleles_jaccard_dist(alleles, counts, kmer_len=5, min_freq=5):
    """
    Given a list of alleles and their observed counts,
    return their mean jaccard distance
    """
    allele_cnt = len(alleles)
    all_kfeats = [make_kmer_sets(_, kmer_len, min_freq) for _ in alleles]
    dist_total = 0
    tot_pairs = 0
    for idx1, idx2 in itertools.combinations(range(allele_cnt), 2):
        dist = jaccard_compare_kmers(*all_kfeats[idx1], *all_kfeats[idx2])
        if dist is None:
            continue
        pair_cnt = counts[idx1] * counts[idx2]
        dist_total += dist * pair_cnt
        tot_pairs += pair_cnt

    return 1 - dist_total / tot_pairs if tot_pairs else None
