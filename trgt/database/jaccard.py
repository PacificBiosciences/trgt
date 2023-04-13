"""
Methods to calculate jaccard distances
"""
import sys
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

def jaccard_compare_seqs(seq1, seq2, kmer_len=5, min_freq=5):
    """
    return the jaccard distance of two sequences
    """
    kmers1, kmers1_freq = make_kmer_sets(seq1, kmer_len, min_freq)
    kmers2, kmers2_freq = make_kmer_sets(seq1, kmer_len, min_freq)

    return jaccard_compare_kmers(kmers1, kmers1_freq, kmers2, kmers2_freq)

def alleles_jaccard_dist(alleles, counts, kmer_len=5, min_freq=5):
    """
    Given a list of alleles and their observed counts,
    return their mean jaccard distance
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

    return 1 - dist_total / tot_pairs if tot_pairs else None

def test():
    import trgt
    import pandas as pd
    data = trgt.load_tdb("/Users/english/code/trgt/test_files/databases/all_vcf.tdb/",
                        lfilters=[("LocusID", "<", 100)])
    a_cnts = trgt.allele_count(data).reset_index().set_index(["LocusID", "allele_number"])
    allele = data['allele'].set_index(["LocusID", "allele_number"])
    a_cnts['sequence'] = allele['sequence']
    result = (a_cnts.reset_index()
               .groupby(['LocusID'])[["sequence", "AC"]]
               .apply(lambda x:
                       alleles_jaccard_dist(x["sequence"].values, x["AC"].values)))
    result.name = "muJI"
    print(pd.concat([data['locus'].set_index("LocusID"), result], axis=1))

if __name__ == '__main__':
    test()
