"""
tests the dna to bytes converter
"""
from itertools import product

import trgt

def test(seq):
    x = trgt.dna_encode(seq)
    y = trgt.dna_decode(x, len(seq))
    # print(seq, x, y)
    assert seq == y, f"{seq} != {y}"

def all_tests():
    test('ATCAGACAGG')
    test('AAAGAGAGA')
    test('CAATACACATACAGACAAGAGTTAG')

    def generate_kmers(chars, k):
        return [''.join(kmer) for kmer in product(chars, repeat=k)]

    for i in range(1, 6):
        for j in generate_kmers(['A', 'T', 'C', 'G'], i):
            test(j)

if __name__ == '__main__':
    all_tests()
