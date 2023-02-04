"""
DNA to bytes encoder/decoder
"""
from io import BytesIO

NUCS = ['A', 'G', 'C', 'T']
def dna_encode(seq):
    """
    Turn a string of DNA to bytes
    """
    def set_bits(seq):
        """
        Turn upto 4 nucleotides to u8
        """
        ret = 0
        for pos, nuc in enumerate(seq):
            pos *= 2
            ret += NUCS.index(nuc) << pos
        return ret

    ret = BytesIO()
    for i in range(0, len(seq), 4):
        mbits = set_bits(seq[i:i+4])
        ret.write(mbits.to_bytes(1, 'big'))
    ret.seek(0)
    return ret.read()

def dna_decode(bstr, m_len):
    """
    Turn bytes to string of DNA
    """
    def miter(b, m_len):
        """
        Iterates bytestring to yield nucleotides
        """
        for i in b:
            for _ in range(min(4, m_len)):
                yield NUCS[i & 3]
                i >>=  2
                m_len -= 1
    return "".join(miter(bstr, m_len))

