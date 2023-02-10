"""
DNA to bytes encoder/decoder
"""
from io import BytesIO

NUCS = {'A':0, 'G':1, 'C':2, 'T':3}
NUCSi = list(NUCS.keys())
def dna_encode(seq):
    """
    Turn a string of DNA to bytes
    """
    ret = BytesIO()
    for i in range(0, len(seq), 4):
        #mbits = sum([NUCS.index(nuc) << pos * 2 for pos,nuc in enumerate(seq[i:i+4])])
        byte = 0
        for pos, nuc in enumerate(seq[i:i+4]):
            byte += NUCS[nuc] << pos * 2
        ret.write(byte.to_bytes(1, 'big'))
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
                yield NUCSi[i & 3]
                i >>=  2
                m_len -= 1
    return "".join(miter(bstr, m_len))
