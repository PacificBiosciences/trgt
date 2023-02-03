"""
DNA to bytes encoder/decoder
"""
from io import BytesIO

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
            pos = pos * 2
            if nuc == 'A':
                ret += 0 << pos
            elif nuc == 'T':
                ret += 1 << pos
            elif nuc == 'C':
                ret += 2 << pos
            elif nuc == 'G':
                ret += 3 << pos
            else:
                raise RuntimeError(f"Unrecognized nucleotide {nuc}")
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
    def check_bits(val):
        """
        Checks first two bits of u8 for nucleotide
        """
        if val & 1:
            if val & 2:
                return 'G'
            else:
                return 'T'
        elif val & 2:
            return 'C'
        else:
            return 'A'

    def miter(b):
        """
        Iterates bytestring to yield nucleotides
        """
        for i in b[:-1]:
            for _ in range(4):
                yield check_bits(i)
                i = i >> 2
        # last byte might be padded
        i = b[-1]
        pad_amt = m_len % 4
        pad_amt = 4 if pad_amt == 0 else pad_amt
        for _ in range(pad_amt):
            yield check_bits(i)
            i = i >> 2
            pad_amt += 1

    return "".join(miter(bstr))

