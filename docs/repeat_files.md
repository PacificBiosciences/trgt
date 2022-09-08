# Repeat definitions

The repeat definition files are BED files containing repeat coordinates and
information about the repeat structure. The first three columns specify the
coordinates of the repeat region. The fourth column contains the following
mandatory fields:

- Repeat identifier (`ID`). The identifier must be globally unique.
- The comma-separated list of repeat motifs (`MOTIFS`) specifying all possible
  motifs that can appear in this region. In particular, each motif present in
  the `STRUC` field must also be present in this field.
- Repeat region structure (`STRUC`). This field describes the overall structure
  of the repeat region. It is composed of one or more tandem repeats, each
  defined by the expression `(<motif>)n`. Tandem repeats can optionally
  separated by interrupting sequences. See the example below.

For example, consider this definition of tandem repeat:

```bash
chr4  3074876  3074966  ID=HTT,MOTIFS=CAG,CCG;STRUC=(CAG)nCAACAG(CCG)n
```

This definition tells us that the repeat region has coordinates
chr4:3074876-3074966. Its identifier is HTT. The region contains two tandem
repeats with motifs CAG and CCG and these tandem repeats are expected to
be separated by a short interrupting sequence CAACAG.
