# Introductory tutorial

In this tutorial, you will learn the basics of TRGT and TRVZ by analyzing a
tiny example dataset included in this repository.

## Prerequisites

- Download the [latest TRGT and TRVZ binaries](https://github.com/PacificBiosciences/trgt/releases)
- Download the [tiny example dataset](https://github.com/PacificBiosciences/trgt/tree/main/example)
- Install recent versions of `samtools` and `bcftools`

## Overview of the input files

Our example dataset consists of a reference genome (`reference.fasta`),
a file with aligned reads (`sample.bam`), and a list of repeats (`repeats.bed`).
The list of repeats is a BED file that, among other information, contains repeat
coordinates, repeat identifiers, and motifs:

```bash
$ cat repeats.bed
chrA   26039021        26039053        ID=chr10_26039021_26039053,STRUC=(TTTG)n
chrA   26041338        26041362        ID=chr10_26041338_26041362,STRUC=(TTG)n
chrA   26041683        26041699        ID=chr10_26041683_26041699,STRUC=(AAAC)n
```

## Genotype repeats

To genotype the repeats, run:

```bash
./trgt example/reference.fasta \
       example/repeats.bed \
       example/sample.bam \
       sample
```

The output consists of files `sample.bcf` and `sample.spanning.bam`. The BCF
file (which is the binary version of a VCF file) contains repeat genotypes while
the BAM file contains pieces of HiFi reads that fully span the repeat sequences.

For example, here is the first entry of the BCF file:

```bash
$ bcftools view --no-header sample.bcf | head -n 1
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample
chrA    3001    .       CAGCAG  CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG,CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG  0       .       TRID=StrA;STRUC=(CAG)n  GT:MC:TSD:RM    1/2:30/40:29:.
```

It says that:

- There is a tandem repeat starting at position 3001 of chrA
- The reference sequence of this repeat is CAGCAG
- This repeat has two non-reference alleles (spanning 30 and 40 CAGs)

## Sort and index the outputs

TRGT outputs are not sorted. So you need to sort and index the BCF:

```bash
bcftools sort -Ob -o sample.sorted.bcf sample.bcf
bcftools index sample.sorted.bcf
```

and then the BAM:

```bash
samtools sort -o sample.spanning.sorted.bam sample.spanning.bam
samtools index sample.spanning.sorted.bam
```

And that's it! The output files are now ready for downstream analysis.

## Visualize a repeat

To visualize the repeat with the identifier "StrA", run:

```bash
./trvz example/reference.fasta \
       example/repeats.bed \
       sample.sorted.bcf \
       sample.spanning.sorted.bam \
       StrA
```

TRVZ outputs two files `StrA.svg` and `StrA.png` that contain the same read
pileup image in two file formats. Note that the SVG file can be directly edited
in vector graphics editing software like [Inkscape](https://inkscape.org/).

The resulting pileup plot shows the sequences of reads spanning each repeat
allele (blue) and the surrounding flanking sequence (green):

![StrA read pileup](figures/StrA.png)
