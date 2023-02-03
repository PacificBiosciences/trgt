1) imports:
trgt run <>
trgt db create <in.vcf> <out.db>
trgt db append <in.vcf> <in.db>
trgt db concat <in.db> <in.db> # This is version 2
trgt db query --chrom --start --end

2) exports:
- Allele counts
	locus - allele - sample count

- Genotype of samples
	sample - locus - allele1 - probably FORMAT stuff
	sample - locus - allele2 - probably FORMAT stuff

- Motif copy numbers across alleles.
	locus - motifid - copy number - number of alleles

0 - CAG - 0 - 10
0 - CAG - 1 - 7
0 - CAG - 2 - 3
0 - CAG - 3 - 29
0 - CAG - 4 - 1



Translate a TRGT vcf into a duckdb

Tables
======

Locus
-----
A region where TRGT was run. This is tied to the input reference bed used by TRGT. 
Assumption: Loci columns are consistently present between different TRGT runs (i.e. are invariant).
	Therefore, when initializing a database, we'll populate these fields, and when appending to a database,
	most 'shared' loci between the two runs will already be present.

columns:
- primary key : UID for the locus
- Chromosome
- Start
- End

Motif
-----
At any given locus, I'm assuming there are potentially multiple INFO/MOTIFS and INFO/STRUC. However, the VCF defines `Number=1`
If there will only ever be a single MOTIF or STRUC per-locus, we can add these columns to the `Locus` table. But, if there are
multiple MOTIF or STRUC, we'll need two other tables.

A three column table describing a MOTIF. If STRUC is always `(MOTIF)n`, this becomes a two column table
- primary key : UID
- MOTIF : the motif sequence
- STRUC : structure of the motif

LocusMotif
----------
A 'joiner' table that would allow us to get the motifs of a locus. Again, this becomes unnecessary if there's only one 
Motif/Struc per-locus
- LocusID : UID from Locus table
- MotifID : UID from Motif table

This table can be removed if we never want to compare loci by motif. An example query would be "Give me every locus with a 
CGG motif"

Sample
------
Simple table of primary key to sample name. This could eventually be expanded to hold metadata of a sample, but
that's outside the scope of TRGT specific DB. For now, it will just be
- primary key : UID
- sample : name of the sample

Allele
------
For every locus, we have alternate alleles. Note that inside of this table, the allele has no relationship to a sample.
This is explained and resolved in the `SampleAlleleProperties` table.

- primary key : every allele will get its own UID.
- LocusID : UID from Locus table
- AL - length of the allele
- AlleleNumber : integer position of the allele (this is equivalent to the GT number e.g. (1/2). 
		The first non-ref allele will have AlleleNumber of 1 (0 is reserved for the REF)
		As a locus collects more alleles, this number increases e.g. if there's 4 Alleles (+REF) at a
		locus, then the `max(AlleleNumber) == 5` 
- Sequence - consensus sequence of the allele

A couple thoughts on reference alleles: We could put a LocusID -> AlleleNumber==0 with the REF sequence in this table.
However, we may not need to keep this as it could be fetched from a reference fasta. 

If the Allele:Sequence column gets too big, I have possible ideas on using variant_key to turn a sequence to a UID. Then
we'd only need to store the UID in this table and we can reference a separate file for lookups of this UID to a specific
sequence.

SampleAlleleProperties
----------------------
Details about properties of an allele for a sample. I'm separating this from Alelle because any given allele could be shared
between samples, but the details (e.g. SD) are dependent on the Sample.

- ALCI : length confidence interval of this sample's allele
- SD : spanning reads supporting the samples's allele

To finish this separation of a description of Allele from the properties of an Allele found in a Sample, I need to
better understand the following FORMAT fields. Are they descriptions of the Allele or the Sample's Allele Properties.
	MC - motif counts of the allele
	MS - motif spans of the allele 
	AP - allele purity
	AM - allele's methlyation - pretty sure this is SampleAlleleProperty

	
