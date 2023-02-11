TRGTdb tutorial
===============


TRGT output vcfs can be collected into a database for easier querying.

The database comprises a set of parquet files:

- locus.pq: The locations of tandem repeats
-- LocusID: identifier for the locus
-- chrom: chromosome of the locus
-- start: 0-based start position of the locus
-- end: 0-based end position of the locus
- allele.pq: Tandem repeat alleles found on a locus
-- LocusID: locus to which the allele belongs
-- allele_number: identifier of the allele on the locus. '0' is the reference allele
-- allele_length: length of the allele
-- sequence: 2bit encoded sequence of the allele
- sample.\*.pq: sample properties of alleles
-- LocusID: identifier for the locus
-- allele_number: identifier of the allele on the locus. '0' is the reference allele
-- spanning_reads: number of spanning reads supporting per allele
-- length_range_lower: allele minimum predicted length
-- length_range_upper: allele maximum predicted length


Creating a database
===================
To create a database, we can use the following command:

```bash
trgt db create -o son.tdb son.vcf.gz
```

The output must end with `.tdb`. We can create a database from multiple vcfs or tdb files.
```bash
trgt db create -o family.tdb son.tdb father.vcf.gz
```

We can also append to an existing database.
```bash
trgt db append --to family.tdb --fr mother.vcf.gz
```

Querying a database
===================
Once we have a database created, we can then perform standard queries. 
- ac       : Locus - allele number - allele count
- as       : Locus - allele_number - sequence
- monref   : Monozygotic reference sites per-sample and overall
- gtmerge  : Collect per-locus genotypes
- metadata : Get table properties e.g. row counts and memory/disk sizes (mb)

By default, results are written as a tsv to stdout. The output can be customized with `-O` and `-o`.
```bash
trgt db query gtmerge family.tdb -O c -o genotypes.csv
```

Python library
==============

TRGT databases can be manipulated with python using the trgt library
```python
import trgt
```

Loading a TRGTdb is simply:
```python
data = trgt.load_tdb("family.tdb")
```

This returns a dictionary with structure
```
  { "locus": pandas.DataFrame,
    "allele": pandas.DataFrame,
    "sample": {
    	"son":pandas.DataFrame, 
    	"father":pandas.DataFrame, 
    	"mother":pandas.DataFrame, 
    }
  }
```

For example:
```python
>>> data['locus'].head()
   LocusID chrom  start    end
0        3  chr1  16682  16774
1        9  chr1  19275  19473
2        4  chr1  20798  20893
3        2  chr1  29714  29822
4        0  chr1  30824  30989
```

Pandas is very helpful for manipulating data. Here's an example query for collecting all loci with at least 3 alleles
```
loci_count = data['allele'].groupby(["LocusID"]).size()
min_3_alleles = loci_count[loci_count >= 3]

data['locus'].set_index("LocusID", inplace=True)
data['locus'].loc[min_3_alleles.index].head()
```

We can also run some of our queries programmatically:
```python
ac = trgt.allele_count(data)
```

When loading a tdb, we can pass filters to pyarrow and only load subsets of data.
```
chr1 = trgt.load_tdb("family.tdb",
		samples=['son', 'mother'],
		lfilters=[("chrom", "=", "chr1"), ("start", ">", 71685300), ("end", "<", 72497289)],
		afilters=[("allele_length", ">=", 1000)])
```
See `help(trgt.load_tdb)` for details on the filters.

To decode allele sequences to a string
```
sequences = data['allele'].apply(trgt.dna_decode_df, axis=1)
```
