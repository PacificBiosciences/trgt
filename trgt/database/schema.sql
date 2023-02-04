CREATE TABLE Locus (
  LocusID INTEGER PRIMARY KEY,
  chrom VARCHAR,
  m_start INTEGER,
  m_end INTEGER,
);

CREATE TABLE Sample (
  SampleID INTEGER PRIMARY KEY,
  name VARCHAR
);

CREATE TABLE Allele (
  LocusID INTEGER,
  allele_number USMALLINT,
  allele_length USMALLINT,
  /* FORMAT/MC, FORMAT/MS, FORMAT/AP */
  sequence BLOB,

  FOREIGN KEY (LocusID) REFERENCES Locus(LocusID)
);

CREATE TABLE SampleAlleleProperties (
  SampleID INTEGER,
  LocusID INTEGER,
  allele_number USMALLINT,
  spanning_reads USMALLINT,

  FOREIGN KEY (SampleID) REFERENCES Sample(SampleID),
  FOREIGN KEY (LocusID) REFERENCES Locus(LocusID),
/*
  confidence_interval INTEGER,  AM 
*/
); 
