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
  AlleleID INTEGER PRIMARY KEY,
  LocusID INTEGER,
  allele_length INTEGER,
  /* FORMAT/MC, FORMAT/MS, FORMAT/AP */
  sequence BLOB,

  FOREIGN KEY (LocusID) REFERENCES Locus(LocusID)
);

CREATE TABLE SampleAlleleProperties (
  SampleID INTEGER,
  AlleleID INTEGER,
  spanning_reads INTEGER,

  FOREIGN KEY (SampleID) REFERENCES Sample(SampleID),
  FOREIGN KEY (AlleleID) REFERENCES Allele(AlleleID)
/*
  confidence_interval INTEGER,  AM 
*/
); 
