"""
Tandem repeat genotyping and visualization from PacBio HiFi data
"""
from trgt.bdna import (
    dna_decode,
    dna_decode_df,
    dna_encode
)

from trgt.database.dbutils import (
    dump_tdb,
    get_tdb_files,
    get_tdb_samplenames,
    load_tdb,
    tdb_consolidate,
    vcf_to_tdb,
)

from trgt.database.query import (
    allele_count,
    monref,
    gtmerge,
)
