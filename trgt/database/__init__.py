""" database commands """
from trgt.database.create import (
    create_main,
    vcf_to_tdb
)

from trgt.database.query import query_main

__all__ = [
    "create_main",
    "query_main"
]
