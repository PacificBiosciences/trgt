import argparse
import truvari
import trgt

def create_main(args):
    """
    Create a new database and fill it with a VCF
    """
    parser = argparse.ArgumentParser(prog="trgt db create", description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("vcf", metavar="VCF", type=str,
                        help="Input TRGT VCF")
    parser.add_argument("dbname", metavar="DB", type=str,
                        help="Output TRGT DB")
    args = parser.parse_args(args)

    truvari.setup_logging()
    trgt.vcf_to_tdb(args.vcf, args.dbname)
    logging.info("Finished")
