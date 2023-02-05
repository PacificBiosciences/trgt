import argparse

import trgt.database as tdb

CMDS = {
    "create": ("Create a TRGTDB", tdb.create_main),
    "query": ("Query a TRGTDB", tdb.query_main)
}

USAGE = "TRGT db commands:\n" + "\n".join([f"    {k:9} {t[0]}" for k,t in CMDS.items()])

def db_main(args):
    """
    Main entrypoint
    """
    parser = argparse.ArgumentParser(prog="trgt db", description=USAGE,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("cmd", metavar="CMD", choices=CMDS.keys(), type=str,
                        help="db command to run")
    parser.add_argument("options", metavar="OPTIONS", nargs=argparse.REMAINDER,
                        help="Options to pass to the command")
    args = parser.parse_args(args)
    CMDS[args.cmd][1](args.options)

