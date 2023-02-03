import argparse

import trgt.database as tdb

CMDS = {
    "create": ("Create a TRGTDB", tdb.create_main),
}

USAGE = "TRGT db:\n" + "\n".join([f"    {k:9} {t[0]}" for k,t in CMDS.items()])

def db_main(args):
    """
    Main entrypoint
    """
    parser = argparse.ArgumentParser(prog="trgt db", description=USAGE,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("cmd", metavar="CMD", choices=CMD.keys(), type=str,
                        help="db command to run")
    parser.add_argument("options", metavar="OPTIONS", nargs=argparse.REMAINDER,
                        help="Options to pass to the command")
    args = parser.parse_args(args)
    CMD[args.cmd][1](args.options)

