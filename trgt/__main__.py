#!/usr/bin/env python3
"""
TRGT main entrypoint
"""
import sys
import argparse

import trgt
from trgt.dbcmds import db_main
from trgt.run import run_main

CMDS = {'run': run_main,
        'db': db_main}

USAGE = f"""\
TRGT v0.0.1 - Tandem repeat genotyping and visualization from PacBio HiFi data

Available commands:
    run         Run TRGT
    db          TRGT database commands
"""

def main():
    """
    Main entrypoint for TRGT
    """
    parser = argparse.ArgumentParser(prog="trgt", description=USAGE,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("cmd", metavar="CMD", choices=CMDS.keys(), type=str, default=None,
                        help="Command to execute")
    parser.add_argument("options", metavar="OPTIONS", nargs=argparse.REMAINDER,
                        help="Options to pass to the command")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit()
    args = parser.parse_args()

    CMDS[args.cmd](args.options)

if __name__ == '__main__':
    main()
