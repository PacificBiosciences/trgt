#!/usr/bin/env python3
"""
TRGT main entrypoint
"""
import sys
import argparse

from trgt.dbcmds import db_main

def run_main(args):
    """
    placeholder
    """
    print(f"Would run trgt executable + {args}")

CMDS = {'run': run_main,
        'db': db_main}

USAGE = """\
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