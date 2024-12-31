# Physkit - A Python toolkit for constants, unit conversions, and equations.
# Copyright (C) 2024 sapphimars
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.


import argparse
import sys
from .conversions import convert_unit
from .config import set_default
from .constants import constants


def main():
    """
    Entry point for the physkit CLI. Parses subcommands to then run the appropriate function.
    """
    parser = argparse.ArgumentParser(
        prog="physkit",
        description="A CLI for the physkit library to display constants or convert units",
    )

    subparsers = parser.add_subparsers(dest="command", help="Sub-command to run")

    # `physkit convert <value> <from_unit> <to_unit>`
    convert_parser = subparsers.add_parser("convert", help="Convert units")
    convert_parser.add_argument("value", type=float, help="Numeric value to convert")
    convert_parser.add_argument("from_unit", type=str, help="Current unit")
    convert_parser.add_argument("to_unit", type=str, help="Target unit")

    # `physkit constant <name> [--system ...]`
    convert_parser = subparsers.add_parser("constant", help="Value of constant")
    convert_parser.add_argument("name", type=str, help="Name of constant to display")
    convert_parser.add_argument(
        "--system",
        type=str,
        default=None,
        help="Optionally set unit system before converting (SI, cgs, etc.)",
    )
    # Parse the CLI args
    args = parser.parse_args()

    # If no subcommand given, print help and exit
    if not args.command:
        parser.print_help()
        sys.exit(1)

    if args.command == "convert":
        run_convert(args)
    elif args.command == "constant":
        run_constant(args)
    else:
        parser.print_help()


def run_convert(args):
    """
    Handle the 'physkit convert' subcommand.
    """

    # Do the conversion
    result = convert_unit(args.value, args.from_unit, args.to_unit)
    print(f"{args.value} {args.from_unit} = {result} {args.to_unit}")


def run_constant(args):
    """
    Handle the 'physkit constant' subcommand.
    """
    if args.system:
        # Change the global default system if specified
        set_default(args.system)

    # Attempt to retrieve the constant by attribute name
    try:
        value = getattr(constants, args.name)
    except AttributeError:
        print(f"No such constant '{args.name}' in the 'constants' class.")
        sys.exit(1)

    val, unit = value
    print(f"{args.name} = {val} {unit}")


# if __name__ == "__main__":
#    main()
