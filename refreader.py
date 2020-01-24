#!/usr/bin/env python3
"""
Demo of majatools.get_hdf_as_array()

python refreader.py --help

example:
    python refreader.py -v refsrs2-L1C_T31TFJ_A017022_20180925T104119-Carpentras.hdf
"""

__author__ = "Jerome Colin"
__license__ = "CC BY"
__version__ = "0.1.4"

import sys
import argparse
from majatools import majatools


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="HDF file")
    parser.add_argument("-v", "--verbose", help="Set verbosity to INFO level", action="store_true")
    args = parser.parse_args()

    majatools.get_hdf_as_array(args.file, verbose=args.verbose)

    sys.exit(0)


if __name__ == "__main__":
    main()
