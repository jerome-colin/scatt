# XML context generator for time series of Maja runs comparison

__author__ = "Jerome Colin"
__license__ = "CC BY"
__version__ = "0.1.0"

import os
import sys
import argparse
from majatools import majatools


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="XML configuration file")
    parser.add_argument("-v", "--verbose", help="Set verbosity to INFO level", action="store_true")

    args = parser.parse_args()

    wrapper(args.config, verbose=args.verbose)

    sys.exit(0)

def wrapper(f_config, verbose=False):

    timeseries = majatools.Timeseries(f_config, verbosity=verbose)

    timeseries.generate()


if __name__ == "__main__":
    main()
