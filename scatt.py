#!/usr/bin/env python3
"""
Simple scatterplot tool to cross-compare MAJA 3 and Maquette outputs

usage: scatt.py [-h] [-v] [-r n] runA runB

positional arguments:
  runA                  XML file describing run A

  runB                  XML file describing run A

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Set verbosity to INFO level
  -r, --resampling      Integer n such that products are resampled to 1/n

Example:

    ./scatt.py M01_MAJA.xml M01_MAQT.xml -r 12

"""

__author__ = "Jerome Colin"
__license__ = "CC BY"
__version__ = "0.1.1"

import os
import sys
import argparse
from majatools import majatools


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("runA", help="XML file describing run A")
    parser.add_argument("runB", help="XML file describing run A")
    parser.add_argument("-v", "--verbose", help="Set verbosity to INFO level", action="store_true")
    parser.add_argument("-r", "--resampling", type=int, \
                        help="Integer n such that products are resampled to 1/n", default=0)

    args = parser.parse_args()

    # Checking resampling argument
    if args.resampling == 0:
        print("WARNING: full resolution processing will fill the RAM, changing to resampling = 6")
        scaling_n = 6
    else:
        scaling_n = args.resampling

    # Creating instance of runs
    if os.path.isfile(args.runA) and os.path.isfile(args.runB):
        A = majatools.Run(args.runA, verbosity=args.verbose)
        B = majatools.Run(args.runB, verbosity=args.verbose)
        if args.verbose:
            print("INFO: Run A is of type %s, Run B is of type %s" % (A.get_type(), B.get_type()))

    else:
        print("ERROR: Couldn't find required input file %s" % args.infile)
        sys.exit(1)

    # Load band mask
    edge_a = A.load_band(name="edge_mask").resample(n=scaling_n)
    edge_b = B.load_band(name="edge_mask").resample(n=scaling_n)

    # Load cloud mask
    clouds_a = A.load_band(name="cloud_mask").resample(n=scaling_n)
    clouds_b = B.load_band(name="cloud_mask").resample(n=scaling_n)

    # Common pure pixels
    common_pure_pixels = clouds_a.band + clouds_b.band + edge_a.band + edge_b.band

    # Scatterplot of AOT
    aot_a_rs = A.load_band(name="aot").resample(n=scaling_n).get_finite(mask=common_pure_pixels)
    aot_b_rs = B.load_band(name="aot").resample(n=scaling_n).get_finite(mask=common_pure_pixels)

    majatools.scatterplot(aot_a_rs.band, aot_b_rs.band, \
                          title=A.context + " " + A.type + " vs " + B.context + " " + B.type, \
                          xt=aot_a_rs.band_name, \
                          yt=aot_b_rs.band_name, \
                          f_savefig=A.context + "_" + aot_a_rs.band_name.replace(" ", "-") \
                                    + "_vs_" \
                                    + B.context + "_" + aot_b_rs.band_name.replace(" ", "-") \
                                    + ".png"
                          )

    sys.exit(0)

if __name__ == "__main__":
    main()
