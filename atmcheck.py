#!/usr/bin/env python3
"""
Simple scatterplot tool to cross-compare SRE = f(TOA,AOT) MAJA 3 and Maquette outputs

usage: atmcheck.py [-h] [-v] [-r] run

positional arguments:
  run                  XML file describing run

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Set verbosity to INFO level
  -r, --resampling      Integer n such that products are resampled to 1/n

Example:

    ./atmcheck.py M01_MAJA.xml -r 12

"""

__author__ = "Jerome Colin"
__license__ = "CC BY"
__version__ = "0.1.3"

import os
import sys
import argparse
from majatools import majatools


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run", help="XML file describing run")
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
    if os.path.isfile(args.run):
        A = majatools.Run(args.run, verbosity=args.verbose)
        if args.verbose:
            print("INFO: Run is of type %s" % (A.get_type()))

    else:
        print("ERROR: you have a typo in one XML file name, please check")
        sys.exit(1)

    # Load band mask
    edge_a = A.load_band(name="edge_mask").resample(n=scaling_n)

    # Load cloud mask
    clouds_a = A.load_band(name="cloud_mask").resample(n=scaling_n)

    # Common pure pixels
    #common_pure_pixels = clouds_a.band + edge_a.band
    common_pure_pixels = edge_a.band

    # Scatterplot of AOT
    aot_a_rs = A.load_band(name="aot").resample(n=scaling_n).get_finite(mask=common_pure_pixels)

    # Scatterplots of SRE
    for s2band in ("B2", "B3", "B4"):
        sre_a_rs = A.load_band(name="sre" + s2band).resample(n=scaling_n).get_finite(mask=common_pure_pixels)
        toa_a_rs = A.load_band(name="toa" + s2band).resample(n=scaling_n).get_finite(mask=common_pure_pixels)

        majatools.atmplot(toa_a_rs.band, sre_a_rs.band, aot_a_rs.band, \
                              title=A.context + " " + A.type, \
                              xt=A.context + " " + toa_a_rs.band_name, \
                              yt=A.context + " " + sre_a_rs.band_name, \
                              f_savefig=A.context + "_" + toa_a_rs.band_name.replace(" ", "-") \
                                        + "_vs_" \
                                        + A.context + "_" + sre_a_rs.band_name.replace(" ", "-") \
                                        + ".png", \
                              mode="sre"
                              )


    sys.exit(0)


if __name__ == "__main__":
    main()
