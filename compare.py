#!/usr/bin/env python3
"""
Compare two runs pixel by pixel on a defined ROI

python compare.py --help

example:
  compare.py MAJA_20180415_NOENV.xml MAJA_20180415_ENV.xml -v -q -s --ulx 660240.0 --lry 4878660.0 --lrx 669240.0 --uly 4887660.0 --report
"""

__author__ = "Jerome Colin"
__license__ = "CC BY"
__version__ = "0.1.4"

import numpy as np
import os
import sys
import argparse
from majatools import majatools


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("runA", help="XML file describing run run_a")
    parser.add_argument("runB", help="XML file describing run run_a")
    parser.add_argument("-v", "--verbose", help="Set verbosity to INFO level", action="store_true")
    parser.add_argument("-r", "--report", help="Get a one line report of cloud free pixels and rmse per band",
                        action="store_true", default=False)
    parser.add_argument("-p", "--plots", help="Get scatterplots as PNG",
                        action="store_true", default=False)
    parser.add_argument("-q", "--quicklook", help="Get a png RGB quicklook",
                        action="store_true", default=False)
    parser.add_argument("-s", "--subset", help="Extract AOI subset", action="store_true", default=False)
    parser.add_argument("--ulx", type=float, help="AOI upper left longitude")
    parser.add_argument("--uly", type=float, help="AOI upper left latitude")
    parser.add_argument("--lrx", type=float, help="AOI lower right longitude")
    parser.add_argument("--lry", type=float, help="AOI lower right latitude")

    args = parser.parse_args()
    subset = args.subset
    verbose = args.verbose
    report = args.report
    quicklook = args.quicklook
    plots = args.plots

    try:
        if subset:
            ulx = args.ulx
            uly = args.uly
            lrx = args.lrx
            lry = args.lry
        else:
            print("WARNING: compare can't yet manage entire products, better use scatt instead. By !")
            sys.exit(0)

        if os.path.isfile(args.runA) and os.path.isfile(args.runB):
            f_run_a = args.runA
            f_run_b = args.runB
        else:
            print("ERROR: you have a typo in one XML file name, please check")
            sys.exit(1)

        compare(f_run_a, f_run_b, verbose, subset, ulx, uly, lrx, lry, report=report, plots=plots, quicklook=quicklook)

    except RuntimeError as e:
        print(e)
        sys.exit(1)

    sys.exit(0)


def compare(f_run_a, f_run_b, verbose, subset, ulx, uly, lrx, lry, report=False, plots=False, quicklook=False):
    # Creating instance of runs
    run_a = majatools.Run(f_run_a, verbosity=verbose)
    run_b = majatools.Run(f_run_b, verbosity=verbose)
    if verbose:
        print("INFO: Run run_a is of type %s, Run run_b is of type %s" % (run_a.get_type(), run_b.get_type()))

    # Load band mask
    edge_a = run_a.load_band(name="edge_mask", subset=subset, ulx=ulx, lry=lry, lrx=lrx, uly=uly)
    edge_b = run_b.load_band(name="edge_mask", subset=subset, ulx=ulx, lry=lry, lrx=lrx, uly=uly)

    # Load cloud mask
    clouds_a = run_a.load_band(name="cloud_mask", subset=subset, ulx=ulx, lry=lry, lrx=lrx, uly=uly)
    clouds_b = run_b.load_band(name="cloud_mask", subset=subset, ulx=ulx, lry=lry, lrx=lrx, uly=uly)

    # Load AOT and vap
    aot_a = run_a.load_band(name="aot", subset=subset, ulx=ulx, lry=lry, lrx=lrx, uly=uly)
    vap_a = run_a.load_band(name="vap", subset=subset, ulx=ulx, lry=lry, lrx=lrx, uly=uly)
    aot = np.nanmean(aot_a.band)
    vap = np.nanmean(vap_a.band)

    # Common pure pixels mask
    common_pure_pixels = clouds_a.band + clouds_b.band + edge_a.band + edge_b.band

    # Bands
    s2bands = ["B2", "B3", "B4", "B8"]

    if quicklook:
        red, green, blue = run_a.get_rgb(subset=subset, ulx=ulx, lry=lry, lrx=lrx, uly=uly)
        majatools.single_quicklook_rgb(red, green, blue, outfile=run_a.context)

    # Reporting
    if report:
        rmses = np.zeros(len(s2bands))

    for n in range(len(s2bands)):
        sre_a_rs = run_a.load_band(name="sre" + s2bands[n], subset=subset, ulx=ulx, lry=lry, lrx=lrx, uly=uly)
        sre_b_rs = run_b.load_band(name="sre" + s2bands[n], subset=subset, ulx=ulx, lry=lry, lrx=lrx, uly=uly)

        cloud_free_ratio, rmse = majatools.single_scatterplot(sre_a_rs, sre_b_rs, common_pure_pixels,
                                                              x_context=run_a.context, y_context=run_b.context,
                                                              mode="sre", png=plots)

        if report:
            rmses[n] = rmse

    if report:
        if cloud_free_ratio >= 0.5:
            if len(s2bands) == 4:
                print("REPORT, %s, %4.2f, %6.4f, %6.4f, %6.4f, %6.4f, %6.4f, %6.4f" % (str(run_a.get_timestamp()), cloud_free_ratio, rmses[0], rmses[1], rmses[2], rmses[3], aot, vap))


if __name__ == "__main__":
    main()
