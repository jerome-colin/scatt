# XML context generator for time series of Maja runs comparison

__author__ = "Jerome Colin"
__license__ = "CC BY"
__version__ = "0.1.0"

import os
import sys
import argparse
from majatools import majatools
from compare import compare


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="XML configuration file")
    parser.add_argument("-v", "--verbose", help="Set verbosity to INFO level", action="store_true")

    args = parser.parse_args()

    wrapper(args.config, verbose=args.verbose)

    sys.exit(0)


def wrapper(f_config, verbose=False, onlyref=False):
    timeseries = majatools.Timeseries(f_config, verbosity=verbose)

    timeseries.generate()

    for product in timeseries.common_product_list:
        context_1 = majatools.Context(timeseries.root_path + timeseries.collection_1_path + product, \
                                      timeseries.type, \
                                      timeseries.context_1, \
                                      timeseries.scale_f_aot, \
                                      timeseries.scale_f_sr, \
                                      timeseries.nodata_aot)

        context_2 = majatools.Context(timeseries.root_path + timeseries.collection_2_path + product, \
                                      timeseries.type, \
                                      timeseries.context_2, \
                                      timeseries.scale_f_aot, \
                                      timeseries.scale_f_sr, \
                                      timeseries.nodata_aot)

        reference = timeseries.find_reference(product)

        if reference is not None:
            print("Got %s as reference " % reference)
            reference_fullpath = timeseries.reference_collection_path + reference

        else:
            reference_fullpath = None

        try:
            compare(context_1, context_2, reference=reference_fullpath, verbose=verbose, subset=True, ulx=timeseries.subset_ulx,
                    uly=timeseries.subset_uly, lrx=timeseries.subset_lrx, lry=timeseries.subset_lry,
                    report=timeseries.report, plots=timeseries.plot, quicklook=timeseries.quicklook)

        except IndexError as e:
            print(e)
            pass


if __name__ == "__main__":
    main()
