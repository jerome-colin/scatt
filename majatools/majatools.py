"""
Majatools PACKAGE

Recommended conda env:
python=3.6
conda install -c conda-forge gdal matplotlib scipy

Note: don't mix gdal packages from base and from forge
"""
__author__ = "Jerome Colin"
__license__ = "CC BY"
__version__ = "0.1.1"

import glob
import sys
import numpy as np
import matplotlib.pyplot as plt

from scipy import stats
from xml.dom import minidom
from osgeo import gdal


class Run:
    """
    Run object described by an XML context file, storing run attributes and providing an load_band method
    """
    def __init__(self, f_config, verbosity = False):
        """
        Initialization of Run
        :param f_config: XML context file passed as argument runA|runB
        :param verbosity: increase verbosity to INFO is True
        """
        try:
            self.xml_config = minidom.parse(f_config)
            self.verbosity = verbosity
            self.path = self.xml_config.getElementsByTagName("path")[0].firstChild.nodeValue
            self.type = self.xml_config.getElementsByTagName("type")[0].firstChild.nodeValue
            self.context = self.xml_config.getElementsByTagName("context")[0].firstChild.nodeValue
            self.scale_f_aot = float(self.xml_config.getElementsByTagName("scale_f_aot")[0].firstChild.nodeValue)
            self.scale_f_sr = float(self.xml_config.getElementsByTagName("scale_f_sr")[0].firstChild.nodeValue)
            self.nodata_aot = float(self.xml_config.getElementsByTagName("nodata_aot")[0].firstChild.nodeValue)

        except IndexError as e:
            print("ERROR: Missing parameter in XML file %s" % f_config)
            print(e)
            sys.exit(1)

    def get_type(self):
        """
        Return Run type
        :return: string of type of the run, typically "maja" or "maqt"
        """
        return self.type

    def load_band(self, name="aot"):
        """
        Get a given image output of Run
        :param name: any string between "aot", "cloud_mask", "edge_mask"
        :return: an Image object of the related band name
        """
        if name == "aot":
            if self.type == "maja":
                f_img = glob.glob(self.path + "*ATB_R1*")[0]
            elif self.type == "maqt":
                f_img = glob.glob(self.path + "ORTHO_SURF_AOT/*10m.tau2")[0]
            else:
                print("ERROR: Unable to find aot product...")
                sys.exit(1)

        if name == "cloud_mask":
            if self.type == "maja":
                f_img = glob.glob(self.path + "MASKS/*CLM_R1.tif")[0]
            elif self.type == "maqt":
                f_img = glob.glob(self.path + "MASQUES/NUAGES/*_10m.nua")[0]
            else:
                print("ERROR: Unable to find aot product...")
                sys.exit(1)

        if name == "edge_mask":
            if self.type == "maja":
                f_img = glob.glob(self.path + "MASKS/*_EDG_R1.tif")[0]
            elif self.type == "maqt":
                f_img = glob.glob(self.path + "MASQUES/BORD/*_10m.bord_final")[0]
            else:
                print("ERROR: Unable to find aot product...")
                sys.exit(1)

        if self.verbosity:
            print("INFO: Found file %s" % f_img.split("/")[-1])

        try:
            data = gdal.Open(f_img)
        except RuntimeError as e:
            print('ERROR: Unable to open ' + f_img)
            print(e)
            sys.exit(1)

        if name == "aot" and self.type == "maja":
            return Image(data.GetRasterBand(2).ReadAsArray(), self.type + " " + name, \
                         self.scale_f_aot, verbosity=self.verbosity)
        elif name == "aot" and self.type == "maqt":
            return Image(data.GetRasterBand(1).ReadAsArray(), self.type + " " + name, \
                         self.scale_f_aot, verbosity=self.verbosity)
        else:
            return Image(data.GetRasterBand(1).ReadAsArray(), self.type + " " + name, verbosity=self.verbosity)


class Image:
    """
    Extends a numpy array with attributes from Run, adds specifics methods
    """
    def __init__(self, band, band_name, scale_f = 1, nodata = -10000, verbosity = False):
        """
        Initialization of Image
        :param band: numpy array
        :param band_name: string name of the band
        :param scale_f: if any, scale factor of the variable in band
        :param nodata: no data value
        :param verbosity: increase verbosity to INFO is True
        """
        self.band = band
        self.band_name = band_name
        self.scale_f = scale_f
        self.nodata = nodata
        self.verbosity = verbosity

        # Changing nodata to NaN
        self.band = self.band.astype('float64')
        search = np.where(self.band == self.nodata)
        self.band[search] = np.nan

        # Apply scale factor
        self.band /= self.scale_f

    def get_finite(self, mask = None):
        """
        :param mask: optional filter to remove any non-zero pixels from self.band
        :return: a vector of finite values of self.band
        """
        if mask is not None:
            search = np.where(mask > 0)
            self.band[search] = np.nan

        return Image(np.extract(np.isfinite(self.band), self.band), self.band_name, verbosity=self.verbosity)

    def resample(self, n=6):
        """
        Downscaling method
        :param n: resampling factor
        :return: a resampled Image object, or self if n = 1
        """
        if n == 1:
            return self
        else:
            band_dim = len(self.band[:, 0])

            if np.mod(band_dim, n) == 0:
                size = int(band_dim / n)
                shift = int(n / 2)

                if self.verbosity:
                    print("INFO: Dimension of resampled_band is now %i x %i pixels" % (size, size))

                subset = np.zeros([size, size])

                for i in np.arange(0, len(subset)):
                    for j in np.arange(0, len(subset)):
                        subset[i, j] = self.band[i * n - shift, j * n - shift]

            else:
                print("Remainder is %i, change scaling..." % np.mod(band_dim, n))

            return Image(subset, self.band_name + " resampled", verbosity=self.verbosity)

    def apply_mask(self):
        pass


def scatterplot(a, b, \
                title = "demo", xt = "x", yt = "y", \
                f_savefig = "demo.png", mode='aot', show = False):
    """

    :param a: numpy array (1D or 2D)
    :param b: numpy array (1D or 2D)
    :param title: string of title
    :param xt: label of x axis
    :param yt: label of y axis
    :param f_savefig: filename to save the figure to
    :param mode: defines axis range, defaults to 'aot'
    :return:
    """
    slope, intercept, r_value, p_value, std_err = stats.linregress(a,b)

    fig, ax = plt.subplots(figsize=(6, 6))
    plt.title("%s (R2 = %6.4f)" % (title, r_value))

    if mode == 'sr':
        ax.plot([0, 0.5], [0, 0.5], 'k--')
    elif mode == 'aot':
        ax.plot([0, 0.6], [0, 0.6], 'k--')

    ax.set_aspect('equal', 'box')
    plt.xlabel(xt)
    plt.ylabel(yt)
    ax.plot(a, b, 'bo', markersize=2)
    plt.savefig(f_savefig, format='png')
    if show == True:
        plt.show()

