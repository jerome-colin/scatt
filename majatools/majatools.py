"""
Majatools PACKAGE

Recommended conda env:
python=3.6
conda install -c conda-forge gdal matplotlib scipy

Note: don't mix gdal packages from base and from forge
"""
__author__ = "Jerome Colin"
__license__ = "CC BY"
__version__ = "0.2.1"

import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
import os

from xml.dom import minidom
from osgeo import gdal


class Aoi:
    """
    Area of interest
    """

    def __init__(self, x_center, y_center, px_range):
        """
        Create an Aoi instance
        :param x_center: central pixel coordinate along the horizontal
        :param y_center: central pixel coordinate along the vertical
        :param px_range: half-height / half-width range of the square Aoi
        """

        self.x_center = x_center
        self.y_center = y_center
        self.px_range = px_range


class Image:
    """
    Extends a numpy array with attributes from Run, adds specifics methods
    """

    def __init__(self, band, band_name, scale_f=1, nodata=-10000, verbosity=False):
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
        self.x_dim, self.y_dim = self.__get_band_size()

        # Changing nodata to NaN
        self.band = self.band.astype('float64')
        search = np.where(self.band == self.nodata)
        self.band[search] = np.nan

        # Apply scale factor
        self.band /= self.scale_f

    def get_finite(self, mask=None):
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

    def subset_aoi(self, aoi):
        """
        Extract square image subset by pixel coordinates
        :param aoi: aoi instance
        :return: an Image instance
        """
        return Image(self.__crop_band_by_aoi(aoi), self.band_name, scale_f=1, nodata=self.nodata)

    def __check_is_inside(self, aoi):
        """
        Returns True if an AOI is within the extent of an image
        :param aoi:
        :return: Boolean
        """
        if (aoi.x_center > self.x_dim) or (aoi.y_center > self.y_dim):
            print("ERROR: central coordinates of AOI outside of the image")
            return False

        if (aoi.x_center - aoi.px_range < 0) or (aoi.x_center + aoi.px_range > self.x_dim):
            print("ERROR: AOI extends outside of the image")
            return False

        if (aoi.y_center - aoi.px_range < 0) or (aoi.y_center + aoi.px_range > self.y_dim):
            print("ERROR: AOI extends outside of the image")
            return False

        return True

    def __crop_band_by_aoi(self, aoi):
        """
        Crops an array according to AOI extent
        :param aoi: an Aoi instance
        :return: an array
        """
        if self.__check_is_inside(aoi):
            return self.band[aoi.y_center - aoi.px_range:aoi.y_center + aoi.px_range, \
                   aoi.x_center - aoi.px_range:aoi.x_center + aoi.px_range]

    def __get_band_size(self):
        """
        Private, returns band dimensions
        :return:
        """
        try:
            return len(self.band[0, :]), len(self.band[:, 0])

        except IndexError:
            return len(self.band), 0

    def apply_mask(self, mask, reverse=False):
        """
        Apply a given mask to self, replacing non-zeros by NaN
        :param mask:
        """
        if reverse:
            search = np.where(mask.band == 0)
            self.band[search] = np.NaN
        else:
            search = np.where(mask.band != 0)
            self.band[search] = np.NaN


class Run:
    """
    Run object described by an XML context file, storing run attributes and providing a load_band method
    """

    def __init__(self, f_config, verbosity=False):
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
            print("ERROR: Mandatory parameter missing in XML file %s" % f_config)
            print(e)
            sys.exit(1)

        # Optional parameters
        try:
            self.dtm_path = self.xml_config.getElementsByTagName("dtm_path")[0].firstChild.nodeValue
            self.water_path = self.xml_config.getElementsByTagName("water_path")[0].firstChild.nodeValue


        except:
            pass

    def get_type(self):
        """
        Return Run type
        :return: string of type of the run, typically "maja" or "maqt"
        """
        return self.type

    def load_band(self, name="aot", subset=False, ulx=0, uly=0, lrx=0, lry=0):
        """
        Get a given image output of Run, optionally a subset by coordinates

        :param name: any string between "aot", "cloud_mask", "edge_mask"
        :param subset: if True use gdal_translate to subset an AOI
        :param uly: upper left y
        :param ulx: upper left x
        :param lrx: lower right x
        :param lry: lower right y
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

        if name[:3] == "sre":
            if self.type == "maja":
                f_img = glob.glob(self.path + "*SRE_" + name[-2:] + "*")[0]
            elif self.type == "maqt":
                f_img = glob.glob(self.path + "ORTHO_SURF_AOT/*10m." + name[-2:])[0]
            else:
                print("ERROR: Unable to find SRE product...")
                sys.exit(1)

        if name[:3] == "toa":
            if self.type == "maja":
                print("ERROR: Not yet implemented...")
            elif self.type == "maqt":
                f_img = glob.glob(self.path + "ORTHO_TOA_ABS/*10m." + name[-2:])[0]
            else:
                print("ERROR: Unable to find TOA product...")
                sys.exit(1)

        if name == "cloud_mask":
            if self.type == "maja":
                f_img = glob.glob(self.path + "MASKS/*CLM_R1.tif")[0]
            elif self.type == "maqt":
                f_img = glob.glob(self.path + "MASQUES/NUAGES/*_10m.nua")[0]
            else:
                print("ERROR: Unable to find cloud_mask product...")
                sys.exit(1)

        if name == "edge_mask":
            if self.type == "maja":
                f_img = glob.glob(self.path + "MASKS/*_EDG_R1.tif")[0]
            elif self.type == "maqt":
                f_img = glob.glob(self.path + "MASQUES/BORD/*_10m.bord_final")[0]
            else:
                print("ERROR: Unable to find edge_mask product...")
                sys.exit(1)

        if name == "water mask":
            try:
                f_img = glob.glob(self.water_path + "*_10m.water")[0]
            except:
                print("WARNING: Unable to find WATER MASK product for %s in %s" % (self.type, self.water_path))

        if name == "dtm":
            try:
                f_img = glob.glob(self.dtm_path + "*10mfloat.mnt")[0]
            except:
                print("WARNING: Unable to find DTM product for %s..." % self.type)

        try:
            data = gdal.Open(f_img)
        except RuntimeError as e:
            print('ERROR: Unable to open ' + f_img)
            print(e)
            sys.exit(1)

        if subset:
            try:
                translate = 'gdal_translate -projwin %s %s %s %s %s %s' % (ulx, uly, lrx, lry, f_img, ".tmp.tif")
                os.system(translate)
                data = gdal.Open(".tmp.tif")
                os.remove(".tmp.tif")

            except RuntimeError as e:
                print('ERROR: Unable to open ' + f_img)
                print(e)
                sys.exit(1)

        if self.verbosity:
            print("INFO: Found file %s" % f_img.split("/")[-1])

        if name == "aot" and self.type == "maja":
            return Image(data.GetRasterBand(2).ReadAsArray(), self.type + " " + name, \
                         scale_f=self.scale_f_aot, verbosity=self.verbosity)
        elif name == "aot" and self.type == "maqt":
            return Image(data.GetRasterBand(1).ReadAsArray(), self.type + " " + name, \
                         scale_f=self.scale_f_aot, verbosity=self.verbosity)
        elif name[:3] == "sre":
            return Image(data.GetRasterBand(1).ReadAsArray(), self.type + " " + name, \
                         scale_f=self.scale_f_sr, verbosity=self.verbosity)
        elif name[:3] == "toa":
            return Image(data.GetRasterBand(1).ReadAsArray(), self.type + " " + name, \
                         scale_f=self.scale_f_sr, verbosity=self.verbosity)
        else:
            return Image(data.GetRasterBand(1).ReadAsArray(), self.type + " " + name, verbosity=self.verbosity)


def diffmap(a, b, mode, with_dtm=False):
    """
    Produce an absolute difference map as image
    :param a: a numpy 2D array
    :param b: b numpy 2D array
    :param mode: either "aot" or "sreXX" where XX refers to a band (eg. B3)
    """

    var_a = a.load_band(name=mode)
    var_b = b.load_band(name=mode)

    edge_a = a.load_band(name="edge_mask")
    edge_b = b.load_band(name="edge_mask")

    if with_dtm:
        try:
            dtm = a.load_band(name="dtm")
            water = a.load_band(name="water mask")
            dtm.apply_mask(water)

        except:
            try:
                dtm = b.load_band(name="dtm")
                water = b.load_band(name="water mask")
                dtm.apply_mask(water)
            except:
                print("WARNING: Unable to find any DTM product or unspecified in XML, skipping option withDTM...")
                with_dtm = False

    var_a.apply_mask(edge_a)
    var_b.apply_mask(edge_b)
    diff = np.abs(var_a.band - var_b.band)

    plt.clf()

    if with_dtm:
        plt.subplots(figsize=(12, 6))
        plt.subplot(121)

    plt.title("Abs(%s %s - %s %s)" % (a.context, a.type, b.context, b.type))
    plt.imshow(diff, interpolation='none')
    plt.colorbar(orientation="horizontal")

    if with_dtm:
        plt.subplot(122)
        plt.imshow(dtm.band, interpolation='none', cmap='terrain')
        plt.colorbar(orientation="horizontal")
        plt.title("MNT")

    plt.savefig("Diffmap_%s_%s-%s_%s.png" % (a.context, a.type, b.context, b.type))


def single_scatterplot(a, b, x_context="A", y_context="B", mode='aot', show=False):
    """
    :param a: sample a, all land and water pixels, numpy array (1D or 2D)
    :param b: sample b, all land and water pixels, numpy array (1D or 2D)
    :param title: string of title
    :param xt: label of x axis
    :param yt: label of y axis
    :param f_savefig: filename to save the figure to
    :param mode: defines axis range, defaults to 'aot'
    :param show: showing plot, defaults to False
    :return: png file
    """

    rmse = np.std(a.band - b.band)

    fig, (ax1) = plt.subplots(1, figsize=(6, 6))
    ax1.set_title("%i cloud-free pixels (rmse = %8.4f)" % (len(a.band), rmse))

    if mode == 'sre':
        ax1.plot([0, 1.0], [0, 1.0], 'k--')
    elif mode == 'aot':
        ax1.plot([0, 0.6], [0, 0.6], 'k--')

    ax1.set_aspect('equal', 'box')
    ax1.set_xlabel(x_context + " " + a.band_name)
    ax1.set_ylabel(y_context + " " + b.band_name)
    ax1.plot(a.band, b.band, 'bo', markersize=2)

    f_savefig = x_context + "_" + a.band_name.replace(" ", "-") + "_vs_" + y_context + "_" + b.band_name.replace(" ",
                                                                                                                 "-") + ".png"

    plt.savefig(f_savefig, format='png')
    if show == True:
        plt.show()


def single_quicklook(r, g, b):
    pass


def scatterplot(a, b, c, d, \
                title="demo", xt="x", yt="y", \
                f_savefig="demo.png", mode='aot', show=False):
    """
    :param a: sample a, all land and water pixels, numpy array (1D or 2D)
    :param b: sample b, all land and water pixels, numpy array (1D or 2D)
    :param c: sample c, only land pixels, numpy array (1D or 2D)
    :param d: sample d, only land pixels, numpy array (1D or 2D)
    :param title: string of title
    :param xt: label of x axis
    :param yt: label of y axis
    :param f_savefig: filename to save the figure to
    :param mode: defines axis range, defaults to 'aot'
    :param show: showing plot, defaults to False
    :return:
    """
    # slope_all, intercept_all, r_value_all, p_value_all, std_err_all = stats.linregress(a, b)
    # slope_ground, intercept_ground, r_value_ground, p_value_ground, std_err_ground = stats.linregress(c, d)

    diff_all = a - b
    std_all = np.sqrt(np.mean(abs(diff_all - diff_all.mean()) ** 2))

    diff_ground = c - d
    std_ground = np.sqrt(np.mean(abs(diff_ground - diff_ground.mean()) ** 2))
    std_ground = np.std(diff_all)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    ax1.set_title("Cloud-free pixels (std dev=%8.4f)" % (std_all))

    if mode == 'sre':
        ax1.plot([0, 1.0], [0, 1.0], 'k--')
    elif mode == 'aot':
        ax1.plot([0, 0.6], [0, 0.6], 'k--')

    ax1.set_aspect('equal', 'box')
    ax1.set_xlabel(xt)
    ax1.set_ylabel(yt)
    ax1.plot(a, b, 'bo', markersize=2)

    ax2.set_title("Cloud-free pixels, land only (std dev=%8.4f)" % (std_ground))
    if mode == 'sre':
        ax2.plot([0, 1.0], [0, 1.0], 'k--')
    elif mode == 'aot':
        ax2.plot([0, 0.6], [0, 0.6], 'k--')

    ax2.set_aspect('equal', 'box')
    ax2.set_xlabel(xt)
    ax2.set_ylabel(yt)
    ax2.plot(c, d, 'go', markersize=2)

    plt.savefig(f_savefig, format='png')
    if show == True:
        plt.show()


def atmplot(toa, rse, aot, \
            title="demo", xt="x", yt="y", \
            f_savefig="demo.png", mode='cloud_free', show=False):
    """
    :param a: sample a, all land and water pixels, numpy array (1D or 2D)
    :param b: sample b, all land and water pixels, numpy array (1D or 2D)
    :param c: sample c, only land pixels, numpy array (1D or 2D)
    :param d: sample d, only land pixels, numpy array (1D or 2D)
    :param title: string of title
    :param xt: label of x axis
    :param yt: label of y axis
    :param f_savefig: filename to save the figure to
    :param mode: defines axis range, defaults to 'aot'
    :param show: showing plot, defaults to False
    :return:
    """
    # slope_all, intercept_all, r_value_all, p_value_all, std_err_all = stats.linregress(a, b)
    # slope_ground, intercept_ground, r_value_ground, p_value_ground, std_err_ground = stats.linregress(c, d)

    fig, ax1 = plt.subplots(1, 1, figsize=(6, 6))

    if mode == "cloud_free":
        ax1.set_title("SRE=f(TOA), cloud-free pixels")
    else:
        ax1.set_title("SRE=f(TOA)")

    max_range = max([aot.max(), rse.max()])
    ax1.plot([0.0, max_range], [0.0, max_range], 'k--')

    ax1.set_aspect('equal', 'box')
    ax1.set_xlabel(xt)
    ax1.set_ylabel(yt)

    print("INFO: min AOT=%6.4f, max AOT=%6.4f" % (aot.min(), aot.max()))

    aot_step = 0.2
    for aot_r in np.arange(0, aot.max(), aot_step):
        aot_s = np.where(np.logical_and(aot >= aot_r, aot < aot_r + aot_step))
        lbl = ("AOT[%4.2f-%4.2f]" % (aot_r, aot_r + aot_step))
        ax1.plot(toa[aot_s], rse[aot_s], '.', markersize=4, label=lbl)
        ax1.legend(loc='upper right')

    plt.savefig(f_savefig, format='png')
    if show == True:
        plt.show()
