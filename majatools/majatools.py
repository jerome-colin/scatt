"""
Majatools PACKAGE

Recommended conda env:
python=3.6
conda install -c conda-forge gdal matplotlib scipy

Note: don't mix gdal packages from base and from forge
"""
__author__ = "Jerome Colin"
__license__ = "CC BY"
__version__ = "0.3.4"

import datetime
import glob
import sys
import numpy as np
from scipy.stats import describe
import matplotlib.pyplot as plt
import os
import re
from PIL import Image as pillow

from xml.dom import minidom
from osgeo import gdal


class Aoi:
    """
    Area of interest
    """

    def __init__(self, x_center, y_center, px_range):
        """Create an Aoi instance from pixel coordinates

        :param x_center: central pixel coordinate along the horizontal
        :param y_center: central pixel coordinate along the vertical
        :param px_range: half-height / half-width range of the square Aoi
        """

        self.x_center = x_center
        self.y_center = y_center
        self.px_range = px_range


class Context:
    """Description of a Maja simulation context used as an alternative to xml context file for Run

    """

    def __init__(self, path, type, context, scale_f_aot, scale_f_sr, nodata_aot, verbosity=False):
        self.verbosity = verbosity
        self.path = path
        self.type = type
        self.context = context
        self.scale_f_aot = scale_f_aot
        self.scale_f_sr = scale_f_sr
        self.nodata_aot = nodata_aot


class Image:
    """
    Extends a numpy array with attributes from Run, adds specifics methods
    """

    def __init__(self, band, band_name, scale_f=1, nodata=-10000, verbosity=False):
        """Initialization of Image

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
        :return: a vector of finite values of Image.band
        """
        if mask is not None:
            search = np.where(mask > 0)
            self.band[search] = np.nan

        return Image(np.extract(np.isfinite(self.band), self.band), self.band_name, verbosity=self.verbosity)

    def resample(self, n=6):
        """Downscaling method

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
        """Extract square image subset by pixel coordinates

        :param aoi: aoi instance
        :return: an Image instance
        """
        return Image(self.__crop_band_by_aoi(aoi), self.band_name, scale_f=1, nodata=self.nodata)

    def __check_is_inside(self, aoi):
        """Returns True if an AOI is within the extent of an image

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
        """Crops an Image band according to AOI extent

        :param aoi: an Aoi instance
        :return: an array
        """
        if self.__check_is_inside(aoi):
            return self.band[aoi.y_center - aoi.px_range:aoi.y_center + aoi.px_range, \
                   aoi.x_center - aoi.px_range:aoi.x_center + aoi.px_range]

    def __get_band_size(self):
        """Private, returns band dimensions

        :return: a tuple of int
        """
        try:
            return len(self.band[0, :]), len(self.band[:, 0])

        except IndexError:
            return len(self.band), 0

    def apply_mask(self, mask, reverse=False):
        """Apply a given mask to self, replacing non-zeros by NaN

        :param mask:
        """
        if reverse:
            search = np.where(mask.band == 0)
            self.band[search] = np.NaN
        else:
            search = np.where(mask.band != 0)
            self.band[search] = np.NaN

    def get_band_uint8(self):
        """Convert a Image.band array to unint8
        :return: an unint8 numpy array
        """
        b_max = np.max(self.band)
        if b_max > 0:
            img = self.band / np.max(self.band) * 256
        else:
            img = self.band * 0

        return img.astype(np.uint8)

    def get_pixel_count(self):
        x, y = self.__get_band_size()
        if y != 0:
            return x * y
        else:
            return x


class Run:
    """
    Run object described by an XML context file, storing run attributes and providing a load_band method

    TODO: replace the XML context by a dict, and implement a Yaml parser outside in a Context object
    """

    def __init__(self, config, verbosity=False):
        """Initialization of Run

        :param config: Context object or XML context file
        :param verbosity: increase verbosity to INFO is True
        """
        if isinstance(config, Context):
            if verbosity:
                print("INFO: Run got a Context object")

            self.verbosity = verbosity
            self.path = config.path
            # Fix missing trailing / in path
            if self.path[-1] != '/':
                self.path += '/'

            self.type = config.type
            self.context = config.context
            self.scale_f_aot = config.scale_f_aot
            self.scale_f_sr = config.scale_f_sr
            self.nodata_aot = config.nodata_aot


        else:
            try:
                self.xml_config = minidom.parse(config)
                self.verbosity = verbosity
                self.path = self.xml_config.getElementsByTagName("path")[0].firstChild.nodeValue
                self.type = self.xml_config.getElementsByTagName("type")[0].firstChild.nodeValue
                self.context = self.xml_config.getElementsByTagName("context")[0].firstChild.nodeValue
                self.scale_f_aot = float(self.xml_config.getElementsByTagName("scale_f_aot")[0].firstChild.nodeValue)
                self.scale_f_sr = float(self.xml_config.getElementsByTagName("scale_f_sr")[0].firstChild.nodeValue)
                self.nodata_aot = float(self.xml_config.getElementsByTagName("nodata_aot")[0].firstChild.nodeValue)

                # Fix missing trailing / in path
                if self.path[-1] != '/':
                    self.path += '/'

            except IndexError as e:
                print("ERROR: Mandatory parameter missing in XML file %s" % config)
                print(e)
                sys.exit(1)

            # Optional parameters
            try:
                self.dtm_path = self.xml_config.getElementsByTagName("dtm_path")[0].firstChild.nodeValue
                self.water_path = self.xml_config.getElementsByTagName("water_path")[0].firstChild.nodeValue

            except:
                pass

    def load_band(self, name="aot", subset=False, ulx=0, uly=0, lrx=0, lry=0):
        """Get a given image output of Run, optionally a subset by coordinates

        :param name: any string between "aot", "cloud_mask", "edge_mask"
        :param subset: if True use gdal_translate to subset an AOI
        :param uly: upper left y
        :param ulx: upper left x
        :param lrx: lower right x
        :param lry: lower right y
        :return: an Image object of the related band name
        """

        # Identify product file in product collection
        f_img = self.__get_f_img(name=name)

        # Extract a Gdal dataset from product file (optionally subset)
        product = get_geodata(f_img, subset=subset, ulx=ulx, uly=uly, lrx=lrx, lry=lry)

        if self.verbosity:
            print("INFO: Found file %s" % f_img.split("/")[-1])

        # Construct Image instance
        if name == "aot" and self.type == "maja":
            return Image(product.GetRasterBand(2).ReadAsArray(), self.type + " " + name, \
                         scale_f=self.scale_f_aot, verbosity=self.verbosity)
        elif name == "aot" and self.type == "maqt":
            return Image(product.GetRasterBand(1).ReadAsArray(), self.type + " " + name, \
                         scale_f=self.scale_f_aot, verbosity=self.verbosity)
        if name == "vap" and self.type == "maja":
            return Image(product.GetRasterBand(1).ReadAsArray(), self.type + " " + name, \
                         scale_f=20, verbosity=self.verbosity)
        elif name == "vap" and self.type == "maqt":
            return Image(product.GetRasterBand(1).ReadAsArray(), self.type + " " + name, \
                         scale_f=self.scale_f_aot, verbosity=self.verbosity)
        elif name[:3] == "sre":
            return Image(product.GetRasterBand(1).ReadAsArray(), self.type + " " + name, \
                         scale_f=self.scale_f_sr, verbosity=self.verbosity)
        elif name[:3] == "toa":
            return Image(product.GetRasterBand(1).ReadAsArray(), self.type + " " + name, \
                         scale_f=self.scale_f_sr, verbosity=self.verbosity)
        else:
            return Image(product.GetRasterBand(1).ReadAsArray(), self.type + " " + name, verbosity=self.verbosity)

    def get_rgb(self, subset=False, ulx=0, uly=0, lrx=0, lry=0):
        """Return B4, B3 and B2 Sentinel 2 SRE bands as a triplet of numpy array for quicklook

        :param subset:
        :param ulx:
        :param uly:
        :param lrx:
        :param lry:
        :return: red, green, blue as numpy array
        """
        red = self.load_band(name="sreB4", subset=subset, ulx=ulx, lry=lry, lrx=lrx, uly=uly)
        green = self.load_band(name="sreB3", subset=subset, ulx=ulx, lry=lry, lrx=lrx, uly=uly)
        blue = self.load_band(name="sreB2", subset=subset, ulx=ulx, lry=lry, lrx=lrx, uly=uly)
        return red, green, blue

    def get_timestamp(self, short=False):
        """Get the timestamp of a given Run object as a string of format %Y%m%d-%H%M%S"

        :return: a string
        """
        if self.type == "maja":
            try:

                if short:
                    pattern = "[0-9]{8}"
                    s = re.findall(pattern, self.path)
                    return s[0]
                else:
                    pattern = "[0-9]{8}\-[0-9]{6}"
                    s = re.findall(pattern, self.path)
                    return datetime.datetime.strptime(s[0], "%Y%m%d-%H%M%S")

            except IndexError:
                print("WARNING: couldn't find any timestamp for %s" % self.context)
                return None
        else:
            pass

    def get_type(self):
        """Return Run type

        :return: string of type of the run, typically "maja" or "maqt"
        """
        return self.type

    def __get_f_img(self, name):
        """Finds the actual file that relates to a variable name for a given product collection

        :param name: variable name
        :return: a file name
        """

        try:
            # Kept for backward compatibility
            if name == "aot":
                if self.type == "maja":
                    f_img = glob.glob(self.path + "*ATB_R1*")[0]
                elif self.type == "maqt":
                    f_img = glob.glob(self.path + "ORTHO_SURF_AOT/*10m.tau2")[0]
                else:
                    print("ERROR: Unable to find aot product...")
                    sys.exit(1)

            if name == "vap":
                if self.type == "maja":
                    f_img = glob.glob(self.path + "*ATB_R1*")[0]
                elif self.type == "maqt":
                    f_img = glob.glob(self.path + "ORTHO_SURF_VAP/*10m.tau2")[0]  # NOT TESTED YET!!
                else:
                    print("ERROR: Unable to find aot product...")
                    sys.exit(1)

            if name == "aot_R1":
                if self.type == "maja":
                    f_img = glob.glob(self.path + "*ATB_R1*")[0]
                elif self.type == "maqt":
                    f_img = glob.glob(self.path + "ORTHO_SURF_AOT/*10m.tau2")[0]
                else:
                    print("ERROR: Unable to find aot product...")
                    sys.exit(1)

            if name == "aot_R2":
                if self.type == "maja":
                    f_img = glob.glob(self.path + "*ATB_R2*")[0]
                elif self.type == "maqt":
                    f_img = glob.glob(self.path + "ORTHO_SURF_AOT/*20m.tau2")[0]
                else:
                    print("ERROR: Unable to find aot product...")
                    sys.exit(1)

            if name[:3] == "sre":
                if self.type == "maja":
                    try:
                        f_img = glob.glob(self.path + "*SRE_" + name[-2:] + "*")[0]
                    except IndexError:
                        if name[-3:] == "B8A":
                            f_img = glob.glob(self.path + "*SRE_" + name[-3:] + "*")[0]
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

            # Kept for backward compatibility
            if name == "cloud_mask":
                if self.type == "maja":
                    f_img = glob.glob(self.path + "MASKS/*CLM_R1.tif")[0]
                elif self.type == "maqt":
                    f_img = glob.glob(self.path + "MASQUES/NUAGES/*_10m.nua")[0]
                else:
                    print("ERROR: Unable to find cloud_mask product...")
                    sys.exit(1)

            if name == "cloud_mask_R1":
                if self.type == "maja":
                    f_img = glob.glob(self.path + "MASKS/*CLM_R1.tif")[0]
                elif self.type == "maqt":
                    f_img = glob.glob(self.path + "MASQUES/NUAGES/*_10m.nua")[0]
                else:
                    print("ERROR: Unable to find cloud_mask product...")
                    sys.exit(1)

            if name == "cloud_mask_R2":
                if self.type == "maja":
                    f_img = glob.glob(self.path + "MASKS/*CLM_R2.tif")[0]
                elif self.type == "maqt":
                    f_img = glob.glob(self.path + "MASQUES/NUAGES/*_20m.nua")[0]
                else:
                    print("ERROR: Unable to find cloud_mask product...")
                    sys.exit(1)

            # Kept for backward compatibility
            if name == "edge_mask":
                if self.type == "maja":
                    f_img = glob.glob(self.path + "MASKS/*_EDG_R1.tif")[0]
                elif self.type == "maqt":
                    f_img = glob.glob(self.path + "MASQUES/BORD/*_10m.bord_final")[0]
                else:
                    print("ERROR: Unable to find edge_mask product...")
                    sys.exit(1)

            if name == "edge_mask_R1":
                if self.type == "maja":
                    f_img = glob.glob(self.path + "MASKS/*_EDG_R1.tif")[0]
                elif self.type == "maqt":
                    f_img = glob.glob(self.path + "MASQUES/BORD/*_10m.bord_final")[0]
                else:
                    print("ERROR: Unable to find edge_mask product...")
                    sys.exit(1)

            if name == "edge_mask_R2":
                if self.type == "maja":
                    f_img = glob.glob(self.path + "MASKS/*_EDG_R2.tif")[0]
                elif self.type == "maqt":
                    f_img = glob.glob(self.path + "MASQUES/BORD/*_20m.bord_final")[0]
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

            return f_img

        except IndexError as e:
            print("ERROR: couldn't load %s from product %s" % (name, self.path))
            print(e)
            sys.exit(1)

    def __get_gdal_dataset(self, f_img, subset=False, ulx=0, uly=0, lrx=0, lry=0):
        """[DEPRECATED] Extract a Gdal object from a product file, optionally a subset from coordinates

        :param f_img: product image file
        :param subset: if True use gdal_translate to subset an AOI
        :param ulx: upper left x
        :param uly: upper left y
        :param lrx: lower right x
        :param lry: lower right y
        :return: a gdal object
        """
        print('WARNING: Run.__get_gdal_dataset is deprecated"')

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
        return data


class Timeseries:
    """
    A collection of Run instances
    """

    def __init__(self, f_config, verbosity=False):
        """Initialization of Timeseries

        :param f_config: XML configuration file
        :param verbosity: increase verbosity to INFO is True
        """
        try:
            self.xml_config = minidom.parse(f_config)
            self.verbosity = verbosity
            self.root_path = self.xml_config.getElementsByTagName("root_path")[0].firstChild.nodeValue
            self.collection_1_path = self.xml_config.getElementsByTagName("collection_1_path")[0].firstChild.nodeValue
            self.collection_2_path = self.xml_config.getElementsByTagName("collection_2_path")[0].firstChild.nodeValue
            self.type = self.xml_config.getElementsByTagName("type")[0].firstChild.nodeValue
            self.context_1 = self.xml_config.getElementsByTagName("context_1")[0].firstChild.nodeValue
            self.context_2 = self.xml_config.getElementsByTagName("context_2")[0].firstChild.nodeValue
            self.scale_f_aot = float(self.xml_config.getElementsByTagName("scale_f_aot")[0].firstChild.nodeValue)
            self.scale_f_sr = float(self.xml_config.getElementsByTagName("scale_f_sr")[0].firstChild.nodeValue)
            self.nodata_aot = float(self.xml_config.getElementsByTagName("nodata_aot")[0].firstChild.nodeValue)
            self.subset_ulx = float(self.xml_config.getElementsByTagName("subset_ulx")[0].firstChild.nodeValue)
            self.subset_uly = float(self.xml_config.getElementsByTagName("subset_uly")[0].firstChild.nodeValue)
            self.subset_lrx = float(self.xml_config.getElementsByTagName("subset_lrx")[0].firstChild.nodeValue)
            self.subset_lry = float(self.xml_config.getElementsByTagName("subset_lry")[0].firstChild.nodeValue)
            self.report = self.__convert_to_bool(self.xml_config.getElementsByTagName("report")[0].firstChild.nodeValue)
            self.plot = self.__convert_to_bool(self.xml_config.getElementsByTagName("plot")[0].firstChild.nodeValue)
            self.quicklook = self.__convert_to_bool(self.xml_config.getElementsByTagName("quicklook")[0].firstChild.nodeValue)
            self.reference_collection_path = self.xml_config.getElementsByTagName("reference_collection_path")[
                0].firstChild.nodeValue
            self.reference_collection_name = self.xml_config.getElementsByTagName("reference_collection_name")[
                0].firstChild.nodeValue
            # Fix missing trailing / in path
            if self.reference_collection_path[-1] != '/':
                self.reference_collection_path += '/'

            if verbosity:
                print("INFO: report set to %r, plot set to %r, quicklook set to %r" % (self.report, self.plot, self.quicklook))

        except IndexError as e:
            print("ERROR: Mandatory parameter missing in XML file %s" % f_config)
            print(e)
            sys.exit(1)

        except ValueError as e:
            print("ERROR: Bad value type in XML file %s" % f_config)
            print(e)
            sys.exit(1)

    def generate(self):
        """Create a collection of XML contexts for a given timeseries collection

        :return: a collection of files
        """
        if self.verbosity:
            print("INFO: gathering products")

        self.__get_product_fullpath_list()

        if self.verbosity:
            print("INFO: reference collection is set to %s" % self.reference_collection_path)

        if self.reference_collection_path is not None:
            self.common_product_reference_list = []
            self.__get_reference_fullpath_list()
            self.__common_reference()

    def find_reference(self, product):
        if self.reference_collection_path is not None:
            try:
                pid = [row[0] for row in self.common_product_reference_list].index(product)
                if self.verbosity:
                    print("INFO: found reference match for %s" % product)

                return self.common_product_reference_list[pid][1]

            except ValueError:
                if self.verbosity:
                    print("INFO: no reference match for %s" % product)

                return None

        else:
            return None

    def __get_product_fullpath_list(self):
        """Generate a list of products available in root_path

        :return: a list of files
        """
        pattern = "SENTINEL2*"
        product_list_1 = [os.path.basename(x) for x in glob.glob(self.root_path + self.collection_1_path + pattern)]
        product_list_2 = [os.path.basename(x) for x in glob.glob(self.root_path + self.collection_2_path + pattern)]

        # Keep common products of both collections in a list
        self.common_product_list = self.__common_member(product_list_1, product_list_2)
        self.common_product_list_len = len(self.common_product_list)

    def __get_reference_fullpath_list(self):

        pattern = "*.hdf"

        self.reference_list = [os.path.basename(x) for x in glob.glob(self.reference_collection_path + pattern)]

        if self.verbosity:
            print("INFO: found %i reference elements in collection" % len(self.reference_list))

    def __convert_to_bool(self, val):
        if val == "True":
            return True
        if val == "False":
            return False
        else:
            return None

    def __common_member(self, a, b):
        """Return common elements of two lists a and b

        :param a: list
        :param b: list
        :return: list
        """
        a_set = set(a)
        b_set = set(b)
        if (a_set & b_set):
            return sorted(list(a_set & b_set))
        else:
            print("WARNING: No common elements in product list")

    def __common_reference(self):
        """Defines a list of common product:reference pairs

        :return: self.common_product_reference_list
        """

        for product in self.common_product_list:
            try:
                pattern = "[0-9]{8}"
                product_timestamp = re.findall(pattern, product)

            except IndexError:
                print("WARNING: couldn't find any timestamp for %s" % self.common_product_list)

            for ref in self.reference_list:
                try:
                    pattern = "[0-9]{8}"
                    reference_timestamp = re.findall(pattern, ref)

                    if product_timestamp == reference_timestamp:
                        if self.verbosity:
                            print("INFO: %s matches %s" % (ref, product))

                        self.common_product_reference_list.append([product, ref])

                except IndexError:
                    print("WARNING: couldn't find any timestamp for %s" % self.common_product_list)

    def __to_context(self):
        pass


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


def get_geodata(f_img, subset=False, ulx=0, uly=0, lrx=0, lry=0):
    """Extract a Gdal object from a product file, optionally a subset from coordinates

    :param f_img: product image file
    :param subset: if True use gdal_translate to subset an AOI
    :param ulx: upper left x
    :param uly: upper left y
    :param lrx: lower right x
    :param lry: lower right y
    :return: a gdal object
    """
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

    return data


def get_hdf_as_array(hdf_file, scale_f=10000, verbose=False):
    """Return and HDF MODIS-like file as array of dimension band:cols:rows

    :param hdf_file: HDF MODIS-like file
    :param scale_f: MODIS reflectance scale factor, defaults to 1/10000
    :param verbose: if True prints stats per band
    """

    # open the dataset
    hdf_ds = gdal.Open(hdf_file, gdal.GA_ReadOnly)
    hdf_subds = hdf_ds.GetSubDatasets()
    print(len(hdf_subds))

    hdf_arr = np.zeros((len(hdf_subds), 900, 900))

    for b in range(len(hdf_subds)):
        hdf_subds_name = hdf_subds[b][0]
        hdf_arr[b, :, :] = gdal.Open(hdf_subds_name, gdal.GA_ReadOnly).ReadAsArray().astype(np.int16) / scale_f

    if verbose:
        for b in range(len(hdf_subds)):
            print(hdf_subds[b][0])
            print(describe(hdf_arr[b].flatten()))

    return hdf_arr


def diffmap(a, b, mode, with_dtm=False):
    """Produce an absolute difference map as image

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


def single_scatterplot(a, b, mask, x_context="A", y_context="B", mode='aot', png=False):
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

    masked_a = a.get_finite(mask=mask)
    masked_b = b.get_finite(mask=mask)

    px_count = a.get_pixel_count()
    mask_count = masked_a.get_pixel_count()

    ratio = mask_count / px_count

    try:
        #rmse = np.std(masked_a.band - masked_b.band)
        rmse = np.sqrt(np.mean((masked_a.band - masked_b.band)**2))


        if png == True:
            fig, (ax1) = plt.subplots(1, figsize=(6, 6))
            ax1.set_title("%3.1f %% cloud-free pixels (rmse = %5.4f)" % (ratio * 100, rmse))

            if mode == 'sre':
                ax1.plot([0, 1.0], [0, 1.0], 'k--')
                ax1.set_xlim(0, 1)
                ax1.set_ylim(0, 1)
            elif mode == 'aot':
                ax1.plot([0, 0.6], [0, 0.6], 'k--')

            ax1.set_aspect('equal', 'box')
            ax1.set_xlabel(x_context + " " + a.band_name)
            ax1.set_ylabel(y_context + " " + b.band_name)
            ax1.plot(masked_a.band, masked_b.band, 'bo', markersize=2)

            f_savefig = x_context + "_" + masked_a.band_name.replace(" ",
                                                                     "-") + "_vs_" + y_context + "_" + masked_b.band_name.replace(
                " ", "-") + ".png"
            plt.savefig(f_savefig, format='png')
            plt.close(fig)

        return ratio, rmse

    except ValueError:
        print("WARNING: got a bad shape in single_scatterplot, had to skip this one : %s" % (a.band_name))
        print("-------: image %s has shape %s" % (a.band_name, str(a.band.shape)))
        print("-------: image %s has shape %s" % (b.band_name, str(b.band.shape)))
        return 0, 0


def single_quicklook_rgb(r, g, b, outfile="toto.png"):
    """Generate an RGB quicklook

    :param r: red Image
    :param g: green Image
    :param b: blue Image
    :param outfile: image file, usually the Run.context without extension
    :return: a PNG file
    """
    array_stack = np.dstack((r.get_band_uint8(), g.get_band_uint8(), b.get_band_uint8()))
    img = pillow.fromarray(array_stack)
    if outfile[-4:] != ".png":
        outfile += ".png"

    img.save(outfile)


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
