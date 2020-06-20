# Defining global variables for model function

import calendar
from datetime import datetime, timedelta, date
import numpy as np
import rasterio
from s3fs.core import S3FileSystem
import os
import boto3
import sys
import yaml
import fiona
import rasterio.mask
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio import shutil as rio_shutil
from rasterio.vrt import WarpedVRT

#TODO - get large extent tif files for testing the RasterManager warping functions


class VegConfig:
    """"""

    # Attributes Here.

    # === Minimalist Beta Version Params ====
    start_year = None
    end_year = None
    start_day = None
    end_day = None
    geoproperties_file = None
    shapefile = None

    # ==== Version 2.0 params ====
    #input/output
    in_root = None
    out_root = None

    precip_settings = {'crs': None, 'cols': None, 'rows': None, 'xres': None, 'yres': None, 'left': None, 'top': None}
    ndvi_settings = {'crs': None, 'cols': None, 'rows': None, 'xres': None, 'yres': None, 'left': None, 'top': None}
    pet_settings = {'crs': None, 'cols': None, 'rows': None, 'xres': None, 'yres': None, 'left': None, 'top': None}
    tmin_settings = {'crs': None, 'cols': None, 'rows': None, 'xres': None, 'yres': None, 'left': None, 'top': None}
    tavg_settings = {'crs': None, 'cols': None, 'rows': None, 'xres': None, 'yres': None, 'left': None, 'top': None}
    tmax_settings = {'crs': None, 'cols': None, 'rows': None, 'xres': None, 'yres': None, 'left': None, 'top': None}

    non_std_inputs = False


    def __init__(self, config_dictionary):
        for key, val in config_dictionary.items():
            setattr(self, key, val)

    def update_config(self, k, v, cfg_path):
        """"""
        setattr(self, k, v)
        with open(cfg_path, 'r') as cfg:
            file_dict = yaml.safe_load(cfg)
        file_dict['{}'.format(k)] = v
        with open(cfg_path, 'w') as w_file:
            w_file.write(yaml.dump(file_dict))
            # # old way
            # append_file.write('{}: {}\n'.format(k, v))
    def update_feature(selfs, k, v, cfg_path):
        """Used to update the config_dict file when a configuration needs to be ADDED or APPENDED to something, like a
        dictionary that gets new entries
        This function does not update the attribute like update_config does. This function is used when an attribute
         has ALREADY been appended or updated in the code but we still need to update the configuration FILE"""
        with open(cfg_path, 'r') as cfg:
            file_dict = yaml.safe_load(cfg)
        # overprint the entries with the new config_dict
        file_dict['{}'.format(k)] = v
        with open(cfg_path, 'w') as w_file:
            w_file.write(yaml.dump(file_dict))



class PathManager:
    """
    This class creates the input paths for the dynamic and static data needed in the model.
    """


    # Todo - configuration non_std_inputs = True, accept parameters such as ndvi_fmt = userndvi_{}_vers1.tif
    #  where {} is the date and date_fmt YYYYmmdd or YYYYdd etc for a filepath. A dictionary can relate user-input
    #  name formatting to model std formats. Right now, non_std_inputs must be false and only 1 type of naming
    #  convention per paramter file is allowed.

    ndvif = None
    pptf = None
    petf = None
    tavgf = None
    tminf = None
    tmaxf = None

    def __init__(self, config):
            self.config = config


    def get_dynamic_data(self, today, settings): # DOY=None, year_doy=None
        """
        This gets dynamic data, such as NDVI, precipitation, temperature, etc..
        :param today: datetime object to retrieve year, month, day, etc. from today
        :param settings: data set characteristics from the configurations file
        :return:
        """

        name_key = 'name_fmt'
        loc_key = 'dir_loc'
        dt_key = 'dt_fmt'
        clim_key = 'climatology'
        doy = today.timetuple().tm_yday

        print('settings', settings)

        if settings[clim_key]:
            # for climatology then we expect a DOY format
            if settings[dt_key] == 'doy':
                dynamic_key = '{:03d}'.format(doy)
            else:
                print('{} is set to climatology but date format from config_dict is {}'.format(settings[name_key],
                                                                                          settings[dt_key]))
                sys.exit(0)
        elif settings[dt_key] == 'YYYYdoy':
            dynamic_key = '{}{:03d}'.format(today.year, doy)
        else:
            print('Hey user, the format of the dt_fmt configuration you gave: {} is not supported at '
                  'this time'.format(settings[dt_key]))
            sys.exit(0)

        fpath = os.path.join(settings[loc_key], settings[name_key].format(dynamic_key))
        return fpath


    def get_static_data(self, settings):
        """
        This gets static data, such as soil data sets.
        :param settings: data set characteristics from the configurations file
        :return:
        """
        fpath = settings['file_loc']
        return fpath


    def s3_delete_local(self, outpath, from_file, bucket, prefix_no_slash):
        """
        This function will move the model outputs from a local folder to a cloud bucket.
        :param from_file: path the the local folder
        :param bucket: name of the cloud bucket = 'dev-et-data'
        :param prefix_no_slash: "folder" in cloud bucket  = 'v1DRB_outputs'
        :return:
        """

        objecta='{}/{}'.format(prefix_no_slash,outpath)
        s3 = boto3.client('s3')
        with open(from_file, "rb") as f:
            s3.upload_fileobj(f, bucket, objecta)
        os.remove(from_file)


class RasterManager:
    """
    This class addresses all the data wrangling needs to get the data sets needed in the model
    to the extent nad resolution chosen by the user.
    The user needs to input a sample shapefile (polygon) and geotiff file that is used to
    get the attributes for processing extent and resolution, etc.
    """

    # TODO - get large extent tif files for testing the RasterManager warping functions
    # TODO - make overlapping of 2 dataset work
    # TODO - rasterio either clips based on a sample raster and shapefile, or based on user defined geo info
    # #TODO - https://rasterio.readthedocs.io/en/latest/topics/masking-by-shapefile.html
    # Todo - test raster manager as a standalone module.

    out_root = None
    # --- geographic info for destination files ---
    crs = None
    cols = None
    rows = None
    xres = None
    yres = None
    # left geo coord 'e.g. Westernmost Longitude"
    left = None
    # top geo coord e.g highest Latitude extent
    top = None
    transform = [xres, 0.0, left, 0.0, yres, top]

    geoproperties_file = None
    # can contain multiple features of interest
    shapefile = None

    temp_folder = None


    def __init__(self, config):
        self.config = config

        self.geoproperties_file = config.geoproperties_file
        self.shapefile = config.shapefile
        self.temp_folder = os.path.join(config.out_root, config.temp_folder)
        if not os.path.exists(self.temp_folder):
            os.makedirs(self.temp_folder)

        if self.geoproperties_file == None or self.shapefile==None:
            print('Assuming the user entered values in the config_dict for boundaries of the AOI not implemented at thsi time')
            sys.exit(0)

    # ----------- create output rasters -----------------
    def output_rasters(self, arr, outdir, outname):
        """
        This function creates geotiff files from the model output arrays.
        """

        outpath = os.path.join(outdir, outname)
        print('the outpath for file {} is {}'.format(outname, outpath))

        # get the geoinfo from sample tiff to output intermediate files
        ds = rasterio.open(self.geoproperties_file)
        band1 = arr
        with rasterio.open(outpath, 'w', driver='GTiff', height=self.rows, width=self.cols,
                           count=1, dtype='float64', crs=self.crs, transform=self.transform) as wrast:
            wrast.write(band1, indexes=1)

        # TODO - Set an AWS Cloud flag in the config_dict file to activate this function or not...
        # delete files created locally and put in bucket
        # PathManager.s3_delete_local(from_file, bucket, prefix_no_slash)


    def set_model_std_grid(self, feat=0):
        """Clips and crops a tiff to the extent of a feature in a shapefile
        :param feat: feat is  the feature id of the shapefile from like a GeoJSON)
        # https://rasterio.readthedocs.io/en/latest/topics/virtual-warping.html
        """
        print(self.shapefile)
        with fiona.open(self.shapefile, 'r') as shapefile:
            # todo - set up an error if user has shapefile with more than one feature.
            # shape = shapefile[0]['geometry']
            shapes = [feature["geometry"] for feature in shapefile]

            for feature in shapefile:
                # matching the FID of the given shapefile from a typical geoJSON (Not Ordered Dict nonsense)
                if feat == feature['id']:
                    shapes = [feature['geometry']]
        # print('This is the shape var:', shapes)

        with rasterio.open(self.geoproperties_file, 'r') as src:
            out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
            out_meta = src.meta
        # once the image is cropped, the image metadata dictionary is updated with the cropped transform and bounds.
        out_meta.update({"driver": "GTiff",
                         "height": out_image.shape[1],
                         "width": out_image.shape[2],
                         "transform": out_transform})

        self.crs = out_meta['crs']
        self.transform = out_meta['transform']
        self.left = self.transform[2]
        self.top = self .transform[5]
        self.cols = out_meta['width']
        self.rows = out_meta['height']
        self.xres = self.transform[0]
        self.yres = self.transform[4]
        # return out_meta

    def normalize_to_std_grid(self, inputs, resamplemethod = 'nearest'):
        """
        Uses rasterio virtual raster to standardize grids of different crs, resolution, boundaries based on  a shapefile geometry feature
        :param inputs: a list of (daily) raster input files for the water balance.
        :param outloc: output locations 'temp' for the virtual files
        :return: list of numpy arrays
        """
        outputs = []
        npy_outputs = []
        if resamplemethod == 'nearest':
            rs = Resampling.nearest
        else:
            print('only nearest neighbor resampling is supported at this time')
            sys.exit(0)

        for i, warpfile in enumerate(inputs):
            print('warpfile', warpfile)
            with rasterio.open(warpfile) as src:
                # create the virtual raster based on the standard rasterio attributes from the sample tiff and shapefile feature.
                with WarpedVRT(src, resampling=rs,
                               crs=self.crs,
                               transform=self.transform,
                               height=self.rows,
                               width=self.cols) as vrt:
                    data = vrt.read()
                    print(type(vrt))
                    # save the file as an enumerated tiff. reopen outside this loop with the outputs list
                    outwarp = os.path.join(self.temp_folder, 'temp_{}.tif'.format(i))
                    rio_shutil.copy(vrt, outwarp, driver='GTiff')
                    outputs.append(outwarp)

        # output each virtual file as a temporary .tif file in a temp folder somewhere in the outputs directory.
        # for each file in the temp directory read in the raster as a numpy array and return the list of numpy arrays
        # from this method for us in the rest of the code.
        for ow in outputs:
            with rasterio.open(ow, 'r') as src:
                arr = src.read(1)
                npy_outputs.append(arr)

        return npy_outputs


class VegET:
    """
    This is the heart of the Veg ET model, a soil-water balance model.
    There are 3 main functions: soil_water, surface_runoff, veg_et which include the model calculations.
    """
    ### Atttributes = None, True False to start out.

    start_year = None
    end_year = None
    start_day = None
    end_day = None

    rf_low_thresh_temp = None
    rf_high_thresh_temp = None
    rf_value = None
    melt_factor = None
    dc_coeff = None
    rf_coeff = None
    k_factor = None
    ndvi_factor = None
    water_factor = None
    bias_corr = None
    alfa_factor = None

    geoproperties_file = None
    outdir = None
    accumulate_mode = None
    # ----- static soil rasters for model run------
    interception_settings = None
    whc_settings = None
    field_capacity_settings = None
    saturation_settings = None
    watermask_settings = None

    # dataset_configurations
    precip_settings = None
    ndvi_settings = None
    pet_settings = None
    tmin_settings = None
    tavg_settings = None
    tmax_settings = None

    def __init__(self, veget_config_path=None):
        if veget_config_path is None:
            print('You need to point to the configure file to get the path')
            sys.exit(1)
        else:
            # this allows for the config_dict to be created from a preexisting file
            self.config_path = veget_config_path
            if os.path.exists(self.config_path):
                with open(self.config_path, 'r') as cfgpath:
                    self.config_dict = yaml.safe_load(cfgpath)
            else:
                print('the path does not exist, check the path you gave VegET()')
                sys.exit(0)
            # create an instance of the VegET model using the configurations from the file.
            self.config = VegConfig(self.config_dict)
            # initialize the classes that manage Raster data and input/output paths to the data
            self.rmanager = RasterManager(config=self.config)
            sgrid = self.rmanager.set_model_std_grid(self)
            self.pmanager = PathManager(config=self.config)

            # set startday and endday
            self.start_day = self.config.start_day
            self.end_day = self.config.end_day
            self.start_year = self.config.start_year
            self.end_year = self.config.end_year
            # set interception settings
            self.interception_settings = self.config.interception_settings
            self.whc_settings = self.config.whc_settings
            self.saturation_settings = self.config.saturation_settings
            self.watermask_settings = self.config.watermask_settings
            self.field_capacity_settings = self.config.field_capacity_settings
            self.ndvi_settings = self.config.ndvi_settings
            self.precip_settings = self.config.precip_settings
            self.pet_settings = self.config.pet_settings
            self.tavg_settings = self.config.tavg_settings
            self.tmin_settings = self.config.tmin_settings
            self.tmax_settings = self.config.tmax_settings

            # set here in init so they can be float
            self.rf_low_thresh_temp = float(self.config.rf_low_thresh_temp)
            self.rf_high_thresh_temp = float(self.config.rf_high_thresh_temp)
            self.rf_value = float(self.config.rf_value)
            self.melt_factor = float(self.config.melt_factor)
            self.dc_coeff = float(self.config.dc_coeff)
            self.rf_coeff = float(self.config.rf_coeff)
            self.k_factor = float(self.config.k_factor)
            self.ndvi_factor = float(self.config.ndvi_factor)
            self.water_factor = float(self.config.water_factor)
            self.bias_corr = float(self.config.bias_corr)
            self.alfa_factor = float(self.config.alfa_factor)


            # accumulation modes
            self.accumulate_mode = self.config.accumulate_mode

            # set the output dir and make it if it doens't exist
            self.outdir = self.config.out_root
            if not os.path.exists(self.outdir):
                os.makedirs(self.outdir)

    def _day_of_year(self, today):
        year = today.year
        print(year, today.month, today.day)
        DOY = '{:03d}'.format(today.timetuple().tm_yday)
        print(f'today is {DOY}')

        return DOY, year

    def _end_of_month(self, day, mon, year):
        #  calendar.monthrange return a tuple
        # (weekday of first day of the month, number of days in month)
        last_day_of_month = calendar.monthrange(year, mon)[1]
        # check if date is last day of month
        if date(year, mon, day) == date(year, mon, last_day_of_month):
            return True

        return False

    def _soil_water(self, i, ppt, interception, tmin, tmax, tavg, melt_factor, rf_high_thresh_temp, rf_low_thresh_temp,
                    yest_swf=None, yest_snwpck=None):

        """
        This function takes precipitation, interception, and temperature data to determine
        how much rain and snow was falling during that day's precipitation event, if there was
        snow, how much snow water equivalent was it, and what is the
        initial/starting soil water balance for the day after adding the water amount (rain/snow).
        :param SWi: initial/starting soil water balance for the day after adding the water amount (rain/snow)
        :param saturation: Soil saturation value (mm)
        :param field_capacity: field capacity value (mm) of soil
        :param whc: Water Holding Capacity (WHC) or available water content (AWC) (mm)
        :param rf_low_thresh_temp : 0 degree Celcius -> temperature below where precipitation is considered snow
        :param rf_high_thresh_temp : 6 degree Celcius -> temperature above where precipitation is considered rain
        :param rf_value : 0.167  slope value
        :param melt_factor : 0.06   from literature
        :return: SWi
        """

        # Check for no data value handling
        print('ppt min', np.min(ppt))
        print('ppt max', np.min(ppt))
        ppt[ppt <= -1] = np.nan
        ppt[ppt == 32767] = np.nan
        print('ppt min', np.min(ppt))
        print('ppt max', np.max(ppt))

        print('tavg min', np.min(tavg))
        tavg[tavg <= -100] = np.nan
        print('tavg min', np.min(tavg))
        tmax[tmax <= -100] = np.nan
        tmin[tmin <= -100] = np.nan

        # Effective precipitation
        effppt = ppt * (1 - (interception / 100.0))
        # Intercepted precipitation
        interception = ppt * (interception / 100.0)
        print('interception min', np.min(interception))

        # Snow pack
        # Usage: Creates a melt rate value based on the relationship between
        # max and min air temperature to determine the snow melt and from there the snow pack extent
        melt_rate = melt_factor * ((tmax * tmax) - (tmax - tmin))
        # initialize the snow melt factor array
        snow_melt_fac = np.zeros(ppt.shape)
        # where avg temp <= high_threshold_temp set to 0, else it is equal to the melt factor rate
        # (tavg <= rf_high_tresh_temp, 0, melt_rate)
        snow_melt_fac[tavg <= rf_high_thresh_temp] = melt_rate[tavg <= rf_high_thresh_temp]
        snow_melt_fac[tavg > rf_high_thresh_temp] = 0

        if i == 0:  # first day of model run to initalize and establish the soil water balance
            # Usage: Creates a fraction value based on average temperature that determines
            # if the incoming precipitation is falling as rain, sleet, or snow.
            # if tavg <= 0, make it 0, else if tavg >= 6, make it 1, else (0.167*(tavg-6))
            rain_frac = np.zeros(ppt.shape)
            rain_frac[tavg <= rf_low_thresh_temp] = 0
            rain_frac[tavg >= rf_high_thresh_temp] = 1
            temp_diff_boolean = (tavg < rf_high_thresh_temp) | (tavg > rf_low_thresh_temp)
            rain_frac[temp_diff_boolean] = self.rf_value * (tavg[temp_diff_boolean] - rf_high_thresh_temp)

            RAIN = rain_frac * effppt
            SWE = np.zeros(ppt.shape)  # inital snowpack raster with only 0 values
            snow_melt = SWE
            SNWpk = np.zeros(ppt.shape)  # inital snowpack raster with only 0 values
            SWi = (self.whc * 0.5) + effppt + snow_melt
            SWi[SWi < 0] = np.nan

            return SWi, SNWpk, RAIN, SWE, snow_melt

        else:
            rain_frac = np.zeros(ppt.shape)  # initialize the rain fraction array
            # Creates a fraction value based on average temperature that determines
            # if the incoming precipitation is falling as rain, sleet, or snow.
            # if tavg <= 0, make it 0, else if tavg >= 6, make it 1, else (0.167*(tavg-6))
            rain_frac[tavg <= rf_low_thresh_temp] = 0
            rain_frac[tavg >= rf_high_thresh_temp] = 1
            temp_diff_boolean = (tavg < rf_high_thresh_temp) & (tavg > rf_low_thresh_temp)
            rain_frac[temp_diff_boolean] = self.rf_value * (tavg[temp_diff_boolean]
            )

            RAIN = rain_frac * effppt
            SWE = (1 - rain_frac) * effppt

            # snow melt
            snow_melt = np.zeros(ppt.shape)
            # snow_melt = if melt_rate <= (SWE + yesterday's snowpack), make it melt_rate, else (SWE + yesterday's snowpack)
            snow_melt_boolean = (melt_rate <= (SWE + yest_snwpck))
            snow_melt[snow_melt_boolean] = melt_rate[snow_melt_boolean]
            snow_melt[~snow_melt_boolean] = SWE[~snow_melt_boolean] + yest_snwpck[~snow_melt_boolean]

            SNWpk = np.zeros(ppt.shape)
            # today's snow pack = yesterdays's snow pack + snow water amount - snowmelt
            SNWpk1 = yest_snwpck + SWE - snow_melt
            # set SNWpk1 to 0 if its > 0
            SNWpk[SNWpk1 < 0] = 0
            SNWpk[SNWpk1 >= 0] = SNWpk1[SNWpk1 >= 0]

            # initial soil water balance = yesterday's soil water balance + Rain + snow melt
            SWi = yest_swf + RAIN + snow_melt

            return SWi, SNWpk, RAIN, SWE, snow_melt

    def _surface_runoff(self, SWi, saturation, field_capacity, whc, rf_coeff, geo_dict=None):
        """
        This function determines the runoff part of the model. Runoff is the total of deep drainage and surface runoff.
        :param SWi: initial/starting soil water balance for the day after adding the water amount (rain/snow)
        :param saturation: Soil saturation value (mm)
        :param field_capacity: field capacity value (mm) of soil
        :param whc: Water Holding Capacity (WHC) or available water content (AWC) (mm)
        :param rf_coeff: determines how much of the total runoff is attributed to the surface runoff only
        :param geo_dict: geo-spatial dictionary with attributes of resolution, extent, etc.
        :return: deep drainage and surface runoff
        """

        saturation[saturation < 0] = np.nan
        field_capacity[field_capacity < 0] = np.nan
        whc[whc < 0] = np.nan

        # total runoff based on water left in soil after SAT-FC
        sat_fc = saturation - field_capacity
        Rf1 = SWi - whc
        # if runoff is < 0, make it 0
        Rf = np.zeros(SWi.shape)
        rf_boolean = (Rf1 >= 0)
        Rf[rf_boolean] = Rf1[rf_boolean]

        # Surface runoff
        SRf = np.zeros(SWi.shape)
        # SRf = if rf <= sat_fc, make it (rf * rf_coeff)(35% of the runoff value), else (rf - sat_fc) + (rf_coeff * sat_fc)
        SRf_boolean = (Rf <= sat_fc)
        SRf[SRf_boolean] = Rf[SRf_boolean] * rf_coeff
        SRf[~SRf_boolean] = (Rf[~SRf_boolean] - sat_fc[~SRf_boolean]) + rf_coeff * sat_fc[~SRf_boolean]
        # Deep Drainage
        # DDrain occurs if SWi > WHC, amount of DDrain is SAT <> WHC with a maximum DDrain of SAT - WHC
        DDrain = Rf - SRf

        return DDrain, SRf

    def _veg_et(self, k_factor, ndvi_factor, water_factor, bias_corr, alfa_factor, watermask, pet, ndvi, SWi):
        """
        :param k_factor: 1.25  -> adjusting value in ETa calculation
        :param ndvi_factor: 0.2 -> adjusting value in ETa calculation
        :param water_factor: 0.7 -> VegET over water is determined as 70% of pet
        :param bias_corr: 0.85 -> bias correction for gridmet pet data
        :param alfa_factor: 1.25 -> adjusting grass pet to alfalfa pet
        :param watermask: raster dataset with all inland water bodies
        :param pet: reference ET raster dataset
        :param ndvi: normalized difference vegetation index (ndvi) raster dataset
        :param SWi: initial/starting soil water balance for the day after adding the water amount (rain/snow)
        :param SWf: final/ending soil water balance for the day after taking out ET
        :return:  eta, SWf
        """

        etasw1 = np.zeros(ndvi.shape)
        etasw3 = np.zeros(ndvi.shape)
        etasw4 = np.zeros(ndvi.shape)
        etasw5 = np.zeros(ndvi.shape)
        etasw = np.zeros(ndvi.shape)
        SWf = np.zeros(ndvi.shape)

        etasw1A = (k_factor * ndvi + ndvi_factor) * (pet * bias_corr)
        etasw1B = (k_factor * ndvi) * (pet * bias_corr)

        # etasw1 = if ndvi > 0.4, make it etasw1A, else etasw1B
        ndvi_boolean = (ndvi > 0.4)
        etasw1[ndvi_boolean] = etasw1A[ndvi_boolean]
        etasw1[~ndvi_boolean] = etasw1B[~ndvi_boolean]

        etasw2 = (SWi / (0.5 * self.whc)) * etasw1

        # etasw3 = if SWi > (0.5 * WHC), make it etasw1, else etasw2
        etasw3_boolean = (etasw3 > SWi)
        etasw3[etasw3_boolean] = etasw1[etasw3_boolean]
        etasw3[~etasw3_boolean] = etasw1[~etasw3_boolean]

        # etasw4 = if etasw3 > SWi, make it SWi, else etasw3
        etasw4_boolean = (etasw3 > SWi)
        etasw4[etasw4_boolean] = SWi[etasw4_boolean]
        etasw4[~etasw4_boolean] = etasw3[~etasw4_boolean]

        # etasw = if etasw4 > WHC, make it WHC, else etasw4
        etasw_boolean = (etasw4 > self.whc)
        etasw5[etasw_boolean] = self.whc[etasw_boolean]
        etasw5[~etasw_boolean] = etasw4[~etasw_boolean]

        # ETa of water bodies = 0.70 * (1.25*0.85*ETo)
        water_var = water_factor * bias_corr * alfa_factor
        print(watermask.shape)

        print(pet.shape)
        etawater_boolean = (watermask == 1)
        print(etawater_boolean.shape)
        # put the final etasw values for no-water regions in the final array
        etasw[~etawater_boolean] = etasw5[~etawater_boolean]
        # if it is a water-region, etasw = (calculated ET of water bodies)
        etasw[etawater_boolean] = pet[etawater_boolean] * water_var
        print(etasw.shape)

        SWf1 = SWi - etasw

        # SWf = if SWi > WHC, make it (WHC - etasw), else (if SWf1 < 0.0, make it 0.0, else SWf1)
        SWf_boolean = (SWi > self.whc)
        SWf_boolean2 = (SWf1 < 0.0)

        SWf[SWf_boolean] = self.whc[SWf_boolean] - etasw[SWf_boolean]
        SWf[SWf_boolean2] = 0
        SWf[~SWf_boolean2] = SWf1[~SWf_boolean2]

        return etasw, SWf, etasw5

    def _run_water_bal(self, i, today, interception, whc, field_capacity, saturation,
                       rf_coeff, k_factor, ndvi_factor, water_factor, bias_corr, alfa_factor, watermask, outdir,
                       yest_snwpck=None, yest_swf=None, geoproperties_file=None, daily_mode=True):
        """Here the water balance functions are combined into the water balance model.
        The needed input datasets are collected from buckets in the cloud, the needed functions executed
        and output datasets set up for daily, monthly, yearly rasters.
        """

        #dynamic inputs to the model
        self.ndvi = self.pmanager.get_dynamic_data(today, self.ndvi_settings)
        self.pet = self.pmanager.get_dynamic_data(today, self.pet_settings)
        self.ppt = self.pmanager.get_dynamic_data(today, self.precip_settings)
        self.tavg = self.pmanager.get_dynamic_data(today, self.tavg_settings)
        self.tmin = self.pmanager.get_dynamic_data(today, self.tmin_settings)
        self.tmax = self.pmanager.get_dynamic_data(today, self.tmax_settings)

        # Call Raster Manager function to standardize all the input dataset.
        dynamic_inpts = [self.ndvi, self.pet, self.ppt, self.tavg, self.tmin, self.tmax]

        # All the variables are now Numpy Arrays!
        self.ndvi, self.pet, self.ppt, self.tavg, self.tmin, self.tmax = \
            self.rmanager.normalize_to_std_grid(inputs=dynamic_inpts, resamplemethod='nearest')

        # ====== Call the functions ======
        # output SWi and SNWpk
        SWi, SNWpk, RAIN, SWE, snow_melt = self._soil_water(i, self.ppt, interception, self.tmin, self.tmax, self.tavg,
                                                            self.melt_factor, self.rf_high_thresh_temp, self.rf_low_thresh_temp,
                                                            yest_swf, yest_snwpck)
        DOY, year = self._day_of_year(today=today)

        SWiout =  f'swi_{year}{DOY}.tif'
        print('swout', SWiout)
        SNWpkout = f'snwpk_{year}{DOY}.tif'
        RAINout =  f'rain_{year}{DOY}.tif'
        SWEout = f'swe_{year}{DOY}.tif'
        snow_meltout =  f'snowmelt_{year}{DOY}.tif'

        if daily_mode:
            self.rmanager.output_rasters(SWi, self.outdir, outname=SWiout)
            self.rmanager.output_rasters(SNWpk, self.outdir, outname=SNWpkout)
            self.rmanager.output_rasters(RAIN, self.outdir, outname=RAINout)
            self.rmanager.output_rasters(SWE, self.outdir, outname=SWEout)
            self.rmanager.output_rasters(snow_melt, self.outdir, outname=snow_meltout)

        # output DDRAIN and SRf
        DDrain, SRf = self._surface_runoff(SWi, saturation=self.saturation, field_capacity=self.field_capacity,
                                           whc=self.whc, rf_coeff=self.rf_coeff)
        DDrainout = f'dd_{year}{DOY}.tif'
        SRfout = f'srf_{year}{DOY}.tif'
        if daily_mode:
            self.rmanager.output_rasters(DDrain, self.outdir, outname=DDrainout)
            self.rmanager.output_rasters(SRf, self.outdir, outname=SRfout)

        # output eta and SWf
        etasw, SWf, etasw5 = self._veg_et(k_factor, ndvi_factor, water_factor, bias_corr, alfa_factor, watermask,
                                          self.pet, self.ndvi, SWi)
        etaswout = f'etasw_{year}{DOY}.tif'
        SWfout = f'swf_{year}{DOY}.tif'
        etasw5out = f'etasw5_{year}{DOY}.tif'
        if daily_mode:
            self.rmanager.output_rasters(etasw, outdir, outname=etaswout)
            self.rmanager.output_rasters(SWf, outdir, outname=SWfout)
            self.rmanager.output_rasters(etasw5, outdir, outname=etasw5out)

        return SWf, SNWpk, etasw, DDrain, SRf

    def run_veg_et(self):
        print(
            '''             _ _            ___  ___  _  _ 
            | | | ___  ___ | __>|_ _|| || |
            | ' |/ ._>/ . || _>  | | |_/|_/
            |__/ \___.\_. ||___> |_| <_><_>
                       <___'                '''
        )

        start_dt = datetime.strptime("{}-{:03d}".format(self.start_year, self.start_day), '%Y-%j')
        print(start_dt)
        end_dt = datetime.strptime("{}-{:03d}".format(self.end_year, self.end_day), '%Y-%j')
        print(end_dt)
        time_interval = end_dt - start_dt
        num_days = time_interval.days
        print(num_days)

        accumulate_mode = self.accumulate_mode

        # initially set output_yearly_arrays and output_monhly array to False and you will change
        # them later depending on what is in the accumulate_mode list
        # todo - set these in config_dict.
        output_monthly_arr = False
        output_yearly_arr = False
        # step daily. It is false if not included by default.
        output_daily_arr = False
        output_daily_arr = True

        # Open static inputs and normalize them to standard numpy arrays

        # static inputs
        self.interception = self.pmanager.get_static_data(self.interception_settings)
        self.whc = self.pmanager.get_static_data(self.whc_settings)
        self.field_capacity = self.pmanager.get_static_data(self.field_capacity_settings)
        self.saturation = self.pmanager.get_static_data(self.saturation_settings)
        self.watermask = self.pmanager.get_static_data(self.watermask_settings)
        # package as a list
        static_inputs = [self.interception, self.whc, self.field_capacity, self.saturation, self.watermask]
        # normalizing.
        self.interception, self.whc, self.field_capacity, self.saturation, self.watermask \
            = self.rmanager.normalize_to_std_grid(inputs=static_inputs,resamplemethod='nearest')


        # set monthly and yearly cumulative arrays (use one of the numpys from the
        # static array that has been normalized):
        model_arr_shape = self.interception.shape
        # A total of six output arrays must be instantiated in case accumulate_mode != None
        # monthly
        et_month_cum_arr = np.zeros(model_arr_shape)
        dd_month_cum_arr = np.zeros(model_arr_shape)
        srf_month_cum_arr = np.zeros(model_arr_shape)
        # yearly
        et_yearly_cum_arr = np.zeros(model_arr_shape)
        dd_yearly_cum_arr = np.zeros(model_arr_shape)
        srf_yearly_cum_arr = np.zeros(model_arr_shape)

        # the soil water fraction and snowpack are none to start out.
        changing_swf = None
        changing_snwpck = None
        for i in range(num_days + 1):
            # so what day is it
            today = start_dt + timedelta(days=i)
            if i == 0:

                swf, snwpck, etasw, DDrain, SRf = self._run_water_bal(i, today, self.interception, self.whc, self.field_capacity,
                                                                      self.saturation, self.rf_coeff, self.k_factor,
                                                                      self.ndvi_factor, self.water_factor, self.bias_corr,
                                                                      self.alfa_factor, self.watermask,
                                                                      outdir=self.outdir, yest_snwpck=None, yest_swf=None,
                                                                      geoproperties_file=self.geoproperties_file, daily_mode=output_daily_arr)
                changing_swf = swf
                changing_snwpck = snwpck

            else:

                # see if today is a day that we need to output a monthly raster
                if 'monthly' in accumulate_mode:
                    d = today.day
                    mo = today.month
                    yr = today.year
                    output_monthly_arr = self._end_of_month(d, mo, yr)

                if 'yearly' in accumulate_mode:
                    # todo - deal with Water Year mode later

                    # this function does calendar years
                    d = today.day
                    mo = today.month

                    if d == 31 and mo == 12:
                        output_yearly_arr = True
                else:
                    output_yearly_arr = False

                print('output monthly is {} and output yearly is {}'.format(output_monthly_arr, output_yearly_arr))

                swf, snwpck, etasw, DDrain, SRf = self._run_water_bal(i, today, self.interception, self.whc,
                                                                      self.field_capacity, self.saturation,
                                                                      self.rf_coeff, self.k_factor, self.ndvi_factor,
                                                                      self.water_factor, self.bias_corr, self.alfa_factor,
                                                                      self.watermask, outdir=self.outdir, yest_snwpck=changing_snwpck,
                                                                      yest_swf=changing_swf, geoproperties_file=self.geoproperties_file,
                                                                      daily_mode=output_daily_arr)

                # monthly
                et_month_cum_arr += etasw
                dd_month_cum_arr += DDrain
                srf_month_cum_arr += SRf
                # yearly
                et_yearly_cum_arr += etasw
                dd_yearly_cum_arr += DDrain
                srf_yearly_cum_arr += SRf

                if output_monthly_arr:
                    # function to create monthly output rasters for each variable
                    self.rmanager.output_rasters(et_month_cum_arr, self.outdir,
                                   'model_outputs/etasw_{}{:02d}.tif'.format(today.year, today.month))
                    self.rmanager.output_rasters(dd_month_cum_arr, self.outdir,
                                   'model_outputs/dd_{}{:02d}.tif'.format(today.year, today.month))
                    self.rmanager.output_rasters(srf_month_cum_arr, self.outdir,
                                   'model_outputs/srf_{}{:02d}.tif'.format(today.year, today.month))

                    # zero-out arrays to start the next month over.
                    et_month_cum_arr = np.zeros(model_arr_shape)
                    dd_month_cum_arr = np.zeros(model_arr_shape)
                    srf_month_cum_arr = np.zeros(model_arr_shape)
                    output_monthly_arr = False

                if output_yearly_arr:
                    # function to create yearly output rasters for each variables
                    self.rmanager.output_rasters(et_yearly_cum_arr, self.outdir, 'model_outputs/etasw_{}.tif'.format(today.year))
                    self.rmanager.output_rasters(dd_yearly_cum_arr, self.outdir, 'model_outputs/dd_{}.tif'.format(today.year))
                    self.rmanager.output_rasters(srf_yearly_cum_arr, self.outdir, 'model_outputs/srf_{}.tif'.format(today.year))

                    # zero-out arrays to start the next year over.
                    et_yearly_cum_arr = np.zeros(model_arr_shape)
                    dd_yearly_cum_arr = np.zeros(model_arr_shape)
                    srf_yearly_cum_arr = np.zeros(model_arr_shape)
                    output_yearly_arr = False

                changing_swf = swf
                changing_snwpck = snwpck

            print('-------------------------------')


class VegETAnalysis:
    """A class that does some standard plotting rendering of the data"""
    pass


        
if __name__ == '__main__':


    print('need to implement a main function')
    # print(datetime.now())
    #
    #
    # # todo - handle leap years, handle bad inputs, need to make command line inputs for all params? returns soils automatically for veg et
    # veg_model = VegET(veget_config_path=r'attributes.yml')
    #
    #
    #
    # # run Veg ET model, Acceptable accumulate modes: 'monthly', 'yearly'
    # veg_model.run_veg_et(start_year, end_year, start_day, end_day,
    #            interception, whc, field_capacity, saturation,
    #            rf_coeff, k_factor, ndvi_factor, water_factor, bias_corr, alfa_factor, watermask,
    #            geo_dict=None, geoproperties_file=geoproperties_file, outdir='',
    #            accumulate_mode=['daily', 'monthly', 'yearly'])
    # print(datetime.now())



