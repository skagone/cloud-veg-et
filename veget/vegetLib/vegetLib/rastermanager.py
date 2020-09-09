import os
import sys
import time
import yaml
import calendar
from datetime import datetime, timedelta, date

from s3fs.core import S3FileSystem
import boto3
import fiona

import pandas as pd
import rasterio.mask
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio import shutil as rio_shutil
from rasterio.vrt import WarpedVRT
from timeit import default_timer as t_now

from .pathmanager import PathManager

from .box_poly import box_create_ugly_proprietary_shapefile_plus_json_from_tile
from .log_logger import log_make_logger
from .optimeister import OptiMeister


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


    def __init__(self, config_dict, shp=None):
        self.optimize = False
        self.log = log_make_logger('COOL_RASTERMANAGER')

        self.config_dict = config_dict

        tile = self.config_dict['tile']

        self.log.info('Tony-s stupid - tile name is - {}'.format(tile))
        
        if 'tile' in tile:
            self.log.info("using scalable tile names {}".format(tile))

            #bucket_name = self.config_dict['out_root'].split('/')[0]
            #today = date.today()
            #print("Current date =", today)
            #date_str=today.strftime("%m_%d_%Y")
            #self.config_dict['out_root'] = bucket_name + '/out/DelawareRiverBasin/Run' + date_str + '/' + tile

            if self.config_dict['optimize']:
                self.optimize = True
                self.opti=OptiMeister(config_dict,shp)

        # self.geoproperties_file = config_dict.geoproperties_file
        # self.shapefile = config_dict.shapefile
        # self.temp_folder = os.path.join(config_dict.out_root, config_dict.temp_folder)

        # self.temp_folder = config_dict['temp_folder']
        self.temp_folder = './' + tile
        self.log.info('temp folder is'.format(self.temp_folder))

        if not os.path.exists(self.temp_folder):
            os.makedirs(self.temp_folder)

        # if the user does not include a shapefile in VegET, a box based on the tile name will be created.
        if shp == None:
            self.shapefile = box_create_ugly_proprietary_shapefile_plus_json_from_tile(self.temp_folder, tile)
        else:
            self.shapefile = shp

        self.geoproperties_file = config_dict['geoproperties_file']

        if self.geoproperties_file == None or self.shapefile==None:
            print('Assuming the user entered values in the config_dict for boundaries of the AOI not implemented at thsi time')
            sys.exit(0)

    # ----------- create output rasters -----------------
    def output_rasters(self, arr, outdir, outname):
        """
        This function creates geotiff files from the model output arrays.
        """
        if self.config_dict['path_mode'] == 'local':
            outpath = os.path.join(outdir, outname)
            print('the outpath for file {} is {}'.format(outname, outpath))

            band1 = arr
            with rasterio.open(outpath, 'w', driver='GTiff', height=self.rows, width=self.cols,
                               count=1, dtype='float64', crs=self.crs, transform=self.transform) as wrast:
                wrast.write(band1, indexes=1)

        elif self.config_dict['path_mode'] == 'aws':
            # later on deleted by s3_delete_local()
            # local_outpath = os.path.join(self.config_dict['temp_folder'], outname)
            local_outname = outname.split('/')[-1]
            local_outpath = os.path.join(self.temp_folder, local_outname)
            self.log.debug('local_outpath {}'.format(local_outpath))

            t0 = t_now()

            band1 = arr
            # write to a temp folder
            with rasterio.open(local_outpath, 'w', driver='GTiff', height=self.rows, width=self.cols,
                               count=1, dtype='float64', crs=self.crs, transform=self.transform) as wrast:
                wrast.write(band1, indexes=1)

            # Buckets are not directories but you can treat them like they are
            # bucket_name = os.path.split(self.config_dict['out_root'])[0]     # dev-et-data
            # bucket_prefix = os.path.split(self.config_dict['out_root'])[-1]  # tile_modelrun1
            bucket_name = self.config_dict['out_root'].split('/')[0]
            bucket_prefix_list = self.config_dict['out_root'].split('/')[1:]
            print(bucket_prefix_list)
            bucket_prefix = '/'.join(bucket_prefix_list)
            print("bucket prefix =", bucket_prefix)
            bucket_filepath = os.path.join(bucket_prefix, outname)   # os.path.join(dev-et-data/tile_modelrun1, outname)

            # uploads to aws bucket with filepath
            self.s3_delete_local(local_file=local_outpath, bucket=bucket_name, bucket_filepath=bucket_filepath)
            t_total = t_now() - t0
            self.log.info("OUTPUT - TIME - {} - {}".format(t_total, bucket_filepath))


    def set_model_std_grid(self, feat=0):
        """Clips and crops a tiff to the extent of a feature in a shapefile
        :param feat: feat is  the feature id of the shapefile from like a GeoJSON)
        # https://rasterio.readthedocs.io/en/latest/topics/virtual-warping.html
        """
        # print(self.shapefile)
        with fiona.open(self.shapefile, 'r') as shapefile:
            # todo - set up an error if user has shapefile with more than one feature. GELP n STEFFI
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
        # TODO - Set Blocksize for sample raster and other useful optimization thingys
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
            # print('warpfile', warpfile)
            with rasterio.open(warpfile) as src:
                # create the virtual raster based on the standard rasterio attributes from the sample tiff and shapefile feature.
                with WarpedVRT(src, resampling=rs,
                               crs=self.crs,
                               transform=self.transform,
                               height=self.rows,
                               width=self.cols) as vrt:
                    data = vrt.read()
                    # print(type(vrt))
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

    def _warp_one(self, warpfile, rs):
        t0 = t_now()
        cnt=10
        while(cnt>0):
            try:
                with rasterio.open(warpfile) as src:
                    # create the virtual raster based on the standard rasterio attributes from the sample tiff and shapefile feature.
                    with WarpedVRT(src, resampling=rs,
                           crs=self.crs,
                           transform=self.transform,
                           height=self.rows,
                           width=self.cols) as vrt:
                        data = vrt.read(1)
                        # print(type(vrt))
                        print("data shape =", data.shape)
                        self.log.info("_warp_one Completed {}".format(warpfile))
                        t_total = t_now() - t0
                        self.log.info("WARP - TIME - {} - {}".format(t_total, warpfile))
                    return data
            except rasterio.errors.RasterioIOError:
                    print("Unexpected error:", sys.exc_info()[0])
                    print('oops',cnt)
                    cnt = cnt - 1
                    time.sleep(4)

    def _warp_inputs(self, inputs, resamplemethod):

        self.log.info("_warp_inputs")
        outputs = []
        npy_outputs = []
        if resamplemethod == 'nearest':
            rs = Resampling.nearest
        else:
            print('only nearest neighbor resampling is supported at this time')
            sys.exit(0)

        for i, warpfile in enumerate(inputs):
            print('warpfile', warpfile)
            if (self.optimize):
                data = self.opti.o_warp_one(warpfile, rs, self.crs, self.transform, self.rows, self.cols)
            else:
                data = self._warp_one(warpfile, rs)
            npy_outputs.append(data)
        return npy_outputs

    def normalize_to_std_grid_fast(self, inputs, resamplemethod='nearest'):
        """
        Uses rasterio virtual raster to standardize grids of different crs, resolution, boundaries based on  a shapefile geometry feature
        :param inputs: a list of (daily) raster input files for the water balance.
        :param outloc: output locations 'temp' for the virtual files
        :return: list of numpy arrays
        """

        npy_outputs = self._warp_inputs(inputs, resamplemethod)

        return npy_outputs



    def s3_delete_local(self, local_file, bucket, bucket_filepath):
        """
        This function will move the model outputs from a local folder to a cloud bucket.
        :param local_file: path the the local geo file
        :param outpath: path of a directory to be created in the cloud bucket
        :param bucket: name of the cloud bucket = 'dev-et-data'
        :param bucket_folder: "folder" in cloud bucket  = 'v1DRB_outputs'
        :return:
        """

        s3 = boto3.client('s3')
        with open(local_file, "rb") as f:
            if 'vsis3' in bucket:
                bucket = bucket.split('/')[-1]
                print(bucket, bucket_filepath)
            s3.upload_fileobj(f, bucket, bucket_filepath)
        os.remove(local_file)



