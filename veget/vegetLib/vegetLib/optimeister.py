from .log_logger import log_make_logger
from timeit import default_timer as t_now

import rasterio.mask
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio import shutil as rio_shutil
from rasterio.vrt import WarpedVRT

import numpy as np
from s3fs.core import S3FileSystem

def _np_load_cloud(bucket_file):

    s3 = S3FileSystem()

    ARY = np.load(s3.open(bucket_file))

    return ARY

def _np_save_cloud(bucket_file, ary):
    # pass
    # with s3.open(bucket + path, “w”) as f:
    # more code needed here with save
    s3 = S3FileSystem()
    np.save(s3.open(bucket_file,"wb"), ary)
    
def _make_npy_cache_name(geotiff_fullname, tile):
    print(geotiff_fullname)
    new_name = '/'.join(geotiff_fullname.split('/')[3:]) + '.npy'
    new_name = 'dev-et-data/cache/' + tile + '/' + new_name
    return new_name


class OptiMeister:
    """
    This class is an experimental optimization class to create cached numpy arrays from input geotiff
    the static inputs are re-read evey year and so perhaps numpy arrays are faster than geotiff warps.
    """
    

    def __init__(self, config_dict, shp=None):
         self.log = log_make_logger('OptiMeister')
         self.tile = config_dict['tile']
         self.cache = {}
    
    def _is_in_cache(self, warpfile):
        if warpfile in self.cache:
            return True
        else:
            return False
        
    
    def _cache_npy(self, warpfile, ary):
        self.log.info('Research Cache this file {}'.format(warpfile))
        cache_name = _make_npy_cache_name(warpfile, self.tile)
        _np_save_cloud(cache_name, ary)
        self.cache[warpfile]=cache_name
        self.log.info('NPY Name {}'.format(cache_name))
    
    def _return_cache_data(self, warpfile):
        self.log.info('Getting Cache this file {}'.format(warpfile))
        cache_name = _make_npy_cache_name(warpfile, self.tile)
        data = _np_load_cloud(cache_name)
        return data
        
            
    def o_warp_one(self, warpfile, rs, crs, transform, rows, cols):
        t0 = t_now()
        
        if self._is_in_cache(warpfile):
            self.log.info('RESEARCH RETRIEVING NPY CACHE ITEM'.format(warpfile))
            data = self._return_cache_data(warpfile)
            return data
        else:    
            cnt=10
            while(cnt>0):
                try:
                    with rasterio.open(warpfile) as src:
                        # create the virtual raster based on the standard rasterio attributes from the sample tiff and shapefile feature.
                        with WarpedVRT(src, resampling=rs,
                               crs=crs,
                               transform=transform,
                               height=rows,
                               width=cols) as vrt:
                            data = vrt.read(1)
                            # print(type(vrt))
                            print("data shape =", data.shape)
                            self.log.info("o_warp_one Completed {}".format(warpfile))
                            t_total = t_now() - t0
                            self.log.info("WARP - TIME - {} - {}".format(t_total, warpfile))
                            if 'NDVI' in warpfile:
                                t0 = t_now()
                                self._cache_npy(warpfile,data)
                                t_total = t_now() - t0
                                self.log.info("Cache_Store - TIME - {} - {}".format(t_total, warpfile))
                        return data
                except rasterio.errors.RasterioIOError:
                        print("Unexpected error:", sys.exc_info()[0])
                        print('oops',cnt)
                        cnt = cnt - 1
                        time.sleep(4)
                    
