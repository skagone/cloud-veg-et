import os
from netCDF4 import Dataset
import sys
import gdal
# from gdalconst import GDT_Float32
import numpy as np
from datetime import datetime as dt
from datetime import timedelta
import rasterio
from rasterio.crs import CRS
from rasterio.windows import Window
print('hello world')

def process_netcdf(path, output_path, projection, name_string, main_var='qav', global_dat=True, hugefile=True, chunkparam=None):
    """
    :param path: path to a netcdf4 file
    :return: et arrays, geotransform and dimensions
    """
    main_dataset = Dataset(path, 'r+', format='NETCDF4')
    print('main dataset variables \n', main_dataset.variables, ' \n =====')
    print('main dataset.vars[time] \n', main_dataset.variables['time'] , '\n ======')
    time_list = np.array(main_dataset.variables['time']).tolist()
    print('time array \n', time_list)
    #below getting units out
    units = main_dataset.variables['time'].units
    print('units\n', units)
    start_date_string = units[-19:-9]
    print('start date string \n', start_date_string)
    start_date = dt.strptime(start_date_string, "%Y-%m-%d")
    print('start date datetime', start_date)
    # for i in time_list
    #     dt.
    date_list = []
    for d in time_list:
        datecount = start_date + timedelta(days=d)
        date_list.append(datecount)

    print('here is your date list \n', date_list)
    year_month = []
    for date in date_list:
        #dt.date_list(year, month, day, minutes, seconds)
        date.timetuple()
        month = str(date.month)
        year = str(date.year)
        year_month.append((year, month))
    print(year_month)
    # get the geotransform information: how big pixels are size of rows and colums
    # print 'ssebop dimensions', main_dataset.dimensions
    lat = main_dataset.dimensions['lat']
    latsze = lat.size
    # print 'lat',lat
    print(f'latitude size {latsze}')
    lon = main_dataset.dimensions['lon']
    lonsze = lon.size
    print(f'longitude size {lonsze}')
    # format of dimensions is (x coords, y coords) todo - maybe reverse this for numpy conventions...
    dimensions = (lonsze, latsze)
    # dimensions = (latsze, lonsze)
    print('dimensions we try to use', dimensions)
    # get pixel width/height from geospatial min and max
    # print 'latmin', main_dataset.geospatial_lat_min
    if global_dat:
        # latmin = -60.0
        # latmax = 90.0
        # # testing
        # latmax = -60
        # latmin = -150
        # todo - this is so effed up!!!
        latmax = -90
        height_pixels = (90.0 + 90.0) / latsze
        lonmin = -180
        lonmax = 180
    else:
        latmin = main_dataset.geospatial_lat_min
        latmax = main_dataset.geospatial_lat_max
        print(f'max n min lat: {latmax}, {latmin}')
        height_pixels = (latmax - latmin) / latsze
        # print 'heightpixels', height_pixels
        lonmin = main_dataset.geospatial_lon_min
        lonmax = main_dataset.geospatial_lon_max
    width_pixels = (lonmax - lonmin)/lonsze
    # print 'width pixels', width_pixels
    transform = (lonmin, width_pixels, 0, latmax, 0, -height_pixels) # todo
    print('the transform:', transform)
    # get the entire et array
    # (z(time),y(lat), x(lon) )
    mainvar_dim = (len(time_list), latsze, lonsze)
    mainvar = main_dataset[main_var]
    # mainvar = np.flipud(mainvar)
    # print('main variable\n', mainvar)
    # print(f'====\n main var dims: {mainvar_dim} \n')
    if not hugefile:
        print('Not implemented')
        sys.exit()
    else:
        for ind, tup in enumerate(year_month):
            chunk_gcd = 32
            outputfilename = "{}{}_{}.tif".format(name_string, 'qav', tup[0])
            # write_raster(array, transform, output_path, outputfilename, dimensions, projection, flip_arr=True)
            outfile = os.path.join(output_path, outputfilename)
            with rasterio.open(outfile, 'w', driver='GTiff', height=latsze, width=lonsze,
                           count=1, dtype='float32', crs=CRS.from_epsg(4326), transform=[width_pixels, 0.0, lonmin, 0.0, height_pixels, latmax]) as wfile:

                width_chunk = 8
                height_chunk = 4
                # so that the chunks are square and width is double the height of
                # this raster we calc width using height chunk
                width_chunksize = int(mainvar_dim[1] / height_chunk)
                height_chunksize = int(mainvar_dim[2] / width_chunk)
                for i in range(width_chunk):
                    for j in range(height_chunk):
                        print(f'doing chunk {(i+1)*(j+1)}')
                        if i == 0:
                            width_indices = (i, i + width_chunksize)
                        elif i < width_chunk:
                            width_indices = (i * width_chunksize, (i + 1) * width_chunksize)
                        if j == 0:
                            height_indices = (j, j + height_chunksize)
                        elif j < height_chunk:
                            height_indices = (j * height_chunksize, (j + 1) * height_chunksize)
                        else:
                            pass
                        array = mainvar[ind, height_indices[0]:height_indices[1],
                                width_indices[0]:width_indices[1]]

                        # === convert from m3/s to mm/year/km2 ===
                        # [(mm -> m depth) * (seconds in year)] / (square meters in a square km)
                        raster_conversion_factor = (1000.0 * (60.0 * 60.0 * 24.0 * 365.0)) / 1000000.0
                        raster_conversion_factor = 1
                        # apply the conversion factor
                        array *= raster_conversion_factor
                        # array = np.flipud(array)
                        # for the Window, the first index of the width and height is the offset
                        print('what goes into the window',
                              (width_indices[0], height_indices[0], array.shape[0], array.shape[1]))
                        wfile.write(array, window=Window(width_indices[0], height_indices[0],
                                                         array.shape[0], array.shape[1]), indexes=1)

def run():
    """
    This is the function where the user specifies where to output net cdf as geotiff, where the netcdf file is, and
    what projection were using...
    :return:
    """
    # path to the net cdf file on our computer
    path = "Z:\Data\Runoff\Global\Flow1K\FLO1K.ts.1960.2015.qav.nc"
    projection = 'EPSG:4326'
    output_path = r'Z:\Data\Runoff\Global\Flow1K\FLO1K_ts_1960_2015_qav_gtiff'
    # what the output geotiffs start with
    name_string = 'GRUN_m3s'
    process_netcdf(path, output_path, projection, name_string)

if __name__ == "__main__":
    run()