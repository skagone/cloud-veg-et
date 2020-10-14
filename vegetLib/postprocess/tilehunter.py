"""A script that can take an input directory and 'hunt' down the
filepaths that make up a full tile, mosaic them, write them temporarily
 and then sample the data"""
import os
import rasterio
import rasterio.mask as msk
import fiona
import numpy as np
import rioxarray
import xarray as xr
from time import time
import pandas as pd


def cum_mm_to_m3(raster_dim_meters, arr):
    # scrap too large datasets
    arr[arr>=3000] = np.nan
    # convert mm to mm* m^2
    arr *= raster_dim_meters ** 2
    # mm to m conversion
    arr /= 1000

    # return m3
    return arr

def xr_build_mosaic_ds(tifs, product='band'):
    start = time()
    my_da_list = []
    topleft_affine = None
    for i, tif in enumerate(tifs):
        with rasterio.open(tif) as src:
            affine = src.transform
            if i == 0:
                topleft_affine=affine
            affine_lst = list(affine)
            print('affine list', affine_lst)
            print('affine e', affine.e)
        da = xr.open_rasterio(tif)
        da = da.squeeze().drop(labels='band')
        da.name = product
        my_da_list.append(da)
        tnow = time()
        elapsed = tnow - start
        print(tif, elapsed)
    ds = xr.merge(my_da_list)
    print('the type ', type(ds))
    print('topleft affine, ', list(topleft_affine))
    print('topleft affine e, ', topleft_affine.e)
    print('the transform tttt', list(ds.rio.transform()))
    ds.rio.write_transform(transform=topleft_affine, inplace=True)
    print('the transform zzzz', list(ds.rio.transform(recalc=False)))
    return ds, topleft_affine

def xr_write_geotiff_from_ds(ds, out_path):
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    out_path = os.path.join(out_path, 'temp.tif')
    print(ds)
    print(f'OUTPUT=={out_path}')
    ds.rio.to_raster(out_path)

def zonal_stats(shape_path, raster_path, parameter,
                year=None, month=None, day=None,
                dict={'cubic_m': [], 'mean': [], "Year": [], "DOY": [], "parameter": [], 'id': []},
                mm_flux=True, raster_dim=None):
    # === Zonal Stats ===
    with fiona.open(shape_path, 'r') as shp:
        features = [feature for feature in shp]
        for i, f in enumerate(features):
            print(f'feat dict {i} \n ', f, '\n')
        # shapes = [feature['geometry'] for feature in shp]
        with rasterio.open(raster_path) as src:
            for f in features:
                src_shp = [f['geometry']]
                watershed_id = f['id']
                # watershed_name = f['properties']['NAME']
                # writer.writerow(["Year", "DOY", "parameter", "zon_mean_forID1", "zon_mean_forID2", "zon_mean_forID3", ...])
                outimage, out_transform = msk.mask(src, src_shp, crop=True)
                ws_mean = np.nanmean(outimage)
                if mm_flux:
                    m3_arr = cum_mm_to_m3(raster_dim, outimage)
                    basin_ro = np.nansum(m3_arr)
                ws_year = year
                ws_doy = day
                dict['mean'].append(ws_mean)
                dict['Year'].append(ws_year)
                dict['DOY'].append(ws_doy)
                dict['cubic_m'].append(basin_ro)
                dict['parameter'].append(parameter)
                dict['id'].append(None)
    return dict

def hunt_tile(dir, var, freq, year=None, month=None, doy=None):
    """returns a list of tile paths that will get mosaicked
    together with xr_build_mosaic_ds"""

    year_format = '{}_{}.tif'
    mo_format = '{}_{}{:02d}.tif'
    day_format = '{}_{}{:03d}.tif'

    if freq == 'yearly':
        search_str = year_format.format(var, year)
    elif freq == 'monthly':
        search_str = mo_format.format(var, year, month)
    elif freq == 'daily':
        search_str = day_format.format(var, year, doy)

    else:
        print('needs to be annual, monthly or daily freq argument')

    tiles_lst = []
    for p, d, files in os.walk(dir):

        for f in files:
            if f == search_str:
                tilepath = os.path.join(p, f)
                tiles_lst.append(tilepath)

    return tiles_lst

def tilestats(var, freq, year, tiles, stat_dict=None, iteration=0):
    print('tiles \n', tiles)



    ds, trans = xr_build_mosaic_ds(tifs=tiles)

    print(type(ds))
    print(trans)

    xr_write_geotiff_from_ds(ds=ds, out_path=tempdir)

    if i == 0:
        sample_data = zonal_stats(shape_path=sample_shape, raster_path=tempfile,
                                  parameter=var, year=year, raster_dim=1000)
    else:
        sample_data = zonal_stats(shape_path=sample_shape, raster_path=tempfile,
                                  parameter=var, year=year, dict=stat_dict, raster_dim=1000)

    return sample_data


if __name__ == "__main__":
    basin = 'MurrayDarling'
    output_root = r'Z:\Projects\VegET_Basins\Murray'
    sample_shape = r'Z:\Projects\VegET_Basins\Murray\shapefiles\Murray_Darling.shp'
    tempdir = os.path.join(output_root, 'temp')
    tempfile = os.path.join(tempdir, 'temp.tif')

    var = 'srf'
    freq = 'yearly'
    yrs_list = [2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010,
                2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019]

    for i, yr in enumerate(yrs_list):
        tiles = hunt_tile(dir=output_root, var=var, freq=freq, year=yr)

        if i == 0:
            stats_dict = tilestats(var=var, freq=freq, year=yr, tiles=tiles)
        else:
            stats_dict = tilestats(var=var, freq=freq, year=yr, tiles=tiles, stat_dict=stats_dict, iteration=i)

    print('stats dict\n', stats_dict)
    for k, v in stats_dict.items():
        print(k)
        print(len(v))

    # output the dictionary
    cols = list(stats_dict.keys())
    # lalala
    nice_dataframe = pd.DataFrame(stats_dict, columns=cols)

    data_out = os.path.join(output_root, 'sampled_files')
    if not os.path.exists(data_out):
        os.mkdir(data_out)

    nice_dataframe.to_csv(path_or_buf=os.path.join(data_out,
                                                   '{}_sampled_{}_{}_{}.csv'.format(basin, var,
                                                                                 stats_dict['Year'][0],
                                                                                 stats_dict['Year'][-1])))