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
            # print('affine list', affine_lst)
            # print('affine e', affine.e)
        da = xr.open_rasterio(tif)
        da = da.squeeze().drop(labels='band')
        da.name = product
        my_da_list.append(da)
        tnow = time()
        elapsed = tnow - start
        # print(tif, elapsed)
    ds = xr.merge(my_da_list)
    # print('the type ', type(ds))
    # print('topleft affine, ', list(topleft_affine))
    # print('topleft affine e, ', topleft_affine.e)
    # print('the transform tttt', list(ds.rio.transform()))
    ds.rio.write_transform(transform=topleft_affine, inplace=True)
    # print('the transform zzzz', list(ds.rio.transform(recalc=False)))
    return ds, topleft_affine

def xr_write_geotiff_from_ds(ds, out_path):
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    out_path = os.path.join(out_path, 'temp.tif')
    # print(ds)
    # print(f'OUTPUT=={out_path}')
    ds.rio.to_raster(out_path)

def zonal_stats(shape_path, raster_path, parameter,
                year=None, month=None, day=None,
                dict={'cubic_m': [], 'mm_basin': [], 'mm_year_mean': [],
                      "Year": [], "DOY": [], "parameter": [], 'id': []},
                mm_flux=True, raster_dim=None, basin_area=None):
    # === Zonal Stats ===
    with fiona.open(shape_path, 'r') as shp:
        features = [feature for feature in shp]
        # for i, f in enumerate(features):
            # print(f'feat dict {i} \n ', f, '\n')
        # shapes = [feature['geometry'] for feature in shp]
        with rasterio.open(raster_path) as src:
            for f in features:
                src_shp = [f['geometry']]
                watershed_id = f['id']
                # watershed_name = f['properties']['NAME']
                # writer.writerow(["Year", "DOY", "parameter", "zon_mean_forID1", "zon_mean_forID2", "zon_mean_forID3", ...])
                outimage, out_transform = msk.mask(src, src_shp, crop=True)
                # scrap too-large datasets
                outimage[outimage >= 3000] = np.nan
                outimage[outimage < 0] = np.nan

                ws_mean = np.nanmean(outimage)
                if mm_flux:
                    m3_arr = cum_mm_to_m3(raster_dim, outimage)
                    basin_ro = np.nansum(m3_arr)
                    m_basin = basin_ro/basin_area
                    mm_basin = m_basin * 1000
                ws_year = year
                ws_doy = day
                # average mm/year across the basin. Should be pretty close to 'mm_year'?
                dict['mm_basin'].append(mm_basin)
                # is the mean the same as the mm_basin?
                dict['mm_year_mean']. append(ws_mean)
                dict['Year'].append(ws_year)
                dict['DOY'].append(ws_doy)
                # the full cubic m of discharge for the basin
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

def tilestats(var, freq, year, tiles, stat_dict=None, basin_area=None, i=None):
    print('tiles \n', tiles)



    ds, trans = xr_build_mosaic_ds(tifs=tiles)

    # print(type(ds))
    # print(trans)

    xr_write_geotiff_from_ds(ds=ds, out_path=tempdir)

    if i == 0:
        sample_data = zonal_stats(shape_path=sample_shape, raster_path=tempfile,
                                  parameter=var, year=year, raster_dim=1000, basin_area=basin_area)
    else:
        sample_data = zonal_stats(shape_path=sample_shape, raster_path=tempfile,
                                  parameter=var, year=year, dict=stat_dict, raster_dim=1000, basin_area=basin_area)

    return sample_data


if __name__ == "__main__":

    # TODO ADD deep drainage to srf to get the total runoff.
    #  Also try eliminating neg values.
    areadict = {'Mississippi': 3228840000000, 'Ganges': 1662310000000,
                'Danube': 801840000000, 'Murray': 1049020000000, 'Nile': 3057390000000,
                'San_Francisco': 621621000000}
    # ================================
    basin = 'Murray'
    shapename = 'Murray_Darling'
    raster_dim_meters = 1000
    # ================================
    stats_dict = {'cubic_m': [], 'mm_basin': [], 'mm_year_mean': [],
                      "Year": [], "DOY": [], "parameter": [], 'id': []}
    basin_area = areadict[basin]
    output_root = r'Z:\Projects\VegET_Basins\{}'.format(basin)
    sample_shape = r'Z:\Projects\VegET_Basins\{}\shapefiles\{}.shp'.format(basin, shapename)
    tempdir = os.path.join(output_root, 'temp')
    if not os.path.exists(tempdir):
        os.mkdir(tempdir)
    tempfile = os.path.join(tempdir, 'temp.tif')

    vars = ['dd', 'srf', 'rain', 'etasw']
    freq = 'yearly'
    yrs_list = [2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010,
                2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019]
    with open(r'D:\Users\gparrish\Desktop\testI.txt', 'w') as wfile:
        pass
    with open(r'D:\Users\gparrish\Desktop\testII.txt', 'w') as wfile:
        pass

    list_of_dicts = []
    for var in vars:
        print(f'var \n {var}')
        sd = stats_dict
        for i, yr in enumerate(yrs_list):
            tiles = hunt_tile(dir=output_root, var=var, freq=freq, year=yr)

            ds, trans = xr_build_mosaic_ds(tifs=tiles)

            xr_write_geotiff_from_ds(ds=ds, out_path=tempdir)

            # === Zonal Stats ===
            with fiona.open(sample_shape, 'r') as shp:
                features = [feature for feature in shp]
                # for i, f in enumerate(features):
                # print(f'feat dict {i} \n ', f, '\n')
                # shapes = [feature['geometry'] for feature in shp]
                with rasterio.open(tempfile) as src:
                    for f in features:
                        src_shp = [f['geometry']]
                        watershed_id = f['id']
                        # watershed_name = f['properties']['NAME']
                        # writer.writerow(["Year", "DOY", "parameter", "zon_mean_forID1", "zon_mean_forID2", "zon_mean_forID3", ...])
                        outimage, out_transform = msk.mask(src, src_shp, crop=True)
                        # scrap too-large datasets
                        outimage[outimage >= 3000] = np.nan
                        outimage[outimage < 0] = np.nan

                        ws_mean = np.nanmean(outimage)
                        # calculating cubic meters of discharge for a basin
                        m3_arr = cum_mm_to_m3(raster_dim_meters, outimage)
                        basin_ro = np.nansum(m3_arr)
                        m_basin = basin_ro / basin_area
                        mm_basin = m_basin * 1000
                        ws_year = yr

                        # average mm/year across the basin. Should be pretty close to 'mm_year'?
                        sd['mm_basin'].append(mm_basin)
                        # is the mean the same as the mm_basin?
                        sd['mm_year_mean'].append(ws_mean)
                        sd['Year'].append(ws_year)
                        # the full cubic m of discharge for the basin
                        sd['cubic_m'].append(basin_ro)
                        sd['parameter'].append(var)

        list_of_dicts.append(sd)

    for sd in list_of_dicts:
        print('=========================\n')
        for k, v in sd.items():
            print(f'===k=== \n {k}\n===v=== \n {v}')
        print('=========================\n')
    # for sd in list_of_dicts:
    #     var = sd['parameter'][0]
    #     for k, v in sd.items():
    #         print(f'===k=== \n {k}\n===v=== \n {v}')
    #
    #     # output the dictionary
    #     cols = list(sd.keys())
    #     # lalala
    #     nice_dataframe = pd.DataFrame(sd, columns=cols)
    #
    #     data_out = os.path.join(output_root, 'sampled_files')
    #     if not os.path.exists(data_out):
    #         os.mkdir(data_out)
    #
    #     nice_dataframe.to_csv(path_or_buf=os.path.join(data_out,
    #                                                    '{}_sampled_{}_{}_{}.csv'.format(basin, var,
    #                                                                                  stats_dict['Year'][0],
    #                                                                                  stats_dict['Year'][-1])))


    # for i, yr in enumerate(yrs_list):
    #     tiles = hunt_tile(dir=output_root, var=vars[i], freq=freq, year=yr)
    #
    #     if i == 0:
    #         stats_dict = tilestats(var=vars[i], freq=freq, year=yr, tiles=tiles, basin_area=basin_area, i=i)
    #     else:
    #         stats_dict = tilestats(var=vars[i], freq=freq, year=yr, tiles=tiles, stat_dict=stats_dict,
    #                                basin_area=basin_area, i=i)
    #
    #     print('stats dict\n', stats_dict)
    #     with open(r'D:\Users\gparrish\Desktop\testII.txt', 'a') as afile:
    #         afile.write(f'{stats_dict}\n')
    #     for k, v in stats_dict.items():
    #         print(k)
    #         print(len(v))
    #
    #     # output the dictionary
    #     cols = list(stats_dict.keys())
    #     # lalala
    #     nice_dataframe = pd.DataFrame(stats_dict, columns=cols)
    #
    #     data_out = os.path.join(output_root, 'sampled_files')
    #     if not os.path.exists(data_out):
    #         os.mkdir(data_out)
    #
    #     nice_dataframe.to_csv(path_or_buf=os.path.join(data_out,
    #                                                    '{}_sampled_{}_{}_{}.csv'.format(basin, vars[i],
    #                                                                                  stats_dict['Year'][0],
    #                                                                                  stats_dict['Year'][-1])))
    #     # reset the stats dict
    #     stats_dict=None