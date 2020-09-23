import os
from time import time
import xarray as xr
import rioxarray
from rasterstats import zonal_stats
import fiona
import rasterio
import rasterio.mask as msk
from affine import Affine
import numpy as np

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
    print(ds)
    print(f'OUTPUT=={out_path}')
    ds.rio.to_raster(out_path)


def zonal_stats(shape_path, raster_path, csv_output, parameter,
                dict={'name': [], 'mean': [], "Year": [], "DOY": [], "parameter": [], 'id': []}):
    # === Zonal Stats ===
    with fiona.open(shape_path, 'r') as shp:

        features = [feature for feature in shp]
        for i, f in enumerate(features):
            print(f'feat dict {i} \n ', f, '\n')
        # shapes = [feature['geometry'] for feature in shp]

        with rasterio.open(mosaic_file) as src:

            for f in features:
                src_shp = [f['geometry']]
                watershed_id = f['id']
                watershed_name = f['properties']['Name']
                parameter = statsvar
                # writer.writerow(["Year", "DOY", "parameter", "zon_mean_forID1", "zon_mean_forID2", "zon_mean_forID3", ...])
                outimage, out_transform = msk.mask(src, src_shp, crop=True)
                ws_mean = np.nanmean(outimage)
                ws_year = filename.split('.')[0][-7:-3]
                ws_doy = filename.split('.')[0][-3:]
                dict['name'].append(watershed_name)
                dict['mean'].append(ws_mean)
                dict['Year'].append(ws_year)
                dict['DOY'].append(ws_doy)
                dict['parameter'].append(parameter)
    return dict


statsvar = 'srf'
filename = f'{statsvar}_2015185.tif'
raster_path_tile6 = r'Z:\Projects\VegET_SSEBopET\CONUS2015_WaterSMART_LocalRun\runconfig6\{}'.format(filename)
raster_path_tile13 = r'Z:\Projects\VegET_SSEBopET\CONUS2015_WaterSMART_LocalRun\runconfig13\{}'.format(filename)
outdir = r'\\IGSKMNCNFS016\watersmartfs1\Projects\DRB\mosaic_test'
if not os.path.exists(outdir):
    os.mkdir(outdir)

shape_path = r'\\IGSKMNCNFS016\watersmartfs1\Projects\DRB\Sanford approach\Small Watersheds 128_Aggregated.shp'



# ===== mosaic stuff ====
mosaic_file = os.path.join(outdir, filename)
# mosaic_file = os.path.join(outdir, 'srf_2015185_unity.tif')
# mosaic_file = r'Z:/Projects/VegET_SSEBopET/CONUS2015_WaterSMART_LocalRun/runconfig6/srf_2015185.tif'
if not os.path.exists(mosaic_file):
    print('mosaicking')
    molist = [raster_path_tile6, raster_path_tile13]
    delaware_ds, topleft_affine = xr_build_mosaic_ds(molist)
    # todo - plotting not working as of yet
    # xr.plot.imshow(delaware_ds)
    xr_write_geotiff_from_ds(delaware_ds['band'], out_path=mosaic_file)




