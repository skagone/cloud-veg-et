import os
import rasterio
from rasterio.mask import mask
import numpy as np
import fiona
import pandas as pd
"""
This script is going to extract the annual PET for a given river basin and give average 
mm of PET per year per basin and output the csv to the Netapp
"""

path = r'\\IGSKMNCNFS016\watersmartfs1\Data\ReferenceET\Global\GDAS_PET'
outpath = r'Z:\Users\Kul\projects\5_green_blue_water_supply_stress_index\runoff\VegET\PET_extractions'
yrs_list = [2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010,
                2011, 2012, 2013, 2014, 2015]
gdas_pet_shape = (181, 360)
#================================
areadict = {'Mississippi': 3228840000000, 'Ganges': 1662310000000,
                'Danube': 801840000000, 'Murray': 1049020000000, 'Nile': 3057390000000,
                'San_Francisco': 621621000000}
basin_list = ['Mississippi', 'Ganges', 'Danube', 'Murray', 'Nile', 'San_Francisco']
shapenames = ['Mississippi', 'Ganges_Brahmaputra', 'Danube', 'Murray_Darling', 'Nile', 'Sao_Francisco']

raster_dim_meters = 1000
data_dictionary = {}
# ================================
for basin, shapename in zip(basin_list, shapenames):
    basin_area = areadict[basin]
    output_root = r'Z:\Projects\VegET_Basins\{}'.format(basin)
    sample_shape = r'Z:\Projects\VegET_Basins\{}\shapefiles\{}.shp'.format(basin, shapename)
    'make a nested dict'
    data_dictionary[basin] = {}
    #================================

    for y in yrs_list:
        print(f'==========\n {y} \n==========')
        ypath = os.path.join(path, f'{y}', 'pet')
        print(sample_shape)
        with fiona.open(sample_shape, 'r') as shp:
            print(shp)
            features = [feature for feature in shp]

            print(os.path.join(ypath, f'gdasPET{y}001.tif'))
            with rasterio.open(os.path.join(ypath, f'gdasPET{y}001.tif'), 'r') as src:
                for f in features:
                    src_shp = [f['geometry']]
                    print(shp, src)
                    arr, out_transform = mask(src, src_shp, crop=True)
                    sample_arr_shp = arr.shape
                    print('sample shape {}'.format(sample_shape))
        cum_arr = np.zeros(sample_arr_shp)
        for f in os.listdir(ypath):

            if f.endswith('.tif'):
                # print('das file?')
                # print(f)
                # === Zonal Stats ===
                with fiona.open(sample_shape, 'r') as shp:
                    features = [feature for feature in shp]

                    print(os.path.join(ypath, f))
                    with rasterio.open(os.path.join(ypath, f), 'r') as src:
                        for f in features:
                            src_shp = [f['geometry']]
                            print(shp, src)
                            arr, out_transform = mask(src, src_shp, crop=True)
                            arr[arr >= 10000] = np.nan
                            arr[arr < 0] = np.nan
                            cum_arr += arr
        basin_mean = np.nanmean(cum_arr)
        data_dictionary[basin][f'{y}'] = basin_mean


print(data_dictionary)

for k, vdict in data_dictionary.items():

    fname = f'{k}_PET.csv'

    fout = os.path.join(outpath, fname)

    with open(fout, 'w') as wfile:
        wfile.write('year,pet_mm\n')
        for ky, val in vdict.items():
            wfile.write(f'{int(ky)},{float(val)}\n')





