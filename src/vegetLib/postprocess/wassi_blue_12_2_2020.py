# ===============================================================================
# Copyright 2019 Gabriel Parrish
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ===============================================================================
import os
import numpy as np
import rasterio
import fiona
import rasterio.mask as msk
# ============= standard library imports ========================
# os.environ['PROJ_LIB'] = r'C:\Users\gparrish\AppData\Local\conda\conda\envs\gdal_env\Library\share\proj'
# os.environ['GDAL_DATA'] = r'C:\Users\gparrish\AppData\Local\conda\conda\envs\gdal_env\Library\share\epsg_csv'

output = r'Z:\Users\Kul\projects\5_green_blue_water_supply_stress_index\WaSSI_green_blue\WaSSI_blue_irrigated'
irrigated_et = r'Z:\Users\Kul\projects\5_green_blue_water_supply_stress_index\et\irrigated'
runoff_file = r'Z:\Users\Kul\projects\5_green_blue_water_supply_stress_index\WaSSI_green_blue\WaSSI_blue_irrigated\runoff.csv'

# 2003-2019
study_years = [2003+i for i in range(17)]
print(study_years)

sample_file = 'et_irrigated_{}.tif'

basins_folder = r'Z:\Users\Kul\projects\5_green_blue_water_supply_stress_index\basins\basins_world_wri'
basins_lst = ['Mississippi', 'Murray_Darling', 'Sao_Francisco'] #'Danube', 'Ganges_Brahmaputra','Nile',

# creating keys with a nul value to hold the cumulative array for 'mean()' calculations
basin_dict = {}
for b in basins_lst:
    basin_dict[b] = None

for ii, y in enumerate(study_years):
    print('study year: ', y)
    fname = sample_file.format(y)
    fpath = os.path.join(irrigated_et, fname)

    for basin in basins_lst:
        print('basin!, ', basin)

        # Todo - Find the appropriate value in the runoff file
        with open(runoff_file, 'r') as rfile:

            for i, line in enumerate(rfile):
                # print(i, rfile)
                # get rid of newline character
                line = line.strip('\n')
                if i == 0:
                    # first line is headers and we don't care
                    pass
                else:
                    line_lst = line.split(',')
                    # get the year, basin and runoff value
                    basin_name = line_lst[1]
                    try:
                        basin_name = basin_name.split('-')[0]
                    except:
                        pass
                    try:
                        basin_name = basin_name.split('_')[0]
                    except:
                        pass
                    try:
                        basin_name = basin_name.split(' ')[0]
                    except:
                        pass
                    basin_year = int(line_lst[0])
                    runoff_val = int(line_lst[2])
                    # print(f'basin name {basin_name}')
                    if y == basin_year and basin_name in basin:
                        # print('match made! \n', basin_name, basin_year, runoff_val)
                        scalar_value = runoff_val

            bpath = os.path.join(basins_folder, f'{basin}.shp')
            # open shapefile to clip
            with fiona.open(bpath, 'r') as shapefile:
                shapes = [feature['geometry'] for feature in shapefile]
                with rasterio.open(fpath) as src:
                    arr, out_transform = msk.mask(src, shapes, crop=True)
                    out_meta = src.meta
                    out_meta.update({"driver": "GTiff",
                                     "height": arr.shape[1],
                                     "width": arr.shape[2],
                                     "transform": out_transform,
                                     'dtype': 'float64'})
                    nd_val = out_meta['nodata'] / float(scalar_value)
                    out_meta.update({'nodata': nd_val})
                    arr = arr / float(scalar_value)

                    # accumulate the array to calculate the mean l8er on.
                    if ii == 0:
                        # set the nodata value from the first nd_val at ii==0
                        basin_dict[f'{basin}_nd'] = nd_val
                        # instantiate the array from None to [basin] = arrarr
                        basin_dict
                        # assign the metadata for the cumulative raster
                        basin_dict[f'{basin}_geo'] = out_meta

                        base_offset = (arr.shape[1], arr.shape[2])
                    else:
                        # set the nd values to the original nd in the array
                        cum_arr = basin_dict[basin]
                        # set to the original nodata
                        cum_arr[cum_arr == nd_val] = basin_dict[f'{basin}_nd']
                        cum_arr[cum_arr == basin_dict[f'{basin}_nd']] = 0.0
                        # that way when the cumulative array is written, the correct nd value is still the same
                        # for each no data pixel.
                        # -- accumulate --
                        cum_arr += arr
                        basin_dict[basin] = cum_arr
                        print('accumulation succesful for shape ', arr.shape, 'basin', basin)

                    # todo - output file with date and basin in the name to output location
                    with rasterio.open(os.path.join(output, f'wassi_blue_{basin}_{y}.tif'), 'w', **out_meta) as dest:
                        dest.write(arr)
# the final number of years is the last value assigned to ii
years_count = ii
# loop through the dictionary and get the mean and write it out.
for basin in basins_lst:
    print('basin!, ', basin)
    cumulative_array = basin_dict[basin]
    mean_arr = cumulative_array / float(years_count)
    metadata = basin_dict[f'{basin}_geo']
    with rasterio.open(os.path.join(output, f'wassi_blue_{basin}_mean_{study_years[0]}_{study_years[-1]}.tif'),
                       'w', **metadata) as dest:
        dest.write(mean_arr)

