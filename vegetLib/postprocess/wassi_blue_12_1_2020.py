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
import rasterio
import fiona
import rasterio.mask as msk
# ============= standard library imports ========================
os.environ['PROJ_LIB'] = r'C:\Users\gparrish\AppData\Local\conda\conda\envs\gdal_env\Library\share\proj'
os.environ['GDAL_DATA'] = r'C:\Users\gparrish\AppData\Local\conda\conda\envs\gdal_env\Library\share\epsg_csv'

output = r'Z:\Users\Kul\projects\5_green_blue_water_supply_stress_index\WaSSI_green_blue\WaSSI_blue_irrigated'
irrigated_et = r'Z:\Users\Kul\projects\5_green_blue_water_supply_stress_index\et\irrigated'
runoff_file = r'Z:\Users\Kul\projects\5_green_blue_water_supply_stress_index\WaSSI_green_blue\WaSSI_blue_irrigated\runoff.csv'

# 2003-2019
study_years = [2003+i for i in range(17)]
print(study_years)

sample_file = 'et_irrigated_{}.tif'

basins_folder = r'Z:\Users\Kul\projects\5_green_blue_water_supply_stress_index\basins\basins_world_wri'
basins_lst = ['Danube', 'Ganges-Brahmaputra', 'Mississippi', 'Murray_Darling', 'Nile', 'Sao_Francisco']


for y in study_years:
    print('study year: ', y)
    fname = sample_file.format(y)
    fpath = os.path.join(irrigated_et, fname)

    for basin in basins_lst:
        print('basin!, ', basin)

        # Todo - Find the appropriate value in the runoff file
        with open(runoff_file, 'r') as rfile:

            # for line in rfile:
            #     print(line)

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
                    basin_year = int(line_lst[0])
                    runoff_val = int(line_lst[2])
                    if y == basin_year and basin_name in basin:
                        print('match made! \n', basin_name, basin_year, runoff_val)
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
                                 "transform": out_transform})

                # todo - divide clipped raster by the runoff scalar
                arr /= scalar_value

                # todo - output file with date and basin in the name to output location
                with rasterio.open(os.path.join(output, f'wassi_blue_{basin}_{y}.tif'), 'w', **out_meta) as dest:
                    dest.write(arr)
