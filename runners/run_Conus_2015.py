from veget.vegetLib.vegetLib.veget import VegET
import os
import rasterio
from matplotlib import pyplot as plt
"""
This runner was made on 8/14/2020 as a demo for G.Senay to show to WaterSMART for year 2015
 for the open-source VegET model using historical temperature data rather than climatological.
 See Config for more info.
"""

config_root = r'Z:\Projects\VegET_SSEBopET\CONUS2015_WaterSMART_LocalRun'
shape_root = r'Z:\Projects\VegET_SSEBopET\shapefiles'

# original run setup as of 8/14/2020
tile_lst = ['CONUS_1.shp', 'CONUS_2.shp', 'CONUS_3.shp', 'CONUS_4.shp', 'CONUS_5.shp', 'CONUS_6.shp', 'CONUS_7.shp',
            'CONUS_8.shp', 'CONUS_9.shp', 'CONUS_10.shp', 'CONUS_11.shp', 'CONUS_12.shp', 'CONUS_13.shp',
            'CONUS_14.shp', 'CONUS_15.shp', 'CONUS_16.shp']
config_list = ['runconfig1', 'runconfig2', 'runconfig3', 'runconfig4', 'runconfig5', 'runconfig6', 'runconfig7',
            'runconfig8', 'runconfig9', 'runconfig10', 'runconfig11', 'runconfig12', 'runconfig13',
            'runconfig14', 'runconfig15', 'runconfig16']

# # ==== VDI was recomposed on 8/21/2020 So tiles 14-16 are reset
# tile_lst = ['CONUS_14.shp', 'CONUS_15.shp', 'CONUS_16.shp']
# config_list = ['runconfig14', 'runconfig15', 'runconfig16']

for conf, tile in zip(config_list, tile_lst):
    tilename = tile.split('.')[0]
    shapefile = os.path.join(shape_root, tile)
    config_path = os.path.join(config_root, conf)

    veggie = VegET(veget_config_path=config_path,
                   tile=tilename,
                   shp=shapefile)
    veggie.run_veg_et()