from veget.vegetLib.vegetLib.veget import VegET
import os

"""
This runner was made on 9/1/2020 from Oct 2014-Dec 2016. This run is being done locally but more runs
 (for the Danube, Ganges-Bramaputra etc) will be run in the cloud if possible...
"""

config_root = r'Z:\Projects\VegET_NileBasin\Nile2016'
shape_root = r'C:\WaterSmart\Users\Gabe\Nile_proj\shapefiles'

tile_lst = ['nile1.shp', 'nile2.shp', 'nile3.shp', 'nile4.shp']
# # testing
# tile_lst = ['nile3.shp']
config_list = ['runconfig1', 'runconfig2', 'runconfig3', 'runconfig4']
# # testing
# config_list = ['runconfig3']

for conf, tile in zip(config_list, tile_lst):
    tilename = tile.split('.')[0]
    shapefile = os.path.join(shape_root, tile)
    config_path = os.path.join(config_root, conf)

    veggie = VegET(veget_config_path=config_path,
                   tile=tilename,
                   shp=shapefile)
    veggie.run_veg_et()