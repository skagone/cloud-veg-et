from src.vegetLib import VegET
import os

"""
This runner was made on 9/21/2020 from Oct 2002-Dec 2019. This run is being done locally but more runs
 (for the Danube, Ganges-Bramaputra etc) will be run in the cloud if possible...
"""

config_root = r'Z:\Projects\VegET_Basins\Murray'
shape_root = r'Z:\Projects\VegET_Basins\Murray\shapefiles'

tile_lst = ['A1.shp', 'A2.shp', 'A3.shp', 'A4.shp', 'A5.shp']
# # testing
# tile_lst = ['N4.shp']
config_list = ['A1', 'A2', 'A3', 'A4', 'A5']
# # testing
# config_list = ['N4']

for conf, tile in zip(config_list, tile_lst):
    tilename = tile.split('.')[0]
    shapefile = os.path.join(shape_root, tile)
    config_path = os.path.join(config_root, conf)

    veggie = VegET(veget_config_path=config_path,
                   tile=tilename,
                   shp=shapefile)
    veggie.run_veg_et()