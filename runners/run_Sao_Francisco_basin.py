from veget.vegetLib.vegetLib.veget import VegET
import os

"""
This runner was made on 9/21/2020 from Oct 2002-Dec 2019. This run is being done locally but more runs
 (for the Danube, Ganges-Bramaputra etc) will be run in the cloud if possible...
"""

config_root = r'C:\WaterSmart\Projects\VegET_Basins\Sao_Francisco'
shape_root = r'C:\WaterSmart\Projects\VegET_Basins\Sao_Francisco'

tile_lst = ['S1.shp', 'S2.shp', 'S3.shp', 'S4.shp', 'S5.shp']
# # testing
# tile_lst = ['N4.shp']
config_list = ['S1', 'S2', 'S3', 'S4', 'S5']
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