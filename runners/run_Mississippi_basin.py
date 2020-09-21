from veget.vegetLib.vegetLib.veget import VegET
import os

"""
This runner was made on 9/21/2020 from Oct 2002-Dec 2019. This run is being done locally but more runs
 (for the Danube, Ganges-Bramaputra etc) will be run in the cloud if possible...
"""

config_root = r'C:\WaterSmart\Projects\VegET_Basins\Mississippi'
shape_root = r'C:\WaterSmart\Projects\VegET_Basins\Mississippi'

tile_lst = ['M1.shp', 'M2.shp', 'M3.shp', 'M4.shp', 'M5.shp', 'M6.shp',
            'M7.shp', 'M8.shp', 'M9.shp', 'M10.shp', 'M11.shp', 'M12.shp']
# testing
# tile_lst = ['M4.shp']
config_list = ['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12']
# testing
# config_list = ['M4']

for conf, tile in zip(config_list, tile_lst):
    tilename = tile.split('.')[0]
    shapefile = os.path.join(shape_root, tile)
    config_path = os.path.join(config_root, conf)

    veggie = VegET(veget_config_path=config_path,
                   tile=tilename,
                   shp=shapefile)
    veggie.run_veg_et()