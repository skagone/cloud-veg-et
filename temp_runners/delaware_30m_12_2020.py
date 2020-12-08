from vegetLib.vegetLib.veget import VegET
import os

"""
This runner was made on 12/7/2020. This run is done as part of testing VegET model
 sensitivity to NDVI spatial resolution.
"""

config_root = r'Z:\Projects\VegET_Basins\Ganges'
shape_root = r'Z:\Projects\VegET_Basins\Ganges\shapefiles'

# tile_lst = ['G1.shp', 'G2.shp', 'G3.shp', 'G4.shp', 'G5.shp', 'G6.shp']
# testing
tile_lst = ['D4.shp']
# config_list = ['G1', 'G2', 'G3', 'G4', 'G5', 'G6']
# testing
config_list = ['D4']

for conf, tile in zip(config_list, tile_lst):
    tilename = tile.split('.')[0]
    shapefile = os.path.join(shape_root, tile)
    config_path = os.path.join(config_root, conf)

    veggie = VegET(veget_config_path=config_path,
                   tile=tilename,
                   shp=shapefile)
    veggie.run_veg_et()