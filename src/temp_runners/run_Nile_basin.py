from src.vegetLib import VegET
import os

"""
This runner was made on 9/21/2020 from Oct 2002-Dec 2019. This run is being done locally but more runs
 (for the Danube, Ganges-Bramaputra etc) will be run in the cloud if possible...
"""

config_root = r'C:\WaterSmart\Projects\VegET_Basins\Nile'
shape_root = r'C:\WaterSmart\Projects\VegET_Basins\Nile\shapefiles'

#tile_lst = ['N1.shp', 'N2.shp', 'N3.shp', 'N4.shp', 'N5.shp', 'N6.shp', 'N7.shp', 'N8.shp', 'N9.shp', 'N10.shp']
# # testing
tile_lst = ['N4.shp']
#config_list = ['N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8', 'N9', 'N10']
# # testing
config_list = ['N4']

for conf, tile in zip(config_list, tile_lst):
    tilename = tile.split('.')[0]
    shapefile = os.path.join(shape_root, tile)
    config_path = os.path.join(config_root, conf)

    veggie = VegET(veget_config_path=config_path,
                   tile=tilename,
                   shp=shapefile)
    print(veggie)
    veggie.run_veg_et()