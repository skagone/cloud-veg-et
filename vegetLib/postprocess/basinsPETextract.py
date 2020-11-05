import os
"""
This script is going to extract the annual PET for a given river basin and give average 
mm of PET per year per basin and output the csv to the Netapp
"""

path = r'\\IGSKMNCNFS016\watersmartfs1\Data\ReferenceET\Global\GDAS_PET\Median'

#
# #================================
# areadict = {'Mississippi': 3228840000000, 'Ganges': 1662310000000,
#                 'Danube': 801840000000, 'Murray': 1049020000000, 'Nile': 3057390000000,
#                 'San_Francisco': 621621000000}
# basin_list = ['Mississippi', 'Ganges', 'Danube', 'Murray', 'Nile', 'San_Francisco']
# shapenames = ['Mississippi', 'Ganges_Brahmaputra', 'Danube', 'Murray_Darling', 'Nile', 'Sao_Francisco']
#
# raster_dim_meters = 1000
# # ================================
# for basin, shapename in zip(basin_list, shapenames):
#     basin_area = areadict[basin]
#     output_root = r'Z:\Projects\VegET_Basins\{}'.format(basin)
#     sample_shape = r'Z:\Projects\VegET_Basins\{}\shapefiles\{}.shp'.format(basin, shapename)
#     # print(sample_shape)
#     #================================