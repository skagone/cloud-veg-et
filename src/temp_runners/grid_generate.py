from src.api_veget.dekade_try01 import GridMeister
'''
This is an example of how to use gridmeister to create shapefiles locally.
'''
# D:\Users\gparrish\PycharmProjects\cloud-veg-et\api_veget\dekade_try01\gridmeister.py

chip_output = r'Z:\Projects\VegET_ndviLSP\aa_30m_run\shapefiles'

exp_extent = (-119.813, 35.859,
              -119.051, 36.471)

# 0.002310233679679207525/250 = x/30
gm = GridMeister(tile_name='centralValley', raster_extent=exp_extent,
                 x_raster_res=0.00027722804, y_raster_res=0.00027722804)


lst = gm.chip_list()
print('resulting chip list \n', lst)
print(f'{len(lst)} many chips')
gm.create_chip_shp(ul_lat=None, ul_lon=None, out_location=chip_output, unit_chip=True)

if lst == None:
    gm.create_chip_shp(ul_lat=None, ul_lon=None, out_location=chip_output, unit_chip=True)
else:
    for l in lst:
        gm.create_chip_shp(l[0], l[-1], out_location=chip_output)
