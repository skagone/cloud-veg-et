import os
import rasterio
from rasterio.mask import mask
import fiona
import numpy as np

def cum_mm_to_m3(raster_dim_meters, arr):
    # convert mm to mm* m^2
    arr *= raster_dim_meters ** 2
    # mm to m conversion
    arr /= 1000

    # return m3
    return arr

# year_climatol_file = r'Z:\Data\Runoff\Global\runoff\cmp_ro.grd'
year_climatol_file = r''

#================================
areadict = {'Mississippi': 3228840000000, 'Ganges': 1662310000000,
                'Danube': 801840000000, 'Murray': 1049020000000, 'Nile': 3057390000000,
                'San_Francisco': 621621000000}
basin_list = ['Mississippi', 'Ganges', 'Danube', 'Murray', 'Nile', 'San_Francisco']
shapenames = ['Mississippi', 'Ganges_Brahmaputra', 'Danube', 'Murray_Darling', 'Nile', 'Sao_Francisco']

raster_dim_meters = 1000
# ================================
for basin, shapename in zip(basin_list, shapenames):
    basin_area = areadict[basin]
    output_root = r'Z:\Projects\VegET_Basins\{}'.format(basin)
    sample_shape = r'Z:\Projects\VegET_Basins\{}\shapefiles\{}.shp'.format(basin, shapename)
    # print(sample_shape)
    #================================

    basin_dict = {}

    # === Zonal Stats ===
    with fiona.open(sample_shape, 'r') as shp:
        features = [feature for feature in shp]
        # for i, f in enumerate(features):
        # print(f'feat dict {i} \n ', f, '\n')
        # shapes = [feature['geometry'] for feature in shp]
        with rasterio.open(year_climatol_file, 'r') as src:
            for f in features:
                src_shp = [f['geometry']]
                watershed_id = f['id']
                # watershed_name = f['properties']['NAME']
                # writer.writerow(["Year", "DOY", "parameter", "zon_mean_forID1", "zon_mean_forID2", "zon_mean_forID3", ...])
                outimage, out_transform = mask(src, src_shp, crop=True)
                # scrap too-large datasets
                outimage[outimage >= 5000] = np.nan
                outimage[outimage < 0] = np.nan

                ws_mean = np.nanmean(outimage)
                # calculating cubic meters of discharge for a basin
                m3_arr = cum_mm_to_m3(raster_dim_meters, outimage)
                basin_ro = np.nansum(m3_arr)

                basin_dict['area'] = basin_area
                basin_dict['mm_runoff'] = ws_mean
                basin_dict['m3_runoff'] = basin_ro
                basin_dict['name'] = basin


    print(basin_dict)