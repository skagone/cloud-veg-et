import os
import fiona
from rasterio.mask import mask
import rasterio
import numpy as np
# import pandas as pd

"""
USGS gridded discharge products: helps us know how well VegET is doing
in Mississippi (the river, not the state).
"""

sample_shape = r'Z:\Projects\VegET_Basins\Mississippi\shapefiles\Mississippi.shp'
usgs_root = r'Z:\Data\Runoff\US_USGSRunoff\yearlyRunoff'
usgs_filenames = ['runoff2001.tif', 'runoff2002.tif', 'runoff2003.tif',
                  'runoff2004.tif', 'runoff2005.tif', 'runoff2006.tif',
                  'runoff2007.tif', 'runoff2008.tif', 'runoff2009.tif',
                  'runoff2010.tif', 'runoff2011.tif', 'runoff2012.tif',
                  'runoff2013.tif', 'runoff2014.tif']
usgs_paths = [os.path.join(usgs_root, i) for i in usgs_filenames]

var_dict = {}
cols = []
# === Zonal Stats ===
with fiona.open(sample_shape, 'r') as shp:
    features = [feature for feature in shp]

    for path, fname in zip(usgs_paths, usgs_filenames):
        rasname = fname.split('.')[0]
        cols.append(rasname)
        with rasterio.open(path) as src:
            for f in features:
                src_shp = [f['geometry']]
                watershed_id = f['id']

                outimage, out_transform = mask(src, src_shp, crop=True)
                # scrap too-large datasets
                outimage[outimage >= 5000] = np.nan
                outimage[outimage < 0] = np.nan

                ws_mean = np.nanmean(outimage)

        var_dict[rasname] = ws_mean


print(var_dict)