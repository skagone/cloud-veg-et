"""Script (with a few mods) taken from:
https://rasterio.readthedocs.io/en/latest/topics/reproject.html
"""
import os
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling

input_fol = r'Z:\Projects\VegET_ndviLSP\daily_ndvi'
output_fol = r'Z:\Projects\VegET_ndviLSP\daily_ndvi_wgs84'
dst_crs = 'EPSG:4326'

for f in os.listdir(input_fol):
    if f.endswith('.tif'):
        fpath = os.path.join(input_fol, f)
        outpath = os.path.join(output_fol, f)

        with rasterio.open(fpath) as src:
            transform, width, height = calculate_default_transform(
                src.crs, dst_crs, src.width, src.height, *src.bounds)
            kwargs = src.meta.copy()
            kwargs.update({
                'crs': dst_crs,
                'transform': transform,
                'width': width,
                'height': height
            })

            with rasterio.open(outpath, 'w', **kwargs) as dst:
                for i in range(1, src.count + 1):
                    reproject(
                        source=rasterio.band(src, i),
                        destination=rasterio.band(dst, i),
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=transform,
                        dst_crs=dst_crs,
                        resampling=Resampling.nearest)