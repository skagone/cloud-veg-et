# ---------------------------------------------------------------------------
# SCRIPT NAME: melt_rate.py
# Created on: May 2019
# VERSION: ArcGIS Pro - Python 3.6
# AUTHOR: Stefanie Kagone
# Project lead: Gabriel Senay
# Review: Naga Manohar Velpuri
# Usage:  Create melt rate parameter used in the VegET model.
# -----------------------------------------------------------------------

import os
import arcpy
from arcpy.sa import *

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = 1

#workDir = r'C:\WaterSmart\Users\Stefanie\Projects\Veg_ET\MENA\MENAdata'
#arcpy.env.cellSize = os.path.join(workDir, 'NDVI_500m_med001_snapgrid_ndfill.tif')
#arcpy.env.snapRaster = os.path.join(workDir, 'NDVI_500m_med001_snapgrid_ndfill.tif')

workDir = r'W:\Projects\Veg_ET\USA_data'
arcpy.env.cellSize = os.path.join(workDir, 'NDVI_2001001_snapgrid.tif')
arcpy.env.snapRaster = os.path.join(workDir, 'NDVI_2001001_snapgrid.tif')

melt_factor = 0.06
#inTmin = r'C:\WaterSmart\Data\Temperature\Global\TA_worldmet\tmin_Cdaily'
inTmin = r'W:\Projects\Veg_ET\USA_data\temperature_gridmet\tmin_gridmet\Med_tmin_1984_2017_C'
arcpy.env.workspace = inTmin
TminList = arcpy.ListRasters()
TminList.sort()

#inTmax = r'C:\WaterSmart\Data\Temperature\Global\TA_worldmet\tmax_Cdaily'
inTmax = r'W:\Projects\Veg_ET\USA_data\temperature_gridmet\tmax_gridmet\Med_tmax_1984_2017_C'
arcpy.env.workspace = inTmax
TmaxList = arcpy.ListRasters()
TmaxList.sort()

for tmin, tmax in zip(TminList, TmaxList):
    print(tmin)
    print(tmax)
    #arcpy.env.extent = r'C:\WaterSmart\Users\Stefanie\Projects\Veg_ET\MENA\MENAdata' \
                       #r'\NDVI_500m_med001_snapgrid_ndfill.tif '
    #arcpy.env.snapRaster = r'C:\WaterSmart\Users\Stefanie\Projects\Veg_ET\MENA\MENAdata' \
                           #r'\NDVI_500m_med001_snapgrid_ndfill.tif'

    arcpy.env.extent = r'W:\Projects\Veg_ET\USA_data\NDVI_2001001_snapgrid.tif'
    arcpy.env.snapRaster = r'W:\Projects\Veg_ET\USA_data\NDVI_2001001_snapgrid.tif'
    tx = Raster(os.path.join(inTmax, tmax))
    tn = Raster(os.path.join(inTmin, tmin))
    a = 0.06 * ((tx * tx) - (tx * tn))
    aname = 'meltrate' + tmin[4:]
    print(aname)
    #a.save(os.path.join(workDir, 'melt_rate', aname))
    output_dir = os.path.join(workDir, 'melt_rate_gridmet')
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    a.save(output_dir +os.sep + aname)

print('done')
