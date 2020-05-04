
# ---------------------------------------------------------------------------
# SCRIPT NAME: rain_fraction.py
# Created on: May 2019
# VERSION: ArcGIS Pro - Python 3.6
# AUTHOR: Stefanie Kagone
# Project lead: Gabriel Senay
# Review: Naga Manohar Velpuri
# Usage:  Create rain fraction parameter used in the VegET model.
# -----------------------------------------------------------------------

import os
import arcpy
from arcpy.sa import *

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = 1

#workDir = r'C:\WaterSmart\Users\Stefanie\Projects\Veg_ET\MENA\MENAdata'
#arcpy.env.cellSize = os.path.join(workDir, 'NDVI_500m_med001_snapgrid_ndfill.tif')
#arcpy.env.snapRaster = os.path.join(workDir, 'NDVI_500m_med001_snapgrid_ndfill.tif')
# arcpy.env.extent = os.path.join(workDir, 'NDVI_500m_med001_snapgrid_ndfill.tif')

workDir = r'W:\Projects\Veg_ET\USA_data'
arcpy.env.cellSize = os.path.join(workDir, 'NDVI_2001001_snapgrid.tif')
arcpy.env.snapRaster = os.path.join(workDir, 'NDVI_2001001_snapgrid.tif')
arcpy.env.extent = os.path.join(workDir, 'NDVI_2001001_snapgrid.tif')

#inTavg = os.path.join(workDir, 'tavg_Cdaily')
inTavg = r'W:\Projects\Veg_ET\USA_data\temperature_gridmet\tavg_gridmet\Med_tavg_1984_2017_C'
arcpy.env.workspace = inTavg
TavgList = arcpy.ListRasters()
TavgList.sort()
print(TavgList)


for tavg in TavgList:
    print(tavg)
    ta = Raster(os.path.join(inTavg, tavg))
    #a = Con(ta <= 0, 0, Con(ta >= 6, 1, 0.167*(ta-6)))
    #a = Con(ta <= 6, 0, Con(ta >= 6, 1))
    a = Con(ta <= 6, 0, Con(ta >= 12, 1, 0.08333 * ta))
    #a = Con(ta <= 279.15, 0, Con(ta >= 285.15, 1, 0.08333 * (ta-273.15)))
    aname = 'rainfrac' + tavg[4:]
    print(aname)
    output_dir = os.path.join(workDir, 'rainfraction_gridmet')
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    a.save(output_dir + os.sep + aname)


print('done')
