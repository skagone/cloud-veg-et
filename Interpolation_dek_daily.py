# ---------------------------------------------------------------------------
# SCRIPT NAME: Interpolation_dek_daily.py
# Created on: May 2019
# VERSION: ArcGIS Pro - Python 3.6
# AUTHOR: Stefanie Kagone
# Project lead: Gabriel Senay
# Review:
# Usage:  Linear Interpolation for temperature data.
# -----------------------------------------------------------------------

# Import system modules
import sys, os
import arcpy
from arcpy.sa import *

print(' imessed up  mikes file')

# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")
arcpy.env.rasterStatistics = "STATISTICS"

#workDir = r'C:\WaterSmart\Users\Stefanie\Projects\Veg_ET\MENA\MENAdata'
workDir = r'W:\Projects\Veg_ET\Global_data\tmin_Cdaily'
#arcpy.env.extent = os.path.join(workDir, 'NDVI_500m_med001_snapgrid_ndfill.tif')
#arcpy.env.snapRaster = os.path.join(workDir, 'NDVI_500m_med001_snapgrid_ndfill.tif')

#year = sys.argv[1]
type = 'tmin' # sys.argv[6]
#arcpy.env.workspace = os.path.join(r'Y:\fewspsnfs2\WaterSmart\Data\Temperature\TA_worldmet', str(type))
arcpy.env.workspace = os.path.join(r'W:\Data\Temperature\Global\TA_worldmet', str(type))
rasterList = arcpy.ListRasters('*.tif')
rasterList.sort()
#print(rasterList)

#outputData1 = os.path.join('C:\WaterSmart\Data\Temperature\MENA\TA_worldmet', str(type)+'_Cdaily')
outputData1 = os.path.join('W:\Projects\Veg_ET\Global_data', str(type)+'_Cdaily')
if not os.path.exists(outputData1):
    os.makedirs(outputData1)

#a = int(sys.argv[1]) # 0
#b = int(sys.argv[2]) # 1

#refDate1 = sys.argv[3] # 1   #011
#refDate2 = sys.argv[4] # 11  #012
#cday = sys.argv[5]
#timeDiff = int(refDate2) - int(refDate1)


refDate1 = sys.argv[1] #011
file_in1 = arcpy.env.workspace +os.sep+ 'tmin'+ refDate1 +'.tif'
refDate2 = sys.argv[2] #012
file_in2 = arcpy.env.workspace +os.sep+ 'tmin'+ refDate2 +'.tif'

cday = sys.argv[3]
#endday = sys.argv[4]
#timeDiff = sys.argv[5]
timeDiff = int(sys.argv[5])


print('Interpolating from day '+ str(refDate1) +' to '+ str(refDate2))
print('Days between the 2 rasters is: '+ str(timeDiff))
#for ras in rasterList:
#    print ras

#ra1 = Raster(rasterList[a]) - 273.15
#ra2 = Raster(rasterList[b]) - 273.15

#ra1 = Raster(file_in1) - 273.15
#ra2 = Raster(file_in2) - 273.15

#slope = (ra2-ra1)/timeDiff

# print('Raster 1 = ' + str(rasterList[a]) + ' and Raster 2 = '+ str(rasterList[b]) +' for slope calculation in Celcius')

for day in range(0,timeDiff):
    print(day)
    slope = (Raster(file_in1) - Raster(file_in2)) / timeDiff
    linInterp = Raster(file_in1) + (slope * day)
    final = Float(linInterp) - 273.15
    c = str(int(cday) + int(day))
    #print(str(type)+ c.zfill(3) + ' = ' +str(rasterList[a]) + ' + (slope * '+ str(day) +')')
    final.save(outputData1+os.sep+ str(type)+c.zfill(3)+'.tif')

    print('-----------next day')
print('-----------dekad done - next one')

