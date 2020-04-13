#!/usr/bin/env python
# coding: utf-8

# ## Veg ET model in Python 3.2

# In[2]:


# ---------------------------------------------------------------------------
# SCRIPT NAME: VegET_batch_continous.py
# Created on: April 2014
# VERSION: ArcGIS 10.2
# AUTHOR: Stefanie Kagone
# Project lead: Gabriel Senay
# Review: Naga Manohar Velpuri
# Usage:  Create runoff and evapotranspiration parameters using the VegET model.
# -----------------------------------------------------------------------


# Module imports - setting up the properties

# In[1]:


import datetime
import os
import sys
import time
import traceback

# Import system modules
import arcpy
from arcpy.sa import *

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = 1

workDir = r'W:\Projects\Veg_ET'
arcpy.env.snapRaster = os.path.join(workDir, 'USA_data', 'NDVI_snapgrid_1km.tif')
arcpy.env.cellSize = os.path.join(workDir, 'USA_data', 'NDVI_snapgrid_1km.tif')


# Start of the code:

# In[ ]:


try:
    start_year = 1984
    o = 274

    con_year = 1984  # Oct 1
    m = 274

    outputDir = os.path.join(workDir, 'modelOutput_1km', 'USA', str(con_year))
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    time_stamp = time.strftime('%Y.%m.%d_%H.%M.%S', (time.localtime(time.time())))
    strFile = os.path.join(outputDir, "Log_%s.txt") % time_stamp
    log_file = open(strFile, "a")

    Dir_NDVI = os.path.join(workDir, 'USA_data', 'NDVI_daily_1_km', 'NDVI_filled')
    Dir_PET = os.path.join(r'W:\Data\ReferenceET\USA\Gridmet\Daily\ETo', 'Median_1984_2017')
    Dir_PPT = os.path.join(workDir, 'USA_data', 'precip_gridmet_tiffs')
    Dir_Tavg = os.path.join(workDir, 'USA_data', 'tavg_Cdaily')
    Dir_MR = os.path.join(workDir, 'USA_data', 'melt_rate')
    Dir_RF = os.path.join(workDir, 'USA_data', 'rainfraction')

    etaDir = os.path.join(outputDir, 'etasw')
    if not os.path.exists(etaDir):
        os.makedirs(etaDir)

    etaMDir = os.path.join(outputDir, 'etasw_monthly')
    if not os.path.exists(etaMDir):
        os.makedirs(etaMDir)

    etaYDir = os.path.join(outputDir, 'etasw_yearly')
    if not os.path.exists(etaYDir):
        os.makedirs(etaYDir)

    effDir = os.path.join(outputDir, 'eff')
    if not os.path.exists(effDir):
        os.makedirs(effDir)

    swiDir = os.path.join(outputDir, 'swi')
    if not os.path.exists(swiDir):
        os.makedirs(swiDir)

    swfDir = os.path.join(outputDir, 'swf')
    if not os.path.exists(swfDir):
        os.makedirs(swfDir)

    rfDir = os.path.join(outputDir, 'rf')
    if not os.path.exists(rfDir):
        os.makedirs(rfDir)

    srfDir = os.path.join(outputDir, 'srf')
    if not os.path.exists(srfDir):
        os.makedirs(srfDir)

    ddDir = os.path.join(outputDir, 'dd')
    if not os.path.exists(ddDir):
        os.makedirs(ddDir)

    intcepDir = os.path.join(outputDir, 'intcep')
    if not os.path.exists(intcepDir):
        os.makedirs(intcepDir)

    pptDir = os.path.join(outputDir, 'ppt')
    if not os.path.exists(pptDir):
        os.makedirs(pptDir)

    rainDir = os.path.join(outputDir, 'rain')
    if not os.path.exists(rainDir):
        os.makedirs(rainDir)

    sweDir = os.path.join(outputDir, 'swe')
    if not os.path.exists(sweDir):
        os.makedirs(sweDir)

    snwpkDir = os.path.join(outputDir, 'snwpk')
    if not os.path.exists(snwpkDir):
        os.makedirs(snwpkDir)

    log_file.write('VegET calculation ' + time_stamp + '\n')
    log_file.write('\n')
    log_file.write('\n')
    log_file.write('Calculated for continous years starting in --> ' + str(con_year) + '\n')
    log_file.write('Output Folder: ' + outputDir + '\n')
    log_file.write('\n')
    log_file.write('Input data: \n')
    log_file.write('NDVI folder: ' + Dir_NDVI + '\n')
    log_file.write('PET folder: ' + Dir_PET + '\n')
    log_file.write('PPT folder: ' + Dir_PPT + '\n')

    Igrid = Raster(os.path.join(workDir, 'USA_data', 'Soil_data', 'Intercept_us_filled.tif'))
    WHCgrid = Raster(os.path.join(workDir, 'USA_data', 'Soil_data', 'whc3_1mwgs_500m.tif'))
    SATgrid = Raster(os.path.join(workDir, 'USA_data', 'Soil_data', 'Soil_saturation_us.tif'))
    FCgrid = Raster(os.path.join(workDir, 'USA_data', 'Soil_data', 'Field_capacity_us.tif'))
    Iswegrid = Raster(os.path.join(workDir, 'USA_data', 'SWE_initial_US.tif'))
    Isnowpackgrid = Raster(os.path.join(workDir, 'USA_data', 'Snowpack_initial_US.tif'))

    log_file.write('Intcep grid: ' + str(Igrid) + '\n')
    log_file.write('WHC grid: ' + str(WHCgrid) + '\n')
    log_file.write('SAT grid: ' + str(SATgrid) + '\n')

    WHC = WHCgrid  # WHC grid is 50% of the original
    SAT = SATgrid  # precomputed Satgrid in mm
    FC = FCgrid
    varA = 1.25
    varB = 0.2
    log_file.write('Variable A: ' + str(varA) + '\n')
    log_file.write('Variable B: ' + str(varB) + '\n')

    log_file.write('\n')
    log_file.write('Start calculation\n')
    log_file.write('\n')

    i = 0  # counting +1
    # m = 274  # counting +1   #1 = Jan1
    #swflist = [r'W:\Projects\Veg_ET\modelOutput_1km\USA\2014\swf\swf14060'] # r'C:\WaterSmart\Users\Stefanie\Projects\Veg_ET\MENA\modelOutput\MENArun_rc035\1983\swf\swf83365']
    #swelist = [r'W:\Projects\Veg_ET\modelOutput_1km\USA\2014\swe\swe14060'] #r'C:\WaterSmart\Users\Stefanie\Projects\Veg_ET\MENA\modelOutput\MENArun_rc035\1983\swe\swe83365']
    #snowpacklist = [r'W:\Projects\Veg_ET\modelOutput_1km\USA\2014\snwpk\snwpk14060'] #r'C:\WaterSmart\Users\Stefanie\Projects\Veg_ET\MENA\modelOutput\MENArun_rc035\1983\snwpk
    # \snwpk83365']

    while con_year <= 2016:
        print('DOY: ' + str(con_year) + ' - ' + str(m))
        log_file.write('DOY: ' + str(con_year) + ' - ' + str(m) + '\n')
        jdate = str(m).zfill(3)
        jjdate = str(con_year)[-2:] + jdate

        # DOY to date
        datea = str(datetime.date(con_year, 1, 1) + datetime.timedelta(int(jdate) - 1))
        print('date: ' + str(datea))
        mon = int(datea.split('-')[1])
        day = int(datea.split('-')[2])

        ndvi = os.path.join(Dir_NDVI, str(con_year), str(con_year) + str(jdate) + '.1_km_16_days_NDVI.tif')
        pet = os.path.join(Dir_PET, 'medianETo' + str(jdate) + '.tif')
        ppt1 = os.path.join(Dir_PPT, str(con_year),'prec_' + str(con_year) + str(jdate) + '.tif')
        tavg = Raster(os.path.join(Dir_Tavg, 'tavg_' + str(jdate) + '.tif'))
        melt_rate = Raster(os.path.join(Dir_MR, 'meltrate_' + str(jdate) + '.tif.tif'))
        rain_frac = Raster(os.path.join(Dir_RF, 'rainfrac_' + str(jdate) + '.tif'))

        arcpy.env.workspace = outputDir
        print('NDVI = ' + ndvi)
        print('PET = ' + pet)
        print('PPT = ' + ppt1)
        print('Tavg = ' + str(tavg))
        print('Melt rate = ' + str(melt_rate))
        print('Rain fraction = ' + str(rain_frac))

        etaNew = os.path.join(etaDir, 'eta' + str(jjdate))  # + '.tif'

        log_file.write('Processing day: ' + datea + '\n')
        log_file.write('Input data for ' + str(jjdate) + ': '
                       + os.path.basename(ndvi) + " ,"
                       + os.path.basename(pet) + " ,"
                       + os.path.basename(ppt1) + " ,"
                       + os.path.basename(str(tavg)))
        log_file.write('Output: ' + etaNew + '\n')

        if jjdate == str(start_year)[-2:] + str(o):   # '83274'
            print('Folder created already')
        else:
            if jjdate == str(con_year)[-2:] + '001':
                outputDir = os.path.join(workDir, 'modelOutput_1km', 'USA', str(con_year))
                if not os.path.exists(outputDir):
                    os.makedirs(outputDir)
                etaDir = os.path.join(outputDir, 'etasw')
                if not os.path.exists(etaDir):
                    os.makedirs(etaDir)
                etaMDir = os.path.join(outputDir, 'etasw_monthly')
                if not os.path.exists(etaMDir):
                    os.makedirs(etaMDir)
                etaYDir = os.path.join(outputDir, 'etasw_yearly')
                if not os.path.exists(etaYDir):
                    os.makedirs(etaYDir)
                effDir = os.path.join(outputDir, 'eff')
                if not os.path.exists(effDir):
                    os.makedirs(effDir)
                swiDir = os.path.join(outputDir, 'swi')
                if not os.path.exists(swiDir):
                    os.makedirs(swiDir)
                swfDir = os.path.join(outputDir, 'swf')
                if not os.path.exists(swfDir):
                    os.makedirs(swfDir)
                rfDir = os.path.join(outputDir, 'rf')
                if not os.path.exists(rfDir):
                    os.makedirs(rfDir)
                srfDir = os.path.join(outputDir, 'srf')
                if not os.path.exists(srfDir):
                    os.makedirs(srfDir)
                ddDir = os.path.join(outputDir, 'dd')
                if not os.path.exists(ddDir):
                    os.makedirs(ddDir)
                intcepDir = os.path.join(outputDir, 'intcep')
                if not os.path.exists(intcepDir):
                    os.makedirs(intcepDir)
                pptDir = os.path.join(outputDir, 'ppt')
                if not os.path.exists(pptDir):
                    os.makedirs(pptDir)
                rainDir = os.path.join(outputDir, 'rain')
                if not os.path.exists(rainDir):
                    os.makedirs(rainDir)
                sweDir = os.path.join(outputDir, 'swe')
                if not os.path.exists(sweDir):
                    os.makedirs(sweDir)
                snwpkDir = os.path.join(outputDir, 'snwpk')
                if not os.path.exists(snwpkDir):
                    os.makedirs(snwpkDir)

        print(jjdate + ' is processing in output folder ' + str(outputDir))

        ppt3 = Raster(ppt1)
        ppt3.save(os.path.join(pptDir, 'ppt' + str(jjdate)))  # + '.tif')
        ppt = os.path.join(pptDir, 'ppt' + str(jjdate))  # + '.tif'

        arcpy.env.snapRaster = os.path.join(workDir, 'USA_data', 'NDVI_snapgrid_1km.tif')
        print('Effective PPT = ' + str(ppt) + '* (1 - (Igrid/100)')
        effppt = ppt * (1 - (Igrid / 100.0))
        arcpy.CopyRaster_management(effppt, os.path.join(effDir, 'effppt' + str(jjdate)), 
                                    "#", "#", "#", "NONE","NONE", "#")
        print('Intercepted PPT = ' + str(ppt) + '* (Igrid/100)')
        intcep = ppt * Float(Igrid / 100.0)
        arcpy.env.compression = 'LZW'
        arcpy.CopyRaster_management(intcep, os.path.join(intcepDir, 'intcep' + str(jjdate)), 
                                    "#", "#", "#", "NONE","NONE", "#")
        print('Created eff Rainfall and Interception')
        log_file.write('Created eff Rainfall and Interception' + '\n')
        log_file.write('\n')

        # Snowpack
        print('Snow component of Precipitation....' + '\n')
        print('Calculating what portion is rain or snow based on average temperature')

        snow_melt_fac = arcpy.sa.Con(tavg <= 6, 0, melt_rate)
        snow_melt_fac.save(os.path.join(snwpkDir, 'swmeltfac' + str(jjdate) + '.tif'))

        if jjdate == str(con_year)[-2:] + str(o):   # '83274'
            print('Initial SWE: ' + str(Iswegrid))
            print('Initial Snowpack: ' + str(Isnowpackgrid))
            print('Initial SWi = (WHC * 0.5) + effppt')

            swflist = []
            swelist = []
            snowpacklist = []

            RAIN = rain_frac * effppt
            arcpy.env.compression = 'LZW'
            arcpy.CopyRaster_management(RAIN, os.path.join(rainDir, 'rain' + str(jjdate)),
                                        "#", "#", "#", "NONE", "NONE", "#", "NONE", "NONE")

            SWE = Iswegrid  # inital snowpack grid with 0 values
            arcpy.env.compression = 'LZW'
            arcpy.CopyRaster_management(SWE, os.path.join(sweDir, 'swe' + str(jjdate)),
                                        "#", "#", "#", "NONE", "NONE", "#", "NONE", "NONE")
            swelist.append(os.path.join(sweDir, 'swe' + str(jjdate)))

            snow_melt = SWE

            SNOW_pack = Isnowpackgrid  # inital snowpack grid with 0 values
            arcpy.env.compression = 'LZW'
            arcpy.CopyRaster_management(SNOW_pack, os.path.join(snwpkDir, 'snwpk' + str(jjdate)),
                                        "#", "#", "#", "NONE", "NONE", "#", "NONE", "NONE")
            snowpacklist.append(os.path.join(snwpkDir, 'snwpk' + str(jjdate)))

            SWi = (WHC * 0.5) + effppt + snow_melt
            SWi.save(os.path.join(swiDir, 'swi' + str(jjdate)))  # + '.tif')

        else:
            if jjdate == str(con_year)[-2:] + '001':

                # add 1 to the year since its at the beginning of the process and otherwise looks for 84366
                if con_year == 1985 or con_year == 1989 or con_year == 1993 or con_year == 1997                         or con_year == 2001 or con_year == 2005 or con_year == 2009 or con_year == 2013:
                    pre_year = int(con_year) - 1
                    pre_name = str(pre_year)[-2:] + '366'
                    swflist = [os.path.join(workDir, 'modelOutput_1km', 'USA', str(pre_year), 'swf', 'swf' + pre_name)]
                    swelist = [os.path.join(workDir, 'modelOutput_1km', 'USA', str(pre_year), 'swe', 'swe' + pre_name)]
                    snowpacklist = [os.path.join(workDir, 'modelOutput_1km', 'USA', str(pre_year), 'snwpk', 'snwpk' + pre_name)]
                else:
                    pre_year = int(con_year) - 1
                    pre_name = str(pre_year)[-2:] +'365'
                    swflist = [os.path.join(workDir, 'modelOutput_1km', 'USA', str(pre_year), 'swf','swf'+pre_name)]
                    swelist = [os.path.join(workDir, 'modelOutput_1km', 'USA', str(pre_year), 'swe', 'swe'+pre_name)]
                    snowpacklist = [os.path.join(workDir, 'modelOutput_1km', 'USA', str(pre_year), 'snwpk', 'snwpk'+pre_name)]

                print('Previous SWE: ' + str(swelist))
                print('Previous SNWpk: ' + str(snowpacklist))
                print('Previous SWi: ' + str(swflist))

            else:
                print('Previous SWE: ' + swelist[-1])
                print('Previous SNWpk: ' + snowpacklist[-1])
                print('Previous SWi: ' + swflist[-1])

            RAIN = rain_frac * effppt
            arcpy.env.compression = 'LZW'
            arcpy.CopyRaster_management(RAIN, os.path.join(rainDir, 'rain' + str(jjdate)), 
                                        "#", "#", "#", "NONE", "NONE", "#", "NONE", "NONE")

            SWE = (1 - rain_frac) * effppt
            SWE.save(os.path.join(sweDir, 'swe' + str(jjdate)))  # + '.tif')
            swelist.append(os.path.join(sweDir, 'swe' + str(jjdate)))

            snow_melt = arcpy.sa.Con(melt_rate <= (SWE + snowpacklist[-1]), melt_rate, (SWE + snowpacklist[-1]))

            SNWpk1 = snowpacklist[-1] + SWE - snow_melt
            SNWpk = arcpy.sa.Con(SNWpk1 < 0, 0, SNWpk1)
            SNWpk.save(os.path.join(snwpkDir, 'snwpk' + str(jjdate)))  # + '.tif')
            snowpacklist.append(os.path.join(snwpkDir, 'snwpk' + str(jjdate)))

            SWi = swflist[-1] + RAIN + snow_melt
            SWi.save(os.path.join(swiDir, 'swi' + str(jjdate)))  # + '.tif')

        print('Created Rainfall fraction, Snow pack, SWE and Soil water')
        log_file.write('Created Rainfall fraction, Snow pack, SWE and Soil water')
        log_file.write('\n')

        # Deep drainage and surface runoff = runoff -> now 2 components
        print('Runoff is sum of deep drainage and surface runoff....' + '\n')
        # log_file('Runoff is sum of deep drainage and surface runoff....' + '\n')
        # WHC = FC - WP
        # SAT_FC = SAT - WHC

        # Deep Drainage
        # DDrain occurs if SWi > WHC, amount of DDrain is SAT <> WHC with a maximum DDrain of SAT - WHC

        sat_fc = SAT - FC
        # sat_whc = SAT - WHC
        rf1 = SWi - WHC

        rf = arcpy.sa.Con(rf1 < 0, 0, rf1)
        dc_coeff = 0.65  # drainage coefficient
        rf_coeff = 1 - dc_coeff  # runoff coefficient
        SRf = arcpy.sa.Con(rf <= sat_fc, rf * rf_coeff, (rf - sat_fc) + rf_coeff * sat_fc)
        DDrain = rf - SRf
        arcpy.env.compression = 'LZW'
        arcpy.CopyRaster_management(DDrain, os.path.join(ddDir, 'dd' + str(jjdate)),
                                    "#", "#", "#", "NONE", "NONE", "#", "NONE", "NONE")
        arcpy.CopyRaster_management(SRf, os.path.join(srfDir, 'srf' + str(jjdate)),
                                    "#", "#", "#", "NONE", "NONE", "#", "NONE", "NONE")

        print('Created Surface Runoff and Deep Drainage components')
        log_file.write('Created Surface Runoff and Deep Drainage components' + '\n')
        log_file.write('\n')

        # Convert from -2000 to 10000 range to 0-200 range (scaled ndvi)
        ndviF1 = Raster(ndvi)
        ndviF = Float(ndviF1)  # / 10000.0
        # ndviF = (ndviF2 * 100)  + 100

        # ETaSW
        log_file.write('ETasw = ' + str(varA) + ' + NDVI + ' + str(varA) + '\n')
        etasw1A = (varA * ndviF + varB) * Raster(pet)
        etasw1B = (varA * ndviF) * Raster(pet)  # * 8)
        etasw1 = arcpy.sa.Con(ndviF > 0.4, etasw1A, etasw1B)
        etasw2 = (SWi / (0.5 * WHC)) * etasw1
        etasw3 = arcpy.sa.Con(SWi > (0.5 * WHC), etasw1, etasw2)
        etasw4 = arcpy.sa.Con(etasw3 > SWi, SWi, etasw3)
        etasw5 = arcpy.sa.Con(etasw4 > WHC, WHC, etasw4)
        etasw = etasw5  # Int(etasw5 + 0.5)
        arcpy.env.compression = 'LZW'
        arcpy.CopyRaster_management(etasw, etaNew, "#", "#", "#", "NONE", "NONE", "#", "NONE", "NONE")
        print('Created ETaSW')
        log_file.write('Created ETaSW' + '\n')
        log_file.write('\n')

        SWf1 = SWi - etasw
        arcpy.env.snapRaster = etasw
        bigswi = (WHC - etasw)
        SWf = arcpy.sa.Con(SWi > WHC, bigswi, arcpy.sa.Con(SWf1 < 0.0, 0.0, SWf1))
        arcpy.env.compression = 'LZW'
        arcpy.CopyRaster_management(SWf, os.path.join(swfDir, 'swf' + str(jjdate)), "#", "#", "#",
                                    "NONE", "NONE", "#", "NONE", "NONE")
        swflist.append(os.path.join(swfDir, 'swf' + str(jjdate)))
        print('Created Soil Water')
        log_file.write('Created Soil Water' + '\n')
        print(swflist)
        if con_year == 1984 or con_year == 1988 or con_year == 1992 or con_year == 1996 
        or con_year == 2000 or con_year == 2004 or con_year == 2008 or con_year == 2012:
            if m == 366:
                m = 0
                con_year += 1
        else:
            if m == 365:
                m = 0
                con_year += 1

        i += 1
        m += 1
        print('------------------------')
        log_file.write('--------------------------------' + '\n')
        del SWi
        del SWf
        del SWf1

    log_file.write("all done...")
    log_file.close()
    os.startfile(strFile)

except Exception as exc:
    tb = sys.exc_info()[2]
    tbinfo = traceback.format_tb(tb)[0]
    pymsg = ('PYTHON ERRORS:\nTraceback Info:\n' + tbinfo +
             '\nError Info:\n       ' + str(exc))
    arcpy.AddError(pymsg)
    print(pymsg)

