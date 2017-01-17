#!/bin/env python
# Created on May 21, 2015
#  by Keith Cherkauer
#
# This script will convert the files extracted from the HWSD and reduced to the 
# majority soil unit per VIC model grid cell using the script SimplifySoilDatabase.py
# into a VIC model Soil Parameter file.  Missing soil paramteres will be computed 
# using the functions from Saxton and Rawls (2006) and compiled into the library
# SoilWaterEquations.py stored in /depot/phig/apps/source/python/lib.  These calculations
# assume the simplier representation of K_theta (the decrease of water conductivity
# as soil water becomes less than saturation as described by Campbell (1976).  This method
# does not make use of residual water, so that value is set to zero for soil parameter
# files generated using this program.
#
# Organic matter should not exceed 8% by weight as solid with such high contents were not 
# used in developing the equations.
# 
# Concerns:
#
# I do not know what to do about estimating the bulk density of organic soil.  Will 
# have to look more closely at the VIC model source code.  KAC
#
# Modifications:
#

import sys, os
import pandas as pd
import numpy as np
import read_arcinfo_files as rwarc
import SoilWaterEquations as swe
import json

#sys.path.append("./Scripts")

from VicPreprocessingLibrary import *

#####
# Setup soil column defaults
# - provide soil layer thicknesses in mm for each defined soil layer (works for HWSD where layers are consistently defined)
# - provide naming format so that each layer can be properly identified
# - provide naming format so that files for all variablecan be found, or created as necessary 
#####

if len(sys.argv) <= 1 or len(sys.argv) > 2:
    # no control file given, exit gracefully
    sys.stderr.write( "\nUsage: %s <filename>\n\n\tThis script was written to calculate soil parameters from fundemental soil properties \n\t(%% sand, bulk density, etc) extracted for each VIC simulation cell from the harmonized \n\tworld soils database (HWSD).\n\n\t<filename> is the name of the control file that allows this program to \n\t\twork with your data.  Use the script CreateCalculatePerSoilDBLayerControlFile.py \n\t\tto build a control file.  \n\n\tExample files should also be in the same directory as the executable script.\n\n" % sys.argv[0] )
    sys.exit()

else:
    # control file provided, parse and get started
    fin = open( sys.argv[1], 'r' )
    ctrlDict = json.load( fin )
    fin.close()
    print ctrlDict
    #sys.exit()

    # adjust a few variables because the default type is problematic
    #ctrlDict['inSoilID'] = ctrlDict['inSoilID'].encode('ascii','ignore')


#numSL = 2 # number of soil layers in the soil database
#SoilLayerThick = [ ( 'SOIL_THICK', 0.300, 'm' ), ( 'SOIL_THICK', 0.700, 'm' ) ] # provide spatial file or default value (mm)
#SoilLayerNameFmt = [ "T_%s", "S_%s" ]
#SoilGridNameFmt = 'SoilsData/VICsetup/SoilProperty_%s.asc'
#outSoilDataTable = 'SoilsData/VICsetup/AllSoilDataTable.csv'
#outSoilAttribTable = 'SoilsData/VICsetup/AllSoilAttribTable.csv'
#RASTER = True # Set to True if input files are ArcInfo ASCII rasters, or False for XYZ tables

    #####
    # Layers needed for calculations
    # - Provide name of each layer as needed to complete ctrlDict['SoilLayerNameFmt'] defined above
    #####
#CellNumFile = 'ldascellnum_simres.txt'
#SimCellMask = 'vic_runcell.txt' # set to 'False' to activate all grid cells
#CellSlopeFile = ( 'SLOPE', 1, 'P' ) # need to calculate in ArcGIS from projected DEM, 'P' is percent, 'D' is degrees

#SandFract = ( 'SAND', np.NaN, 'P' ) # sand fraction of soil, no default, 'P' in percent, 'F' as a fraction
#ClayFract = ( 'CLAY', np.NaN, 'P' ) # clay fraction of soil, no default, 'P' in percent, 'F' as a fraction
#OrganicFract = ( 'OC', np.NaN, 'P' ) # organic matter fraction of soil, no default, 'P' in percent, 'F' as a fraction
#GravelFract = ( 'GRAVEL', np.NaN, 'P' ) # gravel fraction of soil (optional, set to 'None' if not provided), no default, 'P' in percent, 'F' as a fraction
#BulkDensity = ( 'BULK_DENSITY', 1, 'kg/dm3' ) # or should this be the REF_BULK_DENSITY column?
#ImpermLayer = ( 'IL', 1, 'None' ) # not sure of the usefulness of this yet

    #####
    # Define other files that contain information for the soil parameter file that will
    #   be read in and added to output data file, even though they are not layer specific.
    #####
#    OtherFiles = [ { 'colName':'elev', 'name':'dem_simres.txt' }, # VIC cell elevations
#                   { 'colName':'avg_T', 'name':'dummy_grid.txt' }, # Average annual air temp
#                   { 'colName':'off_gmt', 'name':'dummy_grid.txt' }, # time offset from GMT
#                   { 'colName':'annual_prec', 'name':'dummy_grid.txt' }, # average annual precipitation
#                   { 'colName':'July_Tavg', 'name':'dummy_grid.txt' } ] # average annual July air temperature

    #####
    # Get Grid Cell Number
    # - there MUST be a valid grid cell number file provided as it provides information on the 
    #   number of grid cells, and other pertanent header information for creating ARCINFO ASCII
    #   grid files as output.
    #####
    INTflag = True # grid cell number file should be integer
    print ctrlDict['CellNumFile']
    SoilParams, GridInfo, Ncells = GetVariable( ctrlDict['CellNumFile'], np.NaN, 'gridcel', -99, INTflag, ctrlDict['RASTER'] )

    INTflag = False # other data types should be treated as floats by default

    #####
    # Create soil layer thickness files, if not already present
    #####
    SaveNames = [0]*ctrlDict['numSL']
    for idx in range(ctrlDict['numSL']): SaveNames[idx] = {} 
    SoilParams['MaxDepth'] = np.zeros( Ncells )
    for SL in range(ctrlDict['numSL']):
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % ( ctrlDict['SoilLayerThick'][SL][0] )
        print ctrlDict['SoilGridNameFmt'] % ColName, ctrlDict['SoilLayerThick'][SL][1], ColName, Ncells
        tmpDF, tmpInfo, tmpNcells = GetVariable( ctrlDict['SoilGridNameFmt'] % ColName, ctrlDict['SoilLayerThick'][SL][1], ColName, Ncells, INTflag, ctrlDict['RASTER'] )
        # ***** Here would be an opportunity to check the IL for shallow grid cells, but need to decided earlier if such shallow soils should be treated as the majority if deeper soils exist in smaller parts fo the grid cell - KAC ***** 
        SoilParams[ColName] = tmpDF[ColName]
        if tmpNcells == Ncells: WriteDFtoASCIIgrid( SoilParams, GridInfo, ColName, ctrlDict['SoilGridNameFmt'] % ColName )
        SaveNames[SL]['SoilThick'] = ColName
        # calculate maximum depth of soil
        SoilParams['MaxDepth'] = SoilParams['MaxDepth'] + SoilParams[ColName]

    #####
    # Read spatial information for all required soil layer properties 
    #####
    for SL in range(ctrlDict['numSL']):

        # get sand fraction
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % ctrlDict['SandFract'][0]
        tmpDF, tmpInfo, tmpNcells = GetVariable( ctrlDict['SoilGridNameFmt'] % ColName, ctrlDict['SandFract'][1], ColName, Ncells, INTflag, ctrlDict['RASTER'] )
        SoilParams[ColName] = tmpDF[ColName]
        if tmpNcells == Ncells: WriteDFtoASCIIgrid( SoilParams, GridInfo, ColName, ctrlDict['SoilGridNameFmt'] % ColName )
        if ctrlDict['SandFract'][2] == 'P': SoilParams[ColName] = SoilParams[ColName] / 100. # convert percent to fraction
        SaveNames[SL][ctrlDict['SandFract'][0]] = ColName

        # get clay fraction
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % ctrlDict['ClayFract'][0]
        tmpDF, tmpInfo, tmpNcells = GetVariable( ctrlDict['SoilGridNameFmt'] % ColName, ctrlDict['ClayFract'][1], ColName, Ncells, INTflag, ctrlDict['RASTER'] )
        SoilParams[ColName] = tmpDF[ColName]
        if tmpNcells == Ncells: WriteDFtoASCIIgrid( SoilParams, GridInfo, ColName, ctrlDict['SoilGridNameFmt'] % ColName )
        #if tmpNcells == Ncells: WriteDFtoASCIIgrid( SoilParams, GridInfo, ctrlDict['SoilGridNameFmt'] % ColName )
        if ctrlDict['ClayFract'][2] == 'P': SoilParams[ColName] = SoilParams[ColName] / 100. # convert percent to fraction
        SaveNames[SL][ctrlDict['ClayFract'][0]] = ColName

        # get organic fraction
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % ctrlDict['OrganicFract'][0]
        tmpDF, tmpInfo, tmpNcells = GetVariable( ctrlDict['SoilGridNameFmt'] % ColName, ctrlDict['OrganicFract'][1], ColName, Ncells, INTflag, ctrlDict['RASTER'] )
        SoilParams[ColName] = tmpDF[ColName]
        if tmpNcells == Ncells: WriteDFtoASCIIgrid( SoilParams, GridInfo, ColName, ctrlDict['SoilGridNameFmt'] % ColName )
        if ctrlDict['OrganicFract'][2] == 'P': SoilParams[ColName] = SoilParams[ColName] / 100. # convert percent to fraction
        SaveNames[SL][ctrlDict['OrganicFract'][0]] = ColName

        # get bulk density
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % ctrlDict['BulkDensity'][0]
        tmpDF, tmpInfo, tmpNcells = GetVariable( ctrlDict['SoilGridNameFmt'] % ColName, ctrlDict['BulkDensity'][1], ColName, Ncells, INTflag, ctrlDict['RASTER'] )
        SoilParams[ColName] = tmpDF[ColName]
        if tmpNcells == Ncells: WriteDFtoASCIIgrid( SoilParams, GridInfo, ColName, ctrlDict['SoilGridNameFmt'] % ColName )
        if ctrlDict['BulkDensity'][2] == 'kg/dm3': SoilParams[ColName] = SoilParams[ColName] / 0.001 # convert percent to fraction
        else: sys.stderr.write( "WARNING: Do not recognize units given for bulk density, %s, no conversion done.\n" % ctrlDict['BulkDensity'][2] )
        SaveNames[SL][ctrlDict['BulkDensity'][0]] = ColName

    # get impermiable layer depth
    ColName = ctrlDict['ImpermLayer'][0] # this is not a layer specific value, though Dsmax will be
    tmpDF, tmpInfo, tmpNcells = GetVariable( ctrlDict['SoilGridNameFmt'] % ColName, ctrlDict['ImpermLayer'][1], ColName, Ncells, INTflag, ctrlDict['RASTER'] )
    SoilParams[ColName] = tmpDF[ColName]
    if tmpNcells == Ncells: WriteDFtoASCIIgrid( SoilParams, GridInfo, ColName, ctrlDict['SoilGridNameFmt'] % ColName )
    SaveNames[0][ctrlDict['ImpermLayer'][0]] = ColName
    SaveNames[1][ctrlDict['ImpermLayer'][0]] = ColName

    # get slope
    ColName = ctrlDict['CellSlopeFile'][0] # this is not a layer specific value, though Dsmax will be
    tmpDF, tmpInfo, tmpNcells = GetVariable( ctrlDict['SoilGridNameFmt'] % ColName, ctrlDict['CellSlopeFile'][1], ColName, Ncells, INTflag, ctrlDict['RASTER'] )
    SoilParams[ColName] = tmpDF[ColName]
    if tmpNcells == Ncells: WriteDFtoASCIIgrid( SoilParams, GridInfo, ColName, ctrlDict['SoilGridNameFmt'] % ColName )
    if ctrlDict['CellSlopeFile'][2] == 'D': SoilParams[ColName] = np.tan( SoilParams[ColName] * np.pi / 180. ) # convert degrees to fraction
    elif ctrlDict['CellSlopeFile'][2] == 'P': SoilParams[ColName] = SoilParams[ColName] / 100. # convert percent to fraction
    else: sys.stderr.write( "WARNING: Do not recognize units given for slope, %s, no conversion done.\n"  %  ctrlDict['CellSlopeFile'][2] )
    SaveNames[0][ctrlDict['CellSlopeFile'][0]] = ColName
    SaveNames[1][ctrlDict['CellSlopeFile'][0]] = ColName

    #####
    # Estimate parameters for each layer of the soil database
    #####
    for SL in range(ctrlDict['numSL']):

        # estimate wilting point
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % 'Theta_1500'
        SoilParams[ColName] = swe.Est_Theta_1500( SoilParams[SaveNames[SL][ctrlDict['SandFract'][0]]], 
                                                  SoilParams[SaveNames[SL][ctrlDict['ClayFract'][0]]], 
                                                  SoilParams[SaveNames[SL][ctrlDict['OrganicFract'][0]]] )
        SaveNames[SL]['Theta_1500'] = ColName
        #WriteDFtoASCIIgrid( SoilParams, GridInfo, ColName, ctrlDict['SoilGridNameFmt'] % ColName )

        # estimate field capacity
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % 'Theta_33'
        SoilParams[ColName] = swe.Est_Theta_33( SoilParams[SaveNames[SL][ctrlDict['SandFract'][0]]], 
                                                SoilParams[SaveNames[SL][ctrlDict['ClayFract'][0]]], 
                                                SoilParams[SaveNames[SL][ctrlDict['OrganicFract'][0]]] )
        SaveNames[SL]['Theta_33'] = ColName
        #WriteDFtoASCIIgrid( SoilParams, GridInfo, ColName, ctrlDict['SoilGridNameFmt'] % ColName )

        # estimate SAT-33 kPa moisture
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % 'Theta_S33'
        SoilParams[ColName] = swe.Est_Theta_S33( SoilParams[SaveNames[SL][ctrlDict['SandFract'][0]]], 
                                                 SoilParams[SaveNames[SL][ctrlDict['ClayFract'][0]]], 
                                                 SoilParams[SaveNames[SL][ctrlDict['OrganicFract'][0]]] )
        SaveNames[SL]['Theta_S33'] = ColName
        #WriteDFtoASCIIgrid( SoilParams, GridInfo, ColName, ctrlDict['SoilGridNameFmt'] % ColName )

        # estimate saturated moisture
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % 'Theta_S'
        SoilParams[ColName] = swe.Est_Theta_S( SoilParams[SaveNames[SL]['Theta_33']], 
                                               SoilParams[SaveNames[SL]['Theta_S33']], 
                                               SoilParams[SaveNames[SL][ctrlDict['SandFract'][0]]] )
        SaveNames[SL]['Theta_S'] = ColName
        #WriteDFtoASCIIgrid( SoilParams, GridInfo, ColName, ctrlDict['SoilGridNameFmt'] % ColName )

        # estimate tension at air entry (bubbling pressure) in kPa
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % 'Psi_e'
        SoilParams[ColName] = swe.Est_Psi_e( SoilParams[SaveNames[SL][ctrlDict['SandFract'][0]]],
                                             SoilParams[SaveNames[SL][ctrlDict['ClayFract'][0]]], 
                                             SoilParams[SaveNames[SL]['Theta_S33']] )
        SaveNames[SL]['Psi_e'] = ColName
        #WriteDFtoASCIIgrid( SoilParams, GridInfo, ColName, ctrlDict['SoilGridNameFmt'] % ColName )

        # estimate B
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % 'B'
        SoilParams[ColName] = swe.Est_B( SoilParams[SaveNames[SL]['Theta_33']], 
                                         SoilParams[SaveNames[SL]['Theta_1500']] )
        SaveNames[SL]['B'] = ColName
        #WriteDFtoASCIIgrid( SoilParams, GridInfo, ColName, ctrlDict['SoilGridNameFmt'] % ColName )

        # estimate lambda
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % 'Lambda'
        SoilParams[ColName] = swe.Est_Lambda( SoilParams[SaveNames[SL]['B']] )
        SaveNames[SL]['Lambda'] = ColName
        #WriteDFtoASCIIgrid( SoilParams, GridInfo, ColName, ctrlDict['SoilGridNameFmt'] % ColName )

        # estimate Ksat ***** CHECK UNITS *****
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % 'Ksat'
        SoilParams[ColName] = swe.Est_K_sat( SoilParams[SaveNames[SL]['Theta_S']], 
                                             SoilParams[SaveNames[SL]['Theta_33']], 
                                             SoilParams[SaveNames[SL]['Lambda']] )
        SaveNames[SL]['Ksat'] = ColName
        #WriteDFtoASCIIgrid( SoilParams, GridInfo, ColName, ctrlDict['SoilGridNameFmt'] % ColName )

        # estimate Dsmax
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % 'Dsmax'
        SoilParams[ColName] = swe.Est_Dsmax( SoilParams[SaveNames[SL]['Ksat']], 
                                             SoilParams[SaveNames[SL][ctrlDict['CellSlopeFile'][0]]] )
        SaveNames[SL]['Dsmax'] = ColName
        #WriteDFtoASCIIgrid( SoilParams, GridInfo, ColName, ctrlDict['SoilGridNameFmt'] % ColName )

        # estimate Expt
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % 'Expt'
        SoilParams[ColName] = swe.Est_Expt( SoilParams[SaveNames[SL]['Lambda']] )
        SaveNames[SL]['Expt'] = ColName
        #WriteDFtoASCIIgrid( SoilParams, GridInfo, ColName, ctrlDict['SoilGridNameFmt'] % ColName )

        # estimate Bubble
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % 'Bubble'
        SoilParams[ColName] = swe.Est_Bubble( SoilParams[SaveNames[SL]['Lambda']] )
        SaveNames[SL]['Bubble'] = ColName
        #WriteDFtoASCIIgrid( SoilParams, GridInfo, ColName, ctrlDict['SoilGridNameFmt'] % ColName )

        # estimate Wcr_FRACT
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % 'Wcr_FRACT'
        SoilParams[ColName] = swe.Est_Wcr_FRACT( SoilParams[SaveNames[SL]['Theta_33']] )
        SaveNames[SL]['Wcr_FRACT'] = ColName
        #WriteDFtoASCIIgrid( SoilParams, GridInfo, ColName, ctrlDict['SoilGridNameFmt'] % ColName )
        tmpColName = ColName

        # estimate Wpwp_FRACT
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % 'Wpwp_FRACT'
        SoilParams[ColName] = swe.Est_Wpwp_FRACT( SoilParams[SaveNames[SL]['Theta_1500']] )
        SaveNames[SL]['Wpwp_FRACT'] = ColName
        #WriteDFtoASCIIgrid( SoilParams, GridInfo, ColName, ctrlDict['SoilGridNameFmt'] % ColName )

        # check that Wcr is greater than Wpwp, if not then switch them
        for idx in range(len(SoilParams[ColName])):
            if SoilParams[tmpColName][idx] < SoilParams[ColName][idx]:
                tmpValue = SoilParams[tmpColName][idx]
                SoilParams[tmpColName][idx] = SoilParams[ColName][idx]
                SoilParams[ColName][idx] = tmpValue

        #####
        # The following are set to constant values for all soil layers, they are included
        # to support the possability of something better being done in the future :)
        #####

        # set soil density
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % 'soil_density'
        SoilParams[ColName] = pd.DataFrame( [ 2650.0 ]*len(SoilParams) )
        SaveNames[SL]['soil_density'] = ColName

        # set soil density
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % 'soil_dens_org'
        SoilParams[ColName] = pd.DataFrame( [ 1300.0 ]*len(SoilParams) )
        SaveNames[SL]['soil_dens_org'] = ColName

        # set soil moisture diffusion
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % 'phi_s'
        SoilParams[ColName] = pd.DataFrame( [ -999.0 ]*len(SoilParams) )
        SaveNames[SL]['phi_s'] = ColName

        # set residual moisture
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % 'resid_moist'
        SoilParams[ColName] = pd.DataFrame( [ 0.0 ]*len(SoilParams) )
        SaveNames[SL]['resid_moist'] = ColName

        # set quartz content
        ColName = ctrlDict['SoilLayerNameFmt'][SL] % 'quartz'
        SoilParams[ColName] = SoilParams[SaveNames[SL][ctrlDict['SandFract'][0]]].copy()
        SaveNames[SL]['quartz'] = ColName

    # add grid cell active mask
    if ctrlDict['SimCellMask']:
        ColName = 'run_cell'
        tmpDF, tmpInfo, tmpNcells = GetVariable( ctrlDict['SimCellMask'], np.NaN, ColName, Ncells, True, ctrlDict['RASTER'] )
        SoilParams[ColName] = tmpDF[ColName]

    # process other data files required for the soil parameter file
    for idx in range(len(ctrlDict['OtherFiles'])):
        ColName = ctrlDict['OtherFiles'][idx]['colName']
        tmpDF, tmpInfo, tmpNcells = GetVariable( ctrlDict['OtherFiles'][idx]['name'], np.NaN, ColName, Ncells, INTflag, ctrlDict['RASTER'] )
        SoilParams[ColName] = tmpDF[ColName]

    # print summary statistics for soil parameters
    print SoilParams.describe()

    # Write soil properties to output files
    SoilParams.to_csv( ctrlDict['outSoilDataTable'] )
    tmpDF = pd.DataFrame( SaveNames )
    tmpDF.to_csv( ctrlDict['outSoilAttribTable'] )




