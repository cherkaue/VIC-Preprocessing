#!/bin/env python
# Created on July 6, 2015
# by Keith Cherkauer
#
# This script uses ASCII reaster grides of fractional land use coverage
# in IGBP classification to create a standard VIC model VEG_PARAM file.
# Ths VEG_PARAM file will include wetlands / ponded area, but will
# include urban area as short grass.
#
# Things to change:
# - need to read in the ponding fraction file
# - sum up the wetland and open water fractions and check that they are 
#   not bigger than the ponding fraction file
# - identify largest vegetation fraction, and assign it to be the wetland 
#   fraction (so must appear as the first item in the cell vegetation list)
# - extra area can be assigned to the same veg type later in the file
#

import sys, os
import pandas as pd
import numpy as np
import read_arcinfo_files as rwarc

sys.path.append("./Scripts")

from VicPreprocessingLibrary import *

# IGBP remap table
VegRemapDF = pd.read_csv( 'Scripts/IGBP_to_VicPonded.csv' )
# Lake-Wetland parameter file
PondFile = 'WetlandFiles/TiszaWetland.Outlet.0m.txt'

# set to True is using gridded data, otherwise set to False
useRASTER = True

# format string that will be used to create file names for all IGBP ASCII raster files
LURasterFormat = 'LandUseData/landusefract_mcd12q1_2005_igbp_%03i_gcs_ar.asc'
# table (CSV) of default vegetation parameters for each type
VegExtraParamFile = 'Scripts/VicVegParamDefaults.csv'
# file with VIC simulation cell numbers
CellNumFile = 'ClimCellNum_GCS_AR.asc'
# set threshold for minimum vegetation cover to include in the simulations (as fraction)
SimFractThres = 0.05
# set name of vegetation parameter file to be created
outVegParamFile = 'LandUseData/VegParam_Ponding.txt'

# set other variables
VegParamDF = pd.read_csv( VegExtraParamFile, index_col=0 )

# open cell number file
Ncells = -99
IGBPFractDF, GridInfo, Ncells = GetVariable( CellNumFile, np.NaN, 'CellNum', Ncells, RASTER=useRASTER)

# open all IGBP land use raster files
for Remap in range(len(VegRemapDF)):
    tmpDF, GridInfo, Ncells = GetVariable( LURasterFormat % VegRemapDF['inClass'][Remap], np.NaN, VegRemapDF['inClass'][Remap], Ncells, RASTER=useRASTER)
    IGBPFractDF = pd.concat( [ IGBPFractDF, tmpDF[VegRemapDF['inClass'][Remap]] ], axis=1 )
IGBPFractDF = IGBPFractDF[pd.notnull(IGBPFractDF['CellNum'])] # remove any undefined cell numbers
IGBPFractDF['CellNum'] = IGBPFractDF['CellNum'].astype(np.int) # convert cell number to integer for indexing
IGBPFractDF.set_index('CellNum', inplace=True) # set cell number as the index

# reduce IGBP land use classification to VIC standard
OutVegClass = VegRemapDF['outClass'].unique()
NumOutVeg = len(OutVegClass)
VegFractDF = pd.DataFrame( 0, index=IGBPFractDF.index, columns=OutVegClass )
IGBPFractDF.fillna(0, inplace=True) # replace NaNs with 0s so that the math works
for Remap in range(len(VegRemapDF)):
    if VegRemapDF['outClass'][Remap] < 254:
        print "Remapping %i to %i" % ( VegRemapDF['inClass'][Remap], VegRemapDF['outClass'][Remap] )
        VegFractDF[VegRemapDF['outClass'][Remap]] = ( VegFractDF[VegRemapDF['outClass'][Remap]] 
                                                      + IGBPFractDF[VegRemapDF['inClass'][Remap]] )

# remove unassigned land cover
VegFractDF = VegFractDF.drop(255,axis=1)
OutVegClass = VegFractDF.columns
NumOutVeg = len(OutVegClass)

# remove cells without any vegetation
VegFractDF = VegFractDF.dropna(thresh=0)
CellList = VegFractDF.index

# open lake-wetland parameter file and get maximum wetland fraction
fpond = open( PondFile, 'r' )
lines = fpond.readlines()
fpond.close()
tmpCellNum = []
tmpWetFrac = []
while len(lines) > 0:
    tmpStr = lines[0].strip().split()
    tmpCellNum = tmpCellNum + [ int(tmpStr[0]) ]
    del lines[0]
    tmpStr = lines[0].strip().split()
    tmpWetFrac = tmpWetFrac + [ float(tmpStr[1]) ]
    del lines[0]
tmpDict = { 'WetFrac':tmpWetFrac, 'LUid':[0]*len(tmpWetFrac) }
WetFracDF = pd.DataFrame( tmpDict, index=np.asarray(tmpCellNum,dtype=np.int) )

# check that wet fraction (class 20) does not exceed ponding fraction
tmpFillS = pd.Series( 0,index=['WetFrac','LUid'] )
for Cell in VegFractDF.index:
    WaterFrac = VegFractDF.loc[Cell,20]
    tmpSeries = VegFractDF.loc[Cell]
    try:
        tmpSeries = tmpSeries[tmpSeries.index < 20].order(ascending=False)
    except AttributeError, errstr:
        print errstr
        print tmpSeries
        sys.exit()
    MajorLU = tmpSeries.index[0]
    if not Cell in WetFracDF.index:
        # there is no ponded fraction in the current cell
        tmpFillS.name = Cell
        WetFracDF.append(tmpFillS)
    # there is a ponded fraction defined for the current cell
    WetFracDF.set_value(Cell,'LUid',MajorLU)
    # resolve change of LU fractions to fit ponded area
    if WaterFrac > WetFracDF['WetFrac'][Cell]:
        # lake and wetland LU area larger than ponded area fraction
        # - rescale L&WL LU to match ponded area
        # - rescale other LU types to fill
        sys.stderr.write( "Cell %i: Lake and wetland LU area larger than ponded area\n" % Cell )
    else:
        tmpSeries.iloc[0] = tmpSeries.iloc[0] - ( WetFracDF['WetFrac'][Cell] - WaterFrac )
        if tmpSeries.iloc[0] < 0:
            # ponded area larger than L&WL LU plus MajorLU
            # - MajorLU will be fully in ponded area
            # - All other LU rescale to fille rest of cell
            sys.stderr.write( "Cell %i: Lake and wetland LU area larger than ponded area\n" % Cell )
            tmpSeries.iloc[0] = 0
        else:
            # ponded area smaller than L&WL LU plus MajorLU
            # - MajorLU will be fully in ponded area
            # - All other LU rescale to fille rest of cell
            tmpSeries.iloc[0] = tmpSeries.iloc[0] - WetFracDF['WetFrac'][Cell]
    # rescale all other LU types to fit in space
    VegFractDF.loc[Cell,20] = WetFracDF['WetFrac'][Cell] # set ponded area fraction
    VegFractDF.loc[Cell,MajorLU] = tmpSeries.iloc[0] # reset major LU type fract

# check that vegetation cover reaches minimum threshold and rescale for missing fractions
for Cell in VegFractDF.index:
    tmpSeries = VegFractDF.loc[Cell]
    tmpSeries = tmpSeries[tmpSeries.index < 20].order(ascending=False)
    NonWaterFrac = 1. - WetFracDF['WetFrac'][Cell]
    # rescale vegetation fractions so that they total 100% - ponded area
    VegFractSum = tmpSeries.sum()
    if VegFractSum > SimFractThres:
        tmpSeries = tmpSeries / VegFractSum * NonWaterFrac
        # remove vegetation fractions less than defined threshold
        tmpSeries[tmpSeries < SimFractThres] = np.NaN
        # repeat rescale process so that vegetation fractions total to 100% after removals
        VegFractSum = np.ma.masked_invalid(tmpSeries).sum()
        tmpSeries = tmpSeries / VegFractSum * NonWaterFrac
        VegFractDF.loc[Cell] = tmpSeries
    VegFractDF.set_value(Cell,20,WetFracDF['WetFrac'][Cell]) # replace NaNs with ponded fraction to complete file

# create vegetation parameter file
fout = open( outVegParamFile, 'w' )
for Cell in VegFractDF.index:
    # write first line for current cell: <cell number> <number of veg types>
    tmpSeries = VegFractDF.loc[Cell]
    #print "BEFORE:", tmpSeries
    tmpSeries = tmpSeries[tmpSeries.index>0] # exclude bare soil
    #print "AFTER 1", tmpSeries
    tmpSeries = tmpSeries[tmpSeries.index<254] # exclude unused land use types
    #print "AFTER 2", tmpSeries
    tmpSeries = tmpSeries[tmpSeries.notnull()]
    #print "AFTER 3", tmpSeries
    #sys.exit()
    NumVegTypes = len(tmpSeries)
    if NumVegTypes > 0:
        fout.write( "%i\t%i\n" % ( Cell, NumVegTypes ) )
        #print "NumVegType", NumVegTypes
        #print tmpSeries
        # write lines for each vegetation type in the current cell
        # write ponded fraction first
        if 20 in tmpSeries.index:
            if tmpSeries[20] > 0:
                # there is a ponded fraction defined
                #print "Found ponded area"
                vegclass = WetFracDF['LUid'][Cell]
                VegExtraStr = '\t'.join(["%.2f" % number for number in VegParamDF.loc[WetFracDF['LUid'][Cell]]])
                fout.write( "\t%i\t%.5f\t%s\n" % ( WetFracDF['LUid'][Cell], WetFracDF['WetFrac'][Cell], VegExtraStr ) )
        # now write everything else
        tmpSeries = tmpSeries[tmpSeries.index < 20]
        for vegclass in range(NumVegTypes-1):
            VegExtraStr = '\t'.join(["%.2f" % number for number in VegParamDF.loc[tmpSeries.index[vegclass]]])
            fout.write( "\t%i\t%.5f\t%s\n" % ( tmpSeries.index[vegclass], tmpSeries[tmpSeries.index[vegclass]], VegExtraStr ) )

fout.close()

    
