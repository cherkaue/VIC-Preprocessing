#!/bin/env python
# Created on July 6, 2015
# by Keith Cherkauer
#
# This script uses ASCII reaster grides of fractional land use coverage
# in IGBP classification to create a standard VIC model VEG_PARAM file.
# Ths VEG_PARAM file does not include wetlands / ponded area, but will
# include urban area as short grass.
#

import sys, os
import pandas as pd
import numpy as np
import read_arcinfo_files as rwarc

sys.path.append("./Scripts")

from VicPreprocessingLibrary import *

# IGBP remap table
VegRemapDF = pd.read_csv( 'Scripts/IGBP_to_VicNoWater.csv' )

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
outVegParamFile = 'LandUseData/VegParam_NoPonding.txt'

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

# check that vegetation cover reaches minimum threshold and rescale for missing fractions
for Cell in VegFractDF.index:
    tmpArray = VegFractDF.loc[Cell].values
    # rescale vegetation fractions so that they total 100%
    VegFractSum = np.ma.masked_invalid(tmpArray).sum()
    if VegFractSum > SimFractThres:
        tmpArray = tmpArray / VegFractSum
        # remove vegetation fractions less than defined threshold
        tmpArray[tmpArray < SimFractThres] = np.NaN
        # repeat rescale process so that vegetation fractions total to 100% after removals
        VegFractSum = np.ma.masked_invalid(tmpArray).sum()
        tmpArray = tmpArray / VegFractSum
        VegFractDF.loc[Cell] = tmpArray

# create vegetation parameter file
fout = open( outVegParamFile, 'w' )
for Cell in VegFractDF.index:
    # write first line for current cell: <cell number> <number of veg types>
    tmpSeries = VegFractDF.loc[Cell]
    tmpSeries = tmpSeries[tmpSeries.index>0] # exclude bare soil
    tmpSeries = tmpSeries[tmpSeries.index<254] # exclude unused land use types
    tmpSeries = tmpSeries[tmpSeries.notnull()]
    NumVegTypes = len(tmpSeries)
    if NumVegTypes > 0:
        fout.write( "%i\t%i\n" % ( Cell, NumVegTypes ) )
        # write lines for each vegetation type in the current cell
        for vegclass in range(NumVegTypes):
            VegExtraStr = '\t'.join(["%.2f" % number for number in VegParamDF.loc[tmpSeries.index[vegclass]]])
            fout.write( "\t%i\t%.2f\t%s\n" % ( tmpSeries.index[vegclass], tmpSeries[tmpSeries.index[vegclass]], VegExtraStr ) )

fout.close()

    
