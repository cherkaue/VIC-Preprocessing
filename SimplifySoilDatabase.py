#!/bin/env python
# Created on May 15, 2015
#  by Keith Cherkauer
#
# This script was written to simplify the full database for the harmonized world
# soils database (HWSD) for use in setting up the VIC model.  It requires four 
# files to run correctly:
#   1. high resolution file with soil IDs (MUID, MU_GLOBAL, etc) for each element
#   2. high resolution mapping of VIC cell numbers to the same elements used by 
#      the soil file
#   3. a file with cell numbers at VIC resolution
#   4. a comma separated variable (CSV) file of the soil properties.
#
# The first three files can be provided as ArcInfo ASCII raster grid files, or
# XYZ files, where X is 'lat', Y is 'lon', and Z is 'value'.  Note that these names
# must appear in the header line of the input files.  For gridded data the first 
# two files must have the same number of cells, if working with XYZ data then the 
# first two files should have the same number of lines.  The third file should have 
# fewer cells or lines than the first two files.
#
# The HWSD dataset yields more than one soil type per grid cell in many cases (even 
# at the 30 m resolution).  Because the VIC model handles only one soil type per cell 
# the database records must be simplified.  This script simplifies data by finding the 
# majority soil type within the VIC solution cell.
#
# Simplification rules
# - Grid MU_GLOBAL values can map to multiple soil types with "SHARE" indicating 
#   the fraction of the grid cell represented by the current soil, and "SEQ"
#   indicating its rank in the grid cell.
#   - non-soils, so those with no depth defined are removed.
#   - soil types with shallow reporting depths (e.g., less than 100 cm) will 
#     be reduced in priority.
#   - majority soil type will be assigned to each VIC cell.
#   - cells with repeating values will be filtered of fractions less 5%
#   - then characteristics of soils will be compared - if similar than values 
#     will be averaged, otherwise only the majority will be kept to prevent 
#     introduction of unreasonable averaged parameters.
#
# Things to watch out for:
# - the column required for the Soil Unit Mapping Code will vary by application.
#   I have found that SU_CODE90 works for the EU, while SU_CODE74 works for the 
#   US.  Despite what is in the documentation, there is no general SU_CODE column.
# - if changing the soil database, then the check to see if the current MUID 
#   contains an actual soil (inSoilDF.ISSOIL) will have to be updated to point to
#   the correct column.
#
# Modifications
# 2015-10-12 Updated to work with non-gridded data files. Input file must be tab 
#   delmiited, and have heaers 'lat', 'lng' and 'value', where 'value' must be the 
#   column with the unique identifier for the VIC model cell.  Other columns are
#   allowed, but will not be used.  KAC
# 2015-10-12 Improved performance by handling soil unit code entries that are
#   not in fact soils (ISSOIL = 0 in the HWSD).  Processing will report a warning 
#   about the MUID, but otherwise processing will continue.    KAC
# 2016-02-16 Merged with actual finished version from Hungary, so a few minor 
#   mostly to documentation, plus program now dumps data frame to a CSV file
#   if output from program is not supposed to be gridded.    KAC
# 2016-02-18 Updated to make better use of tabular (non-gridded) data as both
#   an input and an output.                                  KAC
# 2016-03-31 modifeid to output non-gridded data into XYZ style files, with one
#   file per parameter.  This is what CalculatePerSoilDBLayer.py is looking for
#   and perhaps makes more sense then writing everything to a single file as it
#   is easier to add files for the missing parameters that to edit the fille 
#   database.                                                KAC  ***** It has not been checked yet !!!! ***
#

import sys
import pandas as pd
import numpy as np
import read_arcinfo_files as rwarc
import json

#####
# USER: Define filenames and important data columns
#####
if len(sys.argv) <= 1 or len(sys.argv) > 2:
    # no control file given, exit gracefully
    sys.stderr.write( "\nUsage: %s <filename>\n\n\tThis script was written to simplify the full database for the harmonized \n\tworld soils database (HWSD) for use in setting up the VIC model.\n\n\t<filename> is the name of the control file that allows this program to \n\t\twork with your data.  Use the script CreateSoilSimplifyControlFile.py \n\t\tto build a control file.  \n\n\tExample files should also be in the same directory as the executable script.\n\n" % sys.argv[0] )
    sys.exit()

else:
    # control file provided, parse and get started
    fin = open( sys.argv[1], 'r' )
    ctrlDict = json.load( fin )
    fin.close()

    # adjust a few variables because the default type is problematic
    ctrlDict['inSoilID'] = ctrlDict['inSoilID'].encode('ascii','ignore')
    ctrlDict['inSoilFract'] = ctrlDict['inSoilFract'].encode('ascii','ignore')

#####
# Read in all data files
#####

# read in soil grid cell numbers
SoilGridDict = rwarc.read_ARCINFO_ASCII_grid( ctrlDict['inSoilGrid'] )

# read VIC grid cell numbers at soil grid cell resolution
CellNumDict = rwarc.read_ARCINFO_ASCII_grid( ctrlDict['inCellNum'] )

# read VIC simulation cell number structure for output
if ctrlDict['outCellNumGrid']:
    # output will be gridded file
    outCellNumDict = rwarc.read_ARCINFO_ASCII_grid( ctrlDict['outCellNum'] )
else:
    # output is not gridded, so must read from a TAB delimited file
    tmpDF = pd.read_table( ctrlDict['outCellNum'], sep=r"\s*", engine='python' ) # headers should be 'lat', 'lng' and 'value' other columns are allowed
    tmpDict = tmpDF.to_dict('records')
    outCellNumDict = { 'Ncells': len(tmpDict), 'cell': tmpDict }

# read in soil database
inSoilDF = pd.read_csv( ctrlDict['inSoilTable'], index_col=ctrlDict['inSoilID'], dtype={'SU_SYM74': object, 'SU_SYM85': object} )
inSoilDF = inSoilDF[inSoilDF.ISSOIL != 0] # remove row if it is not a soil

# check that soil and input cell number grids are same size
if SoilGridDict["Ncells"] != CellNumDict["Ncells"]:
    sys.stderr.write( "ERROR: Soils (%i) and input cell number (%i) grids are not the same resolution, this analysis will stop." % ( SoilGridDict["Ncells"], CellNumDict["Ncells"] ))
    sys.exit()

#####
# Convert from read_arcinfo_file dictionary structure to NumPy arrays
#####

SoilGridArray = np.arange(SoilGridDict["Ncells"])
CellNumArray = np.arange(CellNumDict["Ncells"])
for Idx in range( SoilGridDict["Ncells"] ):
    SoilGridArray[Idx] = SoilGridDict["cell"][Idx]["value"]
    CellNumArray[Idx] = CellNumDict["cell"][Idx]["value"]

#####
# Build series from input soil grid data to control analysis
#####

analysisDF = pd.Series( SoilGridArray, index=CellNumArray )

#####
# Build output dataframe
#####

TmpDFa = pd.DataFrame(np.zeros( (outCellNumDict["Ncells"],), dtype=[(ctrlDict['inSoilID'],'i4'),(ctrlDict['inSoilFract'],'f8'),('NUM_TYPES','i4')] ))
TmpDFb = pd.DataFrame(np.zeros( (outCellNumDict["Ncells"],len(ctrlDict['outParams'])) ), columns=ctrlDict['outParams'] )
outDF = pd.concat( [ TmpDFa, TmpDFb ], axis=1 )

NonSoils = {} # establish dictionary for storing use of non-soils
RemoveCells = [] # array to store VIC simulation cells that do not have valid soil data

#####
# Determine majority soil type for every VIC model grid cell
#####

for cidx in range(outCellNumDict["Ncells"]):
    # Loops through each output simulation cell (gridded or not) and identifies
    # all soil types within the output cell.  Does this by mathcing the high 
    # resolution gridded values for soil MUID and VIC model simulation cell
    # numbers. 

    CurrCell = outCellNumDict["cell"][cidx]["value"] # set to current VIC cell number

    if CurrCell in analysisDF.index:
        # Output cell number is found in the high resolution cell database
        #   - THis is what should happen, if preprocessing was done correctly!

        # Get all soil MUIDs for the current VIC simulation cell
        MUIDsCurrCell = pd.Series(analysisDF.loc[CurrCell])
        # Count the number of MUIDs within the VIC simulation cell
        NumMUIDs = len(MUIDsCurrCell)
        # Get the list of unique MUIDs within the current VIC simulation cell 
        MU_IDs = MUIDsCurrCell.unique()
        # Remove invalid (negative) soil types
        MU_IDs = MU_IDs[MU_IDs>0]
        # make sure that the VIC cell includes valid soil types
        if len(MU_IDs) > 0:
            # Should only get here if there is at least one real soil within 
            # current VIC simulation cell.

            # Get soils information for selected MUIDs
            try: 
                SoilDataCurrCell = inSoilDF.loc[MU_IDs]
                # get soil unit IDs from mixtures within each soil grid cell
                SoilCodeCurrCell = pd.Series( SoilDataCurrCell[ctrlDict['inSoilCode']] ).unique()
                # build Series to hold extent of soils within the VIC model grid cell
                SoilAreas = pd.Series( [0.0]*len(SoilCodeCurrCell),index=SoilCodeCurrCell)
                # loop through all soil grid cells within the VIC grid cell and sum soil unit coverage
                for muid in MU_IDs:
                    if not muid in inSoilDF.index:
                        if not muid in NonSoils.keys():
                            NonSoils["%i" % muid] = 0
                        NonSoils["%i" % muid] = NonSoils["%i" % muid] + 1
                    else:
                        TmpSoilData = inSoilDF.loc[muid]
                        if len(TmpSoilData.shape) == 2:
                            # for MU IDs with more than one soil unit
                            NumCellsCurrMUID = MUIDsCurrCell[MUIDsCurrCell==muid].count()
                            for sinfo in range(len(TmpSoilData.index)):
                                try:
                                    SoilAreas[TmpSoilData.iloc[sinfo][ctrlDict['inSoilCode']]] = SoilAreas[TmpSoilData.iloc[sinfo][ctrlDict['inSoilCode']]] + TmpSoilData.iloc[sinfo][ctrlDict['inSoilFract']] / 100. / NumMUIDs * NumCellsCurrMUID
                                except IndexError, errstr:
                                    print "ERROR: %s" % errstr
                                    print "sinfo", sinfo
                                    print "inSoilFract", ctrlDict['inSoilFract']
                                    print "NumCellsCurrMUID", NumCellsCurrMUID
                                    print "NumMUIDs", NumMUIDs
                                    print TmpSoilData
                                    print TmpSoilData.shape, MU_IDs, len(TmpSoilData.shape)
                                    sys.exit()
                        elif len(TmpSoilData.shape) == 1:
                            # for MU IDs with more than one soil unit
                            NumCellsCurrMUID = MUIDsCurrCell[MUIDsCurrCell==muid].count()
                            try:
                                SoilAreas[TmpSoilData[ctrlDict['inSoilCode']]] = SoilAreas[TmpSoilData[ctrlDict['inSoilCode']]] + TmpSoilData[ctrlDict['inSoilFract']] / 100. / NumMUIDs * NumCellsCurrMUID
                            except IndexError, errstr:
                                print "ERROR: %s" % errstr
                                print "sinfo", sinfo
                                print "inSoilFract", ctrlDict['inSoilFract']
                                print "NumCellsCurrMUID", NumCellsCurrMUID
                                print "NumMUIDs", NumMUIDs
                                print TmpSoilData
                                print TmpSoilData.shape, MU_IDs, len(TmpSoilData.shape)
                                sys.exit()
                        else:
                            # really this is not possible, but better check
                            sys.stderr.write("ERROR: Don't know how you got here, it shouldn't be possible." )
                            sys.exit()

                TmpCode = SoilAreas.idxmax()
                TmpFract = SoilAreas.max()
                TmpNumSoil = len(SoilAreas.index)
                #print CurrCell, SoilAreas.sum(), TmpCode, TmpFract, TmpNumSoil, TmpSoilData.shape, MU_IDs, len(TmpSoilData.shape)

                # store results for current VIC grid cell in the output data frame
                outDF.set_value(cidx,ctrlDict['inSoilID'],TmpCode) # store soil unit ID number
                outDF.set_value(cidx,ctrlDict['inSoilFract'],TmpFract) # store fractional coverage of grid cell
                outDF.set_value(cidx,'NUM_TYPES',TmpNumSoil) # store number of soil units present
                tmpSoilUnitInfo = SoilDataCurrCell[SoilDataCurrCell[ctrlDict['inSoilCode']] == SoilAreas.idxmax()]
                for pidx in ctrlDict['outParams']:
                    outDF.set_value(cidx,pidx,tmpSoilUnitInfo.iloc[0][pidx])

                #SoilCodeCurrCell = pd.Series( SoilDataCurrCell[ctrlDict['inSoilCode']] ).unique()

            except KeyError, errstr:
                # VIC simulation cell does not contain a real soil in this database
                print "Current VIC cell does not contain a real soil, so it will be skipped during processing.  You will need to assess this problem before the model will run."
                print "-->>", outCellNumDict["cell"][cidx]
                RemoveCells.append(CurrCell)


        else:
            sys.stderr.write( "No valid soil types found in cell %i (%f,%f): Cell Number %i\n" % ( cidx, outCellNumDict['cell'][cidx]['lat'], outCellNumDict['cell'][cidx]['lng'], outCellNumDict['cell'][cidx]['value'] ) )
            RemoveCells.append(CurrCell)

    else:
        # Output cell number is NOT found in the high resolution cell database

        # Since the high resolution should be generate from the low resolution, this
        # is a problem related to initial processing, and shoul dbe investigated
        # before continuing to use this program to process soil data!

        sys.stderr.write( "WARNING: Current output cell, %i (%f,%f): Cell Number %i, was not found in the high resoltuion grid of output cells.  That suggests you need to regenerate one or the other to match!  It is possible for this to occur for the HUC cell approach, when a HUC is smaller than the resolution of the soil grid so was not included in the gridding process.  We'll skip the HUC for now, but you should check this before proceeding. \n" % ( cidx, outCellNumDict['cell'][cidx]['lat'], outCellNumDict['cell'][cidx]['lng'], outCellNumDict['cell'][cidx]['value'] ) )
        RemoveCells.append(CurrCell)


#print "outDF statistics:", outDF.describe()

#print "NonSoils:", NonSoils
print "outCellNumGrid", ctrlDict['outCellNumGrid']
#print "tmpDF statistics:", tmpDF.describe()

#####
# Merge latitude and longituce information with output DataFrame for ungridded data
#####
if not ctrlDict['outCellNumGrid']:
    outDF = pd.concat( [ outDF, tmpDF[['lat','lng']] ], axis=1 )

#####
# Need to remove problem grid cells here
#####
print "RemoveCells:", RemoveCells
outDF = outDF.drop(outDF.index[RemoveCells])

#####
# Write majority MUID for each VIC model simulation cell
#####

ListParams = outDF.columns
print outDF.columns
if not ctrlDict['outCellNumGrid']:
    ListParams = ListParams.drop(['lat','lng']) # drop latitude and longitude from index list

for param in ListParams:

    # create output file name with extracted values for current parameter
    outFileName = ctrlDict['outSoilNameFmt'] % param
    print "Writing file %s" % outFileName

    #####
    # Output files will be gridded
    #####

    if ctrlDict['outCellNumGrid']:
        # output will be gridded files
        for idx in range(outCellNumDict["Ncells"]):
            outCellNumDict["cell"][idx]["value"] = outDF.get_value(idx,param) 
        # write output to ASCII gridded ESRI raster file format
        if outDF.get_value(0,param) == 'i': # integer only file
            rwarc.write_ARCINFO_ASCII_grid(outFileName,outCellNumDict,INTflag=1)
        else: # floating point file
            rwarc.write_ARCINFO_ASCII_grid(outFileName,outCellNumDict,INTflag=0)

    #####
    # Output files will be XYZ files
    #####

    else:
        # generate a single output table with all parameters
        outDF.to_csv( outFileName, columns=['lng','lat',param], header=False, sep='\t' )
