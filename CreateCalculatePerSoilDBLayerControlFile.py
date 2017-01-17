#!/bin/env python
# Created on April 6, 2016
#  by Keith Cherkauer
#
# This script prompts the user for information related to their personal use
# of the CalculatePerSoilDBLayer.py script, and uses it to create a JSON
# formatted control file for application of the program.
#
# Modifications:

import json
import sys
import numpy as np

def HandlePrompt( inVar ):
    PromptStr = inVar[0]
    Default = inVar[2]
    inVar = inVar[1]
    while True:
        tmpVar = raw_input("%s --> [%s] = " % ( PromptStr, inVar ) )
        if Default and tmpVar == "":
            return inVar
        elif tmpVar != "":
            return tmpVar
        else:
            print "--> Current variable has no default, must be set!!!  Try again."

# check for command l;ine arguments
if len(sys.argv) <= 1 or len(sys.argv) > 2:
    sys.stderr.write( "\nUsage: %s <filename>\n\n\t<filename> is the name of the control file to be created.\n\n" % sys.argv[0] )
    sys.exit()

# build control file
tmpDict = {}

# setup default values for HWSD
SoilDB = { 'HWSD': { 'numSL': 2, 
                     'SoilLayerThick': [ ( 'SOIL_THICK', 0.300, 'm' ), ( 'SOIL_THICK', 0.700, 'm' ) ], 
                     'SoilLayerNameFmt': [ "T_%s", "S_%s" ], 
                     'SandFract': ( 'SAND', np.NaN, 'P' ), # sand fraction of soil, no default, 'P' in percent, 'F' as a fraction
                     'ClayFract': ( 'CLAY', np.NaN, 'P' ), # clay fraction of soil, no default, 'P' in percent, 'F' as a fraction
                     'OrganicFract': ( 'OC', np.NaN, 'P' ), # organic matter fraction of soil, no default, 'P' in percent, 'F' as a fraction
                     'GravelFract': ( 'GRAVEL', np.NaN, 'P' ), # gravel fraction of soil (optional, set to 'None' if not provided), no default, 'P' in percent, 'F' as a fraction
                     'BulkDensity': ( 'BULK_DENSITY', 1, 'kg/dm3' ), # or should this be the REF_BULK_DENSITY column?
                     'ImpermLayer': ( 'IL', 1, 'None' ) } }

# setup dictionary for user queries, each is a tuple with { Prompt Text, Default Value, Default True/False )
#    as the user is prompted to respond, the tuple is replaced with the final input value.
tmpDict["RASTER"] = ( "Are input soils in ArcGIS ASCII raster format (True/False)?", "False", True )
tmpDict["outSoilDataTable"] = ( "Name of data table to be created.", "No Default: Input data file name", False )
tmpDict["outSoilAttribTable"] = ( "Name of output attribute table.", "No Default: Input attribute file name", False )
tmpDict["SoilGridNameFmt"] = ( "Input directory for output from SimplifySoil files", "SoilsData/VICsetup/", True )
tmpDict["SoilDB"] = ( "Source database for soils data, selected from (%s)" % ','.join(SoilDB.keys()), "HWSD", True )
tmpDict['CellNumFile'] = ( "Name of the file to define simulation cell numbers.", "No Default: Input file name", False )
tmpDict['SimCellMask'] = ( "Name of the file that defines active cells.", "No Default: Use cell number file to run all", False )
tmpDict['CellSlopeFile'] = ( "Name of the file with simulation cell slope in percent (used for setting Dsmax).", "No Default: Input file name", False )
tmpDict['OtherFiles'] = []

# prompt user for other options
print "==========\nThese first options MUST be set by the user:\n"

tmpDict["SoilDB"] = HandlePrompt( tmpDict["SoilDB"] )
if tmpDict['SoilDB'].upper() == "HWSD":
    print "Program setup to use the Harmonized World Soil Database"
    tmpDict.update(SoilDB[tmpDict["SoilDB"]])
    del(tmpDict['SoilDB'])
else:
    print "Requested soil database, %s, is not currently available." % tmpDict['SoilDB'] 
    sys.exit(-1)

# input file names
tmpDict['CellNumFile'] = HandlePrompt( tmpDict["CellNumFile"] )

tmpDict['SimCellMask'] = HandlePrompt( tmpDict["SimCellMask"] )

tmpDict['OtherFiles'] = tmpDict['OtherFiles'] + [ { 'colName':'elev', 'name':HandlePrompt( ("Name of file containing elevations for all simulations cells.", "No Default: Input file name", False ) ) } ]

tmpDict['OtherFiles'] = tmpDict['OtherFiles'] + [ { 'colName':'annual_prec', 'name':HandlePrompt( ("Name of file containing annual average precipitation for all simulations cells.", "No Default: Input file name", False ) ) } ]

tmpDict['CellSlopeFile'] = ( HandlePrompt( tmpDict["CellSlopeFile"] ), 1, 'P' )

# output file names
tmpDict['outSoilDataTable'] = HandlePrompt( tmpDict["outSoilDataTable"] )

tmpDict['outSoilAttribTable'] = HandlePrompt( tmpDict["outSoilAttribTable"] )

tmpVar = HandlePrompt( tmpDict["RASTER"] )
if tmpVar.lower() == "false": tmpDict['RASTER'] = False
else: tmpDict['RASTER'] = True

print "==========\nThe next options will have to be changed, if you did not put the files in the default directory\n"

tmpDict['SoilGridNameFmt'] = HandlePrompt( tmpDict["SoilGridNameFmt"] )
if tmpDict['RASTER']:
    # gridded output filename format
    tmpDict['SoilGridNameFmt'] = tmpDict['SoilGridNameFmt'] + "SoilProperty_\%s.asc"
else:
    # output to single soil table file
    tmpDict['SoilGridNameFmt'] = tmpDict['SoilGridNameFmt'] + "SoilProperty_\%s.txt"

print "==========\nThe following selections are required or optional depending on VIC model simulation options, you must provide a valid file name in the correct format.\n"

tmpDict['OtherFiles'] = tmpDict['OtherFiles'] + [ { 'colName':'avg_T', 'name':HandlePrompt( ("Name of file containing average annual soil temperature for all simulations cells - used to define bottom boundary for FULL ENERGY simulations.", "No Default: Input file name", False ) ) } ]

tmpDict['OtherFiles'] = tmpDict['OtherFiles'] + [ { 'colName':'off_gmt', 'name':HandlePrompt( ("Name of file containing hour offset from GMT.", "No Default: Input file name", False ) ) } ]

tmpDict['OtherFiles'] = tmpDict['OtherFiles'] + [ { 'colName':'July_Tavg', 'name':HandlePrompt( ("Name of file containing average annual July air temperature for all simulations cells - Optionally used to define treeline.", "No Default: Input file name", False ) ) } ]

# use json to write file
f = open( sys.argv[1], 'w' )
json.dump( tmpDict, f )
f.close()
