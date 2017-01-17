#!/bin/env python
# Created on January 19, 2016
#  by Keith Cherkauer
#
# This script prompts the user for information related to their personal use
# of the SimplifySoilDatabase.py script, and uses it to create a JSON
# formatted control file for application of the program.
#
# Modifications:
# March 17, 2016 - Modifications to make it more user friendly and interactive.  KAC

import json
import sys

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
    sys.stderr.write( "\nUsage: %s <filename>\n\n\t<filename> is the name of the control file to be created.\n" % sys.argv[0] )
    sys.exit()

# build control file
tmpDict = {}

# setup default values for HWSD
tmpDict["outCellNum"] = ( "Name of the output cell number file, should be VIC resolution", "No Default: Output Simulation Cells - Simulation Resolution", False )
tmpDict["inSoilFract"] = ( "Column header for fraction or probability of soil type from database", "SHARE", True )
tmpDict["outCellNumGrid"] = ( "The output (and output cell number file) from this program will be a grid file (True/False)?", "False", True )
tmpDict["inSoilTable"] = ( "Name of the CSV file derived from GIS database (exported from Access)", "SoilsData/HWSD_DATA.txt", True )
tmpDict["inCellNum"] = ( "Name of input with VIC model cell numbers, soil file grid resolution", "No Default: Input Simulation Cells - Soil Resolution", False )
tmpDict["inSoilGrid"] = ( "Name of input file with soil MU IDs", "No Default: Input Soil MUIDs - Soil Resolution", False )
tmpDict['inSoilCode'] = ( "Column heading containing identifier for specific soil type within a soil unit - SU_CODE74 used for US, SU_CODE90 used for EU", "SU_CODE74", True )
tmpDict["inSoilID"] = ( "Column header name for the soil ID used to link grid with database", "MU_GLOBAL", True )
tmpDict["outSoilNameFmt"] = ( "Output directory for completed soil files", "SoilsData/VICsetup/", True )
tmpDict['outParams'] = ( "Select Soil Database Type: HWSD", "HWSD", True )

outParams = {}
outParams["HWSD"] = ['T_TEXTURE', 'DRAINAGE', 'REF_DEPTH', 'AWC_CLASS', 'ROOTS', 'IL', 'SWR', 'ADD_PROP', 'T_GRAVEL', 'T_SAND', 'T_SILT', 'T_CLAY', 'T_USDA_TEX_CLASS', 'T_REF_BULK_DENSITY', 'T_BULK_DENSITY', 'T_OC', 'T_PH_H2O', 'T_CEC_CLAY', 'T_CEC_SOIL', 'T_BS', 'T_TEB', 'T_CACO3', 'T_CASO4', 'T_ESP', 'T_ECE', 'S_GRAVEL', 'S_SAND', 'S_SILT', 'S_CLAY', 'S_USDA_TEX_CLASS', 'S_REF_BULK_DENSITY', 'S_BULK_DENSITY', 'S_OC', 'S_PH_H2O', 'S_CEC_CLAY', 'S_CEC_SOIL', 'S_BS', 'S_TEB', 'S_CACO3', 'S_CASO4', 'S_ESP', 'S_ECE']

# prompt user for other options
print "==========\nThese first options MUST be set by the user:\n"

tmpDict["inSoilGrid"] = HandlePrompt( tmpDict["inSoilGrid"] )

tmpDict['inCellNum'] = HandlePrompt( tmpDict["inCellNum"] )

tmpDict['outCellNum'] = HandlePrompt( tmpDict["outCellNum"] )

tmpVar = HandlePrompt( tmpDict["outCellNumGrid"] )
if tmpVar.lower() == "false": tmpDict['outCellNumGrid'] = False
else: tmpDict['outCellNumGrid'] = True

tmpDict['inSoilCode'] = HandlePrompt( tmpDict["inSoilCode"] )

print "==========\nThe next options will have to be changed, if you did not put the files in the default directory\n"

tmpDict["inSoilTable"] = HandlePrompt( tmpDict["inSoilTable"] )

tmpDict['outSoilNameFmt'] = HandlePrompt( tmpDict["outSoilNameFmt"] )

if tmpDict['outCellNumGrid']:
    # gridded output filename format
    tmpDict['outSoilNameFmt'] = tmpDict['outSoilNameFmt'] + "SoilProperty_\%s.asc"
else:
    # output to single soil table file
    tmpDict['outSoilNameFmt'] = tmpDict['outSoilNameFmt'] + "SoilProperty_\%s.txt"

print "==========\nThese options should be correct, if using a pre-defined soil database.\n" 

tmpDict['inSoilID'] = HandlePrompt( tmpDict["inSoilID"] )

tmpDict['inSoilFract'] = HandlePrompt( tmpDict["inSoilFract"] )

# get soil headers to be processed
tmpDict['outParams'] = HandlePrompt( tmpDict["outParams"] )
if tmpDict['outParams'].upper() == "HWSD":
    print "Program setup to use the Harmonized World Soil Database"
    tmpDict['outParams'] = outParams["HWSD"]
else:
    print "Requested soil database, %s, is not currently available." % tmpDict['outParams'] 
    sys.exit()

# use json to write file
f = open( sys.argv[1], 'w' )
json.dump( tmpDict, f )
f.close()
