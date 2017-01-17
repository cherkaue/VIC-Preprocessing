# Created on 2015-15-5
#  by Keith Cherkauer
#
# This script file was written to convert GeoTIFF maps of fractional
# grid cell land use cover into ASCII grids for further analysis outside
# of ArcGIS.
#

import os
import arcpy as ap
import numpy as np
import sys

# set file names and other constants
Workspace = 'Z:\data\Projects\TiszaRiver\ArcHydroAnalysis'
LandUseIDs = 'IGBP_Land_Use_Classes.csv'
OutFilePath = "LandUseData"
InFileFmt = "LandUseFract_MCD12Q1_2005_IGBP_%03i_GCS_AR.tif"
OutFileFmt = "LandUseFract_MCD12Q1_2005_IGBP_%03i_GCS_AR.asc"

# setup working environment
ap.env.workspace = Workspace

#####
# Convert each file to ASCII and save
#####

# read in land use class database
ClassList = np.genfromtxt(LandUseIDs, delimiter=',', names=True)
# for each vegetation class, calculate the fraction area
for Idx in range(len(ClassList)):
    # create filenames
    InFileName = InFileFmt % Idx
    OutFileName = OutFileFmt % Idx
    print "Converting %s to %s" % ( InFileName, OutFileName )
    # check to see if outfile is present, if so delete it
    if os.path.isfile( "%s/%s/%s" % ( Workspace, OutFilePath, OutFileName ) ):
        print "... Deleting exsting output file."
        ap.Delete_management( "%s/%s/%s" % ( Workspace, OutFilePath, OutFileName ))
    # convert the raster file to an ASCII file
    ap.RasterToASCII_conversion("%s\%s\%s" % ( Workspace, OutFilePath, InFileName),
                                "%s\%s\%s" % ( Workspace, OutFilePath, OutFileName) )



