# Created on 2015-14-5
#  by Keith Cherkauer
#
# This script file was written to compute the fraction of each IGBP
# vegetation type in a VIC model grid cell using MODIS MCD12Q1 land cover
# and the 0.1 degree resolution CARPATGRID used by the climte dataset.
#

import os.path
import arcpy as ap
import arcpy.sa as sa
import numpy as np

# set file names and other constants
Workspace = 'Z:\data\Projects\TiszaRiver\ArcHydroAnalysis'
LandUseMap = 'LandUseData/IGBPLandUse2005_GCS_AR.tif'
LandUseIDs = 'IGBP_Land_Use_Classes.csv'
inZoneMap = 'TiszaAnalysis.gdb/ZonalPolygons'
inZoneID = "gridcode"
outCellSize = 0.1
OutFilePath = "LandUseData"
OutFileFmt = "LandUseFract_MCD12Q1_2005_IGBP_%03i_GCS_AR.tif"

# setup working environment
ap.env.workspace = Workspace
ap.env.extent = LandUseMap
ap.env.cellsize = LandUseMap

# Check out the ArcGIS Spatial Analyst extension license
ap.CheckOutExtension("Spatial")

# check for output geodatabase, create if not present
#if not os.path.isfile("%s/%s" % (Workspace,OutGeoDB)):
#    print "Making geodatabase"
#    ap.CreateFileGDB_management(Workspace, OutGeoDB)

#####
# Calculate number of land use grid cells per model grid cell for use
# in calculating fractional coverages 
#####

print "Calculating number of land use cells in model grid cells for fracional calculations"
# replace land use cells with "1"s so that zonal sum does the math for us 
outCon = sa.Con(LandUseMap, 1, "", "VALUE >= 0")
# Execute ZonalStatistics
sumCellsPerZone = sa.ZonalStatistics(inZoneMap, inZoneID, outCon, "SUM", "NODATA")

#####
# Estimate fractional  coverage of all vegetation types, and add resulting grid
# to the Land USe Geodatabase
#####

# read in land use class database
ClassList = np.genfromtxt(LandUseIDs, delimiter=',', names=True)
# for each vegetation class, calculate the fraction area
for Idx in range(len(ClassList)):
    OutFileName = OutFileFmt % Idx
    print "Working on %s" % OutFileName
    # filter land use to current class
    outCon = sa.Con(LandUseMap, 1, "", "VALUE = %i" % Idx)
    # Execute ZonalStatistics
    vegCellsPerZone = sa.ZonalStatistics(inZoneMap, inZoneID, outCon, "SUM", "DATA")
    # calculate fraction of grid cell covered by current land use
    outFract = vegCellsPerZone / sumCellsPerZone
    # resample grid from MODIS to the VIC model resolution
    if os.path.isfile( "%s/%s/%s" % ( Workspace, OutFilePath, OutFileName ) ):
        print "... Deleting exsting file."
        ap.Delete_management( "%s/%s/%s" % ( Workspace, OutFilePath, OutFileName ))
    # resample grid to resolution used by simulation
    ap.Resample_management(outFract, "%s/%s" % ( OutFilePath, OutFileName ), outCellSize, "NEAREST")
    # Save the result - Now completed in previous step
    #outFract.save("%s/%s/%s" % ( Workspace, OutFilePath, OutFileName ))


