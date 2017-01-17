#!/bin/env python
# Created on May 27, 2015
#  by Keith Cherkauer
#
# This script reads the ASCII cell number grid file, and the CSV table developed from
# generating a zonal statistics table in ArcGIS using the high-resolution DEM slope values
# and a polygon layer with the VIC model grid cell outlines.  A polygon is used because
# the DEM and slope raster files are in an equal area projection, not used by the VIC model
# so the polygon file is best at preserving the geographic cell definitions. This was saved 
# from ArcGIS as a DBF database table.  The DBF file was read into Excel and exported as a
# CSV text file. The CSV file can now be read using the pandas read_table() function.  This
# program will replace values of cell number in the ASCII grid file with the mean slope (in 
# percent) for use in calculating Dsmax.  
#

import pandas as pd
import read_arcinfo_files as rwarc

outFile = 'SoilsData/VICsetup/SoilProp_0.1deg_SLOPE.asc'

# read in grid cell file
GridInfo = rwarc.read_ARCINFO_ASCII_grid( 'ClimCellNum_GCS_AR.asc' )
cellDF = pd.DataFrame( GridInfo['cell'] )
cellDF = cellDF.set_index('value')

# read zonal statistics table
ZonalStats = pd.read_csv( '../ArcHydroAnalysis/ZonalSlopeTable.csv', index_col=[0] ) 

# add MEAN slope to cell DataFrame, replacing the original 'value' column
cellDF['value'] = ZonalStats['MEAN']

# Change all NaNs to grid na data value before writing
cellDF = cellDF.where((pd.notnull(cellDF)), GridInfo['NODATA_value'])
 
# convert the cell DataFrame back into a dictionary structure
GridInfo['cell'] = cellDF.to_dict('records')

# write output to a new grid file
rwarc.write_ARCINFO_ASCII_grid( outFile, GridInfo, INTflag=0 )

