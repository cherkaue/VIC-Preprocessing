#!/bin/env python
# Created on July 20, 2015
#  by Keith Cherkauer
#
# This script takes the contents of a Zonal Stats table and produces
# the grid and fracton files required for using the GIS routing model.
#

import os, sys
import pandas as pd
import numpy as np
import read_arcinfo_files as rwarc

# check for command line arguments
if len(sys.argv) != 5:
    sys.stderr.write( "Usage: %s <basin prefix> <zonal CSV stat table> <out path> <cell number grid>\n" % sys.argv[0] )
    sys.exit()

# handle command line arguments
print sys.argv
basinPrefix = sys.argv[1] # 'polgar'
zonalTable = sys.argv[2] # 'RoutingLayers/ZonalTables/zonalstats_polgar.csv'
outPath = sys.argv[3] # '../VicSetup/ROUTING_MODEL/'
CellNumGrid = sys.argv[4] # '../BasinAnalysis/ClimCellNum_GCS_AR.asc'

# read zonal statistics table
rawDF = pd.read_csv( zonalTable, index_col=0 )

# read cell number grid
GridInfo = rwarc.read_ARCINFO_ASCII_grid( CellNumGrid )
if not os.path.isfile( CellNumGrid ):
    # no valid grid cell number file was found
    sys.stderr.write( "Grid cell number file %s does not exist, program cannot function without a valid grid cell file.\n" % ( CellNumGrid ))
    sys.exit()

# valid filename was provided, use contents to create DataFrame from ArcInfo raster
GridInfo = rwarc.read_ARCINFO_ASCII_grid( CellNumGrid )
# record number of cells in current grid extent
Ncells = GridInfo['Ncells']
# Convert cell number file into pandas DataFrame for storing all of the variablesc
cellDF = pd.DataFrame( GridInfo['cell'] )
# Convert value column to integer
cellDF['value'] = cellDF['value'].astype(int)
# set name of column
cellDF = cellDF.rename( columns={'value':'gridcel'} )
# drop rows without a cellnum
cellDF = cellDF[cellDF['gridcel']!=GridInfo['NODATA_value']]
# save cell num list
CellNumList = cellDF.gridcel
# set cell numbers to be DF indices
cellDF.set_index( 'gridcel', inplace=True )

# merge zonal statistics with cell number grid
cellDF['ZonalMean'] = rawDF['MEAN']
cellDF['ZonalStd'] = rawDF['STD']
print 'Quantile', rawDF['COUNT'].quantile(0.25)
cellDF['ZonalFract'] = rawDF['COUNT'] / rawDF['COUNT'].quantile(0.25)
cellDF['ZonalFract'] = cellDF['ZonalFract'].clip(upper=1.0)
cellDF['ZonalMask'] = cellDF['ZonalFract'].clip(lower=1.0)
cellDF['RoutNum'] = cellDF.index
cellDF['VicNum'] = cellDF.index

# write routing files
strFmt = lambda x, y: "fluxes_%.4f_%.4f" % ( x, y )
cellDF['FluxName'] = np.array(['fluxes_']*len(cellDF.index))
for pos in cellDF.index:
    cellDF.set_value(pos,'FluxName',strFmt( cellDF['lat'][pos], cellDF['lng'][pos] ))

# write fraction file
#  Needs columns with <routing model ID> <vic model ID> <flux file name> <fraction>
outFileName = "%s/%s_fract.txt" % ( outPath, basinPrefix )
print "Working on %s" % outFileName
outDF = cellDF[['RoutNum','VicNum','FluxName','ZonalFract']]
outDF = outDF.dropna(how='any')
outDF.to_csv( outFileName, sep='\t', header=False, index=False)
print "Wrote %i lines." % len(outDF.index)

# fill missing values with the No Data value from the grid cell number grid
cellDF.fillna(GridInfo['NODATA_value'],inplace=True)

# mean travel time
outFileName = "%s/%s_t0.asc" % ( outPath, basinPrefix )
print "Working on %s" % outFileName
for idx in CellNumList.index:
    pos = CellNumList.loc[idx]
    if pos in cellDF.index:
        GridInfo["cell"][idx]["value"] = cellDF.at[pos,'ZonalMean']
    else:
        GridInfo["cell"][idx]["value"] = GridInfo['NODATA_value']
rwarc.write_ARCINFO_ASCII_grid(outFileName,GridInfo,INTflag=0)

# std dev travel time
outFileName = "%s/%s_delta.asc" % ( outPath, basinPrefix )
print "Working on %s" % outFileName
for idx in CellNumList.index:
    pos = CellNumList.loc[idx]
    if pos in cellDF.index:
        GridInfo["cell"][idx]["value"] = cellDF.at[pos,'ZonalStd']
    else:
        GridInfo["cell"][idx]["value"] = GridInfo['NODATA_value']
rwarc.write_ARCINFO_ASCII_grid(outFileName,GridInfo,INTflag=0)

# write the mask grid
outFileName = "%s/%s_mask.asc" % ( outPath, basinPrefix )
print "Working on %s" % outFileName
for idx in CellNumList.index:
    pos = CellNumList.loc[idx]
    if pos in cellDF.index:
        GridInfo["cell"][idx]["value"] = cellDF.at[pos,'ZonalMask']
    else:
        GridInfo["cell"][idx]["value"] = GridInfo['NODATA_value']
rwarc.write_ARCINFO_ASCII_grid(outFileName,GridInfo,INTflag=1)

