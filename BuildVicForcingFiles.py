#!/bin/env python
# Created on May 12, 2015
#  by Keith Cherkauer
#
# This script reprocesses daily climate data for the Carpathain Basin, 
# from the CARPATCLIM-EU.org project, into VIC model input forcings that
# will be used to drive the model for the Tisza River project.
#
# INPUT: ASCII column files with YEAR, MONTH, and DAY in the first three 
# columns, followed by a single climate variable (see file name) for each
# of the 5000+ grid cells in the CARPATCLIM-EU domain.  
#
# OUTPUT: VIC model input climate files, with YEAR, MONTH and DAY plus
# columns for precipitation (PREC), maximum daily air temperature (TMAX), 
# minimum daily air temperature (TMIN), 2 meter wind speed (WS2), relative 
# humidity (RH), vapor pressure (PVAP), air pressure (PAIR), and global 
# radiation (RG).
#
# modifications:
# 20150707 Added a '#' to the start of the header line so that the VIC model
#          will skip it in the open_file() routine.   KAC  

import pandas as pd
import numpy as np
import read_arcinfo_files as rwarc

# set intput and output file for daily climate data
BaseInfoFile = "ClimateData/PredtandfilaGrid.dat"
DataFileNameFmt = "ClimateData/CARPATGRID_%s_D.ser"
OutFileNameFmt = "/scratch/conte/c/cherkaue/Workspace/TiszaRiver/ForcingData/data_%.4f_%.4f"

# set output files for required climate statistics
OutAnnualP = 'ClimateData/annual_prec.csv' # name for annual precip file
OutAnnualT = 'ClimateData/annual_Tavg.csv' # name of annual Tavg file
OutJulyT   = 'ClimateData/annual_July_Tavg.csv' # name of annual July_Tavg file

# set which climate variables will be output and in what order
DataTypes = [ "PREC", "TMAX", "TMIN", "WS2", "PVAP", "RH", "PAIR", "RG" ]

# read in base tables for all data files, this includes 'lat' and 'lon' of each cell center
print "Opening %s" % BaseInfoFile
basedata = pd.read_fwf( BaseInfoFile ) # read data from fixed width columns

# read in data variables, each file contains a climate variable for all cells and all days
climdata = {}
ColWidths = [ 4, 3, 3 ] + [ 8 ] * len(basedata) # auto width identification messes up the dates
for Type in DataTypes:
    TmpFileName = DataFileNameFmt % (Type)
    print "Opening %s" % TmpFileName
    climdata[Type] = pd.read_fwf( TmpFileName, parse_dates=[[0,1,2]], index_col=0, widths=ColWidths ) # read data from current variable
    climdata[Type].index.name = 'DATE'

# create annual statistics data frame
numCells = len(basedata)
AnnualStats = basedata[['lat','lon']].copy()
AnnualStats['annualprec'] = np.zeros( numCells )
AnnualStats['avg_T'] = np.zeros( numCells )
AnnualStats['July_Tavg'] = np.zeros( numCells )

# write VIC model forcing files
for GridIdx in range(numCells):
    TmpFileName = OutFileNameFmt % ( basedata.loc[GridIdx]['lat'], basedata.loc[GridIdx]['lon'] )
    print "Working on %s" % TmpFileName
    TmpList = [0] * len(DataTypes)
    for Idx in range(len(DataTypes)):
        TmpList[Idx] = climdata[DataTypes[Idx]]['%i' % (GridIdx+1)]
    dfTmp = pd.concat( TmpList, axis=1, keys=DataTypes )
    fout = open( TmpFileName, 'w' )
    fout.write("#") # adds a comment to the start of the header line, so that the VIC model will skip it
    dfTmp.to_csv( fout, sep='\t', date_format='%Y %m %d' )
    fout.close()

    # compute annual average precipitation
    tmpAnnual = dfTmp['PREC'].resample('A',how='sum')
    AnnualStats.set_value( GridIdx, 'annual_prec', np.mean(tmpAnnual) )

    # compute annual average air temperature
    tmpTavg = ( dfTmp['TMAX'] + dfTmp['TMIN'] ) / 2.
    tmpAnnual = tmpTavg.resample('A',how='mean')
    AnnualStats.set_value( GridIdx, 'avg_T', np.mean(tmpAnnual) )

    # compute annual average July air temperature
    tmpMonths = tmpTavg.resample('M',how='mean')
    listJuly = pd.date_range( '%i-07-31' % tmpMonths.index[0].year, '%i-07-31' % tmpMonths.index[-1].year, freq='12M' )
    tmpAnnual = tmpMonths[tmpMonths.index.isin(listJuly)]
    AnnualStats.set_value( GridIdx, 'July_Tavg', np.mean(tmpAnnual) )

# write annual statistics to XYZ files
AnnualStats[['lon','lat','annual_prec']].to_csv( OutAnnualP )
AnnualStats[['lon','lat','avg_T']].to_csv( OutAnnualT )
AnnualStats[['lon','lat','July_Tavg']].to_csv( OutJulyT )
