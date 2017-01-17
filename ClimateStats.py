#!/bin/env python
# Created on July 19, 2015
#  by Keith Cherkauer
#
# This script is designed to compute annual statistics for climate
# forcing data.  It will prodce files for each VIC simulation cell
# that include a time series of the annual metric of interest.  These
# are passed to a program to calculate the Mann-Kendall tau and slope
# and complete a statitsical test to see if there is a statistically 
# significant trend in the annual metric.
#
# INPUT: VIC model ASCII forcing file creted for the Tisza River project
# OUTPUT: Separate files for each metric calculated
#
# WARNING: I messed up the water year, which should start October 1, not
# September 1, so will have to look back at this in the future.

import sys, os
import pandas as pd
import numpy as np

# check for correct command line syntax
if len(sys.argv) != 3:
    sys.stderr.write( "Usage: %s <VIC forcing file list - with path> <output path>\n" % sys.argv[0] )
    sys.exit()

# handle command line arguments
inList = sys.argv[1]
outPath = sys.argv[2]

# read in contents of forcing file list - should include path as part of all file names
fileList = pd.read_table( inList, names=['fileName'] )

# check that output path exists, otherwise create it
if not os.path.exists(outPath):
    os.makedirs(outPath)

for FileID in fileList.index:

    inFile = fileList.get_value(FileID,'fileName')

    sys.stdout.write( "Processing %i of %i: %s\n" % ( FileID, len(fileList.index), inFile ) )

    inPath, fileName = os.path.split( inFile )

    # open and read forcing file
    rawDF = pd.read_table( inFile, skiprows=1, names=['YEAR','MONTH','DAY','PREC','TMAX','TMIN','WS2','PVAP','RH','PAIR','RG'], delim_whitespace=True, parse_dates=[[0,1,2]], index_col=[0] )

    # compute average daily air temperature
    rawDF['TAVG'] = ( rawDF['TMAX'] + rawDF['TMIN'] ) / 2.

    # calculate annual metrics and write to files

    # clip to annual period
    analysisDF = rawDF['1961-09-01':'2010-08-31']

    # annual PREC sum
    tmpSeries = analysisDF['PREC']
    tmpSeries = tmpSeries.resample('AS-SEP',how='sum')
    outFile = "%s/%s.PREC.Annual.Sum" % ( outPath, fileName )
    tmpSeries.to_csv(outFile,sep='\t',float_format='%.2f')

    ##########
    # PREC
    ##########

    # seasonal PREC sum
    tmpSeries = analysisDF['PREC']
    tmpSeries = tmpSeries.resample('QS-MAR',how='sum')
    outFile = "%s/%s.PREC.Autumn.Sum" % ( outPath, fileName )
    tmpSeries[range(0,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.PREC.Winter.Sum" % ( outPath, fileName )
    tmpSeries[range(1,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.PREC.Spring.Sum" % ( outPath, fileName )
    tmpSeries[range(2,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.PREC.Summer.Sum" % ( outPath, fileName )
    tmpSeries[range(3,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')

    #########
    # TMAX
    #########

    # annual TMAX mean
    tmpSeries = analysisDF['TMAX']
    tmpSeries = tmpSeries.resample('AS-SEP',how='mean')
    outFile = "%s/%s.TMAX.Annual.Mean" % ( outPath, fileName )
    tmpSeries.to_csv(outFile,sep='\t',float_format='%.2f')
    
    # seasonal TMAX mean
    tmpSeries = analysisDF['TMAX']
    tmpSeries = tmpSeries.resample('QS-MAR',how='mean')
    outFile = "%s/%s.TMAX.Autumn.Mean" % ( outPath, fileName )
    tmpSeries[range(0,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TMAX.Winter.Mean" % ( outPath, fileName )
    tmpSeries[range(1,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TMAX.Spring.Mean" % ( outPath, fileName )
    tmpSeries[range(2,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TMAX.Summer.Mean" % ( outPath, fileName )
    tmpSeries[range(3,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')

    # annual TMAX max
    tmpSeries = analysisDF['TMAX']
    tmpSeries = tmpSeries.resample('AS-SEP',how='max')
    outFile = "%s/%s.TMAX.Annual.Max" % ( outPath, fileName )
    tmpSeries.to_csv(outFile,sep='\t',float_format='%.2f')
    
    # seasonal TMAX max
    tmpSeries = analysisDF['TMAX']
    tmpSeries = tmpSeries.resample('QS-MAR',how='max')
    outFile = "%s/%s.TMAX.Autumn.Max" % ( outPath, fileName )
    tmpSeries[range(0,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TMAX.Winter.Max" % ( outPath, fileName )
    tmpSeries[range(1,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TMAX.Spring.Max" % ( outPath, fileName )
    tmpSeries[range(2,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TMAX.Summer.Max" % ( outPath, fileName )
    tmpSeries[range(3,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')

    # annual TMAX min
    tmpSeries = analysisDF['TMAX']
    tmpSeries = tmpSeries.resample('AS-SEP',how='min')
    outFile = "%s/%s.TMAX.Annual.Min" % ( outPath, fileName )
    tmpSeries.to_csv(outFile,sep='\t',float_format='%.2f')
    
    # seasonal TMAX min
    tmpSeries = analysisDF['TMAX']
    tmpSeries = tmpSeries.resample('QS-MAR',how='min')
    outFile = "%s/%s.TMAX.Autumn.Min" % ( outPath, fileName )
    tmpSeries[range(0,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TMAX.Winter.Min" % ( outPath, fileName )
    tmpSeries[range(1,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TMAX.Spring.Min" % ( outPath, fileName )
    tmpSeries[range(2,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TMAX.Summer.Min" % ( outPath, fileName )
    tmpSeries[range(3,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')

    #########
    # TMIN
    #########

    # annual TMIN mean
    tmpSeries = analysisDF['TMIN']
    tmpSeries = tmpSeries.resample('AS-SEP',how='mean')
    outFile = "%s/%s.TMIN.Annual.Mean" % ( outPath, fileName )
    tmpSeries.to_csv(outFile,sep='\t',float_format='%.2f')
    
    # seasonal TMIN mean
    tmpSeries = analysisDF['TMIN']
    tmpSeries = tmpSeries.resample('QS-MAR',how='mean')
    outFile = "%s/%s.TMIN.Autumn.Mean" % ( outPath, fileName )
    tmpSeries[range(0,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TMIN.Winter.Mean" % ( outPath, fileName )
    tmpSeries[range(1,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TMIN.Spring.Mean" % ( outPath, fileName )
    tmpSeries[range(2,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TMIN.Summer.Mean" % ( outPath, fileName )
    tmpSeries[range(3,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')

    # annual TMIN max
    tmpSeries = analysisDF['TMIN']
    tmpSeries = tmpSeries.resample('AS-SEP',how='max')
    outFile = "%s/%s.TMIN.Annual.Max" % ( outPath, fileName )
    tmpSeries.to_csv(outFile,sep='\t',float_format='%.2f')

    # seasonal TMIN max
    tmpSeries = analysisDF['TMIN']
    tmpSeries = tmpSeries.resample('QS-MAR',how='max')
    outFile = "%s/%s.TMIN.Autumn.Max" % ( outPath, fileName )
    tmpSeries[range(0,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TMIN.Winter.Max" % ( outPath, fileName )
    tmpSeries[range(1,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TMIN.Spring.Max" % ( outPath, fileName )
    tmpSeries[range(2,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TMIN.Summer.Max" % ( outPath, fileName )
    tmpSeries[range(3,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')

    # annual TMIN min
    tmpSeries = analysisDF['TMIN']
    tmpSeries = tmpSeries.resample('AS-SEP',how='min')
    outFile = "%s/%s.TMIN.Annual.Min" % ( outPath, fileName )
    tmpSeries.to_csv(outFile,sep='\t',float_format='%.2f')
    
    # seasonal TMIN min
    tmpSeries = analysisDF['TMIN']
    tmpSeries = tmpSeries.resample('QS-MAR',how='min')
    outFile = "%s/%s.TMIN.Autumn.Min" % ( outPath, fileName )
    tmpSeries[range(0,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TMIN.Winter.Min" % ( outPath, fileName )
    tmpSeries[range(1,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TMIN.Spring.Min" % ( outPath, fileName )
    tmpSeries[range(2,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TMIN.Summer.Min" % ( outPath, fileName )
    tmpSeries[range(3,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')

    #########
    # TAVG
    #########

    # annual TAVG mean
    tmpSeries = analysisDF['TAVG']
    tmpSeries = tmpSeries.resample('AS-SEP',how='mean')
    outFile = "%s/%s.TAVG.Annual.Mean" % ( outPath, fileName )
    tmpSeries.to_csv(outFile,sep='\t',float_format='%.2f')
    
    # seasonal TAVG mean
    tmpSeries = analysisDF['TAVG']
    tmpSeries = tmpSeries.resample('QS-MAR',how='mean')
    outFile = "%s/%s.TAVG.Autumn.Mean" % ( outPath, fileName )
    tmpSeries[range(0,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TAVG.Winter.Mean" % ( outPath, fileName )
    tmpSeries[range(1,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TAVG.Spring.Mean" % ( outPath, fileName )
    tmpSeries[range(2,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TAVG.Summer.Mean" % ( outPath, fileName )
    tmpSeries[range(3,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')

    # annual TAVG max
    tmpSeries = analysisDF['TAVG']
    tmpSeries = tmpSeries.resample('AS-SEP',how='max')
    outFile = "%s/%s.TAVG.Annual.Max" % ( outPath, fileName )
    tmpSeries.to_csv(outFile,sep='\t',float_format='%.2f')
    
    # seasonal TAVG max
    tmpSeries = analysisDF['TAVG']
    tmpSeries = tmpSeries.resample('QS-MAR',how='max')
    outFile = "%s/%s.TAVG.Autumn.Max" % ( outPath, fileName )
    tmpSeries[range(0,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TAVG.Winter.Max" % ( outPath, fileName )
    tmpSeries[range(1,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TAVG.Spring.Max" % ( outPath, fileName )
    tmpSeries[range(2,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TAVG.Summer.Max" % ( outPath, fileName )
    tmpSeries[range(3,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')

    # annual TAVG min
    tmpSeries = analysisDF['TAVG']
    tmpSeries = tmpSeries.resample('AS-SEP',how='min')
    outFile = "%s/%s.TAVG.Annual.Min" % ( outPath, fileName )
    tmpSeries.to_csv(outFile,sep='\t',float_format='%.2f')
    
    # seasonal TAVG min
    tmpSeries = analysisDF['TAVG']
    tmpSeries = tmpSeries.resample('QS-MAR',how='min')
    outFile = "%s/%s.TAVG.Autumn.Min" % ( outPath, fileName )
    tmpSeries[range(0,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TAVG.Winter.Min" % ( outPath, fileName )
    tmpSeries[range(1,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TAVG.Spring.Min" % ( outPath, fileName )
    tmpSeries[range(2,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')
    outFile = "%s/%s.TAVG.Summer.Min" % ( outPath, fileName )
    tmpSeries[range(3,len(tmpSeries),4)].to_csv(outFile,sep='\t',float_format='%.2f')




