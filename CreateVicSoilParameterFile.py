#!/bin/env python
# Created on May 28, 2015
#  by Keith Cherkauer
#
# This program uses output from CalculatePerSoilDBLayer.py, which includes
# a table of soil parameter values for all VIC model cell coordinates and 
# a table that provides links between variable names and their specific
# values at different depths in the soil column.  The input file provides 
# soil inforamtion at the native vertical resolution of the soil database.
# This script is designed to generate a soil parameter file for the number
# and position of VIC cell soils.  Most often this program will generate a
# 3 layer file for a traditional VIC model application, or a 12 layer file
# for applications of VIC-CropSyst, but the number of soil layers is user 
# defined.
#
# Modifications:
# 20150707 Modified to convert run_cell and gridcel columns to type INT.
#          This eliminates problems reading the file into VIC. Also added 
#          step to drop rows/cells that are not going to be run.  This
#          cleans up the final file, if a mask has been provided for
#          the run_cell variable.      KAC

'''Definition of VIC model soil file columns (see http://www.hydro.washington.edu/Lettenmaier/Models/VIC/Documentation/SoilParam.shtml
run_cell	N/A	1	1 = Run Grid Cell, 0 = Do Not Run
gridcel	N/A	1	Grid cell number
lat	degrees	1	Latitude of grid cell
lon	degrees	1	Longitude of grid cell
infilt	N/A	1	Variable infiltration curve parameter (binfilt)
Ds	fraction	1	Fraction of Dsmax where non-linear baseflow begins
Dsmax	mm/day	1	Maximum velocity of baseflow
Ws	fraction	1	Fraction of maximum soil moisture where non-linear baseflow occurs
c	N/A	1	Exponent used in baseflow curve, normally set to 2
expt	N/A	Nlayer	Exponent n (=3+2/lambda) in Campbell's eqn for hydraulic conductivity, HBH 5.6 (where lambda = soil pore size distribution parameter). Values should be > 3.0.			
Ksat	mm/day	Nlayer	Saturated hydrologic conductivity
phi_s	mm/mm	Nlayer	Soil moisture diffusion parameter
init_moist	mm	Nlayer	Initial layer moisture content
elev	m	1	Average elevation of grid cell
depth	m	Nlayer	Thickness of each soil moisture layer
avg_T	C	1	Average soil temperature, used as the bottom boundary for soil heat flux solutions
dp	m	1	Soil thermal damping depth (depth at which soil temperature remains constant through the year, ~4 m)
bubble	cm	Nlayer	Bubbling pressure of soil. Values should be > 0.0
quartz	fraction	Nlayer	Quartz content of soil
bulk_density	kg/m3	Nlayer	Bulk density of soil layer
soil_density	kg/m3	Nlayer	Soil particle density, normally 2685 kg/m3
organic	fraction	Nlayer	Fraction of soil layer that is organic
bulk_dens_org	kg/m3	Nlayer	Bulk density of organic portion of soil
soil_dens_org	kg/m3	Nlayer	Soil particle density of organic portion of soil, normally 1300 kg/m3
off_gmt	hours	1	Time zone offset from GMT. This parameter determines how VIC interprets sub-daily time steps relative to the model start date and time.
Wcr_FRACT	fraction	Nlayer	Fractional soil moisture content at the critical point (~70% of field capacity) (fraction of maximum moisture)
Wpwp_FRACT	fraction	Nlayer	Fractional soil moisture content at the wilting point (fraction of maximum moisture)
rough	m	1	Surface roughness of bare soil
snow_rough	m	1	Surface roughness of snowpack
annual_prec	mm	1	Average annual precipitation.
resid_moist	fraction	Nlayer	Soil moisture layer residual moisture.
fs_active	1 or 0	1	If set to 1, then frozen soil algorithm is activated for the grid cell. A 0 indicates that frozen soils are not computed even if soil temperatures fall below 0C.
frost_slope	C	1	Slope of uniform distribution of soil temperature (if SPATIAL_FROST == TRUE in the global parameter file).
max_snow_distrib_slope	m	1	Maximum slope of the snow depth distribution. This is only used if SPATIAL_SNOW == TRUE in the global parameter file.
July_Tavg	C	1	Average July air temperature.
'''

import sys, os
import pandas as pd
import numpy as np

#sys.path.append("./Scripts")

from VicPreprocessingLibrary import *

#########
# USER changes start here
#########

#########
# Define files
#########
SoilDataTable = 'SoilParamBuffalo.txt'
SoilAttribTable = 'SoilParamBuffalo_attrib.txt'
outSoilParamFile = 'VicSetup/SOIL_PARAM.txt'

#########
# Define VIC model soil setup
#########
ModelLayers = [ 0.1, 0.2, 0.7 ] # default layer thicknesses
ModelOptions = [ 'All', 'EB', 'FS', 'Tree' ] # see below for 'options' description
Nlayer = len(ModelLayers) # number of layers
NoData = -999

#########
# Define soil database setup
#########
numSL = 2 # number of soil layers in the soil database
SoilLayerNameFmt = [ "T_%s", "S_%s" ]

#####
# Define names and order of appearance for all soil variables
# - This section should only be changed if the order or number of columns in the VIC
#   model soil parameter file has been changed
# - If the script CalculatePerSoilDBLayer.py or its input file have been changed such
#   that the attribute names are no longer the same as those provided in the dictionary 
#   defined below with the key 'sourceCol'.
# In the dictionary, keys are defined as: 
#     'name' = output soil parameter column name, should include "_%i" for layer 
#              specific values
#     'layers' = 0 (False) not provided for layers, 1 (True) provided for layers 
#     'options' = 'All' for all cases, 'EB' for energy balance, 'FS' for frozen soil, 
#                 'ORG' for organic flag set, 'SpFr' for spatial frost, 'SpSn' for
#                 spatial snow, and 'TL' for compute treeline.
#     'sourceCol' = name of source column in input data file, or np.NaN if not available
#                   in the source file.
#     'default' = a default value if one is suggested, or NaN for no default.
#####

DepthFmt = 'depth_%i' # name format for soil depth values
IMoistFmt = 'init_moist_%i' # name format for initial soil moisture values
SoilParamCols = [ {'name': 'run_cell',               'layers':0, 'options':'All',  'sourceCol':'run_cell',      'default':1 },
                  {'name': 'gridcel',                'layers':0, 'options':'All',  'sourceCol':'gridcel',       'default':False }, 
                  {'name': 'lat',                    'layers':0, 'options':'All',  'sourceCol':'lat',           'default':False }, 
                  {'name': 'lon',                    'layers':0, 'options':'All',  'sourceCol':'lng',           'default':False }, 
                  {'name': 'infilt',                 'layers':0, 'options':'All',  'sourceCol':False,           'default':0.001 }, 
                  {'name': 'Ds',                     'layers':0, 'options':'All',  'sourceCol':False,           'default':0.1 }, 
                  {'name': 'Dsmax',                  'layers':0, 'options':'All',  'sourceCol':'Dsmax',         'default':False }, 
                  {'name': 'Ws',                     'layers':0, 'options':'All',  'sourceCol':False,           'default':0.97 }, 
                  {'name': 'c',                      'layers':0, 'options':'All',  'sourceCol':False,           'default':2 }, 
                  {'name': 'expt_%i',                'layers':1, 'options':'All',  'sourceCol':'Expt',          'default':False }, 
                  {'name': 'Ksat_%i',                'layers':1, 'options':'All',  'sourceCol':'Ksat',          'default':False }, 
                  {'name': 'phi_s_%i',               'layers':1, 'options':'All',  'sourceCol':'phi_s',         'default':NoData }, 
                  {'name': IMoistFmt,                'layers':1, 'options':'All',  'sourceCol':False,           'default':False }, 
                  {'name': 'elev',                   'layers':0, 'options':'All',  'sourceCol':'elev',          'default':False }, 
                  {'name': DepthFmt,                 'layers':1, 'options':'All',  'sourceCol':False,           'default':False }, 
                  {'name': 'avg_T',                  'layers':0, 'options':'EB',   'sourceCol':'avg_T',         'default':False }, 
                  {'name': 'dp',                     'layers':0, 'options':'EB',   'sourceCol':False,           'default':4.0 }, 
                  {'name': 'bubble_%i',              'layers':1, 'options':'FS',   'sourceCol':'Bubble',        'default':False }, 
                  {'name': 'quartz_%i',              'layers':1, 'options':'EB',   'sourceCol':'quartz',        'default':False }, 
                  {'name': 'bulk_density_%i',        'layers':1, 'options':'All',  'sourceCol':'BULK_DENSITY',  'default':False }, 
                  {'name': 'soil_density_%i',        'layers':1, 'options':'All',  'sourceCol':'soil_density',  'default':False }, 
                  {'name': 'organic_%i',             'layers':1, 'options':'ORG',  'sourceCol':'OC',            'default':False }, 
                  {'name': 'bulk_den_org_%i',        'layers':1, 'options':'ORG',  'sourceCol':False,           'default':False }, 
                  {'name': 'soil_dens_org_%i',       'layers':1, 'options':'ORG',  'sourceCol':'soil_dens_org', 'default':False }, 
                  {'name': 'off_gmt',                'layers':0, 'options':'All',  'sourceCol':'off_gmt',       'default':False }, 
                  {'name': 'Wcr_FRACT_%i',           'layers':1, 'options':'All',  'sourceCol':'Wcr_FRACT',     'default':False }, 
                  {'name': 'Wpwp_FRACT_%i',          'layers':1, 'options':'All',  'sourceCol':'Wpwp_FRACT',    'default':False }, 
                  {'name': 'rough',                  'layers':0, 'options':'All',  'sourceCol':False,           'default':0.001 }, 
                  {'name': 'snow_rough',             'layers':0, 'options':'All',  'sourceCol':False,           'default':0.0005 }, 
                  {'name': 'annual_prec',            'layers':0, 'options':'All',  'sourceCol':'annual_prec',   'default':False }, 
                  {'name': 'resid_moist_%i',         'layers':1, 'options':'All',  'sourceCol':'resid_moist',   'default':False }, 
                  {'name': 'fs_active',              'layers':0, 'options':'FS',   'sourceCol':False,           'default':1 }, 
                  {'name': 'frost_slope',            'layers':0, 'options':'SpFr', 'sourceCol':False,           'default':False }, 
                  {'name': 'max_snow_distrib_slope', 'layers':0, 'options':'SpSn', 'sourceCol':False,           'default':False }, 
                  {'name': 'July_Tavg',              'layers':0, 'options':'TL',   'sourceCol':'July_Tavg',     'default':False } ]

########################################################
# NO USER MODIFICATIONS SHOULD BE NECESSARY BELOW HERE #
########################################################

##########
# Create file header
##########
outHeader = []
for idx in range(len(SoilParamCols)):
    if SoilParamCols[idx]['options'] in ModelOptions:
        if SoilParamCols[idx]['layers']:
            for lyr in range(Nlayer):
                outHeader = outHeader + [ SoilParamCols[idx]['name'] % (lyr+1) ]
        else:
            outHeader = outHeader + [ SoilParamCols[idx]['name'] ]

##########
# Read in data table
##########
SoilData = pd.read_csv( SoilDataTable )
tmpInfo = pd.read_csv( SoilAttribTable, index_col=0 )
SoilInfo = tmpInfo.to_dict( 'records' )

Ncells = len(SoilData)

##########
# Create empty DataFrame for soil parameters
##########
SoilParamDF = pd.DataFrame( np.full([Ncells,len(outHeader)],np.nan), index=range(Ncells), columns=outHeader )

#########
# Map soil database layers to VIC model soils layers
# - This will use a nearest neighbor type approach, where the depth of the middle of the 
#   VIC soil layer will be used to determine which soil layer from the database will be 
#   used to set the VIC model soil layer properties.  I have decided that doing an 
#   interpolation between defined soil layers from the database is hard to justify, since
#   average soil perperties are not necessarily what one would find between the soil layers.
#   The relationship is set first assuming all soil layers are the same, then adjusted
#   if it is determined that the soil profile in the database was restrictd by depth.  In 
#   such a case the thickness of soil layers below the surface layer will be reduced 
#   in proportion to the reduction in total soil column depth.
##########
LayerLinkDF = pd.DataFrame( np.zeros((Ncells,Nlayer)), index=range(Ncells), columns=range(Nlayer), dtype='i8' )
ModelLayers = np.array( ModelLayers )
MaxVicDepth = ModelLayers.sum()
print SoilData.keys()
print SoilData['MaxDepth']
MaxDepth = SoilData['MaxDepth'].max()
# find cells with soils depths less than maximum
for cel in range(Ncells):
    if SoilData.get_value(cel,'MaxDepth') != MaxDepth:
        # adjust VIC model soil layer depths for shallow soil
        reducRatio = ( MaxDepth - ModelLayers[0] ) / ( MaxVicDepth - ModelLayers[0] )
    else: reducRatio = 1.
    # calculate depth of VIC layer center
    tmpModelLayers = np.append( ModelLayers[0], ModelLayers[1:]*reducRatio )
    for lyr in range(Nlayer):
        Depth = tmpModelLayers[:(lyr+1)].sum() - tmpModelLayers[lyr]/2. # calculate depth of layer center
        # find corresponding layer in soil database
        if Depth > MaxDepth:
            # this database layer contains VIC soil layer for this grid cell
            LayerLinkDF.set_value(cel,lyr,numSL-1)
        else:
            tmpDepth = 0
            for SL in range(numSL):
                tmpDepth = tmpDepth + SoilData.get_value(cel,SoilInfo[SL]['SoilThick'])
                if Depth <= tmpDepth:
                    # this database layer contains VIC soil layer for this grid cell
                    LayerLinkDF.set_value(cel,lyr,SL)
                    break
	# set VIC cell depth (m)
        SoilParamDF.set_value(cel,DepthFmt % (lyr+1),tmpModelLayers[lyr])
	# set VIC cell initial moisture (mm)
        SoilParamDF.set_value(cel,IMoistFmt % (lyr+1),tmpModelLayers[lyr]*0.7*SoilData.get_value(cel,SoilInfo[LayerLinkDF.get_value(cel,lyr)]['Theta_S'])*1000.)

##########
# Assign layer properties
##########
outCol = 0
for ColInfo in SoilParamCols:
    if ColInfo['options'] in ModelOptions:
        print "Processing %s." % ColInfo['name']
        # selected variable is used with current VIC model option set
        if ColInfo['layers']:
            print "This is a column value"
            # current data column is defined for all model soil layers
            print outHeader[outCol], 'depth' in outHeader[outCol] or 'init_moist' in outHeader[outCol]
            if not ( 'depth' in outHeader[outCol] or 'init_moist' in outHeader[outCol] ):
                print "Saving layer values"
                # setting depth and initial moisture content are special cases, dealt with earlier in the code
                if ColInfo['sourceCol']:
                    print "... from table"
                    # data value is defined in the soil table
                    for lyr in range(Nlayer):
                        for cel in range(Ncells):
                            # get soil database layer
                            SL = int(LayerLinkDF.get_value(cel,lyr))
                            # set output variable from default value
                            SoilParamDF.set_value(cel,outHeader[outCol],SoilData.get_value(cel,SoilInfo[SL][ColInfo['sourceCol']])) 
                        outCol  = outCol + 1
                elif ColInfo['default']:
                    print "... from default"
                    # data value is defined with a default value
                    for lyr in range(Nlayer):
                        for cel in range(Ncells):
                            # get soil database layer
                            SL = int(LayerLinkDF.get_value(cel,lyr))
                            # set output variable from default value
                            SoilParamDF.set_value(cel,outHeader[outCol],ColInfo['default']) 
                        outCol  = outCol + 1
                else:
                    # nither defined - this is an error
                    sys.stderr.write( "ERROR: Current column, %s, has not been assigned a value for current soil layer.\n" % ColInfo['name'] )
                    sys.exit()
            else:
                outCol = outCol + Nlayer # skip over columns for depth and init_moist
        else:
            # current data column is not defined for soil layers
            if ColInfo['sourceCol']:
                # data value is defined in the soil table
                if outHeader[outCol] == 'Dsmax':
                    # Dsmax is handled differently, since it should be the average of the soil layers contributing
                    # Don't know how to handle this case yet, so currently just use the bottom soil value
                    SoilParamDF[outHeader[outCol]] = SoilData[SoilInfo[numSL-1][ColInfo['sourceCol']]]
                else:
                    SoilParamDF[outHeader[outCol]] = SoilData[ColInfo['sourceCol']]
            elif ColInfo['default']:
                # data value is defined with a default value
                SoilParamDF[outHeader[outCol]] = np.full([Ncells],ColInfo['default']) # fill array with default value
            else:
                # nither defined - this is an error
                sys.stderr.write( "ERROR: Current column, %s, has not been assigned a value.\n" % ColInfo['name'] )
                sys.exit()
            outCol = outCol + 1
    else: 
        print "The variable %s not included as option %s not used for current simulation." % (  ColInfo['name'], ColInfo['options'] )

##########
# write soil parameter file
##########

print SoilParamDF.describe()
# replace NaNs with 0 for mask column
MaskCol = SoilParamCols[0]['name']
SoilParamDF[MaskCol] = SoilParamDF[MaskCol].where((pd.notnull(SoilParamDF[MaskCol])), 0)
# drop rows where there is no cell number defined
SoilParamDF = SoilParamDF[pd.notnull(SoilParamDF[SoilParamCols[1]['name']])]
# replace all other NaNs with NoData value
SoilParamDF = SoilParamDF.where((pd.notnull(SoilParamDF)), NoData)
# drop cells that are not set to be run
SoilParamDF = SoilParamDF[SoilParamDF.run_cell != 0]
# convert columns to type INT if required
for ColInfo in SoilParamCols:
    SoilParamDF['run_cell'] = SoilParamDF['run_cell'].astype(int)
    SoilParamDF['gridcel'] = SoilParamDF['gridcel'].astype(int)
    if 'fs_active' in SoilParamDF.columns:
	    SoilParamDF['fs_active'] = SoilParamDF['fs_active'].astype(int)
# write table to a file
fout = open( outSoilParamFile, 'w' )
fout.write("#")
SoilParamDF.to_csv( fout, sep='\t', index=False )
fout.close()
