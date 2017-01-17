#!/bin/env python
# Created on May 21, 2015
#  by Keith Cherkauer
#
# This script will convert the files extracted from the HWSD and reduced to the 
# majority soil unit per VIC model grid cell using the script SimplifySoilDatabase.py
# into a VIC model Soil Parameter file.  Missing soil paramteres will be computed 
# using the functions from Saxton and Rawls (2006) and compiled into the library
# SoilWaterEquations.py stored in /depot/phig/apps/source/python/lib.
# 

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
import read_arcinfo_files as rwarc

#####
# Setup soil column defaults
# - provide soil layer thicknesses in mm for each defined soil layer (works for HWSD where layers are consistently defined)
# - provide naming format so that each layer can be properly identified
# - provide naming format so that files for all variablecan be found, or created as necessary 
#####
SoilLayerDepths = [ 300, 700 ] # in millimeters
SoilLayerNameFmt = [ "T_%s", "S_%s" ]
SoilGridFmt = 'SoilsData/VICsetup/SoilProp_0.1deg_%s.asc'

#####
# Layers needed for calculations
# - Provide name of each layer as needed to complete SoilLayerNameFmt defined above
#####
SandFract = ( 'SAND', 100 ) # sand fraction of soil, 100 = in percent, 1 = as a fraction
ClayFract = ( 'CLAY', 100 ) # clay fraction of soil, 100 = in percent, 1 = as a fraction
OrganicFract = ( 'OC', 100 ) # organic matter fraction of soil, 100 = in percent, 1 = as a fraction
GravelFract = ( 'GRAVEL', 100 ) # gravel fraction of soil (optional, set to 'None' if not provided), 100 = in percent, 1 = as a fraction
BulkDensity = ( 'BULK_DENSITY', 1 ) # or should this be the REF_BULK_DENSITY column?
ImpermLayer = ( 'IL', 1 ) # not sure of the usefulness of this yet
CellSlopes = ( 'None', 1 ) # need to calculate in ArcGIS from projected DEM

#####
# Set file names for required information 
#####
CellNumGrid = 'ClimCellNum_GCS_AR.asc' # cell number grid file

#####
# Output filename
#####
SoilParamFile = 'SoilsData/SOIL_PARAM.asc'

#####
# Define functions
#####
def GetVariable( filename, value, colName, Ncells, INTflag='False', GRIDflag='True' ):
    '''This function will handle the given filename, either by opening and reading the
    contents of the file (if it exists), or creating a temporary array with the given 
    default value, if the filename is set to 'None'.  The column name of the reulting 
    Series is set to the value of colName.  If INTflag is set to tru the returned
    Series is set to an integer, otherwise it will be a floating point value.  A pandas
    Series is returned.'''

    if filename.capitalize() == 'None':
        # no filename has been provided, create variable Series using a constant default
        if value == np.NaN:
            # check that valid default value has been provided
            sys.stderr.write( "ERROR: No valid default value has been set for %s.  If this is the grid cell number file, then a file MUST be provided.\n" % ( filename ))
            sys.exit()
        # set dtype for current Series
        if INTflag == 'True': tmpType = 'i8'
        else: tmpType = 'f8'
        # create 
        tmpDF = pd.DataFrame( [value]*Ncells, dtype=tmpType, columns=[colName] )

    elif not os.path.isfile( filename ):
        # filename was provided, check if grid file exists
        sys.stderr.write( "ERROR: file for %s, %s, not found.\n" % ( colName, filename ))
        sys.exit()

    else:
        # valid filename was provided, use contents to create Series
        GridInfo = rwarc.read_ARCINFO_ASCII_grid( filename )
        if Ncells == np.NaN:
            # record number of cells in current grid extent
            Ncells = GridInfo['Ncells']
        else:
            # check that number of cells in file is the same as that defined
            if Ncells != GridInfo['Ncells']:
                sys.stderr.write( "ERROR: Found %i cells in the file %s, but expected %i.\n" % ( GridInfo['Ncells'], filename, Ncells ) )
                sys.exit()
        # Convert cell number file into pandas DataFrame for storing all of the variables
        tmpDF = pd.DataFrame( GridInfo['cell'] )
        # Convert value column to integers, if required
        if INTflag:
            tmpDF['value'] = tmpDF['value'].astype(int)
        # set name of column
        tmpDF = tmpDF.rename( columns={'value':colName} )

    return( tmpDF, Ncells )

#####
# Get Grid Cell Number
#####
( SoilParams, Ncells ) = GetVariables( CellNumGrid, np.NaN, 'gridcel', np.NaN, INTflag='True' )

#####
# Calculate missing parameters
#####

#####
# Finish building soil parameter table
#####
# set active cell grid file, set to 'None' if all grid cells are active
ActiveCellGrid = 'None'
# if no grid file is given, set all cells to 1
if ActiveCellGrid == 'None':
    # no selection of cells provided, set all to active
    SoilParams['run_cell'] = pd.Series(np.ones((Ncells,), dtype='i8' ), index=SoilParams.index )
elif not os.path.isfile( ActiveCellGrid ):
    # check if grid file exists
    sys.stderr.write( "ERROR: Active cell file, %s, not found.\n" % ( ActiveCellGrid ))
    sys.exit()
else:
    # read cell number grid into storage, since its header will be used for saving intermediate files 
    TmpGrid = rwarc.read_ARCINFO_ASCII_grid( ActiveCellGrid )
    # Convert cell number file into pandas DataFrame for storing all of the variables
    TmpDF = pd.DataFrame( TmpGrid['cell'] )
    # Convert value column to integer
    TmpDF['value'] = TmpDF['value'].astype(int)
    # set name of column
    TmpDF = TmpDF.rename( columns={'value':'run_cell'} )
    # add new column to soil param DataFrame
    SoilParams = pd.concat( [ SoilParams, TmpDF ] )
    # remove temporary files
    del TmpDF
    del TmpGrid

#####
# Get initial infiltion parameter values (typically a calibration parameter)
#####

# set infiltration parameter grid file, set to 'None' to use the starting default of 0.2
InfiltGrid = 'None'
# if no grid file is given, set all cells to 1
if InfiltGrid == 'None':
    # no cell file provided, set all cell to default value
    SoilParams['infilt'] = pd.Series( [0.2]*Ncells, index=SoilParams.index)
elif not os.path.isfile( InfiltGrid ):
    # check if grid file exists
    sys.stderr.write( "ERROR: Active cell file, %s, not found.\n" % ( InfiltGrid ))
    sys.exit()
else:
    # read cell number grid into storage, since its header will be used for saving intermediate files 
    TmpGrid = rwarc.read_ARCINFO_ASCII_grid( InfiltGrid )
    # Convert cell number file into pandas DataFrame for storing all of the variables
    TmpDF = pd.DataFrame( TmpGrid['cell'] )
    # set name of column
    TmpDF = TmpDF.rename( columns={'value':'infilt'} )
    # add new column to soil param DataFrame
    SoilParams = pd.concat( [ SoilParams, TmpDF ] )
    # remove temporary files
    del TmpDF
    del TmpGrid

#####
# Get initial Ds fraction for baseflow (typically a calibration parameter)
#####

# set infiltration parameter grid file, set to 'None' to use the starting default of 0.2
DsGrid = 'None'
# if no grid file is given, set all cells to 1
if DsGrid == 'None':
    # no cell file provided, set all cell to default value
    SoilParams['Ds'] = pd.Series( [0.001]*Ncells, index=SoilParams.index)
elif not os.path.isfile( DsGrid ):
    # check if grid file exists
    sys.stderr.write( "ERROR: Active cell file, %s, not found.\n" % ( DsGrid ))
    sys.exit()
else:
    # read cell number grid into storage, since its header will be used for saving intermediate files 
    TmpGrid = rwarc.read_ARCINFO_ASCII_grid( DsGrid )
    # Convert cell number file into pandas DataFrame for storing all of the variables
    TmpDF = pd.DataFrame( TmpGrid['cell'] )
    # set name of column
    TmpDF = TmpDF.rename( columns={'value':'Ds'} )
    # add new column to soil param DataFrame
    SoilParams = pd.concat( [ SoilParams, TmpDF ] )
    # remove temporary files
    del TmpDF
    del TmpGrid

#####
# Get initial Ws fraction for baseflow (typically a calibration parameter)
#####

# set Ws parameter grid file, set to 'None' to use the starting default of 0.9
WsGrid = 'None'
# if no grid file is given, set all cells to default
if WsGrid == 'None':
    # no cell file provided, set all cell to default value
    SoilParams['Ws'] = pd.Series( [0.9]*Ncells, index=SoilParams.index)
elif not os.path.isfile( WsGrid ):
    # check if grid file exists
    sys.stderr.write( "ERROR: Active cell file, %s, not found.\n" % ( WsGrid ))
    sys.exit()
else:
    # read cell number grid into storage, since its header will be used for saving intermediate files 
    TmpGrid = rwarc.read_ARCINFO_ASCII_grid( WsGrid )
    # Convert cell number file into pandas DataFrame for storing all of the variables
    TmpDF = pd.DataFrame( TmpGrid['cell'] )
    # set name of column
    TmpDF = TmpDF.rename( columns={'value':'Ws'} )
    # add new column to soil param DataFrame
    SoilParams = pd.concat( [ SoilParams, TmpDF ] )
    # remove temporary files
    del TmpDF
    del TmpGrid

#####
# Get initial baseflow curve shape epc fraction for baseflow (typically set to 2)
#####

# set baseflow shape exponent, set to 'None' to use the starting default of 2
cGrid = 'None'
# if no file is given, set all cells to default
if cGrid == 'None':
    # no cell file provided, set all cell to default value
    SoilParams['c'] = pd.Series( [2]*Ncells, index=SoilParams.index)
elif not os.path.isfile( cGrid ):
    # check if grid file exists
    sys.stderr.write( "ERROR: Active cell file, %s, not found.\n" % ( cGrid ))
    sys.exit()
else:
    # read cell number grid into storage, since its header will be used for saving intermediate files 
    TmpGrid = rwarc.read_ARCINFO_ASCII_grid( cGrid )
    # Convert cell number file into pandas DataFrame for storing all of the variables
    TmpDF = pd.DataFrame( TmpGrid['cell'] )
    # set name of column
    TmpDF = TmpDF.rename( columns={'value':'c'} )
    # add new column to soil param DataFrame
    SoilParams = pd.concat( [ SoilParams, TmpDF ] )
    # remove temporary files
    del TmpDF
    del TmpGrid

#####
# Get initial Dsmax (mm/day) for baseflow (typically a calibration parameter)
#####
'''The parameter Dsmax is the maximum velocity of baseflow for each grid cell. This can be estimated using the saturated hydraulic conductivity, Ksat, for each grid cell multiplied by the slope of the grid cell. The values for Ksat can be averaged for the layers for which baseflow will be included. When working in decimal degrees, the elevation data for the basin should be projected to an equal area map projection, in order to have horizontal dimensions in the same units as the vertical dimensions so that the slopes computed in Arc/Info are meaningful values.'''

# set infiltration parameter grid file, set to 'None' to use the starting default of 0.2
DsmaxGrid = 'None'
# if no grid file is given, set all cells to 1
if DsmaxGrid == 'None':
    # no cell file provided, set all cell to default value
    SoilParams['Dsmax'] = pd.Series( [10.0]*Ncells, index=SoilParams.index)
elif not os.path.isfile( DsmaxGrid ):
    # check if grid file exists
    sys.stderr.write( "ERROR: Active cell file, %s, not found.\n" % ( DsmaxGrid ))
    sys.exit()
else:
    # read cell number grid into storage, since its header will be used for saving intermediate files 
    TmpGrid = rwarc.read_ARCINFO_ASCII_grid( DsmaxGrid )
    # Convert cell number file into pandas DataFrame for storing all of the variables
    TmpDF = pd.DataFrame( TmpGrid['cell'] )
    # set name of column
    TmpDF = TmpDF.rename( columns={'value':'Dsmax'} )
    # add new column to soil param DataFrame
    SoilParams = pd.concat( [ SoilParams, TmpDF ] )
    # remove temporary files
    del TmpDF
    del TmpGrid


# expt
# Nlayer
# Exponent (=3+2/lambda
# Ksat - mm/day - Nlayer
# phi_s - mm/mm - Nlayer
# init_moist - mm - Nlayer
# elev - m
# depth - m - Nlayer
# avg_T - deg C
# dp - m
# bubble - cm - Nlayer
# quartz - fraction - Nlayer
# bulk_density - kg/m3 - Nlayer
# soil_density - kg/m3 - Nlayer
# organic - fraction - Nlayer
# bulk_dens_org - kg/m3 - Nlayer
# soil_dens_org - kg/m3 - Nlayer
# off_gmt - hours
# Wcr_FRACT - fraction - Nlayer
# Wpwp_FRACT - fraction - Nlayer
# rough - m
# snow_rough
# annual_prec
# resid_moist
# fs_active
# frost_slope
# max_snow_distrib_slope
# July_Tavg


SoilParams.to_csv( SoilParamFile, sep='\t', float_format='%g' )

