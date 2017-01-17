#!/usr/bin/python
# Created on November 11, 2010
#  by Keith Cherkauer
#
# This script uses a mask file to extract soil and vegetation parameters from
# the full LDAS data files.
# The original c-shell script was simply TOO slow, so I have rewritten it.
#

from read_VIC_files import read_soilparam_file, write_soilparam_file, read_vegparam_file, write_vegparam_file
from read_arcinfo_files import read_ARCINFO_ASCII_grid
import sys, os
from string import atof

#LdasSoilFile = "/home/pasture/b/cherkaue/Data/LDAS/UpdatedLdasSoilParams_20061103.txt"
#LdasVegFile  = "/home/pasture/b/cherkaue/Data/LDAS/ldas_lai.expanded.vegparams"
LdasSoilFile = "/depot/phig/data/VIC_SIMS_DATA/SoilData/LDAS/ldas_us.soil.orig.maod"
LdasVegFile  = "/depot/phig/data/VIC_SIMS_DATA/VegParams/ldas_lai.expanded.vegparams"

if ( len( sys.argv ) < 2 ):

    sys.stderr.write( "\nUsage: %s <mask file> <output dir> [-S <Soil Param File>] [-V <Veg Param File>]\n" % sys.argv[0] )
    sys.stderr.write( "\n\tNOTE: This program uses ASCII column format soil files for input and output.\n" )
    sys.stderr.write( "\tNOTE 2: The mask file can be a list of cell centers or an ArcInfo grid file.\n" )
    sys.stderr.write( "\tNOTE 3: Does not currently handle the snow band file.\n\n" )

else:

    # handle command line options
    MaskFile = sys.argv[1]
    del sys.argv[1]
    OutDir = sys.argv[1]
    del sys.argv[1]
    while ( len( sys.argv ) > 1 ):
        if ( sys.argv[1] == "-s" or sys.argv[1] == "-S" ):
            del sys.argv[1]
            LdasSoilFile = sys.argv[1]
        elif ( sys.argv[1] == "-v" or sys.argv[1] == "-V" ):
            del sys.argv[1]
            LdasVegFile = sys.argv[1]
        else:
            sys.stderr.write( "ERROR - No such option." )
            sys.stderr.write( "\nUsage: %s <mask file> <output dir> [-S <Soil Param File>] [-V <Veg Param File>]\n" % sys.argv[0] )
            sys.stderr.write( "\n\tNOTE: This program uses ASCII column format soil files for input and output.\n" )
            sys.stderr.write( "\tNOTE 2: The mask file can be a list of cell centers (<lat> <lng>) or an ArcInfo grid file.\n" )
            sys.stderr.write( "\tNOTE 3: Does not currently handle the snow band file.\n\n" )
            sys.exit()
        del sys.argv[1]

    # check for output directory and build output file names
    if not os.path.isdir( OutDir ):
        os.makedirs( OutDir )
    OutSoilFile = "%s/SoilParams" % OutDir
    OutVegFile  = "%s/VegParams" % OutDir

    #
    # process mask file
    #
    
    fin = open( MaskFile )
    line = fin.readline()
    fin.close()
    if ( "ncols" in line ):
        # process arcinfo ASCII grid mask file
        MaskInfoTable = read_ARCINFO_ASCII_grid(MaskFile, FLAG="FILTERED", INTflag=1)
        Ncells = MaskInfoTable["Ncells"]
        print "Number of Cells From Mask Grid = ", Ncells
        LatLngList = [0]*Ncells
        for cell in range( MaskInfoTable["Ncells"] ):
            LatLngList[cell] = "%.4f %.4f" % ( MaskInfoTable["cell"][cell]["lat"], MaskInfoTable["cell"][cell]["lng"] )
    else:
        # process list of grid cell centers
        fin = open( MaskFile )
        lines = fin.readlines()
        fin.close()
        Ncells = len( lines )
        print "Number of Cells from Mask File = ", Ncells
        LatLngList = [0]*Ncells
        for cell in range( MaskInfoTable["Ncells"] ):
            TmpVar = lines[cell].strip().split()
            LatLngList[cell] = "%.4f %.4f" % ( atof(TmpVar[0]), atof(TmpVar[1]) )
    print "Ncells from the mask file = %i" % Ncells

    #
    # process soil file
    #
    
    SoilTable = read_soilparam_file( LdasSoilFile, Nlayers=3, SPATIAL_SNOW="FALSE", SPATIAL_FROST="FALSE", EXCESS_ICE="FALSE", JULY_TAVG_SUPPLIED="FALSE" )
    Ncells = SoilTable["Ncells"]
    print "Number of Cells in Soil File = ", Ncells
    CellNumList = []
    DelCells = 0
    print "Ncells from the soil param file = %i" % Ncells
    for cell in range( Ncells-1, -1, -1 ):
        if not ( "%.4f %.4f" % ( SoilTable["cell"][cell]["lat"], SoilTable["cell"][cell]["lng"] ) in LatLngList ):
            # delete cells that do not appear in the mask
            del SoilTable["cell"][cell]
            SoilTable["Ncells"] = SoilTable["Ncells"] - 1
            DelCells = DelCells + 1
        else:
            CellNumList = CellNumList + [ SoilTable["cell"][cell]["CellNum"] ]
    write_soilparam_file( OutSoilFile, SoilTable, Nlayers=3, SPATIAL_SNOW="FALSE", SPATIAL_FROST="FALSE", EXCESS_ICE="FALSE", JULY_TAVG_SUPPLIED="FALSE" )
    print "Deleted %i cells, leaving %i." % ( DelCells, SoilTable["Ncells"] )

    #
    # process vegetation file
    #
    
    VegTable = read_vegparam_file( LdasVegFile, GLOBAL_LAI="TRUE" )
    Ncells = VegTable["Ncells"]
    print "Number of Cells from Vegetation File = ", Ncells
    DelCells = 0
    print "Ncells from the veg param file = %i" % Ncells
    for cellnum in VegTable["cell"].keys():
        if not ( cellnum in CellNumList ):
            # delete cells that do not appear in the mask
            del VegTable["cell"][cellnum]
            VegTable["Ncells"] = VegTable["Ncells"] - 1
            DelCells = DelCells + 1
    write_vegparam_file( OutVegFile, VegTable, GLOBAL_LAI="TRUE" )
    print "Deleted %i cells, leaving %i." % ( DelCells, VegTable["Ncells"] )


