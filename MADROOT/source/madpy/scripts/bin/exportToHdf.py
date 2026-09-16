#!PYTHONEXE

#$Id: exportToHdf.py 3348 2011-04-20 20:17:28Z murco $

usage = """
exportToHdf.py is a script used to convert a Cedar file to hdf5 format.

Required arguments:

    --cedarFilename - full path of existing Madrigal file. 

    --hdf5Filename - full path of hdf5 file to write.

Optional arguments - set these to add layouts, parameters or filters:
    
    --independentSpatialParms - a comma separated list of parameters as mnemonics
        that represent independent spatial variables.  Causes array layout to be added to 
        output Hdf5 file

    --arraySplittingParms - a comma separated list of parameters as mnemonics used to split
        arrays into subarrays.  For example, beamcode would split data with separate beamcodes
        into separate arrays. The number of separate arrays will be up to the product of the number of 
        unique values found for each parameter, with the restirction that combinations with no records will
        not create a separate array.
    
    --extraParameters - These parameters will be added to the output file if 
                        they are not already in the input file. Comma-delimited.
                        
    --filter - Filter argument as in isprint command as string (eg, 'ti.500,2000') Only one allowed.

Example:
    exportToHdf --cedarFilename=/opt/madrigal/experiments/1998/mlh/20jan98/mil20100112.001
                --hdf5Filename=/home/user/data/mil20100112.hdf5
                --independentSpatialParms=range
                --arraySplittingParms=kinst,pl,mdtyp
                --extraParameters=ti,te
                --filter=ti,500,1000
"""

import sys
import os, os.path
import getopt
import traceback
import madrigal.data

# parse command line
arglist = ''
longarglist = ['cedarFilename=',
               'hdf5Filename=',
               'independentSpatialParms=',
               'arraySplittingParms=',
               'extraParameters=',
               'filter=']

optlist, args = getopt.getopt(sys.argv[1:], arglist, longarglist)


# set default values
cedarFilename = None
hdf5Filename = None
independentSpatialParms = []
arraySplittingParms = []
extraParameters = []
filter = None

for opt in optlist:
    if opt[0] == '--cedarFilename':
        cedarFilename = opt[1]
    elif opt[0] == '--hdf5Filename':
        hdf5Filename = opt[1]
    elif opt[0] == '--independentSpatialParms':
        independentSpatialParms = opt[1].split(',')
    elif opt[0] == '--arraySplittingParms':
        arraySplittingParms = opt[1].split(',')
    elif opt[0] == '--extraParameters':
        extraParameters = opt[1].split(',')
    elif opt[0] == '--filter':
        filter = opt[1]
        
    else:
        raise ValueError('Illegal option %s\n%s' % (opt[0], usage))

# check that all required arguments passed in
if cedarFilename == None:
    print('--cedarFilename argument required - must be full path of existing madrigal file')
    print(usage)
    sys.exit(0)

if hdf5Filename == None:
    print('--hdf5Filename argument required - must be full path of hdf5 file to write')
    sys.exit(0)
    
fileObj = madrigal.data.MadrigalFile(cedarFilename)

fileObj.exportToHdf(output = hdf5Filename,
                    independentSpatialParms = independentSpatialParms,
                    arraySplittingParms = arraySplittingParms,
                    extraParameters = extraParameters,
                    filter = filter)

print()
print('The file %s has been converted to hdf5 format: %s' % (cedarFilename, hdf5Filename))
