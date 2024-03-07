#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-09 12:05:00 junaid"
# ================================================================
#
# Created: Muhammad junaid  <mjo147@uregina.ca>
# Copyright (c) trottar & junaid
#

import pandas as pd
import sys, os, time

################################################################################################################################################
'''
User Inputs
'''
ROOTPrefix = sys.argv[1]
runType = sys.argv[2]
runNum = sys.argv[3]
timestmp=sys.argv[4]
column=sys.argv[5]

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
import ltsep as lt 

# Import package for cuts
from ltsep import Root

lt=Root(os.path.realpath(__file__),"efficiency")

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
SIMCPATH=lt.SIMCPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH

################################################################################################################################################

def getTable():
    # Output for luminosity table
    inp_f = UTILPATH+"/scripts/efficiency/OUTPUTS/%s_%s_efficiency_data_%s.csv"  % (ROOTPrefix.replace("replay_",""),runType,timestmp)

    # Converts csv data to dataframe
    try:
        eff_data = pd.read_csv(inp_f)
        print(inp_f)
        #print(eff_data.keys())
    except IOError:
        print("Error: %s does not appear to exist." % inp_f)
        sys.exit(0)
        
    eff_data = eff_data[eff_data['Run_Number'] == float(runNum)]
    return eff_data

def uniqueColumn(column,eff_data):
    for val in column.split(","):
        print('''
    {0}\t{1}
        '''.format(val,eff_data[val].item()))

def main():

    eff_data = getTable()

    if "All" == column:

        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print(eff_data)
    elif "List" == column:
        print(eff_data.keys())
    else:
        uniqueColumn(column,eff_data)

if __name__ == "__main__":
    main()

    
