#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-09 12:10:00 junaid"
# ================================================================
#
# Author:  Muhammad Junaid <mjo147@uregina.ca>
#
# Copyright (c) junaid
#

import pandas as pd
import sys, os, time
import csv

################################################################################################################################################
'''
User Inputs
'''
ROOTPrefix = sys.argv[1]
runType = sys.argv[2]
timestmp=sys.argv[3]
column1=sys.argv[4]
column2=sys.argv[5]

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

def display_columns(csv_file, column1, column2):

    csv_file = UTILPATH+"/scripts/efficiency/OUTPUTS/%s_%s_efficiency_data_%s.csv"  % (ROOTPrefix.replace("replay_",""),runType,timestmp)

    with open(csv_file, 'r', newline='') as file:
        reader = csv.DictReader(file)

        # Print column names
        print("{} --- {}".format(column1, column2))
        # Iterate over rows and print values
        for row in reader:
#            eff = 0.9890
#            if float(row[column2]) < eff:
#               print(row[column1], " ------ " ,row[column2])
	
             print(row[column1], " ------ " ,row[column2])

csv_path = UTILPATH

display_columns(csv_path, column1, column2)
