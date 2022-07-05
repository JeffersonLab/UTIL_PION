#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2022-06-22 06:11:18 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import numpy as np
import pandas as pd
import sys,os

##################################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=1:
    print("!!!!! ERROR !!!!!\n Expected 1 arguments\n Usage is with - RUNNUMBER\n!!!!! ERROR !!!!!")
    sys.exit(1)

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
import ltsep as lt 

# Add this to all files for more dynamic pathing
USER =  lt.SetPath(os.path.realpath(__file__)).getPath("USER") # Grab user info for file finding
HOST = lt.SetPath(os.path.realpath(__file__)).getPath("HOST")
REPLAYPATH = lt.SetPath(os.path.realpath(__file__)).getPath("REPLAYPATH")
UTILPATH = lt.SetPath(os.path.realpath(__file__)).getPath("UTILPATH")
ANATYPE = lt.SetPath(os.path.realpath(__file__)).getPath("ANATYPE")

################################################################################################################################################

# Input params - run number and max number of events
runNum = sys.argv[1]
runList = "runlist_pionLT_2022.csv"

inp_f = UTILPATH+"/"+runList

# Converts csv data to dataframe
try:
    charge_data = pd.read_csv(inp_f)
    #print(charge_data.keys())
except IOError:
    print("Error: %s does not appear to exist." % inp_f)

run_data = charge_data[charge_data["Run"] == int(runNum)]

try:
    run_data["EBeam"].iloc[0]
except IndexError:
    print('''
    ERROR: Invalid run number {0}
           Not found in {1}
    '''.format(runNum,runList))
    sys.exit(0)

kinDict = {
    "EBeam" : run_data["EBeam"].iloc[0],
    "P_HMS" : run_data["P_HMS"].iloc[0],
    "Theta_HMS" : run_data["Theta_HMS"].iloc[0],
    "Theta_SHMS" : run_data["Theta_SHMS"].iloc[0],
    "P_SHMS" : run_data["P_SHMS"].iloc[0],
    "Target" : run_data["Target"].iloc[0],
}

#print(kinDict)

bool_data = [i for i,row in charge_data.iterrows() if kinDict['EBeam'] == row['EBeam'] and kinDict['P_HMS'] == row['P_HMS'] and kinDict['Theta_HMS'] == row['Theta_HMS'] and kinDict['Theta_SHMS'] == row['Theta_SHMS'] and kinDict['P_SHMS'] == row['P_SHMS'] and kinDict['Target'] == row['Target'] and row['Comments'].isdigit()]
#print(bool_data)

charge_data = charge_data.iloc[bool_data]

runVals = charge_data['Run'].values
try:
    PI_TOT = sum([int(val) for val in charge_data['Comments']])
except ValueError:
    print('''
    ERROR: Invalid comment '{0}'
           given for run {1}
    '''.format(run_data['Comments'].values[0],runNum))
    sys.exit(0)
Q_TOT = charge_data['Charge_(mC)'].sum()

print('''

{0}
For runs
{1}
Total Pion count is...
{2}
Total charge is...
{3:.3f} mC
{0}

'''.format("-"*(len(runVals)*6),runVals,PI_TOT,Q_TOT))
