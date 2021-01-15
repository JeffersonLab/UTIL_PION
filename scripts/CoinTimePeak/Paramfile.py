#! /usr/bin/python

# 06/11/20 - Stephen Kay, University of Regina

# Python  script which creates a new timing paramater file based upon a provided list of inputs and a provided existing parameter file

# Import relevant packages
import numpy as np
import root_numpy as rnp
import pandas as pd
import ROOT
import scipy
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import sys, math, os, subprocess
from os import path

sys.path.insert(0, 'python/')
if len(sys.argv)-1!=2:
    print("!!!!! ERROR !!!!!\n Expected 2 arguments\n Usage is with - InputKinList OldParamFile \n!!!!! ERROR !!!!!")
    sys.exit(1)
# Input params - run number and max number of events
InputList = sys.argv[1]
OldParam = sys.argv[2]

USER = subprocess.getstatusoutput("whoami") # Grab user info for file finding
HOST = subprocess.getstatusoutput("hostname")
if ("farm" in HOST[1]):
    REPLAYPATH = "/group/c-kaonlt/USERS/%s/hallc_replay_lt" % USER[1]
elif ("qcd" in HOST[1]):
    REPLAYPATH = "/group/c-kaonlt/USERS/%s/hallc_replay_lt" % USER[1]
elif ("lark.phys.uregina" in HOST[1]):
    REPLAYPATH = "/home/%s/work/JLab/hallc_replay_lt" % USER[1]

# Add more path setting as needed in a similar manner
OUTPATH = "%s/UTIL_KAONLT/scripts/CoinTimePeak/OUTPUT" % REPLAYPATH
PARAMPATH = "%s/UTIL_KAONLT/DB/PARAM" % REPLAYPATH
CUTPATH = "%s/UTIL_KAONLT/DB/CUTS" % REPLAYPATH

print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER[1], HOST[1], REPLAYPATH))
KinList = "%s/UTIL_KAONLT/scripts/CoinTimePeak/Kinematics/%s" % (REPLAYPATH,InputList)
TimingCutFile = "%s/%s.csv" % (PARAMPATH,OldParam)
if (path.exists(KinList) == False or path.isfile(KinList) == False):
    print("!!!!! ERROR !!!!!\n %s does not exist or is not a valid file - check 1st input arg and try again \n!!!!! ERROR !!!!!" % KinList)
    sys.exit(2)
else:
    print("%s used as list of kinematics/settings to grab fit values for" % KinList)
if (path.exists(TimingCutFile) == False or path.isfile(TimingCutFile) == False):
    print("!!!!! ERROR !!!!!\n %s does not exist or is not a valid file - check 2nd input arg and try again \n!!!!! ERROR !!!!!" % TimingCutFile)
    sys.exit(3)
else:
    print("%s being used as old param file to base offsets/windows from" % TimingCutFile)

KinListf = open(KinList)
#KinLineNum = 0 # Uncomment for loop testing if needed
ParamData=[]
FailedParamData=[]
# Open kinematic list file and go through every kinematic in the list
for KinLine in KinListf:
    #KinLineNum += 1 # Uncomment for testing loop if needed
    #print(KinLineNum) # Testing output
    KinLine=KinLine.rstrip();
    KinFile = "%s/%s_Output.csv" % (OUTPATH, KinLine)
    # Check the corresponding output .csv file exists, if it does, open it
    if (path.exists(KinFile) == True and path.isfile(KinFile) == True):
        KinFilef = open(KinFile)
        #KinFileLineNum = 0 #
        for KinFileLine in KinFilef:
            # Check the line is not blank and process it if not
            if (KinFileLine != "\n"):
                RunParamData = [0] * 10 # Initialise a 10 element array for the param data for each run (line) in the KinFile
                # Go through output csv line by line, convert each line to an array
                #KinFileLineNum += 1
                KinFileLine=KinFileLine.rstrip() # Remoe trailing spaces
                KinFileData=KinFileLine.split(",") # Convert csv row to an array
                RunParamData[0]=int(KinFileData[0]) # Run number
                #RunParamData[1]=int(KinFileData[0])+1
                RunParamData[1]=int(KinFileData[0]) # Run start = run end as this is to cover one run at a time
                RunParamData[6]=float(KinFileData[1]) # Pion peak pos
                RunParamData[7]=float(KinFileData[5]) # Kaon peak pos
                RunParamData[8]=float(KinFileData[9]) # Proton peak pos
                # Need to now get the rest of the run param values from the existing param file which we were also given
                TempPar = -1 # To check later
                TimingCutLineNum = 0 # Count Timing cut file line number
                TimingCutFilef = open(TimingCutFile)
                for TimingCutLine in TimingCutFilef:
                    TimingCutLineNum += 1 # Add one to line number at start of each loop
                    if(TimingCutLineNum > 1): # Skip first line
                        TimingCutLine = TimingCutLine.partition('#')[0] # Treat anything after a # as a comment and ignore it
                        TimingCutLine = TimingCutLine.rstrip()
                        OldParamFileArr = TimingCutLine.split(",") # Convert line into an array, anything after a comma is a new entry
                        if(int(RunParamData[0]) in range (int(OldParamFileArr[0]), int(OldParamFileArr[1])+1)): # Check if run number for file is within any of the ranges specified in the cut file
                            TempPar += 2 # If run number is in range, set to non -1 value
                            RunParamData[2]=float(OldParamFileArr[2])
                            RunParamData[3]=float(OldParamFileArr[3])
                            RunParamData[4]=float(OldParamFileArr[4])
                            RunParamData[5]=float(OldParamFileArr[5])
                            RunParamData[9]=float(OldParamFileArr[9])
                # Close loop over timing cut file and close the file
                TimingCutFilef.close()
                # If the run number was in the old param file, our new parameter file entry is now complete, we need to append it to our array for the new param file
                if (TempPar != -1):
                    ParamData.append(RunParamData)
                elif(TempPar == -1): # If the run number wasn't found in the old param file add it to a failed array file
                    FailedParamData.append(RunParamData)
            # This is the else condition of if a line in the output for a kinematic file was emptry
            else:
                print("!!! WARNING !!! - Blank line in %s - check output is OK - !!! WARNING !!!" % KinFile)
        # End loop over output file for a kinematic and close it
        KinFilef.close()
    # If the output csv file for a kinematic doesn't exist or can't be opened, skip it
    else:
        print("!!!!! ERROR !!!!!\n %s does not exist or is not a valid file - Skipping \n!!!!! ERROR !!!!!" % KinFile)
# End loop over kinematic file list and cloe it
KinListf.close()

ParamDataArr=np.array(ParamData, dtype='O')
# Should now have an array of param entries for our file, need to print them to file, first, sort them by run number
ParamDataArr = ParamDataArr[ParamDataArr[:,0].argsort()] # Sort by values in 1st column (starting run number)
# Save to file with appropriate formatting and header
np.savetxt(("%s_TimingParamFile.csv" % InputList), ParamDataArr, fmt="%i,%i,%2.3f,%2.3f,%i,%i,%3.3f,%3.3f,%3.3f,%2.3f", delimiter=",", header='Run_Start,Run_End,Bunch_Spacing,Coin_Offset,nSkip,nWindows,Pion_Prompt_Peak,Kaon_Prompt_Peak,Proton_Prompt_Peak,RF_Offset', comments='')
# If there are any runs that didn't match the fed parameter file, write them as a different csv
if(len(FailedParamData) != 0):
    FailedParamDataArr=np.array(FailedParamData, dtype='O')
    FailedParamDataArr = FailedParamDataArr[FailedParamDataArr[:,0].argsort()] # Sort by values in 1st column (starting run number)
    np.savetxt(("%s_Failed_TimingParamFile.csv" % InputList), FailedParamDataArr, fmt="%i,%i,%2.3f,%2.3f,%i,%i,%3.3f,%3.3f,%3.3f,%2.3f", delimiter=",", header='Run_Start,Run_End,Bunch_Spacing,Coin_Offset,nSkip,nWindows,Pion_Prompt_Peak,Kaon_Prompt_Peak,Proton_Prompt_Peak,RF_Offset', comments='')




