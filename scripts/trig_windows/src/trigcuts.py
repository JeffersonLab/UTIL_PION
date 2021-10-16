#! /usr/bin/python

#
# Description: Script to dynamically set new trigger windows and update the param file with these values
# ================================================================
# Time-stamp: "2021-10-15 03:45:15 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

import uproot as up
import numpy as np
import pandas as pd
import scipy
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import sys, math, os, subprocess

RunType = sys.argv[1]
ROOTPrefix = sys.argv[2]
runNum = sys.argv[3]
MaxEvent=sys.argv[4]

# Add this to all files for more dynamic pathing
USER = subprocess.getstatusoutput("whoami") # Grab user info for file finding
HOST = subprocess.getstatusoutput("hostname")

# Set path depending upon hostname. Change or add more as needed  
if ("farm" in HOST[1]):
    REPLAYPATH="/group/c-pionlt/online_analysis/hallc_replay_lt"
elif ("lark" in HOST[1]):
    REPLAYPATH = "/home/%s/work/JLab/hallc_replay_lt" % USER[1]
elif ("cdaq" in HOST[1]):
    REPLAYPATH = "/home/cdaq/hallc-online/hallc_replay_lt"
elif ("trottar" in HOST[1]):
    REPLAYPATH = "/home/trottar/Analysis/hallc_replay_lt"

# Import package for cuts
sys.path.insert(0, '%s/UTIL_PION/bin/python/' % REPLAYPATH)
import kaonlt as klt

print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER[1], HOST[1], REPLAYPATH))

out_f = "%s/UTIL_PION/scripts/trig_windows/OUTPUTS/trig_data.csv" % REPLAYPATH

# Construct the name of the rootfile based upon the info we provided
OUTPATH = "%s/UTIL_PION/OUTPUT/Analysis/PionLT" % REPLAYPATH        # Output folder location
rootName = "%s/UTIL_PION/ROOTfiles/Analysis/%s/%s_%s_%s.root" % (REPLAYPATH,RunType,ROOTPrefix,runNum,MaxEvent)     # Input file location and variables taking
print ("Attempting to process %s" %(rootName))
if os.path.exists(OUTPATH):
    if os.path.islink(OUTPATH):
        pass
    elif os.path.isdir(OUTPATH):
        pass
    else:
        print ("%s exists but is not a directory or sym link, check your directory/link and try again" % (OUTPATH))
        sys.exit(2)
else:
    print("Output path not found, please make a sym link or directory called OUTPUT in UTIL_PION to store output")
    sys.exit(3)
if os.path.isfile(rootName):
    print ("%s exists, processing" % (rootName))
else:
    print ("%s not found - do you have the correct sym link/folder set up?" % (rootName))
    sys.exit(4)
print("Output path checks out, outputting to %s" % (OUTPATH))

# Open report file to grab prescale values
report = "%s/UTIL_PION/REPORT_OUTPUT/Analysis/%s/%s_%s_%s.report" % (REPLAYPATH,RunType,ROOTPrefix,runNum,MaxEvent)
f = open(report)
psList = ['SW_Ps1_factor','SW_Ps2_factor','SW_Ps3_factor','SW_Ps4_factor','SW_Ps5_factor','SW_Ps6_factor']
   
# Prescale input value (psValue) to its actual DAQ understanding (psActual)
psActual = [-1,1,2,3,5,9,17,33,65,129,257,513,1025,2049,4097,8193,16385,32769]
psValue = [-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]

# Search root file for prescale values, then save as variables
for line in f:
    data = line.split(':')
    for i, obj in enumerate(psList) :
        if (psList[i] in data[0]) : 
            if (i == 0) :  
                ps1_tmp = data[1].strip()
            if (i == 1) : 
                ps2_tmp = data[1].strip()
            if (i == 2) :
                ps3_tmp = data[1].strip()
            if (i == 3) :
                ps4_tmp = data[1].strip()
            if (i == 4) :
                ps5_tmp = data[1].strip()
            if (i == 5) :
                ps6_tmp = data[1].strip()
ps1=int(ps1_tmp)
ps2=int(ps2_tmp)
ps3=int(ps3_tmp)
ps4=int(ps4_tmp)
ps5=int(ps5_tmp)
ps6=int(ps6_tmp)

# Convert the prescale input values to their actual DAQ values
for i,index in enumerate(psActual):
    #psValue
    if (index == ps1) :
        if(index == -1):
            PS1 = 0
        else:
            PS1 = psActual[i]
    if (index == ps2) :
        if(index == -1):
            PS2 = 0
        else:
            PS2 = psActual[i]            
    if (index == ps3) :
        if(index == -1):
            PS3 = 0
        else:
            PS3 = psActual[i]
    if (index == ps4) :
        if(index == -1):
            PS4 = 0
        else:
            PS4 = psActual[i]            
    if (index == ps5) :
        if(index == -1):
            PS5 = 0
        else:
            PS5 = psActual[i]
    if (index == ps6) :
        if(index == -1):
            PS6 = 0
        else:
            PS6 = psActual[i]            
f.close()

print("\nPre-scale values...\nPS1:{0}, PS2:{1}, PS3:{2}, PS4:{3}, PS5:{4}, PS6:{5}\n".format(PS1,PS2,PS3,PS4,PS5,PS6))

# Save only the used prescale triggers to the PS_used list
PS_list = [["PS1",PS1],["PS2",PS2],["PS3",PS3],["PS4",PS4],["PS5",PS5],["PS6",PS6]]
PS_used = []
for val in PS_list:
    if val[1] != 0:
        PS_used.append(val)

# Check if COIN trigger is used by seeing it was saved in the PS_used list
if len(PS_used) > 2:
    PS_names = [PS_used[0][0],PS_used[1][0],PS_used[2][0]]
    SHMS_PS = PS_used[0][1]
    HMS_PS = PS_used[1][1]
    COIN_PS = PS_used[2][1]
else:
    PS_names = [PS_used[0][0],PS_used[1][0]]
    SHMS_PS = PS_used[0][1]
    HMS_PS = PS_used[1][1]

'''
ANALYSIS TREE, T
'''

tree = up.open(rootName)["T"]
branch = klt.pyBranch(tree)

H_bcm_bcm4a_AvgCurrent = tree.array("H.bcm.bcm4a.AvgCurrent")

if PS_names[0] is "PS1":
    T_coin_pTRIG_SHMS_ROC1_tdcTimeRaw = tree.array("T.coin.pTRIG1_ROC1_tdcTimeRaw")
    T_coin_pTRIG_SHMS_ROC2_tdcTimeRaw = tree.array("T.coin.pTRIG1_ROC2_tdcTimeRaw")
    T_coin_pTRIG_SHMS_ROC1_tdcTime = tree.array("T.coin.pTRIG1_ROC1_tdcTime")
    T_coin_pTRIG_SHMS_ROC2_tdcTime = tree.array("T.coin.pTRIG1_ROC2_tdcTime")

if PS_names[0] is "PS2":
    T_coin_pTRIG_SHMS_ROC1_tdcTimeRaw = tree.array("T.coin.pTRIG2_ROC1_tdcTimeRaw")
    T_coin_pTRIG_SHMS_ROC2_tdcTimeRaw = tree.array("T.coin.pTRIG2_ROC2_tdcTimeRaw")
    T_coin_pTRIG_SHMS_ROC1_tdcTime = tree.array("T.coin.pTRIG2_ROC1_tdcTime")
    T_coin_pTRIG_SHMS_ROC2_tdcTime = tree.array("T.coin.pTRIG2_ROC2_tdcTime")

if PS_names[1] is "PS3":
    T_coin_pTRIG_HMS_ROC1_tdcTimeRaw = tree.array("T.coin.pTRIG3_ROC1_tdcTimeRaw")
    T_coin_pTRIG_HMS_ROC2_tdcTimeRaw = tree.array("T.coin.pTRIG3_ROC2_tdcTimeRaw")
    T_coin_pTRIG_HMS_ROC1_tdcTime = tree.array("T.coin.pTRIG3_ROC1_tdcTime")
    T_coin_pTRIG_HMS_ROC2_tdcTime = tree.array("T.coin.pTRIG3_ROC2_tdcTime")

if PS_names[1] is "PS4":
    T_coin_pTRIG_HMS_ROC1_tdcTimeRaw = tree.array("T.coin.pTRIG4_ROC1_tdcTimeRaw")
    T_coin_pTRIG_HMS_ROC2_tdcTimeRaw = tree.array("T.coin.pTRIG4_ROC2_tdcTimeRaw")
    T_coin_pTRIG_HMS_ROC1_tdcTime = tree.array("T.coin.pTRIG4_ROC1_tdcTime")
    T_coin_pTRIG_HMS_ROC2_tdcTime = tree.array("T.coin.pTRIG4_ROC2_tdcTime")

# Check if COIN trigger is used
if len(PS_used) > 2:
    if PS_names[2] is "PS5":
        T_coin_pTRIG_COIN_ROC1_tdcTimeRaw = tree.array("T.coin.pTRIG5_ROC1_tdcTimeRaw")
        T_coin_pTRIG_COIN_ROC2_tdcTimeRaw = tree.array("T.coin.pTRIG5_ROC2_tdcTimeRaw")
        T_coin_pTRIG_COIN_ROC1_tdcTime = tree.array("T.coin.pTRIG5_ROC1_tdcTime")
        T_coin_pTRIG_COIN_ROC2_tdcTime = tree.array("T.coin.pTRIG5_ROC2_tdcTime")
        
    if PS_names[2] is "PS6":
        T_coin_pTRIG_COIN_ROC1_tdcTimeRaw = tree.array("T.coin.pTRIG6_ROC1_tdcTimeRaw")
        T_coin_pTRIG_COIN_ROC2_tdcTimeRaw = tree.array("T.coin.pTRIG6_ROC2_tdcTimeRaw")
        T_coin_pTRIG_COIN_ROC1_tdcTime = tree.array("T.coin.pTRIG6_ROC1_tdcTime")
        T_coin_pTRIG_COIN_ROC2_tdcTime = tree.array("T.coin.pTRIG6_ROC2_tdcTime")

T_coin_pEDTM_tdcTimeRaw = tree.array("T.coin.pEDTM_tdcTimeRaw")
T_coin_pEDTM_tdcTime = tree.array("T.coin.pEDTM_tdcTime")

fout = REPLAYPATH+'/UTIL_PION/DB/CUTS/run_type/lumi.cuts'

# read in cuts file and make dictionary
c = klt.pyPlot(REPLAYPATH)
readDict = c.read_dict(fout,runNum)

def make_cutDict(cut,inputDict=None):
    '''
    This method calls several methods in kaonlt package. It is required to create properly formated
    dictionaries. The evaluation must be in the analysis script because the analysis variables (i.e. the
    leaves of interest) are not defined in the kaonlt package. This makes the system more flexible
    overall, but a bit more cumbersome in the analysis script. Perhaps one day a better solution will be
    implimented.
    '''

    global c

    c = klt.pyPlot(REPLAYPATH,readDict)
    x = c.w_dict(cut)
    print("\n%s" % cut)
    print(x, "\n")
    
    if inputDict == None:
        inputDict = {}
        
    for key,val in readDict.items():
        if key == cut:
            inputDict.update({key : {}})

    for i,val in enumerate(x):
        tmp = x[i]
        if tmp == "":
            continue
        else:
            inputDict[cut].update(eval(tmp))
        
    return inputDict

cutDict = make_cutDict("c_nozero")
c = klt.pyPlot(REPLAYPATH,cutDict)

# Read in the Misc_Parameters.csv cut parameter file which has the trigger window information
inp_f = REPLAYPATH+'/UTIL_PION/DB/PARAM/Misc_Parameters.csv'
try:
    trig_data = pd.read_csv(inp_f)
except IOError:
    print("Error: %s does not appear to exist." % inp_f)
print(trig_data.keys())

def setWindows(runNum):
    '''
    Set the trigger windows by...

    1) Finding the number of events per bin
    2) Based off the threshold for the number of events, the min and max windows are set.  
    '''

    def getBinEdges(branch,numbins,count_thres):
        # finds the number of events per bin (the zero events are cut out beforehand)
        counts, bins = np.histogram(c.add_cut(branch,"c_nozero"),bins=numbins)
        # Finding the bins that are above the set threshold for the number of events
        binVals = [b for c,b in zip(counts,bins) if c > count_thres]
        # Set min and max windows
        minBin = min(binVals)
        maxBin = max(binVals)
        return [minBin, maxBin]

    # Get windows for {SPEC}_ROC1_tdcTimeRaw and pEDTM_tdcTimeRaw
    c_T_coin_pTRIG_HMS_ROC1_tdcTimeRaw = getBinEdges(T_coin_pTRIG_HMS_ROC1_tdcTimeRaw,200,250)
    c_T_coin_pTRIG_SHMS_ROC2_tdcTimeRaw = getBinEdges(T_coin_pTRIG_SHMS_ROC2_tdcTimeRaw,200,250)
    # Check if COIN trigger is used
    if len(PS_used) > 2:
        c_T_coin_pTRIG_COIN_ROC1_tdcTimeRaw = getBinEdges(T_coin_pTRIG_COIN_ROC1_tdcTimeRaw,200,250)
    c_T_coin_pEDTM_tdcTimeRaw = getBinEdges(T_coin_pEDTM_tdcTimeRaw,200,250)

    # Create a dictionary that contains the information that will be uploaded to Misc_Parameters.csv for a particular run
    new_row = {'Run_Start' : "{:.0f}".format(float(runNum)), 'Run_End' : "{:.0f}".format(float(runNum)), 'noedtm' : 0.0, 'edtmLow' : "{:.0f}".format(float(c_T_coin_pEDTM_tdcTimeRaw[0])), 
               'edtmHigh' : "{:.0f}".format(float(c_T_coin_pEDTM_tdcTimeRaw[1])), 'ptrigHMSLow' : "{:.0f}".format(float(c_T_coin_pTRIG_HMS_ROC1_tdcTimeRaw[0])), 
               'ptrigHMSHigh' : "{:.0f}".format(float(c_T_coin_pTRIG_HMS_ROC1_tdcTimeRaw[1])), 'ptrigSHMSLow' : "{:.0f}".format(float(c_T_coin_pTRIG_SHMS_ROC2_tdcTimeRaw[0])), 
               'ptrigSHMSHigh' : "{:.0f}".format(float(c_T_coin_pTRIG_SHMS_ROC2_tdcTimeRaw[1])), 'ptrigCOINLow' : 0.0, 'ptrigCOINHigh' : 10000.0, 'goodstarttime' : 1.0, 'goodscinhit' : 1.0}
    # Check if COIN trigger is used
    if len(PS_used) > 2:
        # Create a dictionary that contains the information that will be uploaded to Misc_Parameters.csv for a particular run
        new_row = {'Run_Start' : "{:.0f}".format(float(runNum)), 'Run_End' : "{:.0f}".format(float(runNum)), 'noedtm' : 0.0, 'edtmLow' : "{:.0f}".format(float(c_T_coin_pEDTM_tdcTimeRaw[0])), 
                   'edtmHigh' : "{:.0f}".format(float(c_T_coin_pEDTM_tdcTimeRaw[1])), 'ptrigHMSLow' : "{:.0f}".format(float(c_T_coin_pTRIG_HMS_ROC1_tdcTimeRaw[0])), 
                   'ptrigHMSHigh' : "{:.0f}".format(float(c_T_coin_pTRIG_HMS_ROC1_tdcTimeRaw[1])), 'ptrigSHMSLow' : "{:.0f}".format(float(c_T_coin_pTRIG_SHMS_ROC2_tdcTimeRaw[0])), 
                   'ptrigSHMSHigh' : "{:.0f}".format(float(c_T_coin_pTRIG_SHMS_ROC2_tdcTimeRaw[1])), 'ptrigCOINLow' : "{:.0f}".format(float(c_T_coin_pTRIG_COIN_ROC1_tdcTimeRaw[0])), 
                   'ptrigCOINHigh' : "{:.0f}".format(float(c_T_coin_pTRIG_COIN_ROC1_tdcTimeRaw[1])), 'goodstarttime' : 1.0, 'goodscinhit' : 1.0}

    return new_row

def reconParam(runNum):
    '''
    Reconstruct /UTIL_PION/DB/PARAM/Misc_Parameters.csv with new window values
    '''
    # Get new windows and set up row to be added to Misc_Parameters.csv
    new_row = setWindows(runNum)

    global trig_data
    print("\nRemoving...\n",trig_data[(trig_data["Run_Start"] <= int(runNum)) & (trig_data["Run_End"] >= int(runNum))],"\n") 

    # Checking if run number is in a row already
    run_row = trig_data[(trig_data["Run_Start"] <= int(runNum)) & (trig_data["Run_End"] >= int(runNum))]

    # Removing row with this run number argument
    run_index = trig_data.index[(trig_data["Run_Start"] <= int(runNum)) & (trig_data["Run_End"] >= int(runNum))].tolist()
    trig_data.drop(run_index, inplace=True)

    # Setting an open window row that will be added to the end of Misc_Parameters.csv. This ensures that the script will run in the future without errors. 
    # This row will not overwrite the windows that are set above
    open_row = {'Run_Start' : 0, 'Run_End' : 99999, 'noedtm' : 0.0, 'edtmLow' : 0.0, 'edtmHigh' : 10000.0, 'ptrigHMSLow' : 0.0, 'ptrigHMSHigh' : 10000.0, 
               'ptrigSHMSLow' : 0.0, 'ptrigSHMSHigh' : 10000.0, 'ptrigCOINLow' : 0.0, 'ptrigCOINHigh' : 10000.0, 'goodstarttime' : 1.0, 'goodscinhit' : 1.0}

    # Add in newly formed rows to dataframe
    trig_data = trig_data.append(new_row,ignore_index=True)
    trig_data = trig_data.append(open_row,ignore_index=True)

    # Update csv with new version of dataframe
    trig_data.to_csv(inp_f, index=False, header=True, mode='w+',)

    print("\n\nNew version of Misc_Parameters.csv...\n",trig_data)
    return trig_data      

def main():

    reconParam(runNum)

if __name__ == '__main__': main()
