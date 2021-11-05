#! /usr/bin/python

#
# Description: Script to dynamically set new trigger windows and update the param file with these values
# ================================================================
# Time-stamp: "2021-11-03 07:27:31 trottar"
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

################################################################################################################################################
'''
User Inputs
'''

RunType = sys.argv[1]
ROOTPrefix = sys.argv[2]
runNum = sys.argv[3]
MaxEvent=sys.argv[4]

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
ANATYPE=lt.SetPath(os.path.realpath(__file__)).getPath("ANATYPE")

################################################################################################################################################

print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))

out_f = UTILPATH+"/scripts/trig_windows/OUTPUTS/trig_data.csv"

################################################################################################################################################
'''
Check that root/output paths and files exist for use
'''

# Construct the name of the rootfile based upon the info we provided
OUTPATH = UTILPATH+"/OUTPUT/Analysis/%sLT" % ANATYPE        # Output folder location
rootName = UTILPATH+"/ROOTfiles/Analysis/Lumi/%s_%s_%s.root" % (ROOTPrefix,runNum,MaxEvent)     # Input file location and variables taking
print ("Attempting to process %s" %(rootName))
lt.SetPath(os.path.realpath(__file__)).checkDir(OUTPATH)
lt.SetPath(os.path.realpath(__file__)).checkFile(rootName)
print("Output path checks out, outputting to %s" % (OUTPATH))

################################################################################################################################################
'''
Grab prescale values and tracking efficiencies from report file
'''

# Open report file to grab prescale values
report = UTILPATH+"/REPORT_OUTPUT/Analysis/%s/%s_%s_%s.report" % (RunType,ROOTPrefix,runNum,MaxEvent)
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

################################################################################################################################################
'''
Define prescale variables
'''

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

################################################################################################################################################

'''
ANALYSIS TREE, T
'''

tree = up.open(rootName)["T"]

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

################################################################################################################################################
'''
Define and set up cuts
'''

fout = UTILPATH+'/DB/CUTS/run_type/lumi.cuts'

cuts = ["c_nozero_edtm","c_nozero_ptrigHMS","c_nozero_ptrigSHMS"]
# Check if COIN trigger is used
if len(PS_used) > 2:
    cuts = ["c_nozero_edtm","c_nozero_ptrigHMS","c_nozero_ptrigSHMS","c_nozero_ptrigCOIN"]

def make_cutDict(cuts,fout,runNum,CURRENT_ENV):
    '''
    This method calls several methods in kaonlt package. It is required to create properly formated
    dictionaries. The evaluation must be in the analysis script because the analysis variables (i.e. the
    leaves of interest) are not defined in the kaonlt package. This makes the system more flexible
    overall, but a bit more cumbersome in the analysis script. Perhaps one day a better solution will be
    implimented.
    '''

    # read in cuts file and make dictionary
    importDict = lt.SetCuts(CURRENT_ENV).importDict(cuts,fout,runNum)
    for i,cut in enumerate(cuts):
        x = lt.SetCuts(CURRENT_ENV,importDict).booleanDict(cut)
        print("\n%s" % cut)
        print(x, "\n")
        if i == 0:
            inputDict = {}
        cutDict = lt.SetCuts(CURRENT_ENV,importDict).readDict(cut,inputDict)
        for j,val in enumerate(x):
            cutDict = lt.SetCuts(CURRENT_ENV,importDict).evalDict(cut,eval(x[j]),cutDict)
    return lt.SetCuts(CURRENT_ENV,cutDict)

c = make_cutDict(cuts,fout,runNum,os.path.realpath(__file__))

################################################################################################################################################
'''
Grabs the misc param file so it can be updated with trigger windows
'''

# Read in the Misc_Parameters.csv cut parameter file which has the trigger window information
inp_f = UTILPATH+'/DB/PARAM/Misc_Parameters.csv'
try:
    trig_data = pd.read_csv(inp_f)
except IOError:
    print("Error: %s does not appear to exist." % inp_f)
print(trig_data.keys())

################################################################################################################################################

def setWindows(runNum):
    '''
    Set the trigger windows by...

    1) Finding the number of events per bin
    2) Based off the threshold for the number of events, the min and max windows are set.  
    '''

    def getBinEdges(branch,cut,numbins):
        # Calculate the geometric mean
        def geo_mean(inp):
            l_inp = np.log(inp[inp != 0])
            return np.exp(l_inp.mean())
        # finds the number of events per bin (the zero events are cut out beforehand)
        counts, bins = np.histogram(c.add_cut(branch,cut),bins=numbins)
        # Set the threshold for the number of events an order of magnitude more than the geometric mean
        thres_count = geo_mean(counts)*10
        # Finding the bins that are above the set threshold for the number of events
        binVals = [b for c,b in zip(counts,bins) if c > thres_count]
        # Set min and max windows
        minBin = min(binVals)
        maxBin = max(binVals)
        return [minBin, maxBin]

    # Get windows for {SPEC}_ROC1_tdcTimeRaw and pEDTM_tdcTimeRaw
    c_T_coin_pTRIG_HMS_ROC1_tdcTimeRaw = getBinEdges(T_coin_pTRIG_HMS_ROC1_tdcTimeRaw,"c_nozero_ptrigHMS",200)
    c_T_coin_pTRIG_SHMS_ROC2_tdcTimeRaw = getBinEdges(T_coin_pTRIG_SHMS_ROC2_tdcTimeRaw,"c_nozero_ptrigSHMS",200)
    # Check if COIN trigger is used
    if len(PS_used) > 2:
        c_T_coin_pTRIG_COIN_ROC1_tdcTimeRaw = getBinEdges(T_coin_pTRIG_COIN_ROC1_tdcTimeRaw,"c_nozero_ptrigCOIN",200)
    c_T_coin_pEDTM_tdcTimeRaw = getBinEdges(T_coin_pEDTM_tdcTimeRaw,"c_nozero_edtm",200)

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

################################################################################################################################################

def reconParam(runNum):
    '''
    Reconstruct /DB/PARAM/Misc_Parameters.csv with new window values
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

################################################################################################################################################

def main():

    reconParam(runNum)

if __name__ == '__main__': main()
