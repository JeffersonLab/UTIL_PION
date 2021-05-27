#! /usr/bin/python

# 26/05/21 - Stephen Kay, University of Regina
# Script to extract the pion and kaon cointime peak from the data and save info as a new rootfile, subsequent script fits the peak to extract the position to examine the stability
# This version is for the HeepCoin data, here we ONLY care about the ep coin time

# Import relevant packages
import uproot as up
import numpy as np
import root_numpy as rnp
import pandas as pd
import root_pandas as rpd
import ROOT
import scipy
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import sys, math, os, subprocess

sys.path.insert(0, 'python/')
# Check the number of arguments provided to the script
if len(sys.argv)-1!=3:
    print("!!!!! ERROR !!!!!\n Expected 3 arguments\n Usage is with - ROOTfilePrefix RunNumber MaxEvents \n!!!!! ERROR !!!!!")
    sys.exit(1)
# Input params - run number and max number of events
ROOTPrefix = sys.argv[1]
runNum = sys.argv[2]
MaxEvent = sys.argv[3]

USER = subprocess.getstatusoutput("whoami") # Grab user info for file finding
HOST = subprocess.getstatusoutput("hostname")
if ("farm" in HOST[1]):
    REPLAYPATH = "/group/c-kaonlt/USERS/%s/hallc_replay_lt" % USER[1]
elif ("qcd" in HOST[1]):
    REPLAYPATH = "/group/c-kaonlt/USERS/%s/hallc_replay_lt" % USER[1]
elif ("phys.uregina" in HOST[1]):
    REPLAYPATH = "/home/%s/work/JLab/hallc_replay_lt" % USER[1]
elif("skynet" in HOST[1]):
    REPLAYPATH = "/home/%s/Work/JLab/hallc_replay_lt" % USER[1]
    
# Add more path setting as needed in a similar manner
OUTPATH = "%s/UTIL_PION/OUTPUT/Analysis/PionLT" % REPLAYPATH
CUTPATH = "%s/UTIL_PION/DB/CUTS" % REPLAYPATH
sys.path.insert(0, '%s/UTIL_PION/bin/python/' % REPLAYPATH)
import kaonlt as klt

print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER[1], HOST[1], REPLAYPATH))
rootName = "%s/UTIL_PION/ROOTfiles/Analysis/PionLT/%s_%s_%s.root" % (REPLAYPATH, ROOTPrefix, runNum, MaxEvent)
if os.path.exists(OUTPATH):
    if os.path.islink(OUTPATH):
        pass
    elif os.path.isdir(OUTPATH):
        pass
    else:
        print ("%s exists but is not a directory or sym link, check your directory/link and try again" % (OUTPATH))
        sys.exit(2)
else:
    print("Output path not found, please make a sym link or directory called OUTPUT in UTIL_PION/scripts/demo to store output")
    sys.exit(3)
print ("Attempting to process %s" %(rootName))
if os.path.isfile(rootName):
    print ("%s exists, attempting to process" % (rootName))
else:
    print ("%s not found - do you have the correct sym link/folder set up?" % (rootName))
    sys.exit(4)
print("Output path checks out, outputting to %s" % (OUTPATH))

# Read stuff from the main event tree
e_tree = up.open(rootName)["T"]
# Timing info
CTime_ePiCoinTime_ROC1 = e_tree.array("CTime.ePiCoinTime_ROC1")
CTime_eKCoinTime_ROC1 = e_tree.array("CTime.eKCoinTime_ROC1")
CTime_epCoinTime_ROC1 = e_tree.array("CTime.epCoinTime_ROC1")
# HMS info
H_gtr_beta = e_tree.array("H.gtr.beta")
H_gtr_xp = e_tree.array("H.gtr.th") # xpfp -> Theta
H_gtr_yp = e_tree.array("H.gtr.ph") # ypfp -> Phi
H_gtr_dp = e_tree.array("H.gtr.dp")
H_cal_etotnorm = e_tree.array("H.cal.etotnorm")
H_cer_npeSum = e_tree.array("H.cer.npeSum")
# SHMS info
P_gtr_beta = e_tree.array("P.gtr.beta")
P_gtr_xp = e_tree.array("P.gtr.th") # xpfp -> Theta
P_gtr_yp = e_tree.array("P.gtr.ph") # ypfp -> Phi
P_gtr_p = e_tree.array("P.gtr.p")
P_gtr_dp = e_tree.array("P.gtr.dp")
P_cal_etotnorm = e_tree.array("P.cal.etotnorm")
P_aero_npeSum = e_tree.array("P.aero.npeSum")
P_hgcer_npeSum = e_tree.array("P.hgcer.npeSum")
P_hgcer_xAtCer = e_tree.array("P.hgcer.xAtCer")
P_hgcer_yAtCer = e_tree.array("P.hgcer.yAtCer")
MMpi = e_tree.array("P.kin.secondary.MMpi")
MMK = e_tree.array("P.kin.secondary.MMK")
MMp = e_tree.array("P.kin.secondary.MMp")
# Relevant branches now stored as NP arrays

r = klt.pyRoot()
fout = '%s/UTIL_PION/DB/CUTS/run_type/coinpeak.cuts' % REPLAYPATH
# read in cuts file and make dictionary
#c = klt.pyPlot(REPLAYPATH,DEBUG=True)
c = klt.pyPlot(REPLAYPATH)
readDict = c.read_dict(fout,runNum)
# This method calls several methods in kaonlt package. It is required to create properly formated
# dictionaries. The evaluation must be in the analysis script because the analysis variables (i.e. the
# leaves of interest) are not defined in the kaonlt package. This makes the system more flexible
# overall, but a bit more cumbersome in the analysis script. Perhaps one day a better solution will be
# implemented.
def make_cutDict(cut,inputDict=None):

    global c

    c = klt.pyPlot(REPLAYPATH,readDict)
    x = c.w_dict(cut)
    print("%s" % cut)
    print("x ", x)
    
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

cutDict = make_cutDict("coin_ep_cut_all")
c = klt.pyPlot(REPLAYPATH,cutDict)

def coin_events(): 
    # Define the array of arrays containing the relevant HMS and SHMS info
    All_Events_Uncut_tmp = [H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_cal_etotnorm, H_cer_npeSum, CTime_ePiCoinTime_ROC1, CTime_eKCoinTime_ROC1, CTime_epCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_cal_etotnorm, P_aero_npeSum, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMpi, MMK, MMp]
    All_Events_Uncut = [(HBeta, Hxp, Hyp, Hdel, HCal, HCer, CTPi, CTK, CTp, pBeta, pxp, pyp, pP, pDel, pCal, pAero, pHGC, pHGCX, pHGCY, mm1, mm2, mm3) for (HBeta, Hxp, Hyp, Hdel, HCal, HCer, CTPi, CTK, CTp, pBeta, pxp, pyp, pP, pDel, pCal, pAero, pHGC, pHGCX, pHGCY, mm1, mm2, mm3) in zip(*All_Events_Uncut_tmp)] 

    # Create array of arrays of protons after cuts, all events, prompt and random
    Proton_Events_All_tmp = []

    # Go over every array in All_Events_Uncut_tmp, append to the other arrays the array after a cut is applied
    for arr in All_Events_Uncut_tmp:
        Proton_Events_All_tmp.append(c.add_cut(arr, "coin_ep_cut_all")) # Apply PID but no cointime cut
      
    Proton_Events_All = [(HBeta, Hxp, Hyp, Hdel, HCal, HCer, CTPi, CTK, CTp, pBeta, pxp, pyp, pP, pDel, pCal, pAero, pHGC, pHGCX, pHGCY, mm1, mm2, mm3) for (HBeta, Hxp, Hyp, Hdel, HCal, HCer, CTPi, CTK, CTp, pBeta, pxp, pyp, pP, pDel, pCal, pAero, pHGC, pHGCX, pHGCY, mm1, mm2, mm3) in zip(*Proton_Events_All_tmp)]

    COIN_EventInfo = {
        "All_Events" : All_Events_Uncut,
        "Protons_All" : Proton_Events_All,
        }

    return COIN_EventInfo

def main():
    COIN_Data = coin_events()

    COIN_All_Data_Header = ["H_gtr_beta","H_gtr_xp","H_gtr_yp","H_gtr_dp","H_cal_etotnorm","H_cer_npeSum","CTime_ePiCoinTime_ROC1","CTime_eKCoinTime_ROC1","CTime_epCoinTime_ROC1","P_gtr_beta","P_gtr_xp","P_gtr_yp","P_gtr_p","P_gtr_dp","P_cal_etotnorm","P_aero_npeSum","P_hgcer_npeSum","P_hgcer_xAtCer","P_hgcer_yAtCer","MMpi","MMK","MMp"]

    data_keys = list(COIN_Data.keys()) # Create a list of all the keys in all dicts added above, each is an array of data
    #print(data_keys)

    for i in range (0, len(data_keys)):
        DFHeader=list(COIN_All_Data_Header)
        if (i == 0):
            pd.DataFrame(COIN_Data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_CTPeak_Data_HeepCoin.root" % (OUTPATH, runNum, MaxEvent), key ="%s" % data_keys[i])
        elif (i != 0):
            pd.DataFrame(COIN_Data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_CTPeak_Data_HeepCoin.root" % (OUTPATH, runNum, MaxEvent), key ="%s" % data_keys[i], mode ='a') 
                    
if __name__ == '__main__':
    main()
