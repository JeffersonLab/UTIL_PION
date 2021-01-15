#! /usr/bin/python

# 15/10/20 - Stephen Kay, University of Regina
# Script to extract the pion and kaon cointime peak from the data and save info as a new rootfile, subsequent script fits the peak to extract the position to examine the stability

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
OUTPATH = "%s/UTIL_KAONLT/scripts/CoinTimePeak/OUTPUT" % REPLAYPATH
CUTPATH = "%s/UTIL_KAONLT/DB/CUTS" % REPLAYPATH
sys.path.insert(0, '%s/UTIL_KAONLT/bin/python/' % REPLAYPATH)
import kaonlt as klt

print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER[1], HOST[1], REPLAYPATH))
#rootName = "%s/UTIL_KAONLT/ROOTfiles/Proton_coin_replay_production_%s_%s.root" % (REPLAYPATH, runNum, MaxEvent)
rootName = "%s/UTIL_KAONLT/ROOTfiles/%s_%s_%s.root" % (REPLAYPATH, ROOTPrefix, runNum, MaxEvent)
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
# Relevant branches now stored as NP arrays

r = klt.pyRoot()
fout = '%s/UTIL_KAONLT/DB/CUTS/run_type/coinpeak.cuts' % REPLAYPATH
# read in cuts file and make dictionary
c = klt.pyPlot(REPLAYPATH)
readDict = c.read_dict(fout,runNum)
# This method calls several methods in kaonlt package. It is required to create properly formated
# dictionaries. The evaluation must be in the analysis script because the analysis variables (i.e. the
# leaves of interest) are not defined in the kaonlt package. This makes the system more flexible
# overall, but a bit more cumbersome in the analysis script. Perhaps one day a better solution will be
# implimented.
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

cutDict = make_cutDict("coin_epi_cut_all")
cutDict = make_cutDict("coin_ek_cut_all", cutDict)
cutDict = make_cutDict("coin_ep_cut_all", cutDict)
cutDict = make_cutDict("coin_epi_cut_prompt", cutDict)
cutDict = make_cutDict("coin_ek_cut_prompt", cutDict)
cutDict = make_cutDict("coin_ep_cut_prompt", cutDict)
cutDict = make_cutDict("coin_epi_cut_peak_only", cutDict)
cutDict = make_cutDict("coin_ek_cut_peak_only", cutDict)
cutDict = make_cutDict("coin_ep_cut_peak_only", cutDict)
c = klt.pyPlot(REPLAYPATH,cutDict)

def coin_events(): 
    # Define the array of arrays containing the relevant HMS and SHMS info
    All_Events_Uncut_tmp = [H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_cal_etotnorm, H_cer_npeSum, CTime_ePiCoinTime_ROC1, CTime_eKCoinTime_ROC1, CTime_epCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_cal_etotnorm, P_aero_npeSum, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer]
    All_Events_Uncut = [(HBeta, Hxp, Hyp, Hdel, HCal, HCer, CTPi, CTK, CTp, KBeta, Kxp, Kyp, KP, KDel, KCal, KAero, KHGC, KHGCX, KHGCY) for (HBeta, Hxp, Hyp, Hdel, HCal, HCer, CTPi, CTK, CTp, KBeta, Kxp, Kyp, KP, KDel, KCal, KAero, KHGC, KHGCX, KHGCY) in zip(*All_Events_Uncut_tmp)] 

    # Create array of arrays of pions and kaons after cuts, all events, prompt and random
    Pion_Events_All_tmp =[]
    Kaon_Events_All_tmp = []
    Proton_Events_All_tmp = []

    # Go over every array in All_Events_Uncut_tmp, append to the other arrays the array after a cut is applied
    for arr in All_Events_Uncut_tmp:
        Pion_Events_All_tmp.append(c.add_cut(arr, "coin_epi_cut_all")) # Apply PID but no cointime cut
        Kaon_Events_All_tmp.append(c.add_cut(arr, "coin_ek_cut_all")) # Apply PID but no cointime cut
        Proton_Events_All_tmp.append(c.add_cut(arr, "coin_ep_cut_all")) # Apply PID but no cointime cut

    Prompt_Events_piCT_Only = np.array(c.add_cut(CTime_ePiCoinTime_ROC1, "coin_epi_cut_peak_only")) # Apply only the prompt CT cut
    Prompt_Pion_Events_piCT_Only = np.array(c.add_cut(CTime_ePiCoinTime_ROC1, "coin_epi_cut_prompt")) # Apply PID and cointime cuts
    Prompt_Events_kCT_Only = np.array(c.add_cut(CTime_eKCoinTime_ROC1, "coin_ek_cut_peak_only")) # Apply only the prompt CT cut
    Prompt_Kaon_Events_kCT_Only = np.array(c.add_cut(CTime_eKCoinTime_ROC1, "coin_ek_cut_prompt")) # Apply PID and cointime cuts
    Prompt_Events_pCT_Only = np.array(c.add_cut(CTime_epCoinTime_ROC1, "coin_ep_cut_peak_only")) # Apply only the prompt CT cut
    Prompt_Proton_Events_pCT_Only = np.array(c.add_cut(CTime_epCoinTime_ROC1, "coin_ep_cut_prompt")) # Apply PID and cointime cuts
      
    Pion_Events_All = [(HBeta, Hxp, Hyp, Hdel, HCal, HCer, CTPi, CTK, CTp, KBeta, Kxp, Kyp, KP, KDel, KCal, KAero, KHGC, KHGCX, KHGCY) for (HBeta, Hxp, Hyp, Hdel, HCal, HCer, CTPi, CTK, CTp, KBeta, Kxp, Kyp, KP, KDel, KCal, KAero, KHGC, KHGCX, KHGCY) in zip(*Pion_Events_All_tmp)] 
    Kaon_Events_All = [(HBeta, Hxp, Hyp, Hdel, HCal, HCer, CTPi, CTK, CTp, KBeta, Kxp, Kyp, KP, KDel, KCal, KAero, KHGC, KHGCX, KHGCY) for (HBeta, Hxp, Hyp, Hdel, HCal, HCer, CTPi, CTK, CTp, KBeta, Kxp, Kyp, KP, KDel, KCal, KAero, KHGC, KHGCX, KHGCY) in zip(*Kaon_Events_All_tmp)]
    Proton_Events_All = [(HBeta, Hxp, Hyp, Hdel, HCal, HCer, CTPi, CTK, CTp, KBeta, Kxp, Kyp, KP, KDel, KCal, KAero, KHGC, KHGCX, KHGCY) for (HBeta, Hxp, Hyp, Hdel, HCal, HCer, CTPi, CTK, CTp, KBeta, Kxp, Kyp, KP, KDel, KCal, KAero, KHGC, KHGCX, KHGCY) in zip(*Proton_Events_All_tmp)]

    COIN_EventInfo = {
        "All_Events" : All_Events_Uncut,
        "Prompt_piCT_Only" : Prompt_Events_piCT_Only,
        "Prompt_kCT_Only" : Prompt_Events_kCT_Only,
        "Prompt_pCT_Only" : Prompt_Events_pCT_Only,
        "Pions_All" : Pion_Events_All,
        "Kaons_All" : Kaon_Events_All,
        "Protons_All" : Proton_Events_All,
        "Prompt_Pions_piCT_Only" : Prompt_Pion_Events_piCT_Only,
        "Prompt_Kaons_kCT_Only" : Prompt_Kaon_Events_kCT_Only,
        "Prompt_Protons_pCT_Only" : Prompt_Proton_Events_pCT_Only,
        }

    return COIN_EventInfo

def main():
    COIN_Data = coin_events()

    COIN_piCT_Only_Header = ["CTime_eKCoinTime_ROC1"]
    COIN_kCT_Only_Header = ["CTime_eKCoinTime_ROC1"]
    COIN_pCT_Only_Header = ["CTime_epCoinTime_ROC1"]
    COIN_All_Data_Header = ["H_gtr_beta","H_gtr_xp","H_gtr_yp","H_gtr_dp","H_cal_etotnorm","H_cer_npeSum","CTime_ePiCoinTime_ROC1","CTime_eKCoinTime_ROC1","CTime_epCoinTime_ROC1","P_gtr_beta","P_gtr_xp","P_gtr_yp","P_gtr_p","P_gtr_dp","P_cal_etotnorm","P_aero_npeSum","P_hgcer_npeSum","P_hgcer_xAtCer","P_hgcer_yAtCer"]

    data_keys = list(COIN_Data.keys()) # Create a list of all the keys in all dicts added above, each is an array of data
    #print(data_keys)

    for i in range (0, len(data_keys)):
        if("piCT" in data_keys[i]):
            DFHeader=list(COIN_piCT_Only_Header)
        elif("kCT" in data_keys[i]):
            DFHeader=list(COIN_kCT_Only_Header)
        elif("pCT" in data_keys[i]):
            DFHeader=list(COIN_pCT_Only_Header)
        else:
            DFHeader=list(COIN_All_Data_Header)
        if (i == 0):
            pd.DataFrame(COIN_Data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_CTPeak_Data.root" % (OUTPATH, runNum, MaxEvent), key ="%s" % data_keys[i])
        elif (i != 0):
            pd.DataFrame(COIN_Data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_CTPeak_Data.root" % (OUTPATH, runNum, MaxEvent), key ="%s" % data_keys[i], mode ='a') 
                    
if __name__ == '__main__':
    main()
