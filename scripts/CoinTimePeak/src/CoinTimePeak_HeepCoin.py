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

# Check the number of arguments provided to the script
if len(sys.argv)-1!=3:
    print("!!!!! ERROR !!!!!\n Expected 3 arguments\n Usage is with - ROOTfilePrefix RunNumber MaxEvents \n!!!!! ERROR !!!!!")
    sys.exit(1)
# Input params - run number and max number of events
ROOTPrefix = sys.argv[1]
runNum = sys.argv[2]
MaxEvent = sys.argv[3]

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

# Add more path setting as needed in a similar manner
OUTPATH = "%s/OUTPUT/Analysis/%sLT" % (UTILPATH,ANATYPE)
CUTPATH = "%s/DB/CUTS" % UTILPATH

################################################################################################################################################
'''
Check that root/output paths and files exist for use
'''

# Construct the name of the rootfile based upon the info we provided
rootName = "%s/ROOTfiles/Analysis/${ANATYPE}LT/%s_%s_%s.root" % (UTILPATH, ROOTPrefix, runNum, MaxEvent)
print ("Attempting to process %s" %(rootName))
lt.SetPath(os.path.realpath(__file__)).checkDir(OUTPATH)
lt.SetPath(os.path.realpath(__file__)).checkFile(rootName)
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

################################################################################################################################################
'''
Define and set up cuts
'''

fout = '%s/DB/CUTS/run_type/coinpeak.cuts' % UTILPATH

cuts = ["coin_ep_cut_all"]

# read in cuts file and make dictionary
def make_cutDict(cuts,fout,runNum,CURRENT_ENV,DEBUG=False):
    '''
    This method calls several methods in kaonlt package. It is required to create properly formated
    dictionaries. The evaluation must be in the analysis script because the analysis variables (i.e. the
    leaves of interest) are not defined in the kaonlt package. This makes the system more flexible
    overall, but a bit more cumbersome in the analysis script. Perhaps one day a better solution will be
    implimented.
    '''

    # read in cuts file and make dictionary
    importDict = lt.SetCuts(CURRENT_ENV).importDict(cuts,fout,runNum,False)
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
