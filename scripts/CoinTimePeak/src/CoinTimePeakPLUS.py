#! /usr/bin/python

# 25/3/21 - Mostly Stephen Kay, University of Regina. Additions by Jacob Murphy
# Script to extract the pion and kaon cointime peak from the data and save info as a new rootfile, subsequent script fits the peak to extract the position to examine the stability
# Additions are to observe effects of new optical matrix for 6.59 GeV/c into HMS
# These mostly include adding in other variables to examine, as well as missing mass

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
    REPLAYPATH = "/group/c-pionlt/USERS/%s/hallc_replay_lt" % USER[1]
elif ("qcd" in HOST[1]):
    REPLAYPATH = "/group/c-pionlt/USERS/%s/hallc_replay_lt" % USER[1]
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
# Relevant branches now stored as NP arrays

#Begin Additions:

ztarvec = e_tree.array("H.tr.vz")
ztar = e_tree.array("H.react.z")

xfp = e_tree.array("H.dc.x_fp")
xpfp = e_tree.array("H.dc.xp_fp")
yfp = e_tree.array("H.dc.y_fp")
ypfp = e_tree.array("H.dc.yp_fp")
ytar = e_tree.array("H.gtr.y")
xtar = e_tree.array("H.gtr.x")

emiss = e_tree.array("P.kin.secondary.emiss")
pmiss = e_tree.array("P.kin.secondary.pmiss")

# Define particle masses for missing mass calculation
Mp = 0.93828
MPi = 0.13957018 
MK = 0.493677
# Calculate missing mass under different particle assumptions
# NOTE, this should be modified for NON different particle assumptions in hcana, this assumes a "kaon" is specified in the kinematics file
MMPi = np.array([math.sqrt(abs(((em+(math.sqrt((MK*MK)+(gtrp*gtrp)))-(math.sqrt((MPi*MPi)+(gtrp*gtrp))))**2)-(pm*pm))) for (em, pm, gtrp) in zip(emiss, pmiss, P_gtr_p)])
MMK = np.array([math.sqrt(abs((em*em)-(pm*pm))) for (em, pm) in zip(emiss, pmiss)])
MMp = np.array([math.sqrt(abs(((em+(math.sqrt((MK*MK)+(gtrp*gtrp)))-(math.sqrt((Mp*Mp)+(gtrp*gtrp))))**2)-(pm*pm))) for (em, pm, gtrp) in zip(emiss, pmiss, P_gtr_p)])



#End Additions

r = klt.pyRoot()
fout = '%s/UTIL_PION/DB/CUTS/run_type/coinpeak.cuts' % REPLAYPATH
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
c = klt.pyPlot(REPLAYPATH,cutDict)

def coin_events(): 
    # Define the array of arrays containing the relevant HMS and SHMS info
    All_Events_Uncut_tmp = [H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_cal_etotnorm, H_cer_npeSum, CTime_ePiCoinTime_ROC1, CTime_eKCoinTime_ROC1, CTime_epCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_cal_etotnorm, P_aero_npeSum, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, ztar, emiss, pmiss, MMPi, MMK, MMp, xfp, xpfp, yfp, ypfp, ytar, xtar]
    All_Events_Uncut = [(HBeta, Hxp, Hyp, Hdel, HCal, HCer, CTPi, CTK, CTp, KBeta, Kxp, Kyp, KP, KDel, KCal, KAero, KHGC, KHGCX, KHGCY, Kzt, Kem, Kpm, mm1, mm2, mm3, Kxfp, Kxpfp, Kyfp, Kypfp, Kytar, Kxtar) for (HBeta, Hxp, Hyp, Hdel, HCal, HCer, CTPi, CTK, CTp, KBeta, Kxp, Kyp, KP, KDel, KCal, KAero, KHGC, KHGCX, KHGCY, Kzt, Kem, Kpm, mm1, mm2, mm3, Kxfp, Kxpfp, Kyfp, Kypfp, Kytar, Kxtar) in zip(*All_Events_Uncut_tmp)] 

    # Create array of arrays of pions and kaons after cuts, all events, prompt and random
    Pion_Events_All_tmp =[]
    Kaon_Events_All_tmp = []
    Proton_Events_All_tmp = []

    # Go over every array in All_Events_Uncut_tmp, append to the other arrays the array after a cut is applied
    for arr in All_Events_Uncut_tmp:
        Pion_Events_All_tmp.append(c.add_cut(arr, "coin_epi_cut_all")) # Apply PID but no cointime cut
        Kaon_Events_All_tmp.append(c.add_cut(arr, "coin_ek_cut_all")) # Apply PID but no cointime cut
        Proton_Events_All_tmp.append(c.add_cut(arr, "coin_ep_cut_all")) # Apply PID but no cointime cut
      
    Pion_Events_All = [(HBeta, Hxp, Hyp, Hdel, HCal, HCer, CTPi, CTK, CTp, KBeta, Kxp, Kyp, KP, KDel, KCal, KAero, KHGC, KHGCX, KHGCY, Kzt, Kem, Kpm, mm1, mm2, mm3, Kxfp, Kxpfp, Kyfp, Kypfp, Kytar, Kxtar) for (HBeta, Hxp, Hyp, Hdel, HCal, HCer, CTPi, CTK, CTp, KBeta, Kxp, Kyp, KP, KDel, KCal, KAero, KHGC, KHGCX, KHGCY, Kzt, Kem, Kpm, mm1, mm2, mm3, Kxfp, Kxpfp, Kyfp, Kypfp, Kytar, Kxtar) in zip(*Pion_Events_All_tmp)] 
    Kaon_Events_All = [(HBeta, Hxp, Hyp, Hdel, HCal, HCer, CTPi, CTK, CTp, KBeta, Kxp, Kyp, KP, KDel, KCal, KAero, KHGC, KHGCX, KHGCY, Kzt, Kem, Kpm, mm1, mm2, mm3, Kxfp, Kxpfp, Kyfp, Kypfp, Kytar, Kxtar) for (HBeta, Hxp, Hyp, Hdel, HCal, HCer, CTPi, CTK, CTp, KBeta, Kxp, Kyp, KP, KDel, KCal, KAero, KHGC, KHGCX, KHGCY, Kzt, Kem, Kpm, mm1, mm2, mm3, Kxfp, Kxpfp, Kyfp, Kypfp, Kytar, Kxtar) in zip(*Kaon_Events_All_tmp)]
    Proton_Events_All = [(HBeta, Hxp, Hyp, Hdel, HCal, HCer, CTPi, CTK, CTp, KBeta, Kxp, Kyp, KP, KDel, KCal, KAero, KHGC, KHGCX, KHGCY, Kzt, Kem, Kpm, mm1, mm2, mm3, Kxfp, Kxpfp, Kyfp, Kypfp, Kytar, Kxtar) for (HBeta, Hxp, Hyp, Hdel, HCal, HCer, CTPi, CTK, CTp, KBeta, Kxp, Kyp, KP, KDel, KCal, KAero, KHGC, KHGCX, KHGCY, Kzt, Kem, Kpm, mm1, mm2, mm3, Kxfp, Kxpfp, Kyfp, Kypfp, Kytar, Kxtar) in zip(*Proton_Events_All_tmp)]

    COIN_EventInfo = {
        "All_Events" : All_Events_Uncut,
        "Pions_All" : Pion_Events_All,
        "Kaons_All" : Kaon_Events_All,
        "Protons_All" : Proton_Events_All,
        }

    return COIN_EventInfo

def main():
    COIN_Data = coin_events()

    COIN_All_Data_Header = ["H_gtr_beta","H_gtr_xp","H_gtr_yp","H_gtr_dp","H_cal_etotnorm","H_cer_npeSum","CTime_ePiCoinTime_ROC1","CTime_eKCoinTime_ROC1","CTime_epCoinTime_ROC1","P_gtr_beta","P_gtr_xp","P_gtr_yp","P_gtr_p","P_gtr_dp","P_cal_etotnorm","P_aero_npeSum","P_hgcer_npeSum","P_hgcer_xAtCer","P_hgcer_yAtCer","ztar","emiss","pmiss","MMPi","MMK","MMp", "xfp", "xpfp", "yfp", "ypfp", "ytar","xtar"]

    data_keys = list(COIN_Data.keys()) # Create a list of all the keys in all dicts added above, each is an array of data
    #print(data_keys)

    for i in range (0, len(data_keys)):
        DFHeader=list(COIN_All_Data_Header)
        if (i == 0):
            pd.DataFrame(COIN_Data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_CTPeak_Data.root" % (OUTPATH, runNum, MaxEvent), key ="%s" % data_keys[i])
        elif (i != 0):
            pd.DataFrame(COIN_Data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_CTPeak_Data.root" % (OUTPATH, runNum, MaxEvent), key ="%s" % data_keys[i], mode ='a') 
                    
if __name__ == '__main__':
    main()
