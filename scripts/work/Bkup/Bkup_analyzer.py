#! /usr/bin/python

# 15/01/21 - Stephen Kay, University of Regina
# 21/06/21 - Muhammad Junaid, University of Regina, Canada

# Python version of the pion analysis script. Now utilises uproot to select event of each type and writes them to a root file
# Intention is to apply PID/selection cutting here and plot in a separate script
# Python should allow for easier reading of databases storing timing offsets e.t.c.
# 27/04/21 - Updated to use new hcana variables, old determinations removed

###################################################################################################################################################

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

##################################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=3:
    print("!!!!! ERROR !!!!!\n Expected 3 arguments\n Usage is with - ROOTfilePrefix RunNumber MaxEvents \n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################

# Input params - run number and max number of events
ROOTPrefix = sys.argv[1]
runNum = sys.argv[2]
MaxEvent = sys.argv[3]

USER = subprocess.getstatusoutput("whoami") # Grab user info for file finding
HOST = subprocess.getstatusoutput("hostname")
if ("farm" in HOST[1]):
#    REPLAYPATH = "/group/c-pionlt/USERS/%s/hallc_replay_lt" % USER[1]
    REPLAYPATH = "/group/c-kaonlt/USERS/%s/hallc_replay_lt" % USER[1]

elif ("qcd" in HOST[1]):
#    REPLAYPATH = "/group/c-pionlt/USERS/%s/hallc_replay_lt" % USER[1]
    REPLAYPATH = "/group/c-kaonlt/USERS/%s/hallc_replay_lt" % USER[1]

elif ("phys.uregina" in HOST[1]):
    REPLAYPATH = "/home/%s/work/JLab/hallc_replay_lt" % USER[1]

elif("skynet" in HOST[1]):
    REPLAYPATH = "/home/%s/Work/JLab/hallc_replay_lt" % USER[1]

################################################################################################################################################

# Add more path setting as needed in a similar manner
OUTPATH = "%s/UTIL_PION/scripts/demo/OUTPUT/Analysis/PionLT" % REPLAYPATH        # Output folder location
CUTPATH = "%s/UTIL_PION/DB/CUTS" % REPLAYPATH
sys.path.insert(0, '%s/UTIL_PION/bin/python/' % REPLAYPATH)

import kaonlt as klt # Import kaonlt module, need the path setting line above prior to importing this

print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER[1], HOST[1], REPLAYPATH))

#################################################################################################################################################

# Construct the name of the rootfile based upon the info we provided
rootName = "%s/ROOTfiles/Analysis/General/%s_%s_%s.root" % (REPLAYPATH, ROOTPrefix, runNum, MaxEvent)     # Input file location and variables taking
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
    print("Output path not found, please make a sym link or directory called OUTPUT in UTIL_PION/scripts/demo to store output")
    sys.exit(3)
print ("Attempting to process %s" %(rootName))
if os.path.isfile(rootName):
    print ("%s exists, attempting to process" % (rootName))
else:
    print ("%s not found - do you have the correct sym link/folder set up?" % (rootName))
    sys.exit(4)
print("Output path checks out, outputting to %s" % (OUTPATH))

###############################################################################################################################################

# Read stuff from the main event tree
e_tree = up.open(rootName)["T"]

# Timing info
CTime_ePiCoinTime_ROC1 = e_tree.array("CTime.ePiCoinTime_ROC1")  #
CTime_eKCoinTime_ROC1 = e_tree.array("CTime.eKCoinTime_ROC1")    #
CTime_epCoinTime_ROC1 = e_tree.array("CTime.epCoinTime_ROC1")    #
#P_RF_tdcTime = e_tree.array("T.coin.pRF_tdcTime")               #
#P_hod_fpHitsTime = e_tree.array("P.hod.fpHitsTime")             #
H_RF_Dist = e_tree.array("RFTime.HMS_RFtimeDist")            #
P_RF_Dist = e_tree.array("RFTime.SHMS_RFtimeDist")           #

# HMS info
H_hod_goodscinhit = e_tree.array("H.hod.goodscinhit")            #
H_hod_goodstarttime = e_tree.array("H.hod.goodstarttime")        #
H_gtr_beta = e_tree.array("H.gtr.beta")                          # Beta is velocity of particle between pairs of hodoscopes
H_gtr_xp = e_tree.array("H.gtr.th")                              # xpfp -> Theta
H_gtr_yp = e_tree.array("H.gtr.ph")                              # ypfp -> Phi
H_gtr_dp = e_tree.array("H.gtr.dp")                              # dp is Delta
H_cal_etotnorm = e_tree.array("H.cal.etotnorm")                  #
H_cal_etottracknorm = e_tree.array("H.cal.etottracknorm")        #
H_cer_npeSum = e_tree.array("H.cer.npeSum")                      #

# SHMS info
P_hod_goodscinhit = e_tree.array("P.hod.goodscinhit")            #
P_hod_goodstarttime = e_tree.array("P.hod.goodstarttime")        #
P_gtr_beta = e_tree.array("P.gtr.beta")                          # Beta is velocity of particle between pairs of hodoscopes
P_gtr_xp = e_tree.array("P.gtr.th")                              # xpfp -> Theta
P_gtr_yp = e_tree.array("P.gtr.ph")                              # ypfp -> Phi
P_gtr_p = e_tree.array("P.gtr.p")                                #
P_gtr_dp = e_tree.array("P.gtr.dp")                              # dp is Delta 
P_cal_etotnorm = e_tree.array("P.cal.etotnorm")                  #
P_cal_etottracknorm = e_tree.array("P.cal.etottracknorm")        #
P_aero_npeSum = e_tree.array("P.aero.npeSum")                    #
P_aero_xAtAero = e_tree.array("P.aero.xAtAero")                  #
P_aero_yAtAero = e_tree.array("P.aero.yAtAero")                  #
P_hgcer_npeSum = e_tree.array("P.hgcer.npeSum")                  #
P_hgcer_xAtCer = e_tree.array("P.hgcer.xAtCer")                  #
P_hgcer_yAtCer = e_tree.array("P.hgcer.yAtCer")                  #

# Kinematic quantitites
Q2 = e_tree.array("H.kin.primary.Q2")                            #
W = e_tree.array("H.kin.primary.W")                              #
epsilon = e_tree.array("H.kin.primary.epsilon")                  #
#ph_q = e_tree.array("P.kin.secondary.ph_xq")                    #
#emiss = e_tree.array("P.kin.secondary.emiss")                   #
#pmiss = e_tree.array("P.kin.secondary.pmiss")                   #
MMpi = e_tree.array("P.kin.secondary.MMpi")                      #
MMK = e_tree.array("P.kin.secondary.MMK")                        #
MMp = e_tree.array("P.kin.secondary.MMp")                        #
#MandelT = e_tree.array("P.kin.secondary.MandelT")               #
#MandelU = e_tree.array("P.kin.secondary.MandelU")               #

# Misc quantities
#fEvtType = e_tree.array("fEvtHdr.fEvtType")                     #
#RFFreq = e_tree.array("MOFC1FREQ")                              #
#RFFreqDiff = e_tree.array("MOFC1DELTA")                         #
#pEDTM = e_tree.array("T.coin.pEDTM_tdcTime")                    #
# Relevant branches now stored as NP arrays

##############################################################################################################################################

r = klt.pyRoot()
#fout = '%s/UTIL_PION/DB/CUTS/run_type/coin_prod.cuts' % REPLAYPATH
fout = '%s/UTIL_PION/DB/CUTS/run_type/analyzer.cuts' % REPLAYPATH

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

#################################################################################################################################################################

cutDict = make_cutDict("coin_epi_cut_all_RF")
cutDict = make_cutDict("coin_epi_cut_prompt_RF", cutDict)
cutDict = make_cutDict("coin_epi_cut_rand_RF", cutDict)
cutDict = make_cutDict("coin_ek_cut_all_RF", cutDict)
cutDict = make_cutDict("coin_ek_cut_prompt_RF", cutDict)
cutDict = make_cutDict("coin_ek_cut_rand_RF", cutDict)
cutDict = make_cutDict("coin_ep_cut_all_RF", cutDict)
cutDict = make_cutDict("coin_ep_cut_prompt_RF", cutDict)
cutDict = make_cutDict("coin_ep_cut_rand_RF", cutDict)

c = klt.pyPlot(REPLAYPATH,cutDict)

####################################################################################################################################################################################################

def coin_pions(): 

    # Define the array of arrays containing the relevant HMS and SHMS info
    NoCut_COIN_Pions = [H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_ePiCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMpi, H_RF_Dist, P_RF_Dist, Q2, W, epsilon]
    Uncut_COIN_Pions = [(H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_ePiCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMpi, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) for (H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_ePiCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMpi, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) in zip(*NoCut_COIN_Pions)] 

    # Create array of arrays of pions after cuts, all events, prompt and random
    Cut_COIN_Pions_tmp = NoCut_COIN_Pions
    Cut_COIN_Pions_all_tmp = []
    Cut_COIN_Pions_prompt_tmp = []
    Cut_COIN_Pions_rand_tmp = []

    for arr in Cut_COIN_Pions_tmp:
        Cut_COIN_Pions_all_tmp.append(c.add_cut(arr, "coin_epi_cut_all_RF"))
        Cut_COIN_Pions_prompt_tmp.append(c.add_cut(arr, "coin_epi_cut_prompt_RF"))
        Cut_COIN_Pions_rand_tmp.append(c.add_cut(arr, "coin_epi_cut_rand_RF"))

    Cut_COIN_Pions_all = [(H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_ePiCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMpi, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) for (H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_ePiCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMpi, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) in zip(*Cut_COIN_Pions_all_tmp)
	if H_hod_goodscinhit == 1
	if H_hod_goodstarttime == 1
	if P_hod_goodscinhit == 1
	if P_hod_goodstarttime == 1
	]

    Cut_COIN_Pions_prompt = [(H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_ePiCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMpi, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) for (H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_ePiCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMpi, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) in zip(*Cut_COIN_Pions_prompt_tmp)
        if H_hod_goodscinhit == 1
        if H_hod_goodstarttime == 1
        if P_hod_goodscinhit == 1
        if P_hod_goodstarttime == 1
        ]

    Cut_COIN_Pions_random = [(H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_ePiCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMpi, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) for (H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_ePiCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMpi, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) in zip(*Cut_COIN_Pions_rand_tmp)
        if H_hod_goodscinhit == 1
        if H_hod_goodstarttime == 1
        if P_hod_goodscinhit == 1
        if P_hod_goodstarttime == 1
        ]

    COIN_Pions = {
        "Uncut_Pion_Events" : Uncut_COIN_Pions,
        "Cut_Pion_Events_All" : Cut_COIN_Pions_all,
        "Cut_Pion_Events_Prompt" : Cut_COIN_Pions_prompt,
        "Cut_Pion_Events_Random" : Cut_COIN_Pions_random,
        }

    return COIN_Pions

########################################################################################################################################################################################

def coin_kaons(): 

    # Define the array of arrays containing the relevant HMS and SHMS info
    NoCut_COIN_Kaons = [H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_eKCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMK, H_RF_Dist, P_RF_Dist, Q2, W, epsilon]
    Uncut_COIN_Kaons = [(H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_eKCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMK, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) for (H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_eKCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMK, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) in zip(*NoCut_COIN_Kaons)] 

    # Create array of arrays of pions after cuts, all events, prompt and random
    Cut_COIN_Kaons_tmp = NoCut_COIN_Kaons
    Cut_COIN_Kaons_all_tmp = []
    Cut_COIN_Kaons_prompt_tmp = []
    Cut_COIN_Kaons_rand_tmp = []

    for arr in Cut_COIN_Kaons_tmp:
        Cut_COIN_Kaons_all_tmp.append(c.add_cut(arr, "coin_ek_cut_all_RF"))
        Cut_COIN_Kaons_prompt_tmp.append(c.add_cut(arr, "coin_ek_cut_prompt_RF"))
        Cut_COIN_Kaons_rand_tmp.append(c.add_cut(arr, "coin_ek_cut_rand_RF"))

    Cut_COIN_Kaons_all = [(H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_eKCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMK, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) for (H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_eKCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMK, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) in zip(*Cut_COIN_Kaons_all_tmp)
        if H_hod_goodscinhit == 1
        if H_hod_goodstarttime == 1
        if P_hod_goodscinhit == 1
        if P_hod_goodstarttime == 1
        ]

    Cut_COIN_Kaons_prompt = [(H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_eKCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMK, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) for (H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_eKCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMK, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) in zip(*Cut_COIN_Kaons_prompt_tmp)
        if H_hod_goodscinhit == 1
        if H_hod_goodstarttime == 1
        if P_hod_goodscinhit == 1
        if P_hod_goodstarttime == 1
        ]

    Cut_COIN_Kaons_random = [(H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_eKCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMK, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) for (H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_eKCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMK, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) in zip(*Cut_COIN_Kaons_rand_tmp)
        if H_hod_goodscinhit == 1
        if H_hod_goodstarttime == 1
        if P_hod_goodscinhit == 1
        if P_hod_goodstarttime == 1
        ]

    COIN_Kaons = {
        "Uncut_Kaon_Events" : Uncut_COIN_Kaons,
        "Cut_Kaon_Events_All" : Cut_COIN_Kaons_all,
        "Cut_Kaon_Events_Prompt" : Cut_COIN_Kaons_prompt,
        "Cut_Kaon_Events_Random" : Cut_COIN_Kaons_random,
        }

    return COIN_Kaons

#################################################################################################################################################################################################

def coin_protons(): 

    # Define the array of arrays containing the relevant HMS and SHMS info
    NoCut_COIN_Protons = [H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_epCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMp, H_RF_Dist, P_RF_Dist, Q2, W, epsilon]
    Uncut_COIN_Protons = [(H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_epCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMp, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) for (H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_epCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMp, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) in zip(*NoCut_COIN_Protons)] 

    # Create array of arrays of pions after cuts, all events, prompt and random
    Cut_COIN_Protons_tmp = NoCut_COIN_Protons
    Cut_COIN_Protons_all_tmp = []
    Cut_COIN_Protons_prompt_tmp =[]
    Cut_COIN_Protons_rand_tmp =[]

    for arr in Cut_COIN_Protons_tmp:
        Cut_COIN_Protons_all_tmp.append(c.add_cut(arr, "coin_ep_cut_all_RF"))
        Cut_COIN_Protons_prompt_tmp.append(c.add_cut(arr, "coin_ep_cut_prompt_RF"))
        Cut_COIN_Protons_rand_tmp.append(c.add_cut(arr, "coin_ep_cut_rand_RF"))

    Cut_COIN_Protons_all = [(H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_epCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMp, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) for (H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_epCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMp, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) in zip(*Cut_COIN_Protons_all_tmp)
        if H_hod_goodscinhit == 1
        if H_hod_goodstarttime == 1
        if P_hod_goodscinhit == 1
        if P_hod_goodstarttime == 1
        ]

    Cut_COIN_Protons_prompt = [(H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_epCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMp, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) for (H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_epCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMp, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) in zip(*Cut_COIN_Protons_prompt_tmp)
        if H_hod_goodscinhit == 1
        if H_hod_goodstarttime == 1
        if P_hod_goodscinhit == 1
        if P_hod_goodstarttime == 1
        ]

    Cut_COIN_Protons_random = [(H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_epCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMp, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) for (H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, CTime_epCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, MMp, H_RF_Dist, P_RF_Dist, Q2, W, epsilon) in zip(*Cut_COIN_Protons_rand_tmp)
        if H_hod_goodscinhit == 1
        if H_hod_goodstarttime == 1
        if P_hod_goodscinhit == 1
        if P_hod_goodstarttime == 1
        ]

    COIN_Protons = {
        "Uncut_Proton_Events" : Uncut_COIN_Protons,
        "Cut_Proton_Events_All" : Cut_COIN_Protons_all,
        "Cut_Proton_Events_Prompt" : Cut_COIN_Protons_prompt,
        "Cut_Proton_Events_Random" : Cut_COIN_Protons_random,
        }

    return COIN_Protons

###################################################################################################################################################################################################

def main():
    COIN_Pion_Data = coin_pions()
    COIN_Kaon_Data = coin_kaons()
    COIN_Proton_Data = coin_protons()

    # This is just the list of branches we use from the initial root file for each dict
    # I don't like re-defining this here as it's very prone to errors if you included (or removed something) earlier but didn't modify it here
    # Should base the branches to include based on some list and just repeat the list here (or call it again directly below)

    COIN_Pion_Data_Header = ["H_gtr_beta","H_gtr_xp","H_gtr_yp","H_gtr_dp","H_hod_goodscinhit","H_hod_goodstarttime","H_cal_etotnorm","H_cal_etottracknorm","H_cer_npeSum","CTime_ePiCoinTime_ROC1","P_gtr_beta","P_gtr_xp","P_gtr_yp","P_gtr_p","P_gtr_dp","P_hod_goodscinhit","P_hod_goodstarttime","P_cal_etotnorm","P_cal_etottracknorm","P_aero_npeSum","P_aero_xAtAero","P_aero_yAtAero","P_hgcer_npeSum","P_hgcer_xAtCer","P_hgcer_yAtCer","MMpi","H_RF_Dist","P_RF_Dist", "Q2", "W", "epsilon"]
    COIN_Kaon_Data_Header = ["H_gtr_beta","H_gtr_xp","H_gtr_yp","H_gtr_dp","H_hod_goodscinhit","H_hod_goodstarttime","H_cal_etotnorm","H_cal_etottracknorm","H_cer_npeSum","CTime_eKCoinTime_ROC1","P_gtr_beta","P_gtr_xp","P_gtr_yp","P_gtr_p","P_gtr_dp","P_hod_goodscinhit","P_hod_goodstarttime","P_cal_etotnorm","P_cal_etottracknorm","P_aero_npeSum","P_aero_xAtAero","P_aero_yAtAero","P_hgcer_npeSum","P_hgcer_xAtCer","P_hgcer_yAtCer","MMK","H_RF_Dist","P_RF_Dist", "Q2", "W", "epsilon"]
    COIN_Proton_Data_Header = ["H_gtr_beta","H_gtr_xp","H_gtr_yp","H_gtr_dp","H_hod_goodscinhit","H_hod_goodstarttime","H_cal_etotnorm","H_cal_etottracknorm","H_cer_npeSum","CTime_epCoinTime_ROC1","P_gtr_beta","P_gtr_xp","P_gtr_yp","P_gtr_p","P_gtr_dp","P_hod_goodscinhit","P_hod_goodstarttime","P_cal_etotnorm","P_cal_etottracknorm","P_aero_npeSum","P_aero_xAtAero","P_aero_yAtAero","P_hgcer_npeSum","P_hgcer_xAtCer","P_hgcer_yAtCer","MMp","H_RF_Dist","P_RF_Dist", "Q2", "W", "epsilon"]

    # Need to create a dict for all the branches we grab
    data = {}

    for d in (COIN_Pion_Data, COIN_Kaon_Data, COIN_Proton_Data): # Convert individual dictionaries into a "dict of dicts"
        data.update(d)
        data_keys = list(data.keys()) # Create a list of all the keys in all dicts added above, each is an array of data

    for i in range (0, len(data_keys)):
        if("Pion" in data_keys[i]):
            DFHeader=list(COIN_Pion_Data_Header)
        elif("Kaon" in data_keys[i]):
            DFHeader=list(COIN_Kaon_Data_Header)
        elif("Proton" in data_keys[i]):
            DFHeader=list(COIN_Proton_Data_Header)
        else:
            continue
            # Uncomment the line below if you want .csv file output, WARNING the files can be very large and take a long time to process!
            #pd.DataFrame(data.get(data_keys[i])).to_csv("%s/%s_%s.csv" % (OUTPATH, data_keys[i], runNum), header=DFHeader, index=False) # Convert array to panda dataframe and write to csv with correct header
        if (i == 0):
            pd.DataFrame(data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_Analysed_Data.root" % (OUTPATH, runNum, MaxEvent), key ="%s" % data_keys[i])
        elif (i != 0):
            pd.DataFrame(data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_Analysed_Data.root" % (OUTPATH, runNum, MaxEvent), key ="%s" % data_keys[i], mode ='a') 
                    
if __name__ == '__main__':
    main()
print ("Processing Complete")
