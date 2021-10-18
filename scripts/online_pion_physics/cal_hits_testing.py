#! /usr/bin/python

# 15/01/21 - Stephen Kay, University of Regina
# 21/06/21 - Edited By - Muhammad Junaid, University of Regina, Canada

# Python version of the pion analysis script. Now utilises uproot to select event of each type and writes them to a root file
# Intention is to apply PID/selection cutting here and plot in a separate script
# Python should allow for easier reading of databases storing timing offsets e.t.c.
# 27/04/21 - Updated to use new hcana variables, old determinations removed
# 05/10/21 - Updated to add in focal plane quantities (Jacob Murphy)
# 12/10/21 - Added in Dipole Exit quantities (Jacob Murphy)
# 13/10/21 - Added in SHMS Calorimeter Block hits (Jacob Murphy)
# 15/10/21 - Added in eplane and earray values to analysis. Will need to replay again to include new branches to get for analysis

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
    REPLAYPATH = "/group/c-pionlt/online_analysis/hallc_replay_lt"
elif ("qcd" in HOST[1]):
    REPLAYPATH = "/group/c-pionlt/USERS/%s/hallc_replay_lt" % USER[1]
elif ("phys.uregina" in HOST[1]):
    REPLAYPATH = "/home/%s/work/JLab/hallc_replay_lt" % USER[1]
elif("skynet" in HOST[1]):
    REPLAYPATH = "/home/%s/Work/JLab/hallc_replay_lt" % USER[1]
elif("cdaq" in HOST[1]):
    REPLAYPATH = "/home/cdaq/hallc-online/hallc_replay_lt"

################################################################################################################################################

# Add more path setting as needed in a similar manner
OUTPATH = "%s/UTIL_PION/OUTPUT/Analysis/PionLT" % REPLAYPATH        # Output folder location
CUTPATH = "%s/UTIL_PION/DB/CUTS" % REPLAYPATH
sys.path.insert(0, '%s/UTIL_PION/bin/python/' % REPLAYPATH)

import kaonlt as klt # Import kaonlt module, need the path setting line above prior to importing this

print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER[1], HOST[1], REPLAYPATH))

#################################################################################################################################################

# Construct the name of the rootfile based upon the info we provided
rootName = "%s/UTIL_PION/ROOTfiles/Analysis/PionLT/%s_%s_%s.root" % (REPLAYPATH, ROOTPrefix, runNum, MaxEvent)     # Input file location and variables taking
#rootName = "/volatile/hallc/c-pionlt/junaid/ROOTfiles/Analysis/PionLT/Pion_coin_replay_production_8076_-1.root" # Hard coded path to a recent file for testing
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

###############################################################################################################################################

# Read stuff from the main event tree
e_tree = up.open(rootName)["T"]

# Timing info
CTime_ePiCoinTime_ROC1 = e_tree.array("CTime.ePiCoinTime_ROC1")  #
CTime_eKCoinTime_ROC1 = e_tree.array("CTime.eKCoinTime_ROC1")    #
CTime_epCoinTime_ROC1 = e_tree.array("CTime.epCoinTime_ROC1")    #
H_RF_Dist = e_tree.array("RFTime.HMS_RFtimeDist")            #
P_RF_Dist = e_tree.array("RFTime.SHMS_RFtimeDist")           #

# HMS info
H_hod_goodscinhit = e_tree.array("H.hod.goodscinhit")            #
H_hod_goodstarttime = e_tree.array("H.hod.goodstarttime")        #
H_gtr_beta = e_tree.array("H.gtr.beta")                          # Beta is velocity of particle between pairs of hodoscopes
H_gtr_xp = e_tree.array("H.gtr.th")                              # xpfp -> Theta
H_gtr_yp = e_tree.array("H.gtr.ph")                              # ypfp -> Phi
H_gtr_dp = e_tree.array("H.gtr.dp")                              # dp is Delta
# JM 05/10/21 Adding in focal plane variables
H_dc_xfp = e_tree.array("H.dc.x_fp")                            # xp is x focal plane, the vertical position in the focal plane
H_dc_xpfp = e_tree.array("H.dc.xp_fp")                          # xpfp is x' focal plane, the vertical angle
H_dc_yfp = e_tree.array("H.dc.y_fp")                            # yp is x focal plane, the vertical position in the focal plane
H_dc_ypfp = e_tree.array("H.dc.yp_fp")                          # ypfp is y' focal plane, the vertical angle
H_cal_etotnorm = e_tree.array("H.cal.etotnorm")                  #
H_cal_etottracknorm = e_tree.array("H.cal.etottracknorm")        #
H_cer_npeSum = e_tree.array("H.cer.npeSum")                      #

# JM 12/10/21 Added in dipole exit variables
H_dc_InsideDipoleExit = e_tree.array("H.dc.InsideDipoleExit")    #
P_dc_InsideDipoleExit = e_tree.array("P.dc.InsideDipoleExit")    #

# SHMS info
P_hod_goodscinhit = e_tree.array("P.hod.goodscinhit")            #
P_hod_goodstarttime = e_tree.array("P.hod.goodstarttime")        #
P_gtr_beta = e_tree.array("P.gtr.beta")                          # Beta is velocity of particle between pairs of hodoscopes
P_gtr_xp = e_tree.array("P.gtr.th")                              # xpfp -> Theta
P_gtr_yp = e_tree.array("P.gtr.ph")                              # ypfp -> Phi
P_gtr_p = e_tree.array("P.gtr.p")                                #
P_gtr_dp = e_tree.array("P.gtr.dp")                              # dp is Delta 
# JM 05/10/21 Added in focal plane variables
P_dc_xfp = e_tree.array("P.dc.x_fp")                             # xp is x focal plane, the vertical position in the focal plane
P_dc_xpfp = e_tree.array("P.dc.xp_fp")                           # xpfp is x' focal plane, the vertical angle
P_dc_yfp = e_tree.array("P.dc.y_fp")                             # yp is x focal plane, the vertical position in the focal plane
P_dc_ypfp = e_tree.array("P.dc.yp_fp")                           # ypfp is y' focal plane, the vertical angle
P_cal_etotnorm = e_tree.array("P.cal.etotnorm")                  #
P_cal_etottracknorm = e_tree.array("P.cal.etottracknorm")        #
P_cal_fly_earray = e_tree.array("P.cal.fly.earray")              #
P_cal_pr_eplane = e_tree.array("P.cal.pr.eplane")                #
# SJDK 13/10/21 - This seems to generate something of type "jagged.array", I suspect it is an indexed array and we can't access it in the same way we access other variables
P_cal_fly_numGoodAdcHits = e_tree.array("P.cal.fly.numGoodAdcHits")# Indexed hits into calorimeter blocks 
# JM 16/10/21 - Attempt at converting hits into calo blocks to numpy array (which works), but it cannot seem to save it back into the root tree, so it skips it in the analysis
# Left uncommented because it seems to only output a skip error for this variable
Cal_Adc_Hits = np.empty([len(P_cal_fly_numGoodAdcHits),224])
for i in range(len(P_cal_fly_numGoodAdcHits)):
    for j in range(224):
        Cal_Adc_Hits[i][j] = P_cal_fly_numGoodAdcHits[i][j]      # This will not save properly and an error outputs in the analysis, but it just skips it so no harm in leaving it in for now until we figure it out 
Cal_pandas = Cal_Adc_Hits
#print(Cal_pandas)
P_aero_npeSum = e_tree.array("P.aero.npeSum")                    #
P_aero_xAtAero = e_tree.array("P.aero.xAtAero")                  #
P_aero_yAtAero = e_tree.array("P.aero.yAtAero")                  #
P_hgcer_npeSum = e_tree.array("P.hgcer.npeSum")                  #
P_hgcer_xAtCer = e_tree.array("P.hgcer.xAtCer")                  #
P_hgcer_yAtCer = e_tree.array("P.hgcer.yAtCer")                  #
P_ngcer_npeSum = e_tree.array("P.ngcer.npeSum")                  #
P_ngcer_xAtCer = e_tree.array("P.ngcer.xAtCer")                  #
P_ngcer_yAtCer = e_tree.array("P.ngcer.yAtCer")                  #

# Kinematic quantitites
Q2 = e_tree.array("H.kin.primary.Q2")                            #
W = e_tree.array("H.kin.primary.W")                              #
epsilon = e_tree.array("H.kin.primary.epsilon")                  #
ph_q = e_tree.array("P.kin.secondary.ph_xq")                     #
MMpi = e_tree.array("P.kin.secondary.MMpi")                      #
MMK = e_tree.array("P.kin.secondary.MMK")                        #
MMp = e_tree.array("P.kin.secondary.MMp")                        #
MandelT = e_tree.array("P.kin.secondary.MandelT")                #
#MandelU = e_tree.array("P.kin.secondary.MandelU")               #

# Misc quantities
#fEvtType = e_tree.array("fEvtHdr.fEvtType")                     #
#RFFreq = e_tree.array("MOFC1FREQ")                              #
#RFFreqDiff = e_tree.array("MOFC1DELTA")                         #
#pEDTM = e_tree.array("T.coin.pEDTM_tdcTime")                    #
# Relevant branches now stored as NP arrays

# Define distances from focal plane (cm)
D_Calo = 292.64
D_Exit = -307.0

# Calculate X and Y Positions along tracks from focal plane

xCalo = np.array([xfp+xpfp*D_Calo for (xfp, xpfp) in zip(P_dc_xfp, P_dc_xpfp)])
yCalo = np.array([yfp+ypfp*D_Calo for (yfp, ypfp) in zip(P_dc_yfp, P_dc_ypfp)])
xExit = np.array([xfp+xpfp*D_Exit for (xfp, xpfp) in zip(P_dc_xfp, P_dc_xpfp)])
yExit = np.array([yfp+ypfp*D_Exit for (yfp, ypfp) in zip(P_dc_yfp, P_dc_ypfp)])

# Unindex Calo Hits


##############################################################################################################################################

# Defining path for cut file
r = klt.pyRoot()
fout = '%s/UTIL_PION/DB/CUTS/run_type/coin_prod.cuts' % REPLAYPATH

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

#################################################################################################################################################################

# defining Cuts
cutDict = make_cutDict("coin_epi_cut_all_noRF")
cutDict = make_cutDict("coin_epi_cut_all_RF", cutDict)
cutDict = make_cutDict("coin_epi_cut_prompt_RF", cutDict)
cutDict = make_cutDict("coin_epi_cut_rand_RF", cutDict)

c = klt.pyPlot(REPLAYPATH,cutDict)

#################################################################################################################################################################

def coin_pions(): 
    # Define the array of arrays containing the relevant HMS and SHMS info
    NoCut_COIN_Pions = [H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_dc_xfp, H_dc_xpfp, H_dc_yfp, H_dc_ypfp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, H_dc_InsideDipoleExit, P_dc_InsideDipoleExit, CTime_ePiCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_dc_xfp, P_dc_xpfp, P_dc_yfp, P_dc_ypfp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_cal_fly_earray, P_cal_pr_eplane, Cal_pandas, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, P_ngcer_npeSum, P_ngcer_xAtCer, P_ngcer_yAtCer, MMpi, H_RF_Dist, P_RF_Dist, xCalo, yCalo, xExit, yExit, Q2, W, epsilon, ph_q, MandelT]
    Uncut_COIN_Pions = [(H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_dc_xfp, H_dc_xpfp, H_dc_yfp, H_dc_ypfp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, H_dc_InsideDipoleExit, P_dc_InsideDipoleExit, CTime_ePiCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_dc_xfp, P_dc_xpfp, P_dc_yfp, P_dc_ypfp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_cal_fly_earray, P_cal_pr_eplane, Cal_pandas, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, P_ngcer_npeSum, P_ngcer_xAtCer, P_ngcer_yAtCer, MMpi, H_RF_Dist, P_RF_Dist, xCalo, yCalo, xExit, yExit, Q2, W, epsilon, ph_q, MandelT) for (H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_dc_xfp, H_dc_xpfp, H_dc_yfp, H_dc_ypfp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, H_dc_InsideDipoleExit, P_dc_InsideDipoleExit, CTime_ePiCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_dc_xfp, P_dc_xpfp, P_dc_yfp, P_dc_ypfp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_cal_fly_earray, P_cal_pr_eplane, Cal_pandas, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, P_ngcer_npeSum, P_ngcer_xAtCer, P_ngcer_yAtCer, MMpi, H_RF_Dist, P_RF_Dist, xCalo, yCalo, xExit, yExit, Q2, W, epsilon, ph_q, MandelT) in zip(*NoCut_COIN_Pions)] 

    # Create array of arrays of pions after cuts, all events, prompt and random
    Cut_COIN_Pions_tmp = NoCut_COIN_Pions
    Cut_COIN_Pions_noRF_tmp = []
    Cut_COIN_Pions_all_tmp = []
    Cut_COIN_Pions_prompt_tmp = []
    Cut_COIN_Pions_rand_tmp = []

    for arr in Cut_COIN_Pions_tmp:
        Cut_COIN_Pions_noRF_tmp.append(c.add_cut(arr, "coin_epi_cut_all_noRF"))
        Cut_COIN_Pions_all_tmp.append(c.add_cut(arr, "coin_epi_cut_all_RF"))
        Cut_COIN_Pions_prompt_tmp.append(c.add_cut(arr, "coin_epi_cut_prompt_RF"))
        Cut_COIN_Pions_rand_tmp.append(c.add_cut(arr, "coin_epi_cut_rand_RF"))
    #print(Cut_COIN_Pions_all_tmp[31].shape)
    test = pd.DataFrame(Cut_COIN_Pions_all_tmp)
    print(type(test[31][31]))
    Cut_COIN_Pions_noRF = [(H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_dc_xfp, H_dc_xpfp, H_dc_yfp, H_dc_ypfp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, H_dc_InsideDipoleExit, P_dc_InsideDipoleExit, CTime_ePiCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_dc_xfp, P_dc_xpfp, P_dc_yfp, P_dc_ypfp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_cal_fly_earray, P_cal_pr_eplane, Cal_pandas, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, P_ngcer_npeSum, P_ngcer_xAtCer, P_ngcer_yAtCer, MMpi, H_RF_Dist, P_RF_Dist, xCalo, yCalo, xExit, yExit, Q2, W, epsilon, ph_q, MandelT) for (H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_dc_xfp, H_dc_xpfp, H_dc_yfp, H_dc_ypfp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, H_dc_InsideDipoleExit, P_dc_InsideDipoleExit, CTime_ePiCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_dc_xfp, P_dc_xpfp, P_dc_yfp, P_dc_ypfp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_cal_fly_earray, P_cal_pr_eplane, Cal_pandas, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, P_ngcer_npeSum, P_ngcer_xAtCer, P_ngcer_yAtCer, MMpi, H_RF_Dist, P_RF_Dist, xCalo, yCalo, xExit, yExit, Q2, W, epsilon, ph_q, MandelT) in zip(*Cut_COIN_Pions_noRF_tmp)
        ]
    Cut_COIN_Pions_all = [(H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_dc_xfp, H_dc_xpfp, H_dc_yfp, H_dc_ypfp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, H_dc_InsideDipoleExit, P_dc_InsideDipoleExit, CTime_ePiCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_dc_xfp, P_dc_xpfp, P_dc_yfp, P_dc_ypfp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_cal_fly_earray, P_cal_pr_eplane, Cal_pandas, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, P_ngcer_npeSum, P_ngcer_xAtCer, P_ngcer_yAtCer, MMpi, H_RF_Dist, P_RF_Dist, xCalo, yCalo, xExit, yExit, Q2, W, epsilon, ph_q, MandelT) for (H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_dc_xfp, H_dc_xpfp, H_dc_yfp, H_dc_ypfp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, H_dc_InsideDipoleExit, P_dc_InsideDipoleExit, CTime_ePiCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_dc_xfp, P_dc_xpfp, P_dc_yfp, P_dc_ypfp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_cal_fly_earray, P_cal_pr_eplane, Cal_pandas, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, P_ngcer_npeSum, P_ngcer_xAtCer, P_ngcer_yAtCer, MMpi, H_RF_Dist, P_RF_Dist, xCalo, yCalo, xExit, yExit, Q2, W, epsilon, ph_q, MandelT) in zip(*Cut_COIN_Pions_all_tmp)
	]
    print(type(pd.DataFrame(Cut_COIN_Pions_all)[31][31]))

    Cut_COIN_Pions_prompt = [(H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_dc_xfp, H_dc_xpfp, H_dc_yfp, H_dc_ypfp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, H_dc_InsideDipoleExit, P_dc_InsideDipoleExit, CTime_ePiCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_dc_xfp, P_dc_xpfp, P_dc_yfp, P_dc_ypfp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_cal_fly_earray, P_cal_pr_eplane, Cal_pandas, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, P_ngcer_npeSum, P_ngcer_xAtCer, P_ngcer_yAtCer, MMpi, H_RF_Dist, P_RF_Dist, xCalo, yCalo, xExit, yExit, Q2, W, epsilon, ph_q, MandelT) for (H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_dc_xfp, H_dc_xpfp, H_dc_yfp, H_dc_ypfp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, H_dc_InsideDipoleExit, P_dc_InsideDipoleExit, CTime_ePiCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_dc_xfp, P_dc_xpfp, P_dc_yfp, P_dc_ypfp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_cal_fly_earray, P_cal_pr_eplane, Cal_pandas, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, P_ngcer_npeSum, P_ngcer_xAtCer, P_ngcer_yAtCer, MMpi, H_RF_Dist, P_RF_Dist, xCalo, yCalo, xExit, yExit, Q2, W, epsilon, ph_q, MandelT) in zip(*Cut_COIN_Pions_prompt_tmp)
        ]

    Cut_COIN_Pions_random = [(H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_dc_xfp, H_dc_xpfp, H_dc_yfp, H_dc_ypfp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, H_dc_InsideDipoleExit, P_dc_InsideDipoleExit, CTime_ePiCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_dc_xfp, P_dc_xpfp, P_dc_yfp, P_dc_ypfp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_cal_fly_earray, P_cal_pr_eplane, Cal_pandas, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, P_ngcer_npeSum, P_ngcer_xAtCer, P_ngcer_yAtCer, MMpi, H_RF_Dist, P_RF_Dist, xCalo, yCalo, xExit, yExit, Q2, W, epsilon, ph_q, MandelT) for (H_gtr_beta, H_gtr_xp, H_gtr_yp, H_gtr_dp, H_dc_xfp, H_dc_xpfp, H_dc_yfp, H_dc_ypfp, H_hod_goodscinhit, H_hod_goodstarttime, H_cal_etotnorm, H_cal_etottracknorm, H_cer_npeSum, H_dc_InsideDipoleExit, P_dc_InsideDipoleExit, CTime_ePiCoinTime_ROC1, P_gtr_beta, P_gtr_xp, P_gtr_yp, P_gtr_p, P_gtr_dp, P_dc_xfp, P_dc_xpfp, P_dc_yfp, P_dc_ypfp, P_hod_goodscinhit, P_hod_goodstarttime, P_cal_etotnorm, P_cal_etottracknorm, P_cal_fly_earray, P_cal_pr_eplane, Cal_pandas, P_aero_npeSum, P_aero_xAtAero, P_aero_yAtAero, P_hgcer_npeSum, P_hgcer_xAtCer, P_hgcer_yAtCer, P_ngcer_npeSum, P_ngcer_xAtCer, P_ngcer_yAtCer, MMpi, H_RF_Dist, P_RF_Dist, xCalo, yCalo, xExit, yExit, Q2, W, epsilon, ph_q, MandelT) in zip(*Cut_COIN_Pions_rand_tmp)
        ]

    COIN_Pions = {
        "Uncut_Pion_Events" : Uncut_COIN_Pions,
        "Cut_Pion_Events_noRF" : Cut_COIN_Pions_noRF,
        "Cut_Pion_Events_All" : Cut_COIN_Pions_all,
        "Cut_Pion_Events_Prompt" : Cut_COIN_Pions_prompt,
        "Cut_Pion_Events_Random" : Cut_COIN_Pions_random,
        }

    return COIN_Pions

#################################################################################################################################################################

def main():
    COIN_Pion_Data = coin_pions()

    # This is just the list of branches we use from the initial root file for each dict
    # I don't like re-defining this here as it's very prone to errors if you included (or removed something) earlier but didn't modify it here
    # Should base the branches to include based on some list and just repeat the list here (or call it again directly below)

    COIN_Pion_Data_Header = ["H_gtr_beta","H_gtr_xp","H_gtr_yp","H_gtr_dp","H_dc_xfp", "H_dc_xpfp", "H_dc_yfp", "H_dc_ypfp","H_hod_goodscinhit","H_hod_goodstarttime","H_cal_etotnorm","H_cal_etottracknorm","H_cer_npeSum","H_dc_InsideDipoleExit","P_dc_InsideDipoleExit","CTime_ePiCoinTime_ROC1","P_gtr_beta","P_gtr_xp","P_gtr_yp","P_gtr_p","P_gtr_dp","P_dc_xfp", "P_dc_xpfp", "P_dc_yfp", "P_dc_ypfp","P_hod_goodscinhit","P_hod_goodstarttime","P_cal_etotnorm","P_cal_etottracknorm","P_cal_fly_earray","P_cal_pr_eplane","Cal_pandas","P_aero_npeSum","P_aero_xAtAero","P_aero_yAtAero","P_hgcer_npeSum","P_hgcer_xAtCer","P_hgcer_yAtCer","P_ngcer_npeSum","P_ngcer_xAtCer","P_ngcer_yAtCer","MMpi","H_RF_Dist","P_RF_Dist","xCalo","yCalo","xExit","yExit","Q2","W","epsilon","ph_q","MandelT"]

    # Need to create a dict for all the branches we grab
    data = {}
    data.update(COIN_Pion_Data)
    print(type(pd.DataFrame(data.get("Cut_Pion_Events_Prompt"))[31][1]))
    data_keys = list(data.keys()) # Create a list of all the keys in all dicts added above, each is an array of data
#    print(data_keys)
 #   print(len(data))
#    print("test")
    for i in range (0, len(data_keys)):
        if("Pion" in data_keys[i]):
            DFHeader=list(COIN_Pion_Data_Header)
        else:
            continue
      #  print(type(COIN_Pion_Data))
     #   print(i)
            # Uncomment the line below if you want .csv file output, WARNING the files can be very large and take a long time to process!
            #pd.DataFrame(data.get(data_keys[i])).to_csv("%s/%s_%s.csv" % (OUTPATH, data_keys[i], runNum), header=DFHeader, index=False) # Convert array to panda dataframe and write to csv with correct header
        if (i == 0):
            pd.DataFrame(data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_Analysed_Data.root" % (OUTPATH, runNum, MaxEvent), key ="%s" % data_keys[i])
        elif (i != 0):
            pd.DataFrame(data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_Analysed_Data.root" % (OUTPATH, runNum, MaxEvent), key ="%s" % data_keys[i], mode ='a') 
                    
if __name__ == '__main__':
    main()
print ("Processing Complete")
