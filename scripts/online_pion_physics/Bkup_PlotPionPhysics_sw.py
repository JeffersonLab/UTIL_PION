#! /usr/bin/python
####################################################################################
# Created - 20/July/21, Author - Muhammad Junaid, University of Regina, Canada
####################################################################################
# Python version of the pion plotting script. Now utilises uproot to select event of each type and writes them to a root file.
# Python should allow for easier reading of databases storing diferent variables.
# This version of script is for shift workers at JLab
# To run this script, execute: python3 scriptname runnumber

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
import array
from ROOT import TCanvas, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TGaxis, TLine, TMath, TPaveText, TArc, TGraphPolar 
from ROOT import kBlack, kBlue, kRed
sys.path.insert(0, 'python/')

##################################################################################################################################################

# Defining some variables here
minrangeuser = 0 # min range for -t vs phi plot
maxrangeuser = 0.7 # max range for -t vs phi plot
minbin = 0.92 # minimum bin for selecting neutrons events in missing mass distribution
maxbin = 0.98 # maximum bin for selecting neutrons events in missing mass distribution

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
    REPLAYPATH = "/group/c-pionlt/USERS/%s/hallc_replay_lt" % USER[1]
elif ("qcd" in HOST[1]):
    REPLAYPATH = "/group/c-pionlt/USERS/%s/hallc_replay_lt" % USER[1]
elif ("cdaq" in HOST[1]):
    REPLAYPATH = "/home/cdaq/hallc-online/hallc_replay_lt"
elif ("phys.uregina" in HOST[1]):
    REPLAYPATH = "/home/%s/work/JLab/hallc_replay_lt" % USER[1]
elif("skynet" in HOST[1]):
    REPLAYPATH = "/home/%s/Work/JLab/hallc_replay_lt" % USER[1]

#################################################################################################################################################

# Add more path setting as needed in a similar manner                                                                                                                                                          
OUTPATH = "%s/UTIL_PION/OUTPUT/Analysis/PionLT" % REPLAYPATH        # Output folder location                                                                                                     
sys.path.insert(0, '%s/UTIL_PION/bin/python/' % REPLAYPATH)
import kaonlt as klt # Import kaonlt module, need the path setting line above prior to importing this                                                                                                         
print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER[1], HOST[1], REPLAYPATH))
Pion_Analysis_Distributions = "%s/%s_%s_sw_Pion_Analysis_Distributions.pdf" % (OUTPATH, runNum, MaxEvent)

#################################################################################################################################################

# Construct the name of the rootfile based upon the info we provided
rootName = "%s/UTIL_PION/OUTPUT/Analysis/PionLT/%s_%s_%s.root" % (REPLAYPATH, runNum, MaxEvent, ROOTPrefix)     # Input file location and variables taking
#rootName = "/home/cdaq/hallc-online/hallc_replay_lt/UTIL_PION/OUTPUT/Analysis/PionLT/8076_-1_Analysed_Data.root" # Hard coded file for testing on cdaq
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
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Section for grabing Prompt/Random selection parameters from PARAM file
PARAMPATH = "%s/UTIL_PION/DB/PARAM" % REPLAYPATH
print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER[1], HOST[1], REPLAYPATH))
TimingCutFile = "%s/Timing_Parameters.csv" % PARAMPATH # This should match the param file actually being used!

TimingCutf = open(TimingCutFile)
try:
    TimingCutFile
except NameError:
    print("!!!!! ERRROR !!!!!\n One (or more) of the cut files not found!\n!!!!! ERRORR !!!!!")
    sys.exit(2)
print("Reading timing cuts from %s" % TimingCutFile)

PromptWindow = [0, 0]
RandomWindows = [0, 0, 0, 0]
linenum = 0 # Count line number we're on
TempPar = -1 # To check later
for line in TimingCutf: # Read all lines in the cut file
    linenum += 1 # Add one to line number at start of loop
    if(linenum > 1): # Skip first line
        line = line.partition('#')[0] # Treat anything after a # as a comment and ignore it
        line = line.rstrip()
        array = line.split(",") # Convert line into an array, anything after a comma is a new entry
        if(int(runNum) in range (int(array[0]), int(array[1])+1)): # Check if run number for file is within any of the ranges specified in the cut file
            TempPar += 2 # If run number is in range, set to non -1 value
            BunchSpacing = float(array[2])
            CoinOffset = float(array[3]) # Coin offset value
            nSkip = float(array[4]) # Number of random windows skipped
            nWindows = float(array[5]) # Total number of random windows
            PromptPeak = float(array[6]) # Pion CT prompt peak positon
TimingCutf.close() # After scanning all lines in file, close file

if(TempPar == -1): # If value is still -1, run number provided din't match any ranges specified so exit
    print("!!!!! ERROR !!!!!\n Run number specified does not fall within a set of runs for which cuts are defined in %s\n!!!!! ERROR !!!!!" % TimingCutFile)
    sys.exit(3)
elif(TempPar > 1):
    print("!!! WARNING!!! Run number was found within the range of two (or more) line entries of %s !!! WARNING !!!" % TimingCutFile)
    print("The last matching entry will be treated as the input, you should ensure this is what you want")

# From our values from the file, reconstruct our windows
PromptWindow[0] = PromptPeak - (BunchSpacing/2) - CoinOffset
PromptWindow[1] = PromptPeak + (BunchSpacing/2) + CoinOffset
RandomWindows[0] = PromptPeak - (BunchSpacing/2) - CoinOffset - (nSkip*BunchSpacing) - ((nWindows/2)*BunchSpacing)
RandomWindows[1] = PromptPeak - (BunchSpacing/2) - CoinOffset - (nSkip*BunchSpacing) 
RandomWindows[2] = PromptPeak + (BunchSpacing/2) + CoinOffset + (nSkip*BunchSpacing) 
RandomWindows[3] = PromptPeak + (BunchSpacing/2) + CoinOffset + (nSkip*BunchSpacing) + ((nWindows/2)*BunchSpacing)

###############################################################################################################################################

# Read stuff from the main event tree
infile = ROOT.TFile.Open(rootName, "READ")
Uncut_Pion_Events_tree = infile.Get("Uncut_Pion_Events")
Cut_Pion_Events_noRF_tree = infile.Get("Cut_Pion_Events_noRF")
Cut_Pion_Events_All_tree = infile.Get("Cut_Pion_Events_All")
Cut_Pion_Events_Prompt_tree = infile.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Random_tree = infile.Get("Cut_Pion_Events_Random")

###############################################################################################################################################

# Defining Histograms for Pions
P_RFTime_Dist_pions_cut_noRF = ROOT.TH1D("P_RFTime_Dist_pions_cut_noRF", "SHMS RFTime; SHMS_RFTime; Counts", 200, 0, 4)

H_beta_pions_uncut = ROOT.TH1D("H_beta_pions_uncut", "HMS #beta; HMS_gtr_#beta; Counts", 200, 0.5, 1.5)
H_xp_pions_uncut = ROOT.TH1D("H_xp_pions_uncut", "HMS x'; HMS_gtr_xp; Counts", 200, -0.2, 0.2)
H_yp_pions_uncut = ROOT.TH1D("H_yp_pions_uncut", "HMS y'; HMS_gtr_yp; Counts", 200, -0.2, 0.2)
H_dp_pions_uncut = ROOT.TH1D("H_dp_pions_uncut", "HMS #delta; HMS_gtr_dp; Counts", 200, -12, 12)
H_hod_goodscinhit_pions_uncut = ROOT.TH1D("H_hod_goodscinhit_pions_uncut", "HMS hod goodscinhit; HMS_hod_goodscinhi; Counts", 200, 0.7, 1.3)
H_hod_goodstarttime_pions_uncut = ROOT.TH1D("H_hod_goodstarttime_pions_uncut", "HMS hod goodstarttime; HMS_hod_goodstarttime; Counts", 200, 0.7, 1.3)
H_cal_etotnorm_pions_uncut = ROOT.TH1D("H_cal_etotnorm_pions_uncut", "HMS cal etotnorm; HMS_cal_etotnorm; Counts", 200, 0.2, 1.8)
H_cal_etottracknorm_pions_uncut = ROOT.TH1D("H_cal_etottracknorm_pions_uncut", "HMS cal etottracknorm; HMS_cal_etottracknorm; Counts", 300, 0.2, 1.8)
H_cer_npe_pions_uncut = ROOT.TH1D("H_cer_npe_pions_uncut", "HMS cer npeSum; HMS_cer_npeSum; Counts", 200, 0, 50)
H_RFTime_pions_uncut = ROOT.TH1D("H_RFTime_pions_uncut", "HMS RFTime; HMS_RFTime; Counts", 200, 0, 4)
P_beta_pions_uncut = ROOT.TH1D("P_beta_pions_uncut", "SHMS #beta; SHMS_gtr_#beta; Counts", 200, 0.5, 1.5)
P_xp_pions_uncut = ROOT.TH1D("P_xp_pions_uncut", "SHMS x'; SHMS_gtr_xp; Counts", 200, -0.2, 0.2)
P_yp_pions_uncut = ROOT.TH1D("P_yp_pions_uncut", "SHMS y'; SHMS_gtr_yp; Counts", 200, -0.2, 0.2)
P_dp_pions_uncut = ROOT.TH1D("P_dp_pions_uncut", "SHMS #delta; SHMS_gtr_dp; Counts", 200, -30, 30)
P_p_pions_uncut = ROOT.TH1D("P_p_pions_uncut", "SHMS p; SHMS_gtr_p; Counts", 200, 4, 8)
P_hod_goodscinhit_pions_uncut = ROOT.TH1D("P_hod_goodscinhit_pions_uncut", "SHMS hod goodscinhit; SHMS_hod_goodscinhit; Counts", 200, 0.7, 1.3)
P_hod_goodstarttime_pions_uncut = ROOT.TH1D("P_hod_goodstarttime_pions_uncut", "SHMS hod goodstarttime; SHMS_hod_goodstarttime; Counts", 200, 0.7, 1.3)
P_cal_etotnorm_pions_uncut = ROOT.TH1D("P_cal_etotnorm_pions_uncut", "SHMS cal etotnorm; SHMS_cal_etotnorm; Counts", 200, 0, 1)
P_cal_etottracknorm_pions_uncut = ROOT.TH1D("P_cal_etottracknorm_pions_uncut", "SHMS cal etottracknorm; SHMS_cal_etottracknorm; Counts", 200, 0, 1.6)
P_hgcer_npe_pions_uncut = ROOT.TH1D("P_hgcer_npe_pions_uncut", "SHMS HGC npeSum; SHMS_hgcer_npeSum; Counts", 200, 0, 50)
P_xhgcer_pions_uncut = ROOT.TH1D("P_xhgcer_pions_uncut", "SHMS HGC xAtCer; SHMS_hgcer_xAtCer; Counts", 200, -60, 60)
P_yhgcer_pions_uncut = ROOT.TH1D("P_yhgcer_pions_uncut", "SHMS HGC yAtCer; SHMS_hgcer_yAtCer; Counts", 200, -50, 50)
P_aero_npe_pions_uncut = ROOT.TH1D("P_aero_npe_pions_uncut", "SHMS aero npeSum; SHMS_aero_npeSum; Counts", 200, 0, 50)
P_xaero_pions_uncut = ROOT.TH1D("P_xacero_pions_uncut", "SHMS aero xAtAero; SHMS_aero_xAtAero; Counts", 200, -60, 60)
P_yaero_pions_uncut = ROOT.TH1D("P_yaero_pions_uncut", "SHMS aero yAtAero; SHMS_aero_yAtAero; Counts", 200, -50, 50)
P_ngcer_npe_pions_uncut = ROOT.TH1D("P_ngcer_npe_pions_uncut", "SHMS NGC npeSum; SHMS_ngcer_npeSum; Counts", 200, 0, 50)
P_xngcer_pions_uncut = ROOT.TH1D("P_xngcer_pions_uncut", "SHMS NGC xAtCer; SHMS_ngcer_xAtCer; Counts", 200, -70, 50)
P_yngcer_pions_uncut = ROOT.TH1D("P_yngcer_pions_uncut", "SHMS NGC yAtCer; SHMS_ngcer_yAtCer; Counts", 200, -50, 50)
P_MMpi_pions_uncut = ROOT.TH1D("P_MMpi_pions_uncut", "MIssing Mass (no cuts); MM_{#pi}; Counts", 200, 0.5, 1.8)
P_RFTime_pions_uncut = ROOT.TH1D("P_RFTime_pions_uncut", "SHMS RFTime; SHMS_RFTime; Counts", 200, 0, 4)
ePiCoinTime_pions_uncut = ROOT.TH1D("ePiCoinTime_pions_uncut", "Electron-Pion CTime (no cuts); e #pi Coin_Time; Counts", 200, -50, 50)
Q2_pions_uncut = ROOT.TH1D("Q2_pions_uncut", "Q2; Q2; Counts", 200, 0, 6)
W_pions_uncut = ROOT.TH1D("W_pions_uncut", "W; W; Counts", 200, 2, 4)
epsilon_pions_uncut = ROOT.TH1D("epsilon_pions_uncut", "epsilon; epsilon; Counts", 200, 0, 0.8)
phiq_pions_uncut = ROOT.TH1D("phiq_pions_uncut", "phiq; #phi; Counts", 200, -10, 10)
t_pions_uncut = ROOT.TH1D("t_pions_uncut", "t; t; Counts", 200, -1.5, 1)

H_beta_pions_cut = ROOT.TH1D("H_beta_pions_cut", "HMS #beta; HMS_gtr_#beta; Counts", 200, 0.5, 1.5)
H_xp_pions_cut = ROOT.TH1D("H_xp_pions_cut", "HMS x'; HMS_gtr_xp; Counts", 200, -0.2, 0.2)
H_yp_pions_cut = ROOT.TH1D("H_yp_pions_cut", "HMS y'; HMS_gtr_yp; Counts", 200, -0.2, 0.2)
H_dp_pions_cut = ROOT.TH1D("H_dp_pions_cut", "HMS #delta; HMS_gtr_dp; Counts", 200, -12, 12)
H_hod_goodscinhit_pions_cut = ROOT.TH1D("H_hod_goodscinhit_pions_cut", "HMS hod goodscinhit; HMS_hod_goodscinhit; Counts", 200, 0.7, 1.3)
H_hod_goodstarttime_pions_cut = ROOT.TH1D("H_hod_goodstarttime_pions_cut", "HMS hod goodstarttime; HMS_hod_goodstarttime; Counts", 200, 0.7, 1.3)
H_cal_etotnorm_pions_cut = ROOT.TH1D("H_cal_etotnorm_pions_cut", "HMS cal etotnorm; HMS_cal_etotnorm; Counts", 200, 0.6, 1.4)
H_cal_etottracknorm_pions_cut = ROOT.TH1D("H_cal_etottracknorm_pions_cut", "HMS cal etottracknorm; HMS_cal_etottracknorm; Counts", 200, 0.6, 1.4)
H_cer_npeSum_pions_cut = ROOT.TH1D("H_cer_npeSum_pions_cut", "HMS cer npeSum; HMS_cer_npeSum; Counts", 200, 0, 50)
H_RFTime_pions_cut = ROOT.TH1D("H_RFTime_pions_cut", "HMS RFTime; HMS_RFTime; Counts", 200, 0, 4)
P_beta_pions_cut = ROOT.TH1D("P_beta_pions_cut", "SHMS #beta; SHMS_gtr_#beta; Counts", 200, 0.5, 1.5)
P_xp_pions_cut = ROOT.TH1D("P_xp_pions_cut", "SHMS x'; SHMS_gtr_xp; Counts", 200, -0.2, 0.2)
P_yp_pions_cut = ROOT.TH1D("P_yp_pions_cut", "SHMS y'; SHMS_gtr_yp; Counts", 200, -0.2, 0.2)
P_dp_pions_cut = ROOT.TH1D("P_dp_pions_cut", "SHMS #delta; SHMS_gtr_dp; Counts", 200, -15, 15)
P_p_pions_cut = ROOT.TH1D("P_p_pions_cut", "SHMS p; SHMS_gtr_p; Counts", 200, 4, 8)
P_hod_goodscinhit_pions_cut = ROOT.TH1D("P_hod_goodscinhit_pions_cut", "SHMS hod goodscinhit; SHMS_hod_goodscinhit; Counts", 200, 0.7, 1.3)
P_hod_goodstarttime_pions_cut = ROOT.TH1D("P_hod_goodstarttime_pions_cut", "SHMS hod goodstarttime; SHMS_hod_goodstarttime; Counts", 200, 0.7, 1.3)
P_cal_etotnorm_pions_cut = ROOT.TH1D("P_cal_etotnorm_pions_cut", "SHMS cal etotnorm; SHMS_cal_etotnorm; Counts", 200, 0, 1.2)
P_cal_etottracknorm_pions_cut = ROOT.TH1D("P_cal_etottracknorm_pions_cut", "SHMS cal etottracknorm; SHMS_cal_etottracknorm; Counts", 200, 0, 1.6)
P_hgcer_npeSum_pions_cut = ROOT.TH1D("P_hgcer_npeSum_pions_cut", "SHMS HGC npeSum; SHMS_hgcer_npeSum; Counts", 200, 0, 50)
P_hgcer_xAtCer_pions_cut = ROOT.TH1D("P_hgcer_xAtCer_pions_cut", "SHMS HGC xAtCer; SHMS_hgcer_xAtCer; Counts", 200, -40, 30)
P_hgcer_yAtCer_pions_cut = ROOT.TH1D("P_hgcer_yAtCer_pions_cut", "SHMS HGC yAtCer; SHMS_hgcer_yAtCer; Counts", 200, -30, 30)
P_aero_npeSum_pions_cut = ROOT.TH1D("P_aero_npeSum_pions_cut", "SHMS aero npeSum; SHMS_aero_npeSum; Counts", 200, 0, 50)
P_aero_xAtAero_pions_cut = ROOT.TH1D("P_acero_xAtAero_pions_cut", "SHMS aero xAtAero; SHMS_aero_xAtAero; Counts", 200, -40, 30)
P_aero_yAtAero_pions_cut = ROOT.TH1D("P_aero_yAtAero_pions_cut", "SHMS aero yAtAero; SHMS_aero_yAtAero; Counts", 200, -30, 30)
P_ngcer_npeSum_pions_cut = ROOT.TH1D("P_ngcer_npeSum_pions_cut", "SHMS NGC npeSum; SHMS_ngcer_npeSum; Counts", 200, 0, 50)
P_ngcer_xAtCer_pions_cut = ROOT.TH1D("P_ngcer_xAtCer_pions_cut", "SHMS NGC xAtCer; SHMS_ngcer_xAtCer; Counts", 200, -40, 30)
P_ngcer_yAtCer_pions_cut = ROOT.TH1D("P_ngcer_yAtCer_pions_cut", "SHMS NGC yAtCer; SHMS_ngcer_yAtCer; Counts", 200, -30, 30)
P_MMpi_pions_cut = ROOT.TH1D("P_MMpi_pions_cut", "Missing Mass (with cuts); MM_{#pi}; Counts", 200, 0.5, 1.8)
P_RFTime_pions_cut = ROOT.TH1D("P_RFTime_pions_cut", "SHMS RFTime; SHMS_RFTime; Counts", 200, 0, 4)
ePiCoinTime_pions_cut = ROOT.TH1D("ePiCoinTime_pions_cut", "Electron-Pion CTime (with cuts); e #pi Coin_Time; Counts", 200, -50, 50)
Q2_pions_cut = ROOT.TH1D("Q2_pions_cut", "Q2; Q2; Counts", 200, 2, 4)
W_pions_cut = ROOT.TH1D("W_pions_cut", "W; W; Counts", 200, 2.2, 4)
epsilon_pions_cut = ROOT.TH1D("epsilon_pions_cut", "epsilon; epsilon; Counts", 200, 0, 0.8)
phiq_pions_cut = ROOT.TH1D("phiq_pions_cut", "phiq; #phi; Counts", 200, -10, 10)
t_pions_cut = ROOT.TH1D("t_pions_cut", "t; t; Counts", 200, -1, 0.5)

P_beta_pions_cut_prompt = ROOT.TH1D("P_beta_pions_cut_prompt", "SHMS beta; SHMS_#beta; Counts", 200, 0.8, 1.2)
P_RFTime_pions_cut_prompt = ROOT.TH1D("P_RFTime_pions_cut_prompt", "SHMS RFTime; SHMS_RFTime; Counts", 200, 0, 4)
ePiCoinTime_pions_cut_prompt = ROOT.TH1D("ePiCoinTime_pions_cut_prompt", "Electron-Pion CTime; e #pi Coin_Time; Counts", 8, -2, 2)
PMMpi_pions_cut_prompt = ROOT.TH1D("P_MMpi_pions_cut_prompt", "Missing Mass; MM_{#pi}; Counts", 200, 0.5, 1.8)

P_beta_pions_cut_random = ROOT.TH1D("P_beta_pions_cut_random", "SHMS beta; SHMS_#beta; Counts", 200, 0.8, 1.2)
P_RFTime_pions_cut_random = ROOT.TH1D("P_RFTime_pions_cut_random", "SHMS RFTime; SHMS_RFTime; Counts", 200, 0, 4)
ePiCoinTime_pions_cut_random = ROOT.TH1D("ePiCoinTime_pions_cut_random", "Electron-Pion CTime; e #pi Coin_Time; Counts", 160, -40, 40)
P_MMpi_pions_cut_random = ROOT.TH1D("P_MMpi_pions_cut_random", "Missing Mass; MM_{#pi}; Counts", 200, 0.5, 1.8)

P_MMpi_pions_cut_random_scaled = ROOT.TH1D("P_MMpi_pions_cut_random_scaled", "Missing Mass; MM_{#pi}; Counts", 200, 0.5, 1.8)
P_MMpi_pions_cut_random_sub = ROOT.TH1D("P_MMpi_pions_cut_random_sub", "Missing Mass Rndm Sub; MM_{#pi}; Counts", 200, 0.5, 1.8)

##############################################################################################################################################

# 2D Histograms for pions
H_cal_etottracknorm_vs_cer_npeSum_pions_uncut = ROOT.TH2D("H_cal_etottracknorm_vs_cer_npeSum_pions_uncut","HMS cal etottracknorm vs HMS cer npeSum (no cut); H_cal_etottracknorm; H_cer_npeSum",100, 0, 2, 100, 0, 40)
P_hgcer_vs_aero_npeSum_pions_uncut = ROOT.TH2D("P_hgcer_vs_aero_npeSum_pions_uncut", "SHMS HGC npeSum vs SHMS Aero npeSum (no cut); SHMS_hgcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
ePiCoinTime_vs_MMpi_pions_uncut = ROOT.TH2D("ePiCoinTime_vs_MMpi_pions_uncut","Electron-Pion CTime vs Missing Mass (no cut); e #pi Coin_Time; MM_{#pi}", 200, -40, 40, 200, 0, 2)
P_hgcer_yAtCer_vs_xAtCer_pions_uncut = ROOT.TH2D("P_hgcer_yAtCer_vs_xAtCer_pions_uncut", "SHMS HGC yAtCer vs SHMS HGC xAtCer (no cut); SHMS_hgcer_yAtCer; SHMS_hgcer_xAtCer", 100, -50, 50, 100, -50, 50)
P_aero_yAtAero_vs_xAtAero_pions_uncut = ROOT.TH2D("P_aero_yAtAero_vs_xAtAero_pions_uncut", "SHMS aero yAtAero vs SHMS aero xAtAero (no cut); SHMS_aero_yAtAero; SHMS_aero_xAtAero", 100, -50, 50, 100, -50, 50)
ePiCoinTime_vs_beta_pions_uncut = ROOT.TH2D("ePiCoinTime_vs_beta_pions_uncut", "Electron-Pion CTime vs SHMS #beta (no cut); e #pi Coin_Time; SHMS_#beta", 200, -40, 40, 200, 0, 2)
P_RFTime_vs_MMpi_pions_uncut = ROOT.TH2D("P_RFTime_vs_MMpi_pions_uncut", "SHMS RFTime vs Missing Mass (no cuts); SHMS_RFTime_Dist; MM_{#pi}", 100, 0, 4, 100, 0, 2)
P_cal_etottracknorm_vs_ngcer_npeSum_pions_uncut = ROOT.TH2D("P_cal_etottracknorm_vs_ngcer_npeSum_pions_uncut", "SHMS cal etottracknorm vs SHMS NGC xAtCer (no cut); SHMS_cal_etottracknorm; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
P_ngcer_yAtCer_vs_xAtCer_pions_uncut = ROOT.TH2D("P_ngcer_yAtCer_vs_xAtCer_pions_uncut", "SHMS NGC yAtCer vs SHMS NGC xAtCer (no cut); SHMS_ngcer_yAtCer; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
P_ngcer_vs_hgcer_npeSum_pions_uncut = ROOT.TH2D("P_ngcer_vs_hgcer_npeSum_pions_uncut", "SHMS NGC npeSum vs SHMS HGC npeSum (no cut); SHMS_ngcer_npeSum; SHMS_hgcer_npeSum", 100, 0, 50, 100, 0, 50)
P_ngcer_vs_aero_npeSum_pions_uncut = ROOT.TH2D("P_ngcer_vs_aero_npeSum_pions_uncut", "SHMS NGC npeSum vs SHMS aero npeSum (no cut); SHMS_ngcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
H_dp_vs_beta_pions_uncut = ROOT.TH2D("H_dp_vs_beta_pions_uncut", "HMS #delta vs HMS #beta (no cut); HMS #delta; HMS_#beta", 200, -12, 12, 200, 0, 2)
P_dp_vs_beta_pions_uncut = ROOT.TH2D("P_dp_vs_beta_pions_uncut", "SHMS #delta vs SHMS #beta (no cut); SHMS #delta; HMS_#beta", 200, -12, 12, 200, 0, 2)
P_MMpi_vs_beta_pions_uncut = ROOT.TH2D("P_MMpi_vs_beta_pions_uncut", "Missing Mass vs SHMS #beta (no cut); MM_{#pi}; SHMS_#beta", 100, 0, 2, 200, 0, 2)

H_cal_etottracknorm_vs_cer_npeSum_pions_cut_all = ROOT.TH2D("H_cal_etottracknorm_vs_cer_npeSum_pions_cut_all","HMS cal etottracknorm vs HMS cer npeSum (with cuts); H_cal_etottracknorm; H_cer_npeSum",100, 0.5, 1.5, 100, 0, 40)
P_hgcer_vs_aero_npeSum_pions_cut_all = ROOT.TH2D("P_hgcer_vs_aero_npeSum_pions_cut_all", "SHMS HGC npeSum vs SHMS aero npeSum (with cuts); SHMS_hgcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
ePiCoinTime_vs_MMpi_pions_cut_all = ROOT.TH2D("ePiCoinTime_vs_MMpi_pions_cut_all","Electron-Pion CTime vs Missing Mass (with PID cuts); e #pi Coin_Time; MM_{#pi}", 100, -2, 2, 100, 0, 2)
P_hgcer_yAtCer_vs_xAtCer_pions_cut_all = ROOT.TH2D("P_hgcer_yAtCer_vs_xAtCer_pions_cut_all", "SHMS HGC yAtCer vs SHMS HGC xAtCer (with cuts); SHMS_hgcer_yAtCer; SHMS_hgcer_xAtCer", 100, -50, 50, 100, -50, 50)
P_aero_yAtAero_vs_xAtAero_pions_cut_all = ROOT.TH2D("P_aero_yAtAero_vs_xAtAero_pions_cut_all", "SHMS aero yAtAero vs SHMS aero xAtAero (with cuts); SHMS_aero_yAtAero; SHMS_aero_xAtAero", 100, -50, 50, 100, -50, 50)
ePiCoinTime_vs_beta_pions_cut_all = ROOT.TH2D("ePiCoinTime_vs_beta_pions_cut_all", "Electron-Pion CTime vs SHMS #beta (with PID cuts); e #pi Coin_Time; SHMS_#beta", 100, -2, 2, 100, 0.6, 1.4)
P_RFTime_vs_MMpi_pions_cut_all = ROOT.TH2D("P_RFTime_vs_MMpi_pions_cut_all", "SHMS RFTime vs Missing Mass (with cuts); SHMS_RFTime_Dist; MM_{#pi}", 100, 0, 4, 100, 0, 2)
P_cal_etottracknorm_vs_ngcer_npeSum_pions_cut_all = ROOT.TH2D("P_cal_etottracknorm_vs_ngcer_npeSum_pions_cut_all", "P cal etottracknorm vs SHMS NGC xAtCer (with cuts); SHMS_cal_etottracknorm; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
P_ngcer_yAtCer_vs_xAtCer_pions_cut_all = ROOT.TH2D("P_ngcer_yAtCer_vs_xAtCer_pions_cut_all", "SHMS NGC yAtCer vs SHMS NGC xAtCer (with cuts); SHMS_ngcer_yAtCer; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
P_ngcer_vs_hgcer_npeSum_pions_cut_all = ROOT.TH2D("P_ngcer_vs_hgcer_npeSum_pions_cut_all", "SHMS NGC npeSum vs SHMS HGC npeSum (with cuts); SHMS_ngcer_npeSum; SHMS_hgcer_npeSum", 100, 0, 50, 100, 0, 50)
P_ngcer_vs_aero_npeSum_pions_cut_all = ROOT.TH2D("P_ngcer_vs_aero_npeSum_pions_cut_all", "SHMS NGC npeSum vs SHMS aero npeSum (with cuts); SHMS_ngcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
H_dp_vs_beta_pions_cut_all = ROOT.TH2D("H_dp_vs_beta_pions_cut_all", "HMS #delta vs HMS #beta (with cut); HMS #delta; HMS_#beta", 200, -12, 12, 200, 0, 2)
P_dp_vs_beta_pions_cut_all = ROOT.TH2D("P_dp_vs_beta_pions_cut_all", "SHMS #delta vs SHMS #beta (with cut); SHMS #delta; SHMS_#beta", 200, -30, 30, 200, 0, 2)
P_MMpi_vs_beta_pions_cut_all = ROOT.TH2D("P_MMpi_vs_beta_pions_cut_all", "Missing Mass vs SHMS #beta (with cut); MM_{#pi}; SHMS_#beta", 100, 0, 2, 200, 0, 2)

MMpi_vs_ePiCoinTime_pions_cut_prompt = ROOT.TH2D("MMpi_vs_ePiCoinTime_pions_cut_prompt","Missing Mass vs Electron-Pion CTime; MM_{#pi}; e #pi Coin_Time",100, 0, 2, 100, -2, 2)
Q2vsW = ROOT.TH2D("Q2vsW", "Q2 vs W; Q2; W", 200, 0.5, 4.5, 200, 2.7, 3.6)
phiqvst = ROOT.TH2D("phiqvst","; #phi ;t", 12, -3.14, 3.14, 24, 0.0, 1.2)

# 11/09/21 - SJDK - Adding some 3D XY NPE plots, need to take projections of these (which I need to figure out how to do in PyRoot. Making these manually from the command line for now.
P_HGC_xy_npe_pions_cut = ROOT.TH3D("P_HGC_xy_npe_pions_cut", "SHMS HGC NPE as fn of yAtCer vs SHMS HGC xAtCer (with cuts); HGC_yAtCer(cm); HGC_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_Aero_xy_npe_pions_cut = ROOT.TH3D("P_Aero_xy_npe_pions_cut", "SHMS Aerogel NPE as fn of yAtCer vs xAtCer (with cuts); Aero_yAtCer(cm); Aero_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_NGC_xy_npe_pions_cut = ROOT.TH3D("P_NGC_xy_npe_pions_cut", "SHMS NGC NPE as fn of yAtCer vs xAtCer (with cuts); NGC_yAtCer(cm); NGC_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
#################################################################################################################################################



# Filling Histograms for Pions
for event in Cut_Pion_Events_noRF_tree:
    P_RFTime_Dist_pions_cut_noRF.Fill(event.P_RF_Dist)

for event in Uncut_Pion_Events_tree:
    #H_gtr_beta_pions_uncut.Fill(event.H_gtr_beta)
    H_gtr_xp_pions_uncut.Fill(event.H_gtr_xp)
    H_gtr_yp_pions_uncut.Fill(event.H_gtr_yp)
    H_gtr_dp_pions_uncut.Fill(event.H_gtr_dp)
    #H_hod_goodscinhit_pions_uncut.Fill(event.H_hod_goodscinhit)
    #H_hod_goodstarttime_pions_uncut.Fill(event.H_hod_goodstarttime)
    #H_cal_etotnorm_pions_uncut.Fill(event.H_cal_etotnorm)
    H_cal_etottracknorm_pions_uncut.Fill(event.H_cal_etottracknorm)
    #H_cer_npeSum_pions_uncut.Fill(event.H_cer_npeSum)
    #H_RFTime_Dist_pions_uncut.Fill(event.H_RF_Dist)
    #P_gtr_beta_pions_uncut.Fill(event.P_gtr_beta)
    P_gtr_xp_pions_uncut.Fill(event.P_gtr_xp)
    P_gtr_yp_pions_uncut.Fill(event.P_gtr_yp)
    P_gtr_dp_pions_uncut.Fill(event.P_gtr_dp)
    #P_gtr_p_pions_uncut.Fill(event.P_gtr_p)
    #P_hod_goodscinhit_pions_uncut.Fill(event.P_hod_goodscinhit)
    #P_hod_goodstarttime_pions_uncut.Fill(event.P_hod_goodstarttime)
    #P_cal_etotnorm_pions_uncut.Fill(event.P_cal_etotnorm)
    P_cal_etottracknorm_pions_uncut.Fill(event.P_cal_etottracknorm)
    P_hgcer_npeSum_pions_uncut.Fill(event.P_hgcer_npeSum)
    #P_hgcer_xAtCer_pions_uncut.Fill(event.P_hgcer_xAtCer)
    #P_hgcer_yAtCer_pions_uncut.Fill(event.P_hgcer_yAtCer)
    P_aero_npeSum_pions_uncut.Fill(event.P_aero_npeSum)
    #P_aero_xAtAero_pions_uncut.Fill(event.P_aero_xAtAero)
    #P_aero_yAtAero_pions_uncut.Fill(event.P_aero_yAtAero)
    P_ngcer_npeSum_pions_uncut.Fill(event.P_ngcer_npeSum)
    #P_ngcer_xAtCer_pions_uncut.Fill(event.P_ngcer_xAtCer)
    #P_ngcer_yAtCer_pions_uncut.Fill(event.P_ngcer_yAtCer)
    P_kin_MMpi_pions_uncut.Fill(event.MMpi)
    P_RFTime_Dist_pions_uncut.Fill(event.P_RF_Dist)
    CTime_ePiCoinTime_ROC1_pions_uncut.Fill(event.CTime_ePiCoinTime_ROC1)
    #Q2_pions_uncut.Fill(event.Q2)
    #W_pions_uncut.Fill(event.W)
    #epsilon_pions_uncut.Fill(event.epsilon)
    #phiq_uncut.Fill(event.ph_q)
    #t_uncut.Fill(-event.MandelT)
    #H_cal_etottracknorm_vs_H_cer_npeSum_pions_uncut.Fill(event.H_cal_etottracknorm, event.H_cer_npeSum)
    P_hgcer_npeSum_vs_aero_npeSum_pions_uncut.Fill(event.P_hgcer_npeSum, event.P_aero_npeSum)
    CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_uncut.Fill(event.CTime_ePiCoinTime_ROC1, event.MMpi)
    P_RFTime_vs_P_kin_MMpi_pions_uncut.Fill(event.P_RF_Dist, event.MMpi)
    #P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_uncut.Fill(event.P_hgcer_yAtCer, event.P_hgcer_xAtCer)
    #P_aero_yAtAero_vs_aero_xAtAero_pions_uncut.Fill(event.P_aero_yAtAero, event.P_aero_xAtAero)
    CTime_ePiCoinTime_ROC1_vs_beta_pions_uncut.Fill(event.CTime_ePiCoinTime_ROC1, event.P_gtr_beta)
    #P_cal_etottracknorm_vs_P_ngcer_npeSum_pions_uncut.Fill(event.P_cal_etottracknorm, event.P_ngcer_npeSum)
    #P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_uncut.Fill(event.P_ngcer_yAtCer, event.P_ngcer_xAtCer)
    P_ngcer_npeSum_vs_hgcer_npeSum_pions_uncut.Fill(event.P_ngcer_npeSum, event.P_hgcer_npeSum)
    P_ngcer_npeSum_vs_aero_npeSum_pions_uncut.Fill(event.P_ngcer_npeSum, event.P_aero_npeSum)
    H_gtr_dp_pions_vs_beta_pions_uncut.Fill(event.H_gtr_dp, event.H_gtr_beta) 
    P_gtr_dp_pions_vs_beta_pions_uncut.Fill(event.P_gtr_dp, event.P_gtr_beta) 
    P_kin_MMpi_pions_vs_beta_pions_uncut.Fill(event.MMpi, event.P_gtr_beta)
    P_HGC_xy_npe_uncut.Fill(event.P_hgcer_yAtCer,event.P_hgcer_xAtCer,event.P_hgcer_npeSum)
    P_Aero_xy_npe_uncut.Fill(event.P_aero_yAtAero,event.P_aero_xAtAero,event.P_aero_npeSum)

for event in Cut_Pion_Events_All_tree:
    #H_gtr_beta_pions_cut_all.Fill(event.H_gtr_beta)
    H_gtr_xp_pions_cut_all.Fill(event.H_gtr_xp)
    H_gtr_yp_pions_cut_all.Fill(event.H_gtr_yp)
    H_gtr_dp_pions_cut_all.Fill(event.H_gtr_dp)
    #H_hod_goodscinhit_pions_cut_all.Fill(event.H_hod_goodscinhit)
    #H_hod_goodstarttime_pions_cut_all.Fill(event.H_hod_goodstarttime)
    #H_cal_etotnorm_pions_cut_all.Fill(event.H_cal_etotnorm)
    H_cal_etottracknorm_pions_cut_all.Fill(event.H_cal_etottracknorm)
    #H_cer_npeSum_pions_cut_all.Fill(event.H_cer_npeSum)
    #H_RFTime_Dist_pions_cut_all.Fill(event.H_RF_Dist)
    #P_gtr_beta_pions_cut_all.Fill(event.P_gtr_beta)
    P_gtr_xp_pions_cut_all.Fill(event.P_gtr_xp)
    P_gtr_yp_pions_cut_all.Fill(event.P_gtr_yp)
    P_gtr_dp_pions_cut_all.Fill(event.P_gtr_dp)
    #P_gtr_p_pions_cut_all.Fill(event.P_gtr_p)
    #P_hod_goodscinhit_pions_cut_all.Fill(event.P_hod_goodscinhit)
    #P_hod_goodstarttime_pions_cut_all.Fill(event.P_hod_goodstarttime)
    #P_cal_etotnorm_pions_cut_all.Fill(event.P_cal_etotnorm)
    P_cal_etottracknorm_pions_cut_all.Fill(event.P_cal_etottracknorm)
    P_hgcer_npeSum_pions_cut_all.Fill(event.P_hgcer_npeSum)
    #P_hgcer_xAtCer_pions_cut_all.Fill(event.P_hgcer_xAtCer)
    #P_hgcer_yAtCer_pions_cut_all.Fill(event.P_hgcer_yAtCer)
    P_aero_npeSum_pions_cut_all.Fill(event.P_aero_npeSum)
    #P_aero_xAtAero_pions_cut_all.Fill(event.P_aero_xAtAero)
    #P_aero_yAtAero_pions_cut_all.Fill(event.P_aero_yAtAero)
    P_ngcer_npeSum_pions_cut_all.Fill(event.P_ngcer_npeSum)
    #P_ngcer_xAtCer_pions_cut_all.Fill(event.P_ngcer_xAtCer)
    #P_ngcer_yAtCer_pions_cut_all.Fill(event.P_ngcer_yAtCer)   
    P_kin_MMpi_pions_cut_all.Fill(event.MMpi)
    P_RFTime_Dist_pions_cut_all.Fill(event.P_RF_Dist)
    CTime_ePiCoinTime_ROC1_pions_cut_all.Fill(event.CTime_ePiCoinTime_ROC1)
    #Q2_pions_cut_all.Fill(event.Q2)
    #W_pions_cut_all.Fill(event.W)
    epsilon_pions_cut_all.Fill(event.epsilon)
    #phiq_cut_all.Fill(event.ph_q)
    #t_cut_all.Fill(-event.MandelT)
    H_cal_etottracknorm_vs_H_cer_npeSum_pions_cut_all.Fill(event.H_cal_etottracknorm, event.H_cer_npeSum)    
    P_hgcer_npeSum_vs_aero_npeSum_pions_cut_all.Fill(event.P_hgcer_npeSum, event.P_aero_npeSum)
    #P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_cut_all.Fill(event.P_hgcer_yAtCer, event.P_hgcer_xAtCer)
    #P_aero_yAtAero_vs_aero_xAtAero_pions_cut_all.Fill(event.P_aero_yAtAero, event.P_aero_xAtAero)
    CTime_ePiCoinTime_ROC1_vs_beta_pions_cut_all.Fill(event.CTime_ePiCoinTime_ROC1, event.P_gtr_beta)
    CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_cut_all.Fill(event.CTime_ePiCoinTime_ROC1, event.MMpi)
    P_RFTime_vs_P_kin_MMpi_pions_cut_all.Fill(event.P_RF_Dist, event.MMpi)    
    Q2vsW.Fill(event.Q2, event.W)
    phiqvst.Fill(event.ph_q, -event.MandelT)
    #P_cal_etottracknorm_vs_P_ngcer_npeSum_pions_cut_all.Fill(event.P_cal_etottracknorm, event.P_ngcer_npeSum)
    #P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_cut_all.Fill(event.P_ngcer_yAtCer, event.P_ngcer_xAtCer)
    P_ngcer_npeSum_vs_hgcer_npeSum_pions_cut_all.Fill(event.P_ngcer_npeSum, event.P_hgcer_npeSum)
    P_ngcer_npeSum_vs_aero_npeSum_pions_cut_all.Fill(event.P_ngcer_npeSum, event.P_aero_npeSum)
    H_gtr_dp_pions_vs_beta_pions_cut_all.Fill(event.H_gtr_dp, event.H_gtr_beta) 
    P_gtr_dp_pions_vs_beta_pions_cut_all.Fill(event.P_gtr_dp, event.P_gtr_beta) 
    P_kin_MMpi_pions_vs_beta_pions_cut_all.Fill(event.MMpi, event.P_gtr_beta) 

for event in Cut_Pion_Events_Prompt_tree:
    #P_gtr_beta_pions_cut_prompt.Fill(event.P_gtr_beta)
    #P_RFTime_Dist_pions_cut_prompt.Fill(event.P_RF_Dist)
    CTime_ePiCoinTime_ROC1_pions_cut_prompt.Fill(event.CTime_ePiCoinTime_ROC1)
    P_kin_MMpi_pions_cut_prompt.Fill(event.MMpi)
    #P_kin_MMpi_vs_CTime_ePiCoinTime_ROC1_pions_cut_prompt.Fill(event.MMpi, event.CTime_ePiCoinTime_ROC1)

for event in Cut_Pion_Events_Random_tree:
    #P_gtr_beta_pions_cut_random.Fill(event.P_gtr_beta)
    #P_RFTime_Dist_pions_cut_random.Fill(event.P_RF_Dist)
    CTime_ePiCoinTime_ROC1_pions_cut_random.Fill(event.CTime_ePiCoinTime_ROC1)
    P_kin_MMpi_pions_cut_random.Fill(event.MMpi)

print("Histograms filled")

##############################################################################################################################################

# Random subtraction from missing mass and coin_Time
for event in Cut_Pion_Events_Random_tree:
    P_kin_MMpi_pions_cut_random_scaled.Fill(event.MMpi)
    P_kin_MMpi_pions_cut_random_scaled.Scale(1.0/nWindows)
P_kin_MMpi_pions_cut_random_sub.Add(P_kin_MMpi_pions_cut_prompt, P_kin_MMpi_pions_cut_random_scaled, 1, -1)

############################################################################################################################################

# Saving histograms in PDF
c1_kin = TCanvas("c1_kin", "Kinematic Distributions", 100, 0, 1000, 900)
c1_kin.Divide(2,2)
c1_kin.cd(1)
Q2vsW.Draw("COLZ")
c1_kin.cd(2)
epsilon_pions_cut_all.Draw()
c1_kin.cd(3)
phiqvst.GetYaxis().SetRangeUser(minrangeuser,maxrangeuser)
phiqvst.Draw("SURF2 POL")
# Section for polar plotting
gStyle.SetPalette(55)
gPad.SetTheta(90)
gPad.SetPhi(180)
tvsphi_title = TPaveText(0.0277092,0.89779,0.096428,0.991854,"NDC")
tvsphi_title.AddText("-t vs #phi")
tvsphi_title.Draw()
ptphizero = TPaveText(0.923951,0.513932,0.993778,0.574551,"NDC")
ptphizero.AddText("#phi = 0")
ptphizero.Draw()
phihalfpi = TLine(0,0,0,0.6)
phihalfpi.SetLineColor(kBlack)
phihalfpi.SetLineWidth(2)
phihalfpi.Draw()
ptphihalfpi = TPaveText(0.417855,0.901876,0.486574,0.996358,"NDC")
ptphihalfpi.AddText("#phi = #frac{#pi}{2}")
ptphihalfpi.Draw()
phipi = TLine(0,0,-0.6,0)
phipi.SetLineColor(kBlack)
phipi.SetLineWidth(2)
phipi.Draw()
ptphipi = TPaveText(0.0277092,0.514217,0.096428,0.572746,"NDC")
ptphipi.AddText("#phi = #pi")
ptphipi.Draw()
phithreepi = TLine(0,0,0,-0.6)
phithreepi.SetLineColor(kBlack)
phithreepi.SetLineWidth(2)
phithreepi.Draw()
ptphithreepi = TPaveText(0.419517,0.00514928,0.487128,0.0996315,"NDC")
ptphithreepi.AddText("#phi = #frac{3#pi}{2}")
ptphithreepi.Draw()
Arc = TArc()
for k in range(0, 7):
     Arc.SetFillStyle(0)
     Arc.SetLineWidth(2)
     # To change the arc radius we have to change number 0.825 in the lower line.
     Arc.DrawArc(0,0,0.825*(k+1)/(10),0.,360.,"same")
tradius = TGaxis(0,0,0.575,0,0,0.7,10,"-+")
tradius.SetLineColor(2)
tradius.SetLabelColor(2)
tradius.Draw()
phizero = TLine(0,0,0.6,0) 
phizero.SetLineColor(kBlack)
phizero.SetLineWidth(2)
phizero.Draw()
# End of polar plotting section
c1_kin.cd(4)
P_kin_MMpi_pions_cut_random_sub.Draw("hist")
# Section for Neutron Peak Events Selection
shadedpeak = P_kin_MMpi_pions_cut_random_sub.Clone()
shadedpeak.SetFillColor(2)
shadedpeak.SetFillStyle(3244)
shadedpeak.GetXaxis().SetRangeUser(minbin, maxbin)
shadedpeak.Draw("samehist")
NeutronEvt = TPaveText(0.58934,0.675,0.95,0.75,"NDC")
BinLow = P_kin_MMpi_pions_cut_random_sub.GetXaxis().FindBin(minbin)
BinHigh = P_kin_MMpi_pions_cut_random_sub.GetXaxis().FindBin(maxbin)
BinIntegral = int(P_kin_MMpi_pions_cut_random_sub.Integral(BinLow, BinHigh))
NeutronEvt.SetLineColor(2)
NeutronEvt.AddText("e #pi n Events: %i" %(BinIntegral))
NeutronEvt.Draw()
# End of Neutron Peak Events Selection Section
c1_kin.Print(Pion_Analysis_Distributions + '(')

c1_acpt = TCanvas("c1_H_kin", "Electron-Pion Acceptance Distributions", 100, 0, 1000, 900)
c1_acpt.Divide(3,2)
c1_acpt.cd(1)
gPad.SetLogy()
H_gtr_xp_pions_uncut.SetLineColor(2)
H_gtr_xp_pions_uncut.Draw()
H_gtr_xp_pions_cut_all.SetLineColor(4)
H_gtr_xp_pions_cut_all.Draw("same")
c1_acpt.cd(2)
gPad.SetLogy()
H_gtr_yp_pions_uncut.SetLineColor(2)
H_gtr_yp_pions_uncut.Draw()
H_gtr_yp_pions_cut_all.SetLineColor(4)
H_gtr_yp_pions_cut_all.Draw("same")
c1_acpt.cd(3)
gPad.SetLogy()
H_gtr_dp_pions_uncut.SetLineColor(2) 
# 12/09/21 - NH - fixed the plotting range issue for HMS Delta
H_gtr_dp_pions_uncut.SetMinimum(0.1*H_gtr_dp_pions_cut_all.GetMinimum()+1) # min of plot should be one order of magnitude below the min bin in cut distrobution
H_gtr_dp_pions_uncut.SetMaximum(10*H_gtr_dp_pions_uncut.GetBinContent(H_gtr_dp_pions_uncut.GetMaximumBin())) # Max of plot should be 1 order of magnitude greater than the max bin in uncut distrobution
H_gtr_dp_pions_cut_all.SetLineColor(4)
H_gtr_dp_pions_uncut.Draw()
H_gtr_dp_pions_cut_all.Draw("same")

# TLegend (x1, y1, x2, y2) 
legend2 = ROOT.TLegend(0.115, 0.8, 0.6, 0.9)
legend2.AddEntry("H_gtr_dp_pions_uncut", "without cuts", "l")
legend2.AddEntry("H_gtr_dp_pions_cut_all", "with cuts (acpt/RF/PID)", "l")
legend2.Draw("same")
c1_acpt.cd(4)
gPad.SetLogy()
P_gtr_xp_pions_uncut.SetLineColor(2)
P_gtr_xp_pions_uncut.Draw()
P_gtr_xp_pions_cut_all.SetLineColor(4)
P_gtr_xp_pions_cut_all.Draw("same")
c1_acpt.cd(5)
gPad.SetLogy()
P_gtr_yp_pions_uncut.SetLineColor(2)
P_gtr_yp_pions_uncut.Draw()
P_gtr_yp_pions_cut_all.SetLineColor(4)
P_gtr_yp_pions_cut_all.Draw("same")
c1_acpt.cd(6)
gPad.SetLogy()
P_gtr_dp_pions_uncut.SetLineColor(2)
P_gtr_dp_pions_uncut.Draw()
P_gtr_dp_pions_cut_all.SetLineColor(4)
P_gtr_dp_pions_cut_all.Draw("same")
c1_acpt.Print(Pion_Analysis_Distributions)

c1_pid = TCanvas("c1_pid", "Electron-Pion CAL Distributions", 100, 0, 1000, 900)
c1_pid.Divide(2,2)
c1_pid.cd(1)
gPad.SetLogy()
H_cal_etottracknorm_pions_uncut.SetLineColor(2)
H_cal_etottracknorm_pions_uncut.Draw()
H_cal_etottracknorm_pions_cut_all.SetLineColor(4)
H_cal_etottracknorm_pions_cut_all.Draw("same")
legend7 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
legend7.AddEntry("H_cal_etottracknorm_pions_uncut", "without cuts", "l")
legend7.AddEntry("H_cal_etottracknorm_pions_cut_all", "with cuts (acpt/RF/PID)", "l")
legend7.Draw("same")
c1_pid.cd(2)
H_cal_etottracknorm_vs_H_cer_npeSum_pions_cut_all.Draw("COLZ")
c1_pid.cd(3)
gPad.SetLogy()
P_cal_etottracknorm_pions_uncut.SetLineColor(2)
P_cal_etottracknorm_pions_uncut.Draw()
P_cal_etottracknorm_pions_cut_all.SetLineColor(4)
P_cal_etottracknorm_pions_cut_all.Draw("same")
legend8 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
legend8.AddEntry("P_cal_etottracknorm_pions_uncut", "without cuts", "l")
legend8.AddEntry("P_cal_etottracknorm_pions_cut_all", "with cuts (acpt/RF/PID)", "l")
legend8.Draw("same")
c1_pid.cd(4)
#
c1_pid.Print(Pion_Analysis_Distributions)

c1_RF = TCanvas("c1_RF", "Electron-Pion RF Distributions", 100, 0, 1000, 900)
c1_RF.Divide(2,2)
c1_RF.cd(1)
P_RFTime_Dist_pions_uncut.SetLineColor(2)
P_RFTime_Dist_pions_uncut.Draw()
P_RFTime_Dist_pions_cut_all.SetLineColor(4)
P_RFTime_Dist_pions_cut_all.Draw("same")
legend9 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
legend9.AddEntry("P_RFTime_Dist_pions_uncut", "without cuts", "l")
legend9.AddEntry("P_RFTime_Dist_pions_cut_all", "with cuts (acpt/RF/PID)", "l")
legend9.Draw("same")
c1_RF.cd(2)
P_RFTime_Dist_pions_cut_noRF.SetLineColor(2)
P_RFTime_Dist_pions_cut_noRF.Draw()
P_RFTime_Dist_pions_cut_all.SetLineColor(4)
P_RFTime_Dist_pions_cut_all.Draw("same")
legend = ROOT.TLegend(0.115, 0.835, 0.415, 0.9)
legend.AddEntry("P_RFTime_Dist_pions_cut_noRF", "noRF_cuts (acpt/PID)", "l")
legend.AddEntry("P_RFTime_Dist_pions_cut_all", "RF_cuts (acpt/PID)", "l")
legend.Draw("same")
c1_RF.cd(3)
P_RFTime_vs_P_kin_MMpi_pions_uncut.Draw("COLZ")
c1_RF.cd(4)
P_RFTime_vs_P_kin_MMpi_pions_cut_all.Draw("COLZ")
c1_RF.Print(Pion_Analysis_Distributions)

c2_pid = TCanvas("c2_pid", "Electron-Pion Aero/HGC/NGC PID Distributions", 100, 0, 1000, 900)
c2_pid.Divide(2,2)
c2_pid.cd(1)
gPad.SetLogy()
P_hgcer_npeSum_pions_uncut.SetLineColor(2)
P_hgcer_npeSum_pions_uncut.Draw()
P_hgcer_npeSum_pions_cut_all.SetLineColor(4)
P_hgcer_npeSum_pions_cut_all.Draw("same")
legend10 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
legend10.AddEntry("P_hgcer_npeSum_pions_uncut", "without cuts", "l")
legend10.AddEntry("P_hgcer_npeSum_pions_cut_all", "with cuts (acpt/RF/PID)", "l")
legend10.Draw("same")
c2_pid.cd(2)
gPad.SetLogy()
P_aero_npeSum_pions_uncut.SetLineColor(2)
P_aero_npeSum_pions_uncut.Draw()
P_aero_npeSum_pions_cut_all.SetLineColor(4)
P_aero_npeSum_pions_cut_all.Draw("same")
legend11 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
legend11.AddEntry("P_aero_npeSum_pions_uncut", "without cuts", "l")
legend11.AddEntry("P_aero_npeSum_pions_cut_all", "with cuts (acpt/RF/PID)", "l")
legend11.Draw("same")
c2_pid.cd(3)
P_ngcer_npeSum_pions_uncut.SetLineColor(2)
P_ngcer_npeSum_pions_uncut.Draw()
P_ngcer_npeSum_pions_cut_all.SetLineColor(4)
P_ngcer_npeSum_pions_cut_all.Draw("same") 
#legend12 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
#legend12.AddEntry("P_ngcer_npeSum_pions_uncut", "without cuts", "l")
#legend12.AddEntry("P_ngcer_npeSum_pions_cut_all", "with cuts (acpt/RF/PID)", "l")
#legend12.Draw("same")
c2_pid.cd(4)
#
c2_pid.Print(Pion_Analysis_Distributions)

c3_pid = TCanvas("c3_pid", "Electron-Pion Aero/HGC/NGC PID Distributions", 100, 0, 1000, 900)
c3_pid.Divide(2,3)
c3_pid.cd(1)
gPad.SetLogz()
P_hgcer_npeSum_vs_aero_npeSum_pions_uncut.Draw("COLZ")
c3_pid.cd(2)
P_hgcer_npeSum_vs_aero_npeSum_pions_cut_all.Draw("COLZ")
c3_pid.cd(3)
P_ngcer_npeSum_vs_hgcer_npeSum_pions_uncut.Draw("COLZ")
c3_pid.cd(4)
P_ngcer_npeSum_vs_hgcer_npeSum_pions_cut_all.Draw("COLZ")
c3_pid.cd(5)
P_ngcer_npeSum_vs_aero_npeSum_pions_uncut.Draw("COLZ")
c3_pid.cd(6)
P_ngcer_npeSum_vs_aero_npeSum_pions_cut_all.Draw("COLZ")
c3_pid.Print(Pion_Analysis_Distributions)

c1_MM = TCanvas("c1_MM", "Electron-Pion CTime/Missing Mass Distributions", 100, 0, 1000, 900)
c1_MM.Divide(3,2)
c1_MM.cd(1)
CTime_ePiCoinTime_ROC1_pions_uncut.SetLineColor(4)
CTime_ePiCoinTime_ROC1_pions_uncut.Draw()
CTime_ePiCoinTime_ROC1_pions_cut_prompt.SetLineColor(6)
CTime_ePiCoinTime_ROC1_pions_cut_prompt.Draw("same")
CTime_ePiCoinTime_ROC1_pions_cut_random.SetLineColor(8)
CTime_ePiCoinTime_ROC1_pions_cut_random.Draw("same")
legend13 = ROOT.TLegend(0.1, 0.815, 0.48, 0.9)
legend13.AddEntry("CTime_ePiCoinTime_ROC1_pions_uncut", "CT_without cuts", "l")
legend13.AddEntry("CTime_ePiCoinTime_ROC1_pions_cut_prompt", "CT_prompt with cuts (acpt/RF/PID)", "l")
legend13.AddEntry("CTime_ePiCoinTime_ROC1_pions_cut_random", "CT_randoms with cuts (acpt/RF/PID)", "l")
legend13.Draw("same")
c1_MM.cd(2)
CTime_ePiCoinTime_ROC1_pions_cut_all.SetLineColor(4)
CTime_ePiCoinTime_ROC1_pions_cut_all.Draw()
CTime_ePiCoinTime_ROC1_pions_cut_prompt.SetLineColor(6)
CTime_ePiCoinTime_ROC1_pions_cut_prompt.Draw("same")
CTime_ePiCoinTime_ROC1_pions_cut_random.SetLineColor(8)
CTime_ePiCoinTime_ROC1_pions_cut_random.Draw("same")
legend14 = ROOT.TLegend(0.1, 0.815, 0.48, 0.9)
legend14.AddEntry("CTime_ePiCoinTime_ROC1_pions_cut_all", "CT_with cuts (acpt/RF/PID)", "l")
legend14.AddEntry("CTime_ePiCoinTime_ROC1_pions_cut_prompt", "CT_prompt with cuts (acpt/RF/PID)", "l")
legend14.AddEntry("CTime_ePiCoinTime_ROC1_pions_cut_random", "CT_randoms with cuts (acpt/RF/PID)", "l")
legend14.Draw("same")
c1_MM.cd(3)
P_kin_MMpi_pions_uncut.Draw()
# Section for Neutron Peak Events Selection
#shadedpeak = P_kin_MMpi_pions_uncut.Clone()
#shadedpeak.SetFillColor(2)
#shadedpeak.SetFillStyle(3244)
#shadedpeak.GetXaxis().SetRangeUser(minbin, maxbin)
#shadedpeak.Draw("samehist")
#NeutronEvt = TPaveText(0.58934,0.675,0.95,0.75,"NDC")
#BinLow = P_kin_MMpi_pions_uncut.GetXaxis().FindBin(minbin)
#BinHigh = P_kin_MMpi_pions_uncut.GetXaxis().FindBin(maxbin)
#BinIntegral = int(P_kin_MMpi_pions_uncut.Integral(BinLow, BinHigh))
#NeutronEvt.SetLineColor(2)
#NeutronEvt.AddText("e #pi n Events: %i" %(BinIntegral))
#NeutronEvt.Draw()
c1_MM.cd(4)
P_kin_MMpi_pions_cut_all.SetLineColor(4)
P_kin_MMpi_pions_cut_all.Draw()
P_kin_MMpi_pions_cut_prompt.SetLineColor(6)
P_kin_MMpi_pions_cut_prompt.Draw("same")
P_kin_MMpi_pions_cut_random.SetLineColor(8)
P_kin_MMpi_pions_cut_random.Draw("same")
legend15 = ROOT.TLegend(0.4, 0.815, 0.78, 0.9)
legend15.AddEntry("P_kin_MMpi_pions_cut_all", "MM with cuts (acpt/RF/PID)", "l")
legend15.AddEntry("P_kin_MMpi_pions_cut_prompt", "MM_prompt with cuts (acpt/RF/PID)", "l")
legend15.AddEntry("P_kin_MMpi_pions_cut_random", "MM_randoms with cuts (acpt/RF/PID)", "l")
legend15.Draw("same")
c1_MM.cd(5)
P_kin_MMpi_pions_vs_beta_pions_uncut.Draw("COLZ")
c1_MM.cd(6)
P_kin_MMpi_pions_vs_beta_pions_cut_all.Draw("COLZ")
c1_MM.Print(Pion_Analysis_Distributions)

c1_CT = TCanvas("c1_CT", "Electron-Pion CTime Distributions", 100, 0, 1000, 900)
c1_CT.Divide(2,2)
c1_CT.cd(1)
CTime_ePiCoinTime_ROC1_vs_beta_pions_uncut.Draw("COLZ")
# TLine (x1, y1, x2, y2)
LowerPrompt1 = TLine(PromptWindow[0],gPad.GetUymin(),PromptWindow[0],2)
LowerPrompt1.SetLineColor(2)
LowerPrompt1.SetLineWidth(2)
LowerPrompt1.Draw("same")
UpperPrompt1 = TLine(PromptWindow[1],gPad.GetUymin(),PromptWindow[1],2)
UpperPrompt1.SetLineColor(2)
UpperPrompt1.SetLineWidth(2)
UpperPrompt1.Draw("same")
LowerRandomL1 = TLine(RandomWindows[0],gPad.GetUymin(),RandomWindows[0],2)
LowerRandomL1.SetLineColor(8)
LowerRandomL1.SetLineWidth(2)
LowerRandomL1.Draw("same")
UpperRandomL1 = TLine(RandomWindows[1],gPad.GetUymin(),RandomWindows[1],2)
UpperRandomL1.SetLineColor(8)
UpperRandomL1.SetLineWidth(2)
UpperRandomL1.Draw("same")
LowerRandomR1 = TLine(RandomWindows[2],gPad.GetUymin(),RandomWindows[2],2)
LowerRandomR1.SetLineColor(8)
LowerRandomR1.SetLineWidth(2)
LowerRandomR1.Draw("same")
UpperRandomR1 = TLine(RandomWindows[3],gPad.GetUymin(),RandomWindows[3],2)
UpperRandomR1.SetLineColor(8)
UpperRandomR1.SetLineWidth(2)
UpperRandomR1.Draw("same")
c1_CT.cd(2)
CTime_ePiCoinTime_ROC1_vs_beta_pions_cut_all.Draw("COLZ")
c1_CT.cd(3)
CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_uncut.Draw("COLZ")
LowerPrompt2 = TLine(PromptWindow[0],gPad.GetUymin(),PromptWindow[0],2)
LowerPrompt2.SetLineColor(2)
LowerPrompt2.SetLineWidth(2)
LowerPrompt2.Draw("same")
UpperPrompt2 = TLine(PromptWindow[1],gPad.GetUymin(),PromptWindow[1],2)
UpperPrompt2.SetLineColor(2)
UpperPrompt2.SetLineWidth(2)
UpperPrompt2.Draw("same")
LowerRandomL2 = TLine(RandomWindows[0],gPad.GetUymin(),RandomWindows[0],2)
LowerRandomL2.SetLineColor(8)
LowerRandomL2.SetLineWidth(2)
LowerRandomL2.Draw("same")
UpperRandomL2 = TLine(RandomWindows[1],gPad.GetUymin(),RandomWindows[1],2)
UpperRandomL2.SetLineColor(8)
UpperRandomL2.SetLineWidth(2)
UpperRandomL2.Draw("same")
LowerRandomR2 = TLine(RandomWindows[2],gPad.GetUymin(),RandomWindows[2],2)
LowerRandomR2.SetLineColor(8)
LowerRandomR2.SetLineWidth(2)
LowerRandomR2.Draw("same")
UpperRandomR2 = TLine(RandomWindows[3],gPad.GetUymin(),RandomWindows[3],2)
UpperRandomR2.SetLineColor(8)
UpperRandomR2.SetLineWidth(2)
UpperRandomR2.Draw("same")
c1_CT.cd(4)
CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_cut_all.Draw("COLZ")
c1_CT.Print(Pion_Analysis_Distributions)

c1_delta = TCanvas("c1_delta", "Delta Debugging", 100, 0, 1000, 900)
c1_delta.Divide(2,2)
c1_delta.cd(1)
H_gtr_dp_pions_vs_beta_pions_uncut.Draw("COLZ")
c1_delta.cd(2)
H_gtr_dp_pions_vs_beta_pions_cut_all.Draw("COLZ")
c1_delta.cd(3)
P_gtr_dp_pions_vs_beta_pions_uncut.Draw("COLZ")
c1_delta.cd(4)
P_gtr_dp_pions_vs_beta_pions_cut_all.Draw("COLZ")
c1_delta.Print(Pion_Analysis_Distributions + ')')

#############################################################################################################################################

# Making directories in output file
outHistFile = ROOT.TFile.Open("%s/%s_%s_Output_Data.root" % (OUTPATH, runNum, MaxEvent), "RECREATE")                                                                                                    
d_Uncut_Pion_Events = outHistFile.mkdir("Uncut_Pion_Events")
d_Cut_Pion_Events_All = outHistFile.mkdir("Cut_Pion_Events_All")
d_Cut_Pion_Events_Prompt = outHistFile.mkdir("Cut_Pion_Events_Prompt")
d_Cut_Pion_Events_Random = outHistFile.mkdir("Cut_Pion_Events_Random")

# Writing Histograms in output root file
P_HGC_xy_npe_uncut.Write()                                           
P_Aero_xy_npe_uncut.Write();

d_Uncut_Pion_Events.cd()
H_gtr_beta_pions_uncut.Write()
H_gtr_xp_pions_uncut.Write()
H_gtr_yp_pions_uncut.Write()
H_gtr_dp_pions_uncut.Write()
H_hod_goodscinhit_pions_uncut.Write()
H_hod_goodstarttime_pions_uncut.Write()
H_cal_etotnorm_pions_uncut.Write()
H_cal_etottracknorm_pions_uncut.Write()
H_cer_npeSum_pions_uncut.Write()
H_RFTime_Dist_pions_uncut.Write()
P_gtr_beta_pions_uncut.Write()
P_gtr_xp_pions_uncut.Write()
P_gtr_yp_pions_uncut.Write()
P_gtr_dp_pions_uncut.Write()
P_gtr_p_pions_uncut.Write()
P_hod_goodscinhit_pions_uncut.Write()
P_hod_goodstarttime_pions_uncut.Write()
P_cal_etotnorm_pions_uncut.Write()
P_cal_etottracknorm_pions_uncut.Write()
P_hgcer_npeSum_pions_uncut.Write()
P_hgcer_xAtCer_pions_uncut.Write()
P_hgcer_yAtCer_pions_uncut.Write()
P_aero_npeSum_pions_uncut.Write()
P_aero_xAtAero_pions_uncut.Write()
P_aero_yAtAero_pions_uncut.Write()
P_ngcer_npeSum_pions_uncut.Write()
P_ngcer_xAtCer_pions_uncut.Write()
P_ngcer_yAtCer_pions_uncut.Write() 
P_kin_MMpi_pions_uncut.Write()
P_RFTime_Dist_pions_uncut.Write()
CTime_ePiCoinTime_ROC1_pions_uncut.Write()
Q2_pions_uncut.Write()
W_pions_uncut.Write()
epsilon_pions_uncut.Write()
phiq_uncut.Write()
t_uncut.Write()
H_cal_etottracknorm_vs_H_cer_npeSum_pions_uncut.Write()
P_hgcer_npeSum_vs_aero_npeSum_pions_uncut.Write()
CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_uncut.Write()
P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_uncut.Write()
P_aero_yAtAero_vs_aero_xAtAero_pions_uncut.Write()
CTime_ePiCoinTime_ROC1_vs_beta_pions_uncut.Write()
P_RFTime_vs_P_kin_MMpi_pions_uncut.Write()
P_cal_etottracknorm_vs_P_ngcer_npeSum_pions_uncut.Write()
P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_uncut.Write()
P_ngcer_npeSum_vs_hgcer_npeSum_pions_uncut.Write()
P_ngcer_npeSum_vs_aero_npeSum_pions_uncut.Write()
H_gtr_dp_pions_vs_beta_pions_uncut.Write()
P_gtr_dp_pions_vs_beta_pions_uncut.Write()
P_kin_MMpi_pions_vs_beta_pions_uncut.Write()

d_Cut_Pion_Events_All.cd()
H_gtr_beta_pions_cut_all.Write()
H_gtr_xp_pions_cut_all.Write()
H_gtr_yp_pions_cut_all.Write()
H_gtr_dp_pions_cut_all.Write()
H_hod_goodscinhit_pions_cut_all.Write()
H_hod_goodstarttime_pions_cut_all.Write()
H_cal_etotnorm_pions_cut_all.Write()
H_cal_etottracknorm_pions_cut_all.Write()
H_cer_npeSum_pions_cut_all.Write()
H_RFTime_Dist_pions_cut_all.Write()
P_gtr_beta_pions_cut_all.Write()
P_gtr_xp_pions_cut_all.Write()
P_gtr_yp_pions_cut_all.Write()
P_gtr_dp_pions_cut_all.Write()
P_gtr_p_pions_cut_all.Write()
P_hod_goodscinhit_pions_cut_all.Write()
P_hod_goodstarttime_pions_cut_all.Write()
P_cal_etotnorm_pions_cut_all.Write()
P_cal_etottracknorm_pions_cut_all.Write()
P_hgcer_npeSum_pions_cut_all.Write()
P_hgcer_xAtCer_pions_cut_all.Write()
P_hgcer_yAtCer_pions_cut_all.Write()
P_aero_npeSum_pions_cut_all.Write()
P_aero_xAtAero_pions_cut_all.Write()
P_aero_yAtAero_pions_cut_all.Write()
P_ngcer_npeSum_pions_cut_all.Write()
P_ngcer_xAtCer_pions_cut_all.Write()
P_ngcer_yAtCer_pions_cut_all.Write()
P_kin_MMpi_pions_cut_all.Write()
P_RFTime_Dist_pions_cut_all.Write()
CTime_ePiCoinTime_ROC1_pions_cut_all.Write()
Q2_pions_cut_all.Write()
W_pions_cut_all.Write()
epsilon_pions_cut_all.Write()
phiq_cut_all.Write()
t_cut_all.Write()
Q2vsW.Write()
phiqvst.Write()
H_cal_etottracknorm_vs_H_cer_npeSum_pions_cut_all.Write()
P_hgcer_npeSum_vs_aero_npeSum_pions_cut_all.Write()
CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_cut_all.Write()
P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_cut_all.Write()
P_aero_yAtAero_vs_aero_xAtAero_pions_cut_all.Write()
CTime_ePiCoinTime_ROC1_vs_beta_pions_cut_all.Write()
P_RFTime_vs_P_kin_MMpi_pions_cut_all.Write()
P_kin_MMpi_pions_cut_random_sub.Write()
P_cal_etottracknorm_vs_P_ngcer_npeSum_pions_cut_all.Write()
P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_cut_all.Write()
P_ngcer_npeSum_vs_hgcer_npeSum_pions_cut_all.Write()
P_ngcer_npeSum_vs_aero_npeSum_pions_cut_all.Write()
H_gtr_dp_pions_vs_beta_pions_cut_all.Write()
P_gtr_dp_pions_vs_beta_pions_cut_all.Write()
P_kin_MMpi_pions_vs_beta_pions_cut_all.Write()

d_Cut_Pion_Events_Prompt.cd()
P_gtr_beta_pions_cut_prompt.Write()
P_RFTime_Dist_pions_cut_prompt.Write()
CTime_ePiCoinTime_ROC1_pions_cut_prompt.Write()
P_kin_MMpi_pions_cut_prompt.Write()
P_kin_MMpi_vs_CTime_ePiCoinTime_ROC1_pions_cut_prompt.Write()

d_Cut_Pion_Events_Random.cd()
P_gtr_beta_pions_cut_random.Write()
P_RFTime_Dist_pions_cut_random.Write()
CTime_ePiCoinTime_ROC1_pions_cut_random.Write()
P_kin_MMpi_pions_cut_random.Write() 

outHistFile.Close()
infile.Close() 
print ("Processing Complete")
print("!!!!!!!!\n %i pi-n events \n!!!!!!!!" % BinIntegral)
