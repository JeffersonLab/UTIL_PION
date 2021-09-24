#! /usr/bin/python
###########################################################################################################################
# Created - 20/July/21, Author - Muhammad Junaid (mjo147@uregina.ca), University of Regina, Canada (Copyright (c) junaid) #
###########################################################################################################################
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
maxrangeuser = 1.2 # max range for -t vs phi plot
minbin = 0.92 # minimum bin for selecting neutrons events in missing mass distribution
maxbin = 0.98 # maximum bin for selecting neutrons events in missing mass distribution

##################################################################################################################################################

# Check the number of arguments provided to the script
FilenameOverride=False # SJDK 21/09/21 - Added a secret 4th argument so that a full kinematic can be processed, the run number is needed for cuts, but for the kinematic analysis, the filename is not by run number
if len(sys.argv)-1!=3:
    if len(sys.argv)-1!=4:
        print("!!!!! ERROR !!!!!\n Expected 3 arguments\n Usage is with - ROOTfileSuffix RunNumber MaxEvents \n!!!!! ERROR !!!!!")
        sys.exit(1)
    else:
        print ("!!!!! Running with secret 4th argument - FilenameOverride - Taking file name to process as stated EXACTLY in 4th arg !!!!!")
        FilenameOverride=sys.argv[4] # If 4 arguments provided, set the FilenameOverride value to be arg 4

##################################################################################################################################################

# Input params - run number and max number of events
ROOTSuffix = sys.argv[1]
runNum = sys.argv[2]
MaxEvent = sys.argv[3]

USER = subprocess.getstatusoutput("whoami") # Grab user info for file finding
HOST = subprocess.getstatusoutput("hostname")

if ("farm" in HOST[1]):
    REPLAYPATH = "/group/c-pionlt/online_analysis/hallc_replay_lt"
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
if (FilenameOverride == False): # Standard running condition, construct file name from run number and max events e.t.c.
    rootName = "%s/UTIL_PION/OUTPUT/Analysis/PionLT/%s_%s_%s.root" % (REPLAYPATH, runNum, MaxEvent, ROOTSuffix)     # Input file location and variables taking
elif (FilenameOverride != False): # Special condition, with 4th arg, use 4th arg as file name
    rootName = "%s/UTIL_PION/OUTPUT/Analysis/PionLT/%s" % (REPLAYPATH, FilenameOverride)
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
P_RFTime_pions_cut_noRF = ROOT.TH1D("P_RFTime_pions_cut_noRF", "SHMS RFTime; SHMS_RFTime; Counts", 200, 0, 4)

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
H_cer_npe_pions_cut = ROOT.TH1D("H_cer_npe_pions_cut", "HMS cer npeSum; HMS_cer_npeSum; Counts", 200, 0, 50)
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
P_hgcer_npe_pions_cut = ROOT.TH1D("P_hgcer_npe_pions_cut", "SHMS HGC npeSum; SHMS_hgcer_npeSum; Counts", 200, 0, 50)
P_xhgcer_pions_cut = ROOT.TH1D("P_xhgcer_pions_cut", "SHMS HGC xAtCer; SHMS_hgcer_xAtCer; Counts", 200, -40, 30)
P_yhgcer_pions_cut = ROOT.TH1D("P_yhgcer_pions_cut", "SHMS HGC yAtCer; SHMS_hgcer_yAtCer; Counts", 200, -30, 30)
P_aero_npe_pions_cut = ROOT.TH1D("P_aero_npe_pions_cut", "SHMS aero npeSum; SHMS_aero_npeSum; Counts", 200, 0, 50)
P_xaero_pions_cut = ROOT.TH1D("P_xaero_pions_cut", "SHMS aero xAtAero; SHMS_aero_xAtAero; Counts", 200, -40, 30)
P_yaero_pions_cut = ROOT.TH1D("P_yaero_pions_cut", "SHMS aero yAtAero; SHMS_aero_yAtAero; Counts", 200, -30, 30)
P_ngcer_npe_pions_cut = ROOT.TH1D("P_ngcer_npe_pions_cut", "SHMS NGC npeSum; SHMS_ngcer_npeSum; Counts", 200, 0, 50)
P_xngcer_pions_cut = ROOT.TH1D("P_xngcer_pions_cut", "SHMS NGC xAtCer; SHMS_ngcer_xAtCer; Counts", 200, -40, 30)
P_yngcer_pions_cut = ROOT.TH1D("P_yngcer_pions_cut", "SHMS NGC yAtCer; SHMS_ngcer_yAtCer; Counts", 200, -30, 30)
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
P_MMpi_pions_cut_prompt = ROOT.TH1D("P_MMpi_pions_cut_prompt", "Missing Mass; MM_{#pi}; Counts", 200, 0.5, 1.8)

P_beta_pions_cut_randm = ROOT.TH1D("P_beta_pions_cut_randm", "SHMS beta; SHMS_#beta; Counts", 200, 0.8, 1.2)
P_RFTime_pions_cut_randm = ROOT.TH1D("P_RFTime_pions_cut_randm", "SHMS RFTime; SHMS_RFTime; Counts", 200, 0, 4)
ePiCoinTime_pions_cut_randm = ROOT.TH1D("ePiCoinTime_pions_cut_randm", "Electron-Pion CTime; e #pi Coin_Time; Counts", 160, -40, 40)
P_MMpi_pions_cut_randm = ROOT.TH1D("P_MMpi_pions_cut_randm", "Missing Mass; MM_{#pi}; Counts", 200, 0.5, 1.8)

P_MMpi_pions_cut_randm_scaled = ROOT.TH1D("P_MMpi_pions_cut_randm_scaled", "Missing Mass; MM_{#pi}; Counts", 200, 0.5, 1.8)
P_MMpi_pions_cut_randm_sub = ROOT.TH1D("P_MMpi_pions_cut_randm_sub", "Missing Mass Rndm Sub; MM_{#pi}; Counts", 200, 0.5, 1.8)

##############################################################################################################################################

# 2D Histograms for pions
H_cal_etottracknorm_vs_cer_npe_pions_uncut = ROOT.TH2D("H_cal_etottracknorm_vs_cer_npe_pions_uncut","HMS cal etottracknorm vs HMS cer npeSum (no cut); H_cal_etottracknorm; H_cer_npeSum",100, 0, 2, 100, 0, 40)
P_hgcer_vs_aero_npe_pions_uncut = ROOT.TH2D("P_hgcer_vs_aero_npe_pions_uncut", "SHMS HGC npeSum vs SHMS Aero npeSum (no cut); SHMS_hgcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
ePiCoinTime_vs_MMpi_pions_uncut = ROOT.TH2D("ePiCoinTime_vs_MMpi_pions_uncut","Electron-Pion CTime vs Missing Mass (no cut); e #pi Coin_Time; MM_{#pi}", 200, -40, 40, 200, 0, 2)
P_hgcer_yx_pions_uncut = ROOT.TH2D("P_hgcer_yx_pions_uncut", "SHMS HGC yAtCer vs SHMS HGC xAtCer (no cut); SHMS_hgcer_yAtCer; SHMS_hgcer_xAtCer", 100, -50, 50, 100, -50, 50)
P_aero_yx_pions_uncut = ROOT.TH2D("P_aero_yx_pions_uncut", "SHMS aero yAtAero vs SHMS aero xAtAero (no cut); SHMS_aero_yAtAero; SHMS_aero_xAtAero", 100, -50, 50, 100, -50, 50)
ePiCoinTime_vs_beta_pions_uncut = ROOT.TH2D("ePiCoinTime_vs_beta_pions_uncut", "Electron-Pion CTime vs SHMS #beta (no cut); e #pi Coin_Time; SHMS_#beta", 200, -40, 40, 200, 0, 2)
P_RFTime_vs_MMpi_pions_uncut = ROOT.TH2D("P_RFTime_vs_MMpi_pions_uncut", "SHMS RFTime vs Missing Mass (no cuts); SHMS_RFTime_Dist; MM_{#pi}", 100, 0, 4, 100, 0, 2)
P_cal_etottracknorm_vs_ngcer_npe_pions_uncut = ROOT.TH2D("P_cal_etottracknorm_vs_ngcer_npe_pions_uncut", "SHMS cal etottracknorm vs SHMS NGC xAtCer (no cut); SHMS_cal_etottracknorm; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
P_ngcer_yx_pions_uncut = ROOT.TH2D("P_ngcer_yx_pions_uncut", "SHMS NGC yAtCer vs SHMS NGC xAtCer (no cut); SHMS_ngcer_yAtCer; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
P_ngcer_vs_hgcer_npe_pions_uncut = ROOT.TH2D("P_ngcer_vs_hgcer_npe_pions_uncut", "SHMS NGC npeSum vs SHMS HGC npeSum (no cut); SHMS_ngcer_npeSum; SHMS_hgcer_npeSum", 100, 0, 50, 100, 0, 50)
P_ngcer_vs_aero_npe_pions_uncut = ROOT.TH2D("P_ngcer_vs_aero_npe_pions_uncut", "SHMS NGC npeSum vs SHMS aero npeSum (no cut); SHMS_ngcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
H_dp_vs_beta_pions_uncut = ROOT.TH2D("H_dp_vs_beta_pions_uncut", "HMS #delta vs HMS #beta (no cut); HMS #delta; HMS_#beta", 200, -12, 12, 200, 0, 2)
P_dp_vs_beta_pions_uncut = ROOT.TH2D("P_dp_vs_beta_pions_uncut", "SHMS #delta vs SHMS #beta (no cut); SHMS #delta; HMS_#beta", 200, -12, 12, 200, 0, 2)
P_MMpi_vs_beta_pions_uncut = ROOT.TH2D("P_MMpi_vs_beta_pions_uncut", "Missing Mass vs SHMS #beta (no cut); MM_{#pi}; SHMS_#beta", 100, 0, 2, 200, 0, 2)

H_cal_etottracknorm_vs_cer_npe_pions_cut = ROOT.TH2D("H_cal_etottracknorm_vs_cer_npe_pions_cut","HMS cal etottracknorm vs HMS cer npeSum (with cuts); H_cal_etottracknorm; H_cer_npeSum",100, 0.5, 1.5, 100, 0, 40)
P_hgcer_vs_aero_npe_pions_cut = ROOT.TH2D("P_hgcer_vs_aero_npe_pions_cut", "SHMS HGC npeSum vs SHMS aero npeSum (with cuts); SHMS_hgcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
ePiCoinTime_vs_MMpi_pions_cut = ROOT.TH2D("ePiCoinTime_vs_MMpi_pions_cut","Electron-Pion CTime vs Missing Mass (with PID cuts); e #pi Coin_Time; MM_{#pi}", 100, -2, 2, 100, 0, 2)
P_hgcer_yx_pions_cut = ROOT.TH2D("P_hgcer_yx_pions_cut", "SHMS HGC yAtCer vs SHMS HGC xAtCer (with cuts); SHMS_hgcer_yAtCer; SHMS_hgcer_xAtCer", 100, -50, 50, 100, -50, 50)
P_aero_yx_pions_cut = ROOT.TH2D("P_aero_yx_pions_cut", "SHMS aero yAtAero vs SHMS aero xAtAero (with cuts); SHMS_aero_yAtAero; SHMS_aero_xAtAero", 100, -50, 50, 100, -50, 50)
ePiCoinTime_vs_beta_pions_cut = ROOT.TH2D("ePiCoinTime_vs_beta_pions_cut", "Electron-Pion CTime vs SHMS #beta (with PID cuts); e #pi Coin_Time; SHMS_#beta", 100, -2, 2, 100, 0.6, 1.4)
P_RFTime_vs_MMpi_pions_cut = ROOT.TH2D("P_RFTime_vs_MMpi_pions_cut", "SHMS RFTime vs Missing Mass (with cuts); SHMS_RFTime_Dist; MM_{#pi}", 100, 0, 4, 100, 0, 2)
P_cal_etottracknorm_vs_ngcer_npe_pions_cut = ROOT.TH2D("P_cal_etottracknorm_vs_ngcer_npe_pions_cut", "P cal etottracknorm vs SHMS NGC xAtCer (with cuts); SHMS_cal_etottracknorm; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
P_ngcer_yx_pions_cut = ROOT.TH2D("P_ngcer_yx_pions_cut", "SHMS NGC yAtCer vs SHMS NGC xAtCer (with cuts); SHMS_ngcer_yAtCer; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
P_ngcer_vs_hgcer_npe_pions_cut = ROOT.TH2D("P_ngcer_vs_hgcer_npe_pions_cut", "SHMS NGC npeSum vs SHMS HGC npeSum (with cuts); SHMS_ngcer_npeSum; SHMS_hgcer_npeSum", 100, 0, 50, 100, 0, 50)
P_ngcer_vs_aero_npe_pions_cut = ROOT.TH2D("P_ngcer_vs_aero_npe_pions_cut", "SHMS NGC npeSum vs SHMS aero npeSum (with cuts); SHMS_ngcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
H_dp_vs_beta_pions_cut = ROOT.TH2D("H_dp_vs_beta_pions_cut", "HMS #delta vs HMS #beta (with cut); HMS #delta; HMS_#beta", 200, -12, 12, 200, 0, 2)
P_dp_vs_beta_pions_cut = ROOT.TH2D("P_dp_vs_beta_pions_cut", "SHMS #delta vs SHMS #beta (with cut); SHMS #delta; SHMS_#beta", 200, -30, 30, 200, 0, 2)
P_MMpi_vs_beta_pions_cut = ROOT.TH2D("P_MMpi_vs_beta_pions_cut", "Missing Mass vs SHMS #beta (with cut); MM_{#pi}; SHMS_#beta", 100, 0, 2, 200, 0, 2)

MMpi_vs_ePiCoinTime_pions_cut_prompt = ROOT.TH2D("MMpi_vs_ePiCoinTime_pions_cut_prompt","Missing Mass vs Electron-Pion CTime; MM_{#pi}; e #pi Coin_Time",100, 0, 2, 100, -2, 2)
Q2vsW_pions_cut = ROOT.TH2D("Q2vsW_pions_cut", "Q2 vs W; Q2; W", 200, 3.0, 8.0, 200, 2.7, 3.6)
phiqvst_pions_cut = ROOT.TH2D("phiqvst_pions_cut","; #phi ;t", 12, -3.14, 3.14, 24, 0.0, 1.2)

# 11/09/21 - SJDK - Adding some 3D XY NPE plots, need to take projections of these (which I need to figure out how to do in PyRoot. Making these manually from the command line for now.
P_HGC_xy_npe_pions_uncut = ROOT.TH3D("P_HGC_xy_npe_pions_uncut", "SHMS HGC NPE as fn of yAtCer vs SHMS HGC xAtCer (no cuts); HGC_yAtCer(cm); HGC_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_Aero_xy_npe_pions_uncut = ROOT.TH3D("P_Aero_xy_npe_pions_uncut", "SHMS Aerogel NPE as fn of yAtCer vs xAtCer (no cuts); Aero_yAtCer(cm); Aero_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_NGC_xy_npe_pions_uncut = ROOT.TH3D("P_NGC_xy_npe_pions_uncut", "SHMS NGC NPE as fn of yAtCer vs xAtCer (no cuts); NGC_yAtCer(cm); NGC_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_HGC_xy_npe_pions_cut = ROOT.TH3D("P_HGC_xy_npe_pions_cut", "SHMS HGC NPE as fn of yAtCer vs SHMS HGC xAtCer (with cuts); HGC_yAtCer(cm); HGC_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_Aero_xy_npe_pions_cut = ROOT.TH3D("P_Aero_xy_npe_pions_cut", "SHMS Aerogel NPE as fn of yAtCer vs xAtCer (with cuts); Aero_yAtCer(cm); Aero_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_NGC_xy_npe_pions_cut = ROOT.TH3D("P_NGC_xy_npe_pions_cut", "SHMS NGC NPE as fn of yAtCer vs xAtCer (with cuts); NGC_yAtCer(cm); NGC_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
#################################################################################################################################################

# Filling Histograms for Pions
for event in Cut_Pion_Events_noRF_tree:
    P_RFTime_pions_cut_noRF.Fill(event.P_RF_Dist)

for event in Uncut_Pion_Events_tree:
#    H_beta_pions_uncut.Fill(event.H_gtr_beta)
    H_xp_pions_uncut.Fill(event.H_gtr_xp)
    H_yp_pions_uncut.Fill(event.H_gtr_yp)
    H_dp_pions_uncut.Fill(event.H_gtr_dp)
#    H_hod_goodscinhit_pions_uncut.Fill(event.H_hod_goodscinhit)
#    H_hod_goodstarttime_pions_uncut.Fill(event.H_hod_goodstarttime)
#    H_cal_etotnorm_pions_uncut.Fill(event.H_cal_etotnorm)
    H_cal_etottracknorm_pions_uncut.Fill(event.H_cal_etottracknorm)
#    H_cer_npe_pions_uncut.Fill(event.H_cer_npeSum)
#    H_RFTime_pions_uncut.Fill(event.H_RF_Dist)
#    P_beta_pions_uncut.Fill(event.P_gtr_beta)
    P_xp_pions_uncut.Fill(event.P_gtr_xp)
    P_yp_pions_uncut.Fill(event.P_gtr_yp)
    P_dp_pions_uncut.Fill(event.P_gtr_dp)
#    P_p_pions_uncut.Fill(event.P_gtr_p)
#    P_hod_goodscinhit_pions_uncut.Fill(event.P_hod_goodscinhit)
#    P_hod_goodstarttime_pions_uncut.Fill(event.P_hod_goodstarttime)
#    P_cal_etotnorm_pions_uncut.Fill(event.P_cal_etotnorm)
    P_cal_etottracknorm_pions_uncut.Fill(event.P_cal_etottracknorm)
    P_hgcer_npe_pions_uncut.Fill(event.P_hgcer_npeSum)
#    P_xhgcer_pions_uncut.Fill(event.P_hgcer_xAtCer)
#    P_yhgcer_pions_uncut.Fill(event.P_hgcer_yAtCer)
    P_aero_npe_pions_uncut.Fill(event.P_aero_npeSum)
#    P_xaero_pions_uncut.Fill(event.P_aero_xAtAero)
#    P_yaero_pions_uncut.Fill(event.P_aero_yAtAero)
    P_ngcer_npe_pions_uncut.Fill(event.P_ngcer_npeSum)
#    P_xngcer_pions_uncut.Fill(event.P_ngcer_xAtCer)
#    P_yngcer_pions_uncut.Fill(event.P_ngcer_yAtCer)
    P_MMpi_pions_uncut.Fill(event.MMpi)
    P_RFTime_pions_uncut.Fill(event.P_RF_Dist)
    ePiCoinTime_pions_uncut.Fill(event.CTime_ePiCoinTime_ROC1)
#    Q2_pions_uncut.Fill(event.Q2)
#    W_pions_uncut.Fill(event.W)
#    epsilon_pions_uncut.Fill(event.epsilon)
#    phiq_pions_uncut.Fill(event.ph_q)
#    t_pions_uncut.Fill(-event.MandelT)
#    H_cal_etottracknorm_vs_cer_npe_pions_uncut.Fill(event.H_cal_etottracknorm, event.H_cer_npeSum)
    P_hgcer_vs_aero_npe_pions_uncut.Fill(event.P_hgcer_npeSum, event.P_aero_npeSum)
    ePiCoinTime_vs_MMpi_pions_uncut.Fill(event.CTime_ePiCoinTime_ROC1, event.MMpi)
    P_RFTime_vs_MMpi_pions_uncut.Fill(event.P_RF_Dist, event.MMpi)
#    P_hgcer_yx_pions_uncut.Fill(event.P_hgcer_yAtCer, event.P_hgcer_xAtCer)
#    P_aero_yx_pions_uncut.Fill(event.P_aero_yAtAero, event.P_aero_xAtAero)
    ePiCoinTime_vs_beta_pions_uncut.Fill(event.CTime_ePiCoinTime_ROC1, event.P_gtr_beta)
    P_cal_etottracknorm_vs_ngcer_npe_pions_uncut.Fill(event.P_cal_etottracknorm, event.P_ngcer_npeSum)
#    P_ngcer_yx_pions_uncut.Fill(event.P_ngcer_yAtCer, event.P_ngcer_xAtCer)
    P_ngcer_vs_hgcer_npe_pions_uncut.Fill(event.P_ngcer_npeSum, event.P_hgcer_npeSum)
    P_ngcer_vs_aero_npe_pions_uncut.Fill(event.P_ngcer_npeSum, event.P_aero_npeSum)
    H_dp_vs_beta_pions_uncut.Fill(event.H_gtr_dp, event.H_gtr_beta) 
    P_dp_vs_beta_pions_uncut.Fill(event.P_gtr_dp, event.P_gtr_beta) 
    P_MMpi_vs_beta_pions_uncut.Fill(event.MMpi, event.P_gtr_beta)
    P_HGC_xy_npe_pions_uncut.Fill(event.P_hgcer_yAtCer,event.P_hgcer_xAtCer,event.P_hgcer_npeSum)
    P_Aero_xy_npe_pions_uncut.Fill(event.P_aero_yAtAero,event.P_aero_xAtAero,event.P_aero_npeSum)
    P_NGC_xy_npe_pions_uncut.Fill(event.P_ngcer_yAtCer,event.P_ngcer_xAtCer,event.P_ngcer_npeSum)

for event in Cut_Pion_Events_All_tree:
#    H_beta_pions_cut.Fill(event.H_gtr_beta)
    H_xp_pions_cut.Fill(event.H_gtr_xp)
    H_yp_pions_cut.Fill(event.H_gtr_yp)
    H_dp_pions_cut.Fill(event.H_gtr_dp)
#    H_hod_goodscinhit_pions_cut.Fill(event.H_hod_goodscinhit)
#    H_hod_goodstarttime_pions_cut.Fill(event.H_hod_goodstarttime)
#    H_cal_etotnorm_pions_cut.Fill(event.H_cal_etotnorm)
    H_cal_etottracknorm_pions_cut.Fill(event.H_cal_etottracknorm)
#    H_cer_npe_pions_cut.Fill(event.H_cer_npeSum)
#    H_RFTime_pions_cut.Fill(event.H_RF_Dist)
#    P_beta_pions_cut.Fill(event.P_gtr_beta)
    P_xp_pions_cut.Fill(event.P_gtr_xp)
    P_yp_pions_cut.Fill(event.P_gtr_yp)
    P_dp_pions_cut.Fill(event.P_gtr_dp)
#    P_p_pions_cut.Fill(event.P_gtr_p)
#    P_hod_goodscinhit_pions_cut.Fill(event.P_hod_goodscinhit)
#    P_hod_goodstarttime_pions_cut.Fill(event.P_hod_goodstarttime)
#    P_cal_etotnorm_pions_cut.Fill(event.P_cal_etotnorm)
    P_cal_etottracknorm_pions_cut.Fill(event.P_cal_etottracknorm)
    P_hgcer_npe_pions_cut.Fill(event.P_hgcer_npeSum)
#    P_xhgcer_pions_cut.Fill(event.P_hgcer_xAtCer)
#    P_yhgcer_pions_cut.Fill(event.P_hgcer_yAtCer)
    P_aero_npe_pions_cut.Fill(event.P_aero_npeSum)
#    P_xaero_pions_cut.Fill(event.P_aero_xAtAero)
#    P_yaero_pions_cut.Fill(event.P_aero_yAtAero)
    P_ngcer_npe_pions_cut.Fill(event.P_ngcer_npeSum)
#    P_xngcer_pions_cut.Fill(event.P_ngcer_xAtCer)
#    P_yngcer_pions_cut.Fill(event.P_ngcer_yAtCer)   
    P_MMpi_pions_cut.Fill(event.MMpi)
    P_RFTime_pions_cut.Fill(event.P_RF_Dist)
    ePiCoinTime_pions_cut.Fill(event.CTime_ePiCoinTime_ROC1)
#    Q2_pions_cut.Fill(event.Q2)
#    W_pions_cut.Fill(event.W)
    epsilon_pions_cut.Fill(event.epsilon)
#    phiq_pions_cut.Fill(event.ph_q)
#    t_pions_cut.Fill(-event.MandelT)
    H_cal_etottracknorm_vs_cer_npe_pions_cut.Fill(event.H_cal_etottracknorm, event.H_cer_npeSum)    
    P_hgcer_vs_aero_npe_pions_cut.Fill(event.P_hgcer_npeSum, event.P_aero_npeSum)
#    P_hgcer_yx_pions_cut.Fill(event.P_hgcer_yAtCer, event.P_hgcer_xAtCer)
#    P_aero_yx_pions_cut.Fill(event.P_aero_yAtAero, event.P_aero_xAtAero)
    ePiCoinTime_vs_beta_pions_cut.Fill(event.CTime_ePiCoinTime_ROC1, event.P_gtr_beta)
    ePiCoinTime_vs_MMpi_pions_cut.Fill(event.CTime_ePiCoinTime_ROC1, event.MMpi)
    P_RFTime_vs_MMpi_pions_cut.Fill(event.P_RF_Dist, event.MMpi)    
    Q2vsW_pions_cut.Fill(event.Q2, event.W)
    phiqvst_pions_cut.Fill(event.ph_q, -event.MandelT)
    P_cal_etottracknorm_vs_ngcer_npe_pions_cut.Fill(event.P_cal_etottracknorm, event.P_ngcer_npeSum)
#    P_ngcer_yx_pions_cut.Fill(event.P_ngcer_yAtCer, event.P_ngcer_xAtCer)
    P_ngcer_vs_hgcer_npe_pions_cut.Fill(event.P_ngcer_npeSum, event.P_hgcer_npeSum)
    P_ngcer_vs_aero_npe_pions_cut.Fill(event.P_ngcer_npeSum, event.P_aero_npeSum)
    H_dp_vs_beta_pions_cut.Fill(event.H_gtr_dp, event.H_gtr_beta) 
    P_dp_vs_beta_pions_cut.Fill(event.P_gtr_dp, event.P_gtr_beta) 
    P_MMpi_vs_beta_pions_cut.Fill(event.MMpi, event.P_gtr_beta) 
    P_HGC_xy_npe_pions_cut.Fill(event.P_hgcer_yAtCer,event.P_hgcer_xAtCer,event.P_hgcer_npeSum)
    P_Aero_xy_npe_pions_cut.Fill(event.P_aero_yAtAero,event.P_aero_xAtAero,event.P_aero_npeSum)
    P_NGC_xy_npe_pions_cut.Fill(event.P_ngcer_yAtCer,event.P_ngcer_xAtCer,event.P_ngcer_npeSum)

for event in Cut_Pion_Events_Prompt_tree:
#    P_beta_pions_cut_prompt.Fill(event.P_gtr_beta)
#    P_RFTime_pions_cut_prompt.Fill(event.P_RF_Dist)
    ePiCoinTime_pions_cut_prompt.Fill(event.CTime_ePiCoinTime_ROC1)
    P_MMpi_pions_cut_prompt.Fill(event.MMpi)
#    MMpi_vs_ePiCoinTime_pions_cut_prompt.Fill(event.MMpi, event.CTime_ePiCoinTime_ROC1)

for event in Cut_Pion_Events_Random_tree:
#    P_beta_pions_cut_randm.Fill(event.P_gtr_beta)
#    P_RFTime_pions_cut_randm.Fill(event.P_RF_Dist)
    ePiCoinTime_pions_cut_randm.Fill(event.CTime_ePiCoinTime_ROC1)
    P_MMpi_pions_cut_randm.Fill(event.MMpi)
7
print("Histograms filled")

##############################################################################################################################################

# Random subtraction from missing mass and coin_Time
for event in Cut_Pion_Events_Random_tree:
    P_MMpi_pions_cut_randm_scaled.Fill(event.MMpi)
    P_MMpi_pions_cut_randm_scaled.Scale(1.0/nWindows)
P_MMpi_pions_cut_randm_sub.Add(P_MMpi_pions_cut_prompt, P_MMpi_pions_cut_randm_scaled, 1, -1)

############################################################################################################################################

# Saving histograms in PDF
c1_kin = TCanvas("c1_kin", "Kinematic Distributions", 100, 0, 1000, 900)
c1_kin.Divide(2,2)
c1_kin.cd(1)
Q2vsW_pions_cut.Draw("COLZ")
c1_kin.cd(2)
epsilon_pions_cut.Draw()
c1_kin.cd(3)
phiqvst_pions_cut.GetYaxis().SetRangeUser(minrangeuser,maxrangeuser)
phiqvst_pions_cut.Draw("SURF2 POL")
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
tradius = TGaxis(0,0,0.575,0,minrangeuser,maxrangeuser,10,"-+")
tradius.SetLineColor(2)
tradius.SetLabelColor(2)
tradius.Draw()
phizero = TLine(0,0,0.6,0) 
phizero.SetLineColor(kBlack)
phizero.SetLineWidth(2)
phizero.Draw()
# End of polar plotting section
c1_kin.cd(4)
P_MMpi_pions_cut_randm_sub.Draw("hist")
# Section for Neutron Peak Events Selection
shadedpeak = P_MMpi_pions_cut_randm_sub.Clone()
shadedpeak.SetFillColor(2)
shadedpeak.SetFillStyle(3244)
shadedpeak.GetXaxis().SetRangeUser(minbin, maxbin)
shadedpeak.Draw("samehist")
NeutronEvt = TPaveText(0.58934,0.675,0.95,0.75,"NDC")
BinLow = P_MMpi_pions_cut_randm_sub.GetXaxis().FindBin(minbin)
BinHigh = P_MMpi_pions_cut_randm_sub.GetXaxis().FindBin(maxbin)
BinIntegral = int(P_MMpi_pions_cut_randm_sub.Integral(BinLow, BinHigh))
NeutronEvt.SetLineColor(2)
NeutronEvt.AddText("e #pi n Events: %i" %(BinIntegral))
NeutronEvt.Draw()
# End of Neutron Peak Events Selection Section
c1_kin.Print(Pion_Analysis_Distributions + '(')

c1_acpt = TCanvas("c1_H_kin", "Electron-Pion Acceptance Distributions", 100, 0, 1000, 900)
c1_acpt.Divide(3,2)
c1_acpt.cd(1)
gPad.SetLogy()
H_xp_pions_uncut.SetLineColor(2)
H_xp_pions_uncut.Draw()
H_xp_pions_cut.SetLineColor(4)
H_xp_pions_cut.Draw("same")
c1_acpt.cd(2)
gPad.SetLogy()
H_yp_pions_uncut.SetLineColor(2)
H_yp_pions_uncut.Draw()
H_yp_pions_cut.SetLineColor(4)
H_yp_pions_cut.Draw("same")
c1_acpt.cd(3)
gPad.SetLogy() 
# 12/09/21 - NH - fixed the plotting range issue for HMS Delta
H_dp_pions_uncut.SetMinimum(0.1*H_dp_pions_cut.GetMinimum()+1) # min of plot should be one order of magnitude below the min bin in cut distribution
H_dp_pions_uncut.SetMaximum(10*H_dp_pions_uncut.GetBinContent(H_dp_pions_uncut.GetMaximumBin())) # Max of plot should be 1 order of magnitude greater than the max bin in uncut distribution
H_dp_pions_uncut.SetLineColor(2)
H_dp_pions_cut.SetLineColor(4)
H_dp_pions_uncut.Draw()
H_dp_pions_cut.Draw("same")
# TLegend (x1, y1, x2, y2) 
legend2 = ROOT.TLegend(0.115, 0.8, 0.6, 0.9)
legend2.AddEntry("H_dp_pions_uncut", "without cuts", "l")
legend2.AddEntry("H_dp_pions_cut", "with cuts (acpt/RF/PID)", "l")
legend2.Draw("same")
c1_acpt.cd(4)
gPad.SetLogy()
P_xp_pions_uncut.SetLineColor(2)
P_xp_pions_uncut.Draw()
P_xp_pions_cut.SetLineColor(4)
P_xp_pions_cut.Draw("same")
c1_acpt.cd(5)
gPad.SetLogy()
P_yp_pions_uncut.SetLineColor(2)
P_yp_pions_uncut.Draw()
P_yp_pions_cut.SetLineColor(4)
P_yp_pions_cut.Draw("same")
c1_acpt.cd(6)
gPad.SetLogy()
P_dp_pions_uncut.SetMinimum(0.1*P_dp_pions_cut.GetMinimum()+1) # SJDK 18/09/21 - Implemented same fixed as used above for HMS
P_dp_pions_uncut.SetMaximum(10*P_dp_pions_uncut.GetBinContent(P_dp_pions_uncut.GetMaximumBin()))
P_dp_pions_uncut.SetLineColor(2)
P_dp_pions_uncut.Draw()
P_dp_pions_cut.SetLineColor(4)
P_dp_pions_cut.Draw("same")
c1_acpt.Print(Pion_Analysis_Distributions)

c1_pid = TCanvas("c1_pid", "Electron-Pion CAL Distributions", 100, 0, 1000, 900)
c1_pid.Divide(2,2)
c1_pid.cd(1)
gPad.SetLogy()
H_cal_etottracknorm_pions_uncut.SetLineColor(2)
H_cal_etottracknorm_pions_uncut.Draw()
H_cal_etottracknorm_pions_cut.SetLineColor(4)
H_cal_etottracknorm_pions_cut.Draw("same")
legend7 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
legend7.AddEntry("H_cal_etottracknorm_pions_uncut", "without cuts", "l")
legend7.AddEntry("H_cal_etottracknorm_pions_cut", "with cuts (acpt/RF/PID)", "l")
legend7.Draw("same")
c1_pid.cd(2)
H_cal_etottracknorm_vs_cer_npe_pions_cut.Draw("COLZ")
c1_pid.cd(3)
gPad.SetLogy()
P_cal_etottracknorm_pions_uncut.SetLineColor(2)
P_cal_etottracknorm_pions_uncut.Draw()
P_cal_etottracknorm_pions_cut.SetLineColor(4)
P_cal_etottracknorm_pions_cut.Draw("same")
legend8 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
legend8.AddEntry("P_cal_etottracknorm_pions_uncut", "without cuts", "l")
legend8.AddEntry("P_cal_etottracknorm_pions_cut", "with cuts (acpt/RF/PID)", "l")
legend8.Draw("same")
c1_pid.cd(4)
#
c1_pid.Print(Pion_Analysis_Distributions)

c1_RF = TCanvas("c1_RF", "Electron-Pion RF Distributions", 100, 0, 1000, 900)
c1_RF.Divide(2,2)
c1_RF.cd(1)
P_RFTime_pions_uncut.SetLineColor(2)
P_RFTime_pions_uncut.Draw()
P_RFTime_pions_cut.SetLineColor(4)
P_RFTime_pions_cut.Draw("same")
legend9 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
legend9.AddEntry("P_RFTime_pions_uncut", "without cuts", "l")
legend9.AddEntry("P_RFTime_pions_cut", "with cuts (acpt/RF/PID)", "l")
legend9.Draw("same")
c1_RF.cd(2)
P_RFTime_pions_cut_noRF.SetLineColor(2)
P_RFTime_pions_cut_noRF.Draw()
P_RFTime_pions_cut.SetLineColor(4)
P_RFTime_pions_cut.Draw("same")
legend = ROOT.TLegend(0.115, 0.835, 0.415, 0.9)
legend.AddEntry("P_RFTime_pions_cut_noRF", "noRF_cuts (acpt/PID)", "l")
legend.AddEntry("P_RFTime_pions_cut", "RF_cuts (acpt/PID)", "l")
legend.Draw("same")
c1_RF.cd(3)
P_RFTime_vs_MMpi_pions_uncut.Draw("COLZ")
c1_RF.cd(4)
P_RFTime_vs_MMpi_pions_cut.Draw("COLZ")
c1_RF.Print(Pion_Analysis_Distributions)

c2_pid = TCanvas("c2_pid", "Electron-Pion Aero/HGC/NGC PID Distributions", 100, 0, 1000, 900)
c2_pid.Divide(2,2)
c2_pid.cd(1)
gPad.SetLogy()
P_hgcer_npe_pions_uncut.SetLineColor(2)
P_hgcer_npe_pions_uncut.Draw()
P_hgcer_npe_pions_cut.SetLineColor(4)
P_hgcer_npe_pions_cut.Draw("same")
legend10 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
legend10.AddEntry("P_hgcer_npe_pions_uncut", "without cuts", "l")
legend10.AddEntry("P_hgcer_npe_pions_cut", "with cuts (acpt/RF/PID)", "l")
legend10.Draw("same")
c2_pid.cd(2)
gPad.SetLogy()
P_aero_npe_pions_uncut.SetLineColor(2)
P_aero_npe_pions_uncut.Draw()
P_aero_npe_pions_cut.SetLineColor(4)
P_aero_npe_pions_cut.Draw("same")
legend11 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
legend11.AddEntry("P_aero_npe_pions_uncut", "without cuts", "l")
legend11.AddEntry("P_aero_npe_pions_cut", "with cuts (acpt/RF/PID)", "l")
legend11.Draw("same")
c2_pid.cd(3)
P_ngcer_npe_pions_uncut.SetLineColor(2)
P_ngcer_npe_pions_uncut.Draw()
P_ngcer_npe_pions_cut.SetLineColor(4)
P_ngcer_npe_pions_cut.Draw("same") 
#legend12 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
#legend12.AddEntry("P_ngcer_npe_pions_uncut", "without cuts", "l")
#legend12.AddEntry("P_ngcer_npe_pions_cut", "with cuts (acpt/RF/PID)", "l")
#legend12.Draw("same")
c2_pid.cd(4)
#
c2_pid.Print(Pion_Analysis_Distributions)

c3_pid = TCanvas("c3_pid", "Electron-Pion Aero/HGC/NGC PID Distributions", 100, 0, 1000, 900)
c3_pid.Divide(2,3)
c3_pid.cd(1)
gPad.SetLogz()
P_hgcer_vs_aero_npe_pions_uncut.Draw("COLZ")
c3_pid.cd(2)
P_hgcer_vs_aero_npe_pions_cut.Draw("COLZ")
c3_pid.cd(3)
P_ngcer_vs_hgcer_npe_pions_uncut.Draw("COLZ")
c3_pid.cd(4)
P_ngcer_vs_hgcer_npe_pions_cut.Draw("COLZ")
c3_pid.cd(5)
P_ngcer_vs_aero_npe_pions_uncut.Draw("COLZ")
c3_pid.cd(6)
P_ngcer_vs_aero_npe_pions_cut.Draw("COLZ")
c3_pid.Print(Pion_Analysis_Distributions)

c1_MM = TCanvas("c1_MM", "Electron-Pion CTime/Missing Mass Distributions", 100, 0, 1000, 900)
c1_MM.Divide(3,2)
c1_MM.cd(1)
ePiCoinTime_pions_uncut.SetLineColor(4)
ePiCoinTime_pions_uncut.Draw()
ePiCoinTime_pions_cut_prompt.SetLineColor(6)
ePiCoinTime_pions_cut_prompt.Draw("same")
ePiCoinTime_pions_cut_randm.SetLineColor(8)
ePiCoinTime_pions_cut_randm.Draw("same")
legend13 = ROOT.TLegend(0.1, 0.815, 0.48, 0.9)
legend13.AddEntry("ePiCoinTime_pions_uncut", "CT_without cuts", "l")
legend13.AddEntry("ePiCoinTime_pions_cut_prompt", "CT_prompt with cuts (acpt/RF/PID)", "l")
legend13.AddEntry("ePiCoinTime_pions_cut_randm", "CT_randoms with cuts (acpt/RF/PID)", "l")
legend13.Draw("same")
c1_MM.cd(2)
ePiCoinTime_pions_cut.SetLineColor(4)
ePiCoinTime_pions_cut.Draw()
ePiCoinTime_pions_cut_prompt.SetLineColor(6)
ePiCoinTime_pions_cut_prompt.Draw("same")
ePiCoinTime_pions_cut_randm.SetLineColor(8)
ePiCoinTime_pions_cut_randm.Draw("same")
legend14 = ROOT.TLegend(0.1, 0.815, 0.48, 0.9)
legend14.AddEntry("ePiCoinTime_pions_cut", "CT_with cuts (acpt/RF/PID)", "l")
legend14.AddEntry("ePiCoinTime_pions_cut_prompt", "CT_prompt with cuts (acpt/RF/PID)", "l")
legend14.AddEntry("ePiCoinTime_pions_cut_randm", "CT_randoms with cuts (acpt/RF/PID)", "l")
legend14.Draw("same")
c1_MM.cd(3)
P_MMpi_pions_uncut.Draw()
# Section for Neutron Peak Events Selection
#shadedpeak = P_MMpi_pions_uncut.Clone()
#shadedpeak.SetFillColor(2)
#shadedpeak.SetFillStyle(3244)
#shadedpeak.GetXaxis().SetRangeUser(minbin, maxbin)
#shadedpeak.Draw("samehist")
#NeutronEvt = TPaveText(0.58934,0.675,0.95,0.75,"NDC")
#BinLow = P_MMpi_pions_uncut.GetXaxis().FindBin(minbin)
#BinHigh = P_MMpi_pions_uncut.GetXaxis().FindBin(maxbin)
#BinIntegral = int(P_MMpi_pions_uncut.Integral(BinLow, BinHigh))
#NeutronEvt.SetLineColor(2)
#NeutronEvt.AddText("e #pi n Events: %i" %(BinIntegral))
#NeutronEvt.Draw()
c1_MM.cd(4)
P_MMpi_pions_cut.SetLineColor(4)
P_MMpi_pions_cut.Draw()
P_MMpi_pions_cut_prompt.SetLineColor(6)
P_MMpi_pions_cut_prompt.Draw("same")
P_MMpi_pions_cut_randm.SetLineColor(8)
P_MMpi_pions_cut_randm.Draw("same")
legend15 = ROOT.TLegend(0.4, 0.815, 0.78, 0.9)
legend15.AddEntry("P_MMpi_pions_cut", "MM with cuts (acpt/RF/PID)", "l")
legend15.AddEntry("P_MMpi_pions_cut_prompt", "MM_prompt with cuts (acpt/RF/PID)", "l")
legend15.AddEntry("P_MMpi_pions_cut_randm", "MM_randoms with cuts (acpt/RF/PID)", "l")
legend15.Draw("same")
c1_MM.cd(5)
P_MMpi_vs_beta_pions_uncut.Draw("COLZ")
c1_MM.cd(6)
P_MMpi_vs_beta_pions_cut.Draw("COLZ")
c1_MM.Print(Pion_Analysis_Distributions)

c1_CT = TCanvas("c1_CT", "Electron-Pion CTime Distributions", 100, 0, 1000, 900)
c1_CT.Divide(2,2)
c1_CT.cd(1)
ePiCoinTime_vs_beta_pions_uncut.Draw("COLZ")
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
ePiCoinTime_vs_beta_pions_cut.Draw("COLZ")
c1_CT.cd(3)
ePiCoinTime_vs_MMpi_pions_uncut.Draw("COLZ")
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
ePiCoinTime_vs_MMpi_pions_cut.Draw("COLZ")
c1_CT.Print(Pion_Analysis_Distributions)

c1_delta = TCanvas("c1_delta", "Delta Debugging", 100, 0, 1000, 900)
c1_delta.Divide(2,2)
c1_delta.cd(1)
H_dp_vs_beta_pions_uncut.Draw("COLZ")
c1_delta.cd(2)
H_dp_vs_beta_pions_cut.Draw("COLZ")
c1_delta.cd(3)
P_dp_vs_beta_pions_uncut.Draw("COLZ")
c1_delta.cd(4)
P_dp_vs_beta_pions_cut.Draw("COLZ")
c1_delta.Print(Pion_Analysis_Distributions)

c1_proj = TCanvas("c1_proj", "HGC/NGC/Aero XY Projection", 100, 0, 1000, 900)
c1_proj.Divide(2,3)
c1_proj.cd(1)
HGC_proj_yx_uncut = ROOT.TProfile2D(P_HGC_xy_npe_pions_uncut.Project3DProfile("yx"))
HGC_proj_yx_uncut.Draw("COLZ")
c1_proj.cd(2)
HGC_proj_yx_cut = ROOT.TProfile2D(P_HGC_xy_npe_pions_cut.Project3DProfile("yx"))
HGC_proj_yx_cut.Draw("COLZ")
c1_proj.cd(3)
NGC_proj_yx_uncut = ROOT.TProfile2D(P_NGC_xy_npe_pions_uncut.Project3DProfile("yx"))
NGC_proj_yx_uncut.Draw("COLZ")
c1_proj.cd(4)
NGC_proj_yx_cut = ROOT.TProfile2D(P_NGC_xy_npe_pions_cut.Project3DProfile("yx"))
NGC_proj_yx_cut.Draw("COLZ")
c1_proj.cd(5)
Aero_proj_yx_uncut = ROOT.TProfile2D(P_Aero_xy_npe_pions_uncut.Project3DProfile("yx"))
Aero_proj_yx_uncut.Draw("COLZ")
c1_proj.cd(6)
Aero_proj_yx_cut = ROOT.TProfile2D(P_Aero_xy_npe_pions_cut.Project3DProfile("yx"))
Aero_proj_yx_cut.Draw("COLZ")
c1_proj.Print(Pion_Analysis_Distributions + ')')

#############################################################################################################################################

# Making directories in output file
outHistFile = ROOT.TFile.Open("%s/%s_%s_Output_Data.root" % (OUTPATH, runNum, MaxEvent), "RECREATE")                                                                                                    
d_Uncut_Pion_Events = outHistFile.mkdir("Uncut_Pion_Events")
d_Cut_Pion_Events_All = outHistFile.mkdir("Cut_Pion_Events_All")
d_Cut_Pion_Events_Prompt = outHistFile.mkdir("Cut_Pion_Events_Prompt")
d_Cut_Pion_Events_Random = outHistFile.mkdir("Cut_Pion_Events_Random")

# Writing Histograms in output root file
d_Uncut_Pion_Events.cd()
H_beta_pions_uncut.Write()
H_xp_pions_uncut.Write()
H_yp_pions_uncut.Write()
H_dp_pions_uncut.Write()
H_hod_goodscinhit_pions_uncut.Write()
H_hod_goodstarttime_pions_uncut.Write()
H_cal_etotnorm_pions_uncut.Write()
H_cal_etottracknorm_pions_uncut.Write()
H_cer_npe_pions_uncut.Write()
H_RFTime_pions_uncut.Write()
P_beta_pions_uncut.Write()
P_xp_pions_uncut.Write()
P_yp_pions_uncut.Write()
P_dp_pions_uncut.Write()
P_p_pions_uncut.Write()
P_hod_goodscinhit_pions_uncut.Write()
P_hod_goodstarttime_pions_uncut.Write()
P_cal_etotnorm_pions_uncut.Write()
P_cal_etottracknorm_pions_uncut.Write()
P_hgcer_npe_pions_uncut.Write()
P_xhgcer_pions_uncut.Write()
P_yhgcer_pions_uncut.Write()
P_aero_npe_pions_uncut.Write()
P_xaero_pions_uncut.Write()
P_yaero_pions_uncut.Write()
P_ngcer_npe_pions_uncut.Write()
P_xngcer_pions_uncut.Write()
P_yngcer_pions_uncut.Write() 
P_MMpi_pions_uncut.Write()
P_RFTime_pions_uncut.Write()
ePiCoinTime_pions_uncut.Write()
Q2_pions_uncut.Write()
W_pions_uncut.Write()
epsilon_pions_uncut.Write()
phiq_pions_uncut.Write()
t_pions_uncut.Write()
H_cal_etottracknorm_vs_cer_npe_pions_uncut.Write()
P_hgcer_vs_aero_npe_pions_uncut.Write()
ePiCoinTime_vs_MMpi_pions_uncut.Write()
P_hgcer_yx_pions_uncut.Write()
P_aero_yx_pions_uncut.Write()
ePiCoinTime_vs_beta_pions_uncut.Write()
P_RFTime_vs_MMpi_pions_uncut.Write()
P_cal_etottracknorm_vs_ngcer_npe_pions_uncut.Write()
P_ngcer_yx_pions_uncut.Write()
P_ngcer_vs_hgcer_npe_pions_uncut.Write()
P_ngcer_vs_aero_npe_pions_uncut.Write()
H_dp_vs_beta_pions_uncut.Write()
P_dp_vs_beta_pions_uncut.Write()
P_MMpi_vs_beta_pions_uncut.Write()
HGC_proj_yx_uncut.Write()
NGC_proj_yx_uncut.Write()
Aero_proj_yx_uncut.Write()
P_HGC_xy_npe_pions_uncut.Write()
P_Aero_xy_npe_pions_uncut.Write()
P_NGC_xy_npe_pions_uncut.Write()

d_Cut_Pion_Events_All.cd()
H_beta_pions_cut.Write()
H_xp_pions_cut.Write()
H_yp_pions_cut.Write()
H_dp_pions_cut.Write()
H_hod_goodscinhit_pions_cut.Write()
H_hod_goodstarttime_pions_cut.Write()
H_cal_etotnorm_pions_cut.Write()
H_cal_etottracknorm_pions_cut.Write()
H_cer_npe_pions_cut.Write()
H_RFTime_pions_cut.Write()
P_beta_pions_cut.Write()
P_xp_pions_cut.Write()
P_yp_pions_cut.Write()
P_dp_pions_cut.Write()
P_p_pions_cut.Write()
P_hod_goodscinhit_pions_cut.Write()
P_hod_goodstarttime_pions_cut.Write()
P_cal_etotnorm_pions_cut.Write()
P_cal_etottracknorm_pions_cut.Write()
P_hgcer_npe_pions_cut.Write()
P_xhgcer_pions_cut.Write()
P_yhgcer_pions_cut.Write()
P_aero_npe_pions_cut.Write()
P_xaero_pions_cut.Write()
P_yaero_pions_cut.Write()
P_ngcer_npe_pions_cut.Write()
P_xngcer_pions_cut.Write()
P_yngcer_pions_cut.Write()
P_MMpi_pions_cut.Write()
P_RFTime_pions_cut.Write()
ePiCoinTime_pions_cut.Write()
Q2_pions_cut.Write()
W_pions_cut.Write()
epsilon_pions_cut.Write()
phiq_pions_cut.Write()
t_pions_cut.Write()
Q2vsW_pions_cut.Write()
phiqvst_pions_cut.Write()
H_cal_etottracknorm_vs_cer_npe_pions_cut.Write()
P_hgcer_vs_aero_npe_pions_cut.Write()
ePiCoinTime_vs_MMpi_pions_cut.Write()
P_hgcer_yx_pions_cut.Write()
P_aero_yx_pions_cut.Write()
ePiCoinTime_vs_beta_pions_cut.Write()
P_RFTime_vs_MMpi_pions_cut.Write()
P_MMpi_pions_cut_randm_sub.Write()
P_cal_etottracknorm_vs_ngcer_npe_pions_cut.Write()
P_ngcer_yx_pions_cut.Write()
P_ngcer_vs_hgcer_npe_pions_cut.Write()
P_ngcer_vs_aero_npe_pions_cut.Write()
H_dp_vs_beta_pions_cut.Write()
P_dp_vs_beta_pions_cut.Write()
P_MMpi_vs_beta_pions_cut.Write()
P_HGC_xy_npe_pions_cut.Write()
P_Aero_xy_npe_pions_cut.Write()
P_NGC_xy_npe_pions_cut.Write()
HGC_proj_yx_cut.Write()
NGC_proj_yx_cut.Write()
Aero_proj_yx_cut.Write()

d_Cut_Pion_Events_Prompt.cd()
P_beta_pions_cut_prompt.Write()
P_RFTime_pions_cut_prompt.Write()
ePiCoinTime_pions_cut_prompt.Write()
P_MMpi_pions_cut_prompt.Write()
MMpi_vs_ePiCoinTime_pions_cut_prompt.Write()

d_Cut_Pion_Events_Random.cd()
P_beta_pions_cut_randm.Write()
P_RFTime_pions_cut_randm.Write()
ePiCoinTime_pions_cut_randm.Write()
P_MMpi_pions_cut_randm.Write() 

outHistFile.Close()
infile.Close() 
print ("Processing Complete")
print("!!!!!!!!\n %i pi-n events \n!!!!!!!!" % BinIntegral)
