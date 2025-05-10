#! /usr/bin/python
#
# Description:
# ================================================================
# Time-stamp: "2024-03-15 01:29:19 junaid"
# ================================================================
#
# Author:  Muhammad Junaid III <mjo147@uregina.ca>
#
# Copyright (c) junaid
#
###################################################################################################################################################

# Import relevant packages
import uproot
import uproot as up
import numpy as np

np.bool = bool
np.float = float

import root_numpy as rnp
import pandas as pd
import root_pandas as rpd
import ROOT
import scipy
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import sys, math, os, subprocess
import array
import csv
from ROOT import TCanvas, TList, TPaveLabel, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TLegend, TGaxis, TLine, TMath, TLatex, TPaveText, TArc, TGraphPolar, TText, TString
from ROOT import kBlack, kCyan, kRed, kGreen, kMagenta, kBlue
from functools import reduce
import math as ma
import re # Regexp package - for string manipulation

##################################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=7:
    print("!!!!! ERROR !!!!!\n Expected 8 arguments\n Usage is with - ROOTfileSuffixs Beam Energy MaxEvents RunList CVSFile\n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################

# Defining some constants here
minbin = 0.90 # minimum bin for selecting neutrons events in missing mass distribution
maxbin = 0.98 # maximum bin for selecting neutrons events in missing mass distribution

##################################################################################################################################################

# Input params - run number and max number of events
BEAM_ENERGY = sys.argv[1]
Q2 = sys.argv[2]
W = sys.argv[3]
ptheta = sys.argv[4]
DATA_Suffix = sys.argv[5]
MaxEvent = sys.argv[6]
runNum = sys.argv[7]
################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
from ltsep import Root

lt=Root(os.path.realpath(__file__), "Plot_Prod")

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH

#################################################################################################################################################

# Output PDF File Name
print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))
Pion_Analysis_Distributions = "%s/%s_%s_%s_%s_%s_%s_Pion_PID_Analysis_Distributions.pdf" % (OUTPATH, BEAM_ENERGY, Q2, W, ptheta, DATA_Suffix, MaxEvent)

# Input file location and variables taking
rootFile_DATA = "%s/%s_%s_%s_%s_%s_%s.root" % (OUTPATH, BEAM_ENERGY, Q2, W, ptheta, DATA_Suffix, MaxEvent)

###############################################################################################################################################

# Section for grabing Prompt/Random selection parameters from PARAM file
PARAMPATH = "%s/DB/PARAM" % UTILPATH
print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))
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

##############################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Read stuff from the main event tree

infile_DATA = ROOT.TFile.Open(rootFile_DATA, "READ")

Uncut_Pion_Events_Data_tree = infile_DATA.Get("Uncut_Pion_Events")
Cut_Pion_Events_Accpt_Data_tree = infile_DATA.Get("Cut_Pion_Events_Accpt")
Cut_Pion_Events_Prompt_Data_tree = infile_DATA.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Random_Data_tree = infile_DATA.Get("Cut_Pion_Events_Random")
nEntries_TBRANCH_DATA  = Cut_Pion_Events_Accpt_Data_tree.GetEntries()

###################################################################################################################################################
nbins = 200

# Defining Histograms for Pions
# Uncut Data Histograms
H_gtr_beta_pions_data_uncut = ROOT.TH1D("H_gtr_beta_pions_data_uncut", "HMS #beta; HMS_gtr_#beta; Counts", nbins, 0.8, 1.2)
H_cal_etottracknorm_pions_data_uncut = ROOT.TH1D("H_cal_etottracknorm_pions_data_uncut", "HMS cal etottracknorm (uncut); HMS_cal_etottracknorm; Counts", nbins, 0.0, 1.8)
H_cer_npeSum_pions_data_uncut = ROOT.TH1D("H_cer_npeSum_pions_data_uncut", "HMS cer npeSum (uncut); HMS_cer_npeSum; Counts", nbins, 0, 50)
H_RFTime_Dist_pions_data_uncut = ROOT.TH1D("H_RFTime_Dist_pions_data_uncut", "HMS RFTime (uncut); HMS_RFTime; Counts", nbins, 0, 4)
P_gtr_beta_pions_data_uncut = ROOT.TH1D("P_gtr_beta_pions_data_uncut", "SHMS #beta (uncut); SHMS_gtr_#beta; Counts", nbins, 0.0, 1.4)
P_gtr_dp_pions_data_uncut = ROOT.TH1D("P_gtr_dp_protons_pions_data_uncut", "SHMS #delta (uncut); SHMS_gtr_dp; Counts", nbins, -30, 30)
P_cal_etottracknorm_pions_data_uncut = ROOT.TH1D("P_cal_etottracknorm_pions_data_uncut", "SHMS cal etottracknorm (uncut); SHMS_cal_etottracknorm; Counts", nbins, 0.0, 1.8)
P_hgcer_npeSum_pions_data_uncut = ROOT.TH1D("P_hgcer_npeSum_pions_data_uncut", "SHMS HGC npeSum (uncut); SHMS_hgcer_npeSum; Counts", nbins, 0, 50)
P_hgcer_xAtCer_pions_data_uncut = ROOT.TH1D("P_hgcer_xAtCer_pions_data_uncut", "SHMS HGC xAtCer (uncut); SHMS_hgcer_xAtCer; Counts", nbins, -60, 60)
P_hgcer_yAtCer_pions_data_uncut = ROOT.TH1D("P_hgcer_yAtCer_pions_data_uncut", "SHMS HGC yAtCer (uncut); SHMS_hgcer_yAtCer; Counts", nbins, -60, 60)
P_ngcer_npeSum_pions_data_uncut = ROOT.TH1D("P_ngcer_npeSum_pions_data_uncut", "SHMS NGC npeSum (uncut); SHMS_ngcer_npeSum; Counts", nbins, 0, 50)
P_ngcer_xAtCer_pions_data_uncut = ROOT.TH1D("P_ngcer_xAtCer_pions_data_uncut", "SHMS NGC xAtCer (uncut); SHMS_ngcer_xAtCer; Counts", nbins, -60, 60)
P_ngcer_yAtCer_pions_data_uncut = ROOT.TH1D("P_ngcer_yAtCer_pions_data_uncut", "SHMS NGC yAtCer (uncut); SHMS_ngcer_yAtCer; Counts", nbins, -60, 60)
P_aero_npeSum_pions_data_uncut = ROOT.TH1D("P_aero_npeSum_pions_data_uncut", "SHMS aero npeSum (uncut); SHMS_aero_npeSum; Counts", nbins, 0, 50)
P_aero_xAtAero_pions_data_uncut = ROOT.TH1D("P_acero_xAtAero_pions_data_uncut", "SHMS aero xAtAero (uncut); SHMS_aero_xAtAero; Counts", nbins, -60, 60)
P_aero_yAtAero_pions_data_uncut = ROOT.TH1D("P_aero_yAtAero_pions_data_uncut", "SHMS aero yAtAero (uncut); SHMS_aero_yAtAero; Counts", nbins, -60, 60)
P_kin_MMpi_pions_data_uncut = ROOT.TH1D("P_kin_MMpi_pions_data_uncut", "MIssing Mass data (uncut); MM_{\pi}; Counts", nbins, 0, 2.0)
P_RFTime_Dist_pions_data_uncut = ROOT.TH1D("P_RFTime_Dist_pions_data_uncut", "SHMS RFTime (uncut); SHMS_RFTime; Counts", nbins, 0, 4)
CTime_ePiCoinTime_ROC1_pions_data_uncut = ROOT.TH1D("CTime_ePiCoinTime_ROC1_pions_data_uncut", "Electron-Pion CTime (uncut); e pi Coin_Time; Counts", nbins, -50, 50)

# Acceptance Cut Data Histograms
H_gtr_beta_pions_data_accpt_cut_all = ROOT.TH1D("H_gtr_beta_pions_data_accpt_cut_all", "HMS #beta (accpt_cut); HMS_gtr_#beta; Counts", nbins, 0.0, 1.2)
H_cal_etottracknorm_pions_data_accpt_cut_all = ROOT.TH1D("H_cal_etottracknorm_pions_data_accpt_cut_all", "HMS cal etottracknorm (accpt_cut); HMS_cal_etottracknorm; Counts", nbins, 0.0, 1.8)
H_cer_npeSum_pions_data_accpt_cut_all = ROOT.TH1D("H_cer_npeSum_pions_data_accpt_cut_all", "HMS cer npeSum (accpt_cut); HMS_cer_npeSum; Counts", nbins, 0, 50)
H_RFTime_Dist_pions_data_accpt_cut_all = ROOT.TH1D("H_RFTime_Dist_pions_data_accpt_cut_all", "HMS RFTime (accpt_cut); HMS_RFTime; Counts", nbins, 0, 4)
P_gtr_beta_pions_data_accpt_cut_all = ROOT.TH1D("P_gtr_beta_pions_data_accpt_cut_all", "SHMS #beta (accpt_cut); SHMS_gtr_#beta; Counts", nbins, 0.0, 1.4)
P_gtr_dp_pions_data_accpt_cut_all = ROOT.TH1D("P_gtr_dp_protons_pions_data_accpt_cut_all", "SHMS #delta (accpt_cut); SHMS_gtr_dp; Counts", nbins, -30, 30)
P_cal_etottracknorm_pions_data_accpt_cut_all = ROOT.TH1D("P_cal_etottracknorm_pions_data_accpt_cut_all", "SHMS cal etottracknorm (accpt_cut_all); SHMS_cal_etottracknorm; Counts", nbins, 0, 1.8)
P_hgcer_npeSum_pions_data_accpt_cut_all = ROOT.TH1D("P_hgcer_npeSum_pions_data_accpt_cut_all", "SHMS HGC npeSum (accpt_cut); SHMS_hgcer_npeSum; Counts", nbins, 0, 50)
P_hgcer_xAtCer_pions_data_accpt_cut_all = ROOT.TH1D("P_hgcer_xAtCer_pions_data_accpt_cut_all", "SHMS HGC xAtCer (accpt_cut); SHMS_hgcer_xAtCer; Counts", nbins, -60, 60)
P_hgcer_yAtCer_pions_data_accpt_cut_all = ROOT.TH1D("P_hgcer_yAtCer_pions_data_accpt_cut_all", "SHMS HGC yAtCer (accpt_cut); SHMS_hgcer_yAtCer; Counts", nbins, -60, 60)
P_ngcer_npeSum_pions_data_accpt_cut_all = ROOT.TH1D("P_ngcer_npeSum_pions_data_accpt_cut_all", "SHMS NGC npeSum (accpt_cut); SHMS_ngcer_npeSum; Counts", nbins, 0, 50)
P_ngcer_xAtCer_pions_data_accpt_cut_all = ROOT.TH1D("P_ngcer_xAtCer_pions_data_accpt_cut_all", "SHMS NGC xAtCer (accpt_cut); SHMS_ngcer_xAtCer; Counts", nbins, -60, 60)
P_ngcer_yAtCer_pions_data_accpt_cut_all = ROOT.TH1D("P_ngcer_yAtCer_pions_data_accpt_cut_all", "SHMS NGC yAtCer (accpt_cut); SHMS_ngcer_yAtCer; Counts", nbins, -60, 60)
P_aero_npeSum_pions_data_accpt_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_accpt_cut_all", "SHMS aero npeSum (accpt_cut); SHMS_aero_npeSum; Counts", nbins, 0, 50)
P_aero_xAtAero_pions_data_accpt_cut_all = ROOT.TH1D("P_acero_xAtAero_pions_data_accpt_cut_all", "SHMS aero xAtAero (accpt_cut); SHMS_aero_xAtAero; Counts", nbins, -60, 60)
P_aero_yAtAero_pions_data_accpt_cut_all = ROOT.TH1D("P_aero_yAtAero_pions_data_accpt_cut_all", "SHMS aero yAtAero (accpt_cut); SHMS_aero_yAtAero; Counts", nbins, -60, 60)
P_kin_MMpi_pions_data_accpt_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_accpt_cut_all", "MIssing Mass data (accpt_cut); MM_{\pi}; Counts", nbins, 0, 2.0)
P_RFTime_Dist_pions_data_accpt_cut_all = ROOT.TH1D("P_RFTime_Dist_pions_data_accpt_cut_all", "SHMS RFTime (accpt_cut); SHMS_RFTime; Counts", nbins, 0, 4)
CTime_ePiCoinTime_ROC1_pions_data_accpt_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC1_pions_data_accpt_cut_all", "Electron-Pion CTime (accpt_cut); e pi Coin_Time; Counts", nbins, -50, 50)

# Prompt Cut Data Histograms
H_gtr_beta_pions_data_prompt_cut_all = ROOT.TH1D("H_gtr_beta_pions_data_prompt_cut_all", "HMS #beta (prompt+accpt_cut); HMS_gtr_#beta; Counts", nbins, 0.0, 1.2)
H_cal_etottracknorm_pions_data_prompt_cut_all = ROOT.TH1D("H_cal_etottracknorm_pions_data_prompt_cut_all", "HMS cal etottracknorm (prompt+accpt_cut); HMS_cal_etottracknorm; Counts", nbins, 0.0, 1.8)
H_cer_npeSum_pions_data_prompt_cut_all = ROOT.TH1D("H_cer_npeSum_pions_data_prompt_cut_all", "HMS cer npeSum (prompt+accpt_cut); HMS_cer_npeSum; Counts", nbins, 0, 50)
H_RFTime_Dist_pions_data_prompt_cut_all = ROOT.TH1D("H_RFTime_Dist_pions_data_prompt_cut_all", "HMS RFTime (prompt+accpt_cut); HMS_RFTime; Counts", nbins, 0, 4)
P_gtr_beta_pions_data_prompt_cut_all = ROOT.TH1D("P_gtr_beta_pions_data_prompt_cut_all", "SHMS #beta (prompt+accpt_cut); SHMS_gtr_#beta; Counts", nbins, 0.0, 1.4)
P_gtr_dp_pions_data_prompt_cut_all = ROOT.TH1D("P_gtr_dp_protons_pions_data_prompt_cut_all", "SHMS #delta (prompt+accpt_cut); SHMS_gtr_dp; Counts", nbins, -30, 30)
P_cal_etottracknorm_pions_data_prompt_cut_all = ROOT.TH1D("P_cal_etottracknorm_pions_data_prompt_cut_all", "SHMS cal etottracknorm (prompt+accpt_cut); SHMS_cal_etottracknorm; Counts", nbins, 0, 1.8)
P_hgcer_npeSum_pions_data_prompt_cut_all = ROOT.TH1D("P_hgcer_npeSum_pions_data_prompt_cut_all", "SHMS HGC npeSum (prompt+accpt_cut); SHMS_hgcer_npeSum; Counts", nbins, 0, 50)
P_hgcer_xAtCer_pions_data_prompt_cut_all = ROOT.TH1D("P_hgcer_xAtCer_pions_data_prompt_cut_all", "SHMS HGC xAtCer (prompt+accpt_cut); SHMS_hgcer_xAtCer; Counts", nbins, -60, 60)
P_hgcer_yAtCer_pions_data_prompt_cut_all = ROOT.TH1D("P_hgcer_yAtCer_pions_data_prompt_cut_all", "SHMS HGC yAtCer (prompt+accpt_cut); SHMS_hgcer_yAtCer; Counts", nbins, -60, 60)
P_ngcer_npeSum_pions_data_prompt_cut_all = ROOT.TH1D("P_ngcer_npeSum_pions_data_prompt_cut_all", "SHMS NGC npeSum (prompt+accpt_cut); SHMS_ngcer_npeSum; Counts", nbins, 0, 50)
P_ngcer_xAtCer_pions_data_prompt_cut_all = ROOT.TH1D("P_ngcer_xAtCer_pions_data_prompt_cut_all", "SHMS NGC xAtCer (prompt+accpt_cut); SHMS_ngcer_xAtCer; Counts", nbins, -60, 60)
P_ngcer_yAtCer_pions_data_prompt_cut_all = ROOT.TH1D("P_ngcer_yAtCer_pions_data_prompt_cut_all", "SHMS NGC yAtCer (prompt+accpt_cut); SHMS_ngcer_yAtCer; Counts", nbins, -60, 60)
P_aero_npeSum_pions_data_prompt_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_prompt_cut_all", "SHMS aero npeSum (prompt+accpt_cut); SHMS_aero_npeSum; Counts", nbins, 0, 50)
P_aero_xAtAero_pions_data_prompt_cut_all = ROOT.TH1D("P_acero_xAtAero_pions_data_prompt_cut_all", "SHMS aero xAtAero (prompt+accpt_cut); SHMS_aero_xAtAero; Counts", nbins, -60, 60)
P_aero_yAtAero_pions_data_prompt_cut_all = ROOT.TH1D("P_aero_yAtAero_pions_data_prompt_cut_all", "SHMS aero yAtAero (prompt+accpt_cut); SHMS_aero_yAtAero; Counts", nbins, -60, 60)
P_kin_MMpi_pions_data_prompt_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_prompt_cut_all", "MIssing Mass data (prompt+accpt_cut); MM_{\pi}; Counts", nbins, 0.0, 2.0)
P_RFTime_Dist_pions_data_prompt_cut_all = ROOT.TH1D("P_RFTime_Dist_pions_data_prompt_cut_all", "SHMS RFTime (prompt+accpt_cut); SHMS_RFTime; Counts", nbins, 0, 4)
CTime_ePiCoinTime_ROC1_pions_data_prompt_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC1_pions_data_prompt_cut_all", "Electron-Pion CTime (prompt+accpt_cut); e pi Coin_Time; Counts", nbins, -50, 50)

# Random Cut Data Histograms
H_gtr_beta_pions_data_random_cut_all = ROOT.TH1D("H_gtr_beta_pions_data_random_cut_all", "HMS #beta (random+accpt_cut); HMS_gtr_#beta; Counts", nbins, 0.0, 1.2)
H_cal_etottracknorm_pions_data_random_cut_all = ROOT.TH1D("H_cal_etottracknorm_pions_data_random_cut_all", "HMS cal etottracknorm (random+accpt_cut); HMS_cal_etottracknorm; Counts", nbins, 0.0, 1.8)
H_cer_npeSum_pions_data_random_cut_all = ROOT.TH1D("H_cer_npeSum_pions_data_random_cut_all", "HMS cer npeSum (random+accpt_cut); HMS_cer_npeSum; Counts", nbins, 0, 50)
H_RFTime_Dist_pions_data_random_cut_all = ROOT.TH1D("H_RFTime_Dist_pions_data_random_cut_all", "HMS RFTime (random+accpt_cut); HMS_RFTime; Counts", nbins, 0, 4)
P_gtr_beta_pions_data_random_cut_all = ROOT.TH1D("P_gtr_beta_pions_data_random_cut_all", "SHMS #beta (random+accpt_cut); SHMS_gtr_#beta; Counts", nbins, 0.0, 1.4)
P_gtr_dp_pions_data_random_cut_all = ROOT.TH1D("P_gtr_dp_protons_pions_data_random_cut_all", "SHMS #delta (random+accpt_cut); SHMS_gtr_dp; Counts", nbins, -30, 30)
P_cal_etottracknorm_pions_data_random_cut_all = ROOT.TH1D("P_cal_etottracknorm_pions_data_random_cut_all", "SHMS cal etottracknorm (random+accpt_cut); SHMS_cal_etottracknorm; Counts", nbins, 0, 1.8)
P_hgcer_npeSum_pions_data_random_cut_all = ROOT.TH1D("P_hgcer_npeSum_pions_data_random_cut_all", "SHMS HGC npeSum (random+accpt_cut); SHMS_hgcer_npeSum; Counts", nbins, 0, 50)
P_hgcer_xAtCer_pions_data_random_cut_all = ROOT.TH1D("P_hgcer_xAtCer_pions_data_random_cut_all", "SHMS HGC xAtCer (random+accpt_cut); SHMS_hgcer_xAtCer; Counts", nbins, -60, 60)
P_hgcer_yAtCer_pions_data_random_cut_all = ROOT.TH1D("P_hgcer_yAtCer_pions_data_random_cut_all", "SHMS HGC yAtCer (random+accpt_cut); SHMS_hgcer_yAtCer; Counts", nbins, -60, 60)
P_ngcer_npeSum_pions_data_random_cut_all = ROOT.TH1D("P_ngcer_npeSum_pions_data_random_cut_all", "SHMS NGC npeSum (random+accpt_cut); SHMS_ngcer_npeSum; Counts", nbins, 0, 50)
P_ngcer_xAtCer_pions_data_random_cut_all = ROOT.TH1D("P_ngcer_xAtCer_pions_data_random_cut_all", "SHMS NGC xAtCer (random+accpt_cut); SHMS_ngcer_xAtCer; Counts", nbins, -60, 60)
P_ngcer_yAtCer_pions_data_random_cut_all = ROOT.TH1D("P_ngcer_yAtCer_pions_data_random_cut_all", "SHMS NGC yAtCer (random+accpt_cut); SHMS_ngcer_yAtCer; Counts", nbins, -60, 60)
P_aero_npeSum_pions_data_random_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_random_cut_all", "SHMS aero npeSum (random+accpt_cut); SHMS_aero_npeSum; Counts", nbins, 0, 50)
P_aero_xAtAero_pions_data_random_cut_all = ROOT.TH1D("P_acero_xAtAero_pions_data_random_cut_all", "SHMS aero xAtAero (random+accpt_cut); SHMS_aero_xAtAero; Counts", nbins, -60, 60)
P_aero_yAtAero_pions_data_random_cut_all = ROOT.TH1D("P_aero_yAtAero_pions_data_random_cut_all", "SHMS aero yAtAero (random+accpt_cut); SHMS_aero_yAtAero; Counts", nbins, -60, 60)
P_kin_MMpi_pions_data_random_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_random_cut_all", "MIssing Mass data (random+accpt_cut); MM_{\pi}; Counts", nbins, 0, 2.0)
P_RFTime_Dist_pions_data_random_cut_all = ROOT.TH1D("P_RFTime_Dist_pions_data_random_cut_all", "SHMS RFTime (random+accpt_cut); SHMS_RFTime; Counts", nbins, 0, 4)
CTime_ePiCoinTime_ROC1_pions_data_random_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC1_pions_data_random_cut_all", "Electron-Pion CTime (random+accpt_cut); e pi Coin_Time; Counts", nbins, -50, 50)
CTime_ePiCoinTime_ROC1_pions_data_random_unsub_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC1_pions_data_random_unsub_cut_all", "Electron-Pion CTime (random+accpt_cut); e pi Coin_Time; Counts", nbins, -50, 50)

# Cut All Data Histograms
H_gtr_beta_pions_data_cut_all = ROOT.TH1D("H_gtr_beta_pions_data_cut_all", "HMS #beta (all_cut); HMS_gtr_#beta; Counts", nbins, 0.0, 1.2)
H_cal_etottracknorm_pions_data_cut_all = ROOT.TH1D("H_cal_etottracknorm_pions_data_cut_all", "HMS cal etottracknorm (all_cut); HMS_cal_etottracknorm; Counts", nbins, 0.0, 1.8)
H_cer_npeSum_pions_data_cut_all = ROOT.TH1D("H_cer_npeSum_pions_data_cut_all", "HMS cer npeSum (all_cut); HMS_cer_npeSum; Counts", nbins, 0, 50)
H_RFTime_Dist_pions_data_cut_all = ROOT.TH1D("H_RFTime_Dist_pions_data_cut_all", "HMS RFTime (all_cut); HMS_RFTime; Counts", nbins, 0, 4)
P_gtr_beta_pions_data_cut_all = ROOT.TH1D("P_gtr_beta_pions_data_cut_all", "SHMS #beta (all_cut); SHMS_gtr_#beta; Counts", nbins, 0.0, 1.4)
P_gtr_dp_pions_data_cut_all = ROOT.TH1D("P_gtr_dp_protons_pions_data_cut_all", "SHMS #delta (cut_all); SHMS_gtr_dp; Counts", nbins, -30, 30)
P_cal_etottracknorm_pions_data_cut_all = ROOT.TH1D("P_cal_etottracknorm_pions_data_cut_all", "SHMS cal etottracknorm (cut_all); SHMS_cal_etottracknorm; Counts", nbins, 0, 1.8)
P_hgcer_npeSum_pions_data_cut_all = ROOT.TH1D("P_hgcer_npeSum_pions_data_cut_all", "SHMS HGC npeSum (all_cut); SHMS_hgcer_npeSum; Counts", nbins, 0, 50)
P_hgcer_xAtCer_pions_data_cut_all = ROOT.TH1D("P_hgcer_xAtCer_pions_data_cut_all", "SHMS HGC xAtCer (all_cut); SHMS_hgcer_xAtCer; Counts", nbins, -60, 60)
P_hgcer_yAtCer_pions_data_cut_all = ROOT.TH1D("P_hgcer_yAtCer_pions_data_cut_all", "SHMS HGC yAtCer (all_cut); SHMS_hgcer_yAtCer; Counts", nbins, -60, 60)
P_ngcer_npeSum_pions_data_cut_all = ROOT.TH1D("P_ngcer_npeSum_pions_data_cut_all", "SHMS NGC npeSum (all_cut); SHMS_ngcer_npeSum; Counts", nbins, 0, 50)
P_ngcer_xAtCer_pions_data_cut_all = ROOT.TH1D("P_ngcer_xAtCer_pions_data_cut_all", "SHMS NGC xAtCer (all_cut); SHMS_ngcer_xAtCer; Counts", nbins, -60, 60)
P_ngcer_yAtCer_pions_data_cut_all = ROOT.TH1D("P_ngcer_yAtCer_pions_data_cut_all", "SHMS NGC yAtCer (all_cut); SHMS_ngcer_yAtCer; Counts", nbins, -60, 60)
P_aero_npeSum_pions_data_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_cut_all", "SHMS aero npeSum (all_cut); SHMS_aero_npeSum; Counts", nbins, 0, 50)
P_aero_xAtAero_pions_data_cut_all = ROOT.TH1D("P_acero_xAtAero_pions_data_cut_all", "SHMS aero xAtAero (all_cut); SHMS_aero_xAtAero; Counts", nbins, -60, 60)
P_aero_yAtAero_pions_data_cut_all = ROOT.TH1D("P_aero_yAtAero_pions_data_cut_all", "SHMS aero yAtAero (all_cut); SHMS_aero_yAtAero; Counts", nbins, -60, 60)
P_kin_MMpi_pions_data_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_cut_all", "MIssing Mass data (all_cut); MM_{\pi}; Counts", nbins, 0, 2.0)
P_RFTime_Dist_pions_data_cut_all = ROOT.TH1D("P_RFTime_Dist_pions_data_cut_all", "SHMS RFTime (all_cut); SHMS_RFTime; Counts", nbins, 0, 4)
CTime_ePiCoinTime_ROC1_pions_data_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC1_pions_data_cut_all", "Electron-Pion CTime (all_cut); e pi Coin_Time; Counts", nbins, -50, 50)

# Histogram for PID Test
P_kin_MMpi_pions_data_prompt_aerocut_all = ROOT.TH1D("P_kin_MMpi_pions_data_prompt_aerocut_all", "MIssing Mass data (prompt+accpt_cut+aero_cut); MM_{\pi}; Counts", nbins, 0.0, 2.0)
P_kin_MMpi_pions_data_prompt_RFcut_all = ROOT.TH1D("P_kin_MMpi_pions_data_prompt_RFcut_all", "MIssing Mass data (prompt+accpt_cut+RF_cut); MM_{\pi}; Counts", nbins, 0.0, 2.0)
P_kin_MMpi_pions_data_random_aerocut_all = ROOT.TH1D("P_kin_MMpi_pions_data_random_aerocut_all", "MIssing Mass data (random+accpt_cut+aero_cut); MM_{\pi}; Counts", nbins, 0, 2.0)
P_kin_MMpi_pions_data_random_RFcut_all = ROOT.TH1D("P_kin_MMpi_pions_data_random_RFcut_all", "MIssing Mass data (random+accpt_cut+RF_cut); MM_{\pi}; Counts", nbins, 0, 2.0)
P_kin_MMpi_pions_data_aerocut_all = ROOT.TH1D("P_kin_MMpi_pions_data_aerocut_all", "MIssing Mass data (CT+accpt+aero_cut); MM_{\pi}; Counts", nbins, 0, 2.0)
P_kin_MMpi_pions_data_RFcut_all = ROOT.TH1D("P_kin_MMpi_pions_data_RFcut_all", "MIssing Mass data (CT+accpt+RF_cut); MM_{\pi}; Counts", nbins, 0, 2.0)
CTime_ePiCoinTime_ROC1_pions_data_pid_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC1_pions_data_pid_cut_all", "Electron-Pion CTime (accpt_PID cut); e pi Coin_Time; Counts", nbins, -50, 50)

# 2D Histograms
# Uncut Histograms
P_RFTime_vs_gtr_dp_pions_uncut = ROOT.TH2D("P_RFTime_vs_gtr_dp_pions_uncut","P_RFTime_Dist vs P_gtr_dp (uncut); P_RFTime_Dist; P_gtr_dp", 200, 0, 4, 200, -30, 30)
H_cal_etottracknorm_vs_cer_npeSum_pions_uncut = ROOT.TH2D("H_cal_etottracknorm_vs_cer_npeSum_pions_uncut","HMS cal etottracknorm vs HMS cer npeSum (uncut); H_cal_etottracknorm; H_cer_npeSum",100, 0, 2, 100, 0, 50)
P_hgcer_vs_aero_npe_pions_uncut = ROOT.TH2D("P_hgcer_vs_aero_npe_pions_uncut", "SHMS HGC npeSum vs SHMS Aero npeSum (uncut); SHMS_hgcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
P_ngcer_vs_hgcer_npe_pions_uncut = ROOT.TH2D("P_ngcer_vs_hgcer_npe_pions_uncut", "SHMS NGC npeSum vs SHMS HGC npeSum (uncut); SHMS_ngcer_npeSum; SHMS_hgcer_npeSum", 100, 0, 50, 100, 0, 50)
P_ngcer_vs_aero_npe_pions_uncut = ROOT.TH2D("P_ngcer_vs_aero_npe_pions_uncut", "SHMS NGC npeSum vs SHMS aero npeSum (uncut); SHMS_ngcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_uncut = ROOT.TH2D("P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_uncut", "SHMS HGC yAtCer vs SHMS HGC xAtCer (uncut); SHMS_hgcer_yAtCer; SHMS_hgcer_xAtCer", 100, -50, 50, 100, -50, 50)
P_aero_yAtAero_vs_aero_xAtAero_pions_uncut = ROOT.TH2D("P_aero_yAtAero_vs_aero_xAtAero_pions_uncut", "SHMS aero yAtAero vs SHMS aero xAtAero (uncut); SHMS_aero_yAtAero; SHMS_aero_xAtAero", 100, -50, 50, 100, -50, 50)
P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_uncut = ROOT.TH2D("P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_uncut", "SHMS NGC yAtCer vs SHMS NGC xAtCer (uncut); SHMS_ngcer_yAtCer; SHMS_ngcer_xAtCer", 100, -50, 50, 100, -50, 50)
P_cal_etottracknorm_vs_ngcer_npe_pions_uncut = ROOT.TH2D("P_cal_etottracknorm_vs_ngcer_npe_pions_uncut", "SHMS cal etottracknorm vs SHMS NGC xAtCer (uncut); SHMS_cal_etottracknorm; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
CTime_ePiCoinTime_vs_MMpi_pions_uncut = ROOT.TH2D("CTime_ePiCoinTime_vs_MMpi_pions_uncut","Electron-Pion CTime vs Missing Mass (uncut); e #pi Coin_Time; MM_{#pi}", 300, -30, 30, 100, 0, 2)
CTime_ePiCoinTime_vs_beta_pions_uncut = ROOT.TH2D("CTime_ePiCoinTime_vs_beta_pions_uncut", "Electron-Pion CTime vs SHMS #beta (uncut); e #pi Coin_Time; SHMS_#beta", 120, -30, 30, 200, 0, 2)
P_RFTime_vs_MMpi_pions_uncut = ROOT.TH2D("P_RFTime_vs_MMpi_pions_uncut", "SHMS RFTime vs Missing Mass (uncut); SHMS_RFTime_Dist; MM_{#pi}", 100, 0, 4, 100, 0, 2)
CTime_ePiCoinTime_vs_RFTime_pions_uncut = ROOT.TH2D("CTime_ePiCoinTime_vs_RFTime_pions_uncut", "Electron-Pion CTime vs SHMS RFTime (uncut); e #pi Coin_Time; SHMS_RFTime_Dist", 300, -30, 30, 200, 0, 4)

# Acceptance Cut Histograms
P_RFTime_vs_gtr_dp_pions_accpt_cut_all = ROOT.TH2D("P_RFTime_vs_gtr_dp_pions_accpt_cut_all","P_RFTime_Dist vs P_gtr_dp (accpt_cut); P_RFTime_Dist; P_gtr_dp", 200, 0, 4, 200, -30, 30)
H_cal_etottracknorm_vs_cer_npeSum_pions_accpt_cut_all = ROOT.TH2D("H_cal_etottracknorm_vs_cer_npeSum_pions_accpt_cut_all","HMS cal etottracknorm vs HMS cer npeSum (accpt_cut); H_cal_etottracknorm; H_cer_npeSum",100, 0, 2, 100, 0, 40)
P_hgcer_vs_aero_npe_pions_accpt_cut_all = ROOT.TH2D("P_hgcer_vs_aero_npe_pions_accpt_cut_all", "SHMS HGC npeSum vs SHMS Aero npeSum (accpt_cut); SHMS_hgcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
P_ngcer_vs_hgcer_npe_pions_accpt_cut_all = ROOT.TH2D("P_ngcer_vs_hgcer_npe_pions_accpt_cut_all", "SHMS NGC npeSum vs SHMS HGC npeSum (accpt_cut); SHMS_ngcer_npeSum; SHMS_hgcer_npeSum", 100, 0, 50, 100, 0, 50)
P_ngcer_vs_aero_npe_pions_accpt_cut_all = ROOT.TH2D("P_ngcer_vs_aero_npe_pions_accpt_cut_all", "SHMS NGC npeSum vs SHMS aero npeSum (accpt_cut); SHMS_ngcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_accpt_cut_all = ROOT.TH2D("P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_accpt_cut_all", "SHMS HGC yAtCer vs SHMS HGC xAtCer (accpt_cut); SHMS_hgcer_yAtCer; SHMS_hgcer_xAtCer", 100, -50, 50, 100, -50, 50)
P_aero_yAtAero_vs_aero_xAtAero_pions_accpt_cut_all = ROOT.TH2D("P_aero_yAtAero_vs_aero_xAtAero_pions_accpt_cut_all", "SHMS aero yAtAero vs SHMS aero xAtAero (accpt_cut); SHMS_aero_yAtAero; SHMS_aero_xAtAero", 100, -50, 50, 100, -50, 50)
P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_accpt_cut_all = ROOT.TH2D("P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_accpt_cut_all", "SHMS NGC yAtCer vs SHMS NGC xAtCer (accpt_cut); SHMS_ngcer_yAtCer; SHMS_ngcer_xAtCer", 100, -50, 50, 100, -50, 50)
P_cal_etottracknorm_vs_ngcer_npe_pions_accpt_cut_all = ROOT.TH2D("P_cal_etottracknorm_vs_ngcer_npe_pions_accpt_cut_all", "SHMS cal etottracknorm vs SHMS NGC xAtCer (accpt_cut); SHMS_cal_etottracknorm; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
CTime_ePiCoinTime_vs_MMpi_pions_accpt_cut_all = ROOT.TH2D("CTime_ePiCoinTime_vs_MMpi_pions_accpt_cut_all","Electron-Pion CTime vs Missing Mass (accpt_cut); e #pi Coin_Time; MM_{#pi}", 300, -30, 30, 100, 0, 2)
CTime_ePiCoinTime_vs_beta_pions_accpt_cut_all = ROOT.TH2D("CTime_ePiCoinTime_vs_beta_pions_accpt_cut_all", "Electron-Pion CTime vs SHMS #beta (accpt_cut); e #pi Coin_Time; SHMS_#beta", 120, -30, 30, 200, 0, 2)
P_RFTime_vs_MMpi_pions_accpt_cut_all = ROOT.TH2D("P_RFTime_vs_MMpi_pions_accpt_cut_all", "SHMS RFTime vs Missing Mass (accpt_cut); SHMS_RFTime_Dist; MM_{#pi}", 100, 0, 4, 100, 0, 2)
CTime_ePiCoinTime_vs_RFTime_pions_accpt_cut_all = ROOT.TH2D("CTime_ePiCoinTime_vs_RFTime_pions_accpt_cut_all", "Electron-Pion CTime vs SHMS RFTime (accpt_cut); e #pi Coin_Time; SHMS_RFTime_Dist", 300, -30, 30, 200, 0, 4)

# Prompt + Acceptance Cut Histograms
P_RFTime_vs_gtr_dp_pions_prompt_cut_all = ROOT.TH2D("P_RFTime_vs_gtr_dp_pions_prompt_cut_all","P_RFTime_Dist vs P_gtr_dp (prompt+accpt_cut); P_RFTime_Dist; P_gtr_dp", 200, 0, 4, 200, -30, 30)
H_cal_etottracknorm_vs_cer_npeSum_pions_prompt_cut_all = ROOT.TH2D("H_cal_etottracknorm_vs_cer_npeSum_pions_prompt_cut_all","HMS cal etottracknorm vs HMS cer npeSum (prompt+accpt_cut); H_cal_etottracknorm; H_cer_npeSum",100, 0, 2, 100, 0, 40)
P_hgcer_vs_aero_npe_pions_prompt_cut_all = ROOT.TH2D("P_hgcer_vs_aero_npe_pions_prompt_cut_all", "SHMS HGC npeSum vs SHMS Aero npeSum (prompt+accpt_cut); SHMS_hgcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
P_ngcer_vs_hgcer_npe_pions_prompt_cut_all = ROOT.TH2D("P_ngcer_vs_hgcer_npe_pions_prompt_cut_all", "SHMS NGC npeSum vs SHMS HGC npeSum (prompt+accpt_cut); SHMS_ngcer_npeSum; SHMS_hgcer_npeSum", 100, 0, 50, 100, 0, 50)
P_ngcer_vs_aero_npe_pions_prompt_cut_all = ROOT.TH2D("P_ngcer_vs_aero_npe_pions_prompt_cut_all", "SHMS NGC npeSum vs SHMS aero npeSum (prompt+accpt_cut); SHMS_ngcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_prompt_cut_all = ROOT.TH2D("P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_prompt_cut_all", "SHMS HGC yAtCer vs SHMS HGC xAtCer (prompt+accpt_cut); SHMS_hgcer_yAtCer; SHMS_hgcer_xAtCer", 100, -50, 50, 100, -50, 50)
P_aero_yAtAero_vs_aero_xAtAero_pions_prompt_cut_all = ROOT.TH2D("P_aero_yAtAero_vs_aero_xAtAero_pions_prompt_cut_all", "SHMS aero yAtAero vs SHMS aero xAtAero (prompt+accpt_cut); SHMS_aero_yAtAero; SHMS_aero_xAtAero", 100, -50, 50, 100, -50, 50)
P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_prompt_cut_all = ROOT.TH2D("P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_prompt_cut_all", "SHMS NGC yAtCer vs SHMS NGC xAtCer (prompt+accpt_cut); SHMS_ngcer_yAtCer; SHMS_ngcer_xAtCer", 100, -50, 50, 100, -50, 50)
P_cal_etottracknorm_vs_ngcer_npe_pions_prompt_cut_all = ROOT.TH2D("P_cal_etottracknorm_vs_ngcer_npe_pions_prompt_cut_all", "SHMS cal etottracknorm vs SHMS NGC xAtCer (prompt+accpt_cut); SHMS_cal_etottracknorm; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
CTime_ePiCoinTime_vs_MMpi_pions_prompt_cut_all = ROOT.TH2D("CTime_ePiCoinTime_vs_MMpi_pions_prompt_cut_all","Electron-Pion CTime vs Missing Mass (prompt+accpt_cut); e #pi Coin_Time; MM_{#pi}", 300, -30, 30, 200, 0, 2)
CTime_ePiCoinTime_vs_beta_pions_prompt_cut_all = ROOT.TH2D("CTime_ePiCoinTime_vs_beta_pions_prompt_cut_all", "Electron-Pion CTime vs SHMS #beta (prompt+accpt_cut); e #pi Coin_Time; SHMS_#beta", 120, -30, 30, 200, 0, 2)
P_RFTime_vs_MMpi_pions_prompt_cut_all = ROOT.TH2D("P_RFTime_vs_MMpi_pions_prompt_cut_all", "SHMS RFTime vs Missing Mass (prompt+accpt_cut); SHMS_RFTime_Dist; MM_{#pi}", 100, 0, 4, 100, 0, 2)
CTime_ePiCoinTime_vs_RFTime_pions_prompt_cut_all = ROOT.TH2D("CTime_ePiCoinTime_vs_RFTime_pions_prompt_cut_all", "Electron-Pion CTime vs SHMS RFTime (prompt+accpt_cut); e #pi Coin_Time; SHMS_RFTime_Dist", 300, -30, 30, 100, 0, 4)

# Acceptance + Random Cut Histograms
P_RFTime_vs_gtr_dp_pions_random_cut_all = ROOT.TH2D("P_RFTime_vs_gtr_dp_pions_random_cut_all","P_RFTime_Dist vs P_gtr_dp (random+accpt_cut); P_RFTime_Dist; P_gtr_dp", 200, 0, 4, 200, -30, 30)
H_cal_etottracknorm_vs_cer_npeSum_pions_random_cut_all = ROOT.TH2D("H_cal_etottracknorm_vs_cer_npeSum_pions_random_cut_all","HMS cal etottracknorm vs HMS cer npeSum (random+accpt_cut); H_cal_etottracknorm; H_cer_npeSum",100, 0, 2, 100, 0, 40)
P_hgcer_vs_aero_npe_pions_random_cut_all = ROOT.TH2D("P_hgcer_vs_aero_npe_pions_random_cut_all", "SHMS HGC npeSum vs SHMS Aero npeSum (random+accpt_cut); SHMS_hgcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
P_ngcer_vs_hgcer_npe_pions_random_cut_all = ROOT.TH2D("P_ngcer_vs_hgcer_npe_pions_random_cut_all", "SHMS NGC npeSum vs SHMS HGC npeSum (random+accpt_cut); SHMS_ngcer_npeSum; SHMS_hgcer_npeSum", 100, 0, 50, 100, 0, 50)
P_ngcer_vs_aero_npe_pions_random_cut_all = ROOT.TH2D("P_ngcer_vs_aero_npe_pions_random_cut_all", "SHMS NGC npeSum vs SHMS aero npeSum (random+accpt_cut); SHMS_ngcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_random_cut_all = ROOT.TH2D("P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_random_cut_all", "SHMS HGC yAtCer vs SHMS HGC xAtCer (random+accpt_cut); SHMS_hgcer_yAtCer; SHMS_hgcer_xAtCer", 100, -50, 50, 100, -50, 50)
P_aero_yAtAero_vs_aero_xAtAero_pions_random_cut_all = ROOT.TH2D("P_aero_yAtAero_vs_aero_xAtAero_pions_random_cut_all", "SHMS aero yAtAero vs SHMS aero xAtAero (random+accpt_cut); SHMS_aero_yAtAero; SHMS_aero_xAtAero", 100, -50, 50, 100, -50, 50)
P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_random_cut_all = ROOT.TH2D("P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_random_cut_all", "SHMS NGC yAtCer vs SHMS NGC xAtCer (random+accpt_cut); SHMS_ngcer_yAtCer; SHMS_ngcer_xAtCer", 100, -50, 50, 100, -50, 50)
P_cal_etottracknorm_vs_ngcer_npe_pions_random_cut_all = ROOT.TH2D("P_cal_etottracknorm_vs_ngcer_npe_pions_random_cut_all", "SHMS cal etottracknorm vs SHMS NGC xAtCer (random+accpt_cut); SHMS_cal_etottracknorm; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
CTime_ePiCoinTime_vs_MMpi_pions_random_cut_all = ROOT.TH2D("CTime_ePiCoinTime_vs_MMpi_pions_random_cut_all","Electron-Pion CTime vs Missing Mass (random+accpt_cut); e #pi Coin_Time; MM_{#pi}", 300, -30, 30, 200, 0, 2)
CTime_ePiCoinTime_vs_beta_pions_random_cut_all = ROOT.TH2D("CTime_ePiCoinTime_vs_beta_pions_random_cut_all", "Electron-Pion CTime vs SHMS #beta (random+accpt_cut); e #pi Coin_Time; SHMS_#beta", 120, -30, 30, 200, 0, 2)
P_RFTime_vs_MMpi_pions_random_cut_all = ROOT.TH2D("P_RFTime_vs_MMpi_pions_random_cut_all", "SHMS RFTime vs Missing Mass (random+accpt_cut); SHMS_RFTime_Dist; MM_{#pi}", 100, 0, 4, 100, 0, 2)
CTime_ePiCoinTime_vs_RFTime_pions_random_cut_all = ROOT.TH2D("CTime_ePiCoinTime_vs_RFTime_pions_random_cut_all", "Electron-Pion CTime vs SHMS RFTime (random+accpt_cut); e #pi Coin_Time; SHMS_RFTime_Dist", 300, -30, 30, 100, 0, 4)

# All Cuts Histograms
P_RFTime_vs_gtr_dp_pions_cut_all = ROOT.TH2D("P_RFTime_vs_gtr_dp_pions_cut_all","P_RFTime_Dist vs P_gtr_dp (all_cut); P_RFTime_Dist; P_gtr_dp", 200, 0, 4, 200, -30, 30)
H_cal_etottracknorm_vs_cer_npeSum_pions_cut_all = ROOT.TH2D("H_cal_etottracknorm_vs_cer_npeSum_pions_cut_all","HMS cal etottracknorm vs HMS cer npeSum (all_cut); H_cal_etottracknorm; H_cer_npeSum",100, 0, 2, 100, 0, 40)
P_hgcer_vs_aero_npe_pions_cut_all = ROOT.TH2D("P_hgcer_vs_aero_npe_pions_cut_all", "SHMS HGC npeSum vs SHMS Aero npeSum (all_cut); SHMS_hgcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
P_ngcer_vs_hgcer_npe_pions_cut_all = ROOT.TH2D("P_ngcer_vs_hgcer_npe_pions_cut_all", "SHMS NGC npeSum vs SHMS HGC npeSum (all_cut); SHMS_ngcer_npeSum; SHMS_hgcer_npeSum", 100, 0, 50, 100, 0, 50)
P_ngcer_vs_aero_npe_pions_cut_all = ROOT.TH2D("P_ngcer_vs_aero_npe_pions_cut_all", "SHMS NGC npeSum vs SHMS aero npeSum (all_cut); SHMS_ngcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_cut_all = ROOT.TH2D("P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_cut_all", "SHMS HGC yAtCer vs SHMS HGC xAtCer (all_cut); SHMS_hgcer_yAtCer; SHMS_hgcer_xAtCer", 100, -50, 50, 100, -50, 50)
P_aero_yAtAero_vs_aero_xAtAero_pions_cut_all = ROOT.TH2D("P_aero_yAtAero_vs_aero_xAtAero_pions_cut_all", "SHMS aero yAtAero vs SHMS aero xAtAero (all_cut); SHMS_aero_yAtAero; SHMS_aero_xAtAero", 100, -50, 50, 100, -50, 50)
P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_cut_all = ROOT.TH2D("P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_cut_all", "SHMS NGC yAtCer vs SHMS NGC xAtCer (all_cut); SHMS_ngcer_yAtCer; SHMS_ngcer_xAtCer", 100, -50, 50, 100, -50, 50)
P_cal_etottracknorm_vs_ngcer_npe_pions_cut_all = ROOT.TH2D("P_cal_etottracknorm_vs_ngcer_npe_pions_cut_all", "SHMS cal etottracknorm vs SHMS NGC xAtCer (all_cut); SHMS_cal_etottracknorm; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
CTime_ePiCoinTime_vs_MMpi_pions_cut_all = ROOT.TH2D("CTime_ePiCoinTime_vs_MMpi_pions_cut_all","Electron-Pion CTime vs Missing Mass (all_cut); e #pi Coin_Time; MM_{#pi}", 300, -30, 30, 200, 0.0, 2.0)
CTime_ePiCoinTime_vs_beta_pions_cut_all = ROOT.TH2D("CTime_ePiCoinTime_vs_beta_pions_cut_all", "Electron-Pion CTime vs SHMS #beta (all_cut); e #pi Coin_Time; SHMS_#beta", 120, -30, 30, 200, 0, 2)
P_RFTime_vs_MMpi_pions_cut_all = ROOT.TH2D("P_RFTime_vs_MMpi_pions_cut_all", "SHMS RFTime vs Missing Mass (all_cut); SHMS_RFTime_Dist; MM_{#pi}", 100, 0, 4, 100, 0, 2)
CTime_ePiCoinTime_vs_RFTime_pions_cut_all = ROOT.TH2D("CTime_ePiCoinTime_vs_RFTime_pions_cut_all", "Electron-Pion CTime vs SHMS RFTime (all_cut); e #pi Coin_Time; SHMS_RFTime_Dist", 300, -30, 30, 100, 0, 4)

# 3D Histograms
P_HGC_xy_npe_pions_uncut = ROOT.TH3D("P_HGC_xy_npe_pions_uncut", "SHMS HGC NPE as fn of yAtCer vs SHMS HGC xAtCer (no cuts); HGC_yAtCer(cm); HGC_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_Aero_xy_npe_pions_uncut = ROOT.TH3D("P_Aero_xy_npe_pions_uncut", "SHMS Aerogel NPE as fn of yAtCer vs xAtCer (no cuts); Aero_yAtCer(cm); Aero_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_NGC_xy_npe_pions_uncut = ROOT.TH3D("P_NGC_xy_npe_pions_uncut", "SHMS NGC NPE as fn of yAtCer vs xAtCer (no cuts); NGC_yAtCer(cm); NGC_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_HGC_xy_npe_pions_accpt_cut_all = ROOT.TH3D("P_HGC_xy_npe_pions_accpt_cut_all", "SHMS HGC NPE as fn of yAtCer vs SHMS HGC xAtCer (accpt_cut); HGC_yAtCer(cm); HGC_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_Aero_xy_npe_pions_accpt_cut_all = ROOT.TH3D("P_Aero_xy_npe_pions_accpt_cut_all", "SHMS Aerogel NPE as fn of yAtCer vs xAtCer (accpt_cut); Aero_yAtCer(cm); Aero_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_NGC_xy_npe_pions_accpt_cut_all = ROOT.TH3D("P_NGC_xy_npe_pions_accpt_cut_all", "SHMS NGC NPE as fn of yAtCer vs xAtCer (accpt_cut); NGC_yAtCer(cm); NGC_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_HGC_xy_npe_pions_prompt_cut_all = ROOT.TH3D("P_HGC_xy_npe_pions_prompt_cut_all", "SHMS HGC NPE as fn of yAtCer vs SHMS HGC xAtCer (prompt+accpt_cut); HGC_yAtCer(cm); HGC_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_Aero_xy_npe_pions_prompt_cut_all = ROOT.TH3D("P_Aero_xy_npe_pions_prompt_cut_all", "SHMS Aerogel NPE as fn of yAtCer vs xAtCer (prompt+accpt_cut); Aero_yAtCer(cm); Aero_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_NGC_xy_npe_pions_prompt_cut_all = ROOT.TH3D("P_NGC_xy_npe_pions_prompt_cut_all", "SHMS NGC NPE as fn of yAtCer vs xAtCer (prompt+accpt_cut); NGC_yAtCer(cm); NGC_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_HGC_xy_npe_pions_random_cut_all = ROOT.TH3D("P_HGC_xy_npe_pions_random_cut_all", "SHMS HGC NPE as fn of yAtCer vs SHMS HGC xAtCer (random+accpt_cut); HGC_yAtCer(cm); HGC_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_Aero_xy_npe_pions_random_cut_all = ROOT.TH3D("P_Aero_xy_npe_pions_random_cut_all", "SHMS Aerogel NPE as fn of yAtCer vs xAtCer (random+accpt_cut); Aero_yAtCer(cm); Aero_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_NGC_xy_npe_pions_random_cut_all = ROOT.TH3D("P_NGC_xy_npe_pions_random_cut_all", "SHMS NGC NPE as fn of yAtCer vs xAtCer (random+accpt_cut); NGC_yAtCer(cm); NGC_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_HGC_xy_npe_pions_cut_all = ROOT.TH3D("P_HGC_xy_npe_pions_cut_all", "SHMS HGC NPE as fn of yAtCer vs SHMS HGC xAtCer (with cuts); HGC_yAtCer(cm); HGC_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_Aero_xy_npe_pions_cut_all = ROOT.TH3D("P_Aero_xy_npe_pions_cut_all", "SHMS Aerogel NPE as fn of yAtCer vs xAtCer (with cuts); Aero_yAtCer(cm); Aero_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_NGC_xy_npe_pions_cut_all = ROOT.TH3D("P_NGC_xy_npe_pions_cut_all", "SHMS NGC NPE as fn of yAtCer vs xAtCer (with cuts); NGC_yAtCer(cm); NGC_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)

#################################################################################################################################################

# PID Cut Values
H_cer_npeSum_cut_value = 1.5
H_cal_etottracknorm_cut_value = 0.7
P_aero_npeSum_cut_value = 2.5
#P_hgcer_npeSum_cut_value = 1.5
#P_ngcer_npeSum_cut_value = 0.5
#P_RF_Dist_low_cut_value = 1.2
#P_RF_Dist_high_cut_value = 3.4
P_RF_Dist_low_cut_value = 0.0
P_RF_Dist_high_cut_value = 1.6

# Fill Uncut Hitograms
for event in Uncut_Pion_Events_Data_tree:
    H_gtr_beta_pions_data_uncut.Fill(event.H_gtr_beta)
    H_cal_etottracknorm_pions_data_uncut.Fill(event.H_cal_etottracknorm)
    H_cer_npeSum_pions_data_uncut.Fill(event.H_cer_npeSum)
    H_RFTime_Dist_pions_data_uncut.Fill(event.H_RF_Dist)
    P_gtr_beta_pions_data_uncut.Fill(event.P_gtr_beta)
    P_gtr_dp_pions_data_uncut.Fill(event.P_gtr_dp)
    P_cal_etottracknorm_pions_data_uncut.Fill(event.P_cal_etottracknorm)
    P_hgcer_npeSum_pions_data_uncut.Fill(event.P_hgcer_npeSum)
    P_hgcer_xAtCer_pions_data_uncut.Fill(event.P_hgcer_xAtCer)
    P_hgcer_yAtCer_pions_data_uncut.Fill(event.P_hgcer_yAtCer)
    P_ngcer_npeSum_pions_data_uncut.Fill(event.P_ngcer_npeSum)
    P_ngcer_xAtCer_pions_data_uncut.Fill(event.P_ngcer_xAtCer)
    P_ngcer_yAtCer_pions_data_uncut.Fill(event.P_ngcer_yAtCer)
    P_aero_npeSum_pions_data_uncut.Fill(event.P_aero_npeSum)
    P_aero_xAtAero_pions_data_uncut.Fill(event.P_aero_xAtAero)
    P_aero_yAtAero_pions_data_uncut.Fill(event.P_aero_yAtAero)
    P_kin_MMpi_pions_data_uncut.Fill(event.MMpi)
    P_RFTime_Dist_pions_data_uncut.Fill(event.P_RF_Dist)
    CTime_ePiCoinTime_ROC1_pions_data_uncut.Fill(event.CTime_ePiCoinTime_ROC1 + 0.25)
    H_cal_etottracknorm_vs_cer_npeSum_pions_uncut.Fill(event.H_cal_etottracknorm, event.H_cer_npeSum)
    P_hgcer_vs_aero_npe_pions_uncut.Fill(event.P_hgcer_npeSum, event.P_aero_npeSum)
    P_ngcer_vs_hgcer_npe_pions_uncut.Fill(event.P_ngcer_npeSum, event.P_hgcer_npeSum)
    P_ngcer_vs_aero_npe_pions_uncut.Fill(event.P_ngcer_npeSum, event.P_aero_npeSum)
    P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_uncut.Fill(event.P_hgcer_yAtCer, event.P_hgcer_xAtCer)
    P_aero_yAtAero_vs_aero_xAtAero_pions_uncut.Fill(event.P_aero_yAtAero, event.P_aero_xAtAero)
    P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_uncut.Fill(event.P_ngcer_yAtCer, event.P_ngcer_xAtCer)
    P_cal_etottracknorm_vs_ngcer_npe_pions_uncut.Fill(event.P_cal_etottracknorm, event.P_ngcer_npeSum)
    CTime_ePiCoinTime_vs_MMpi_pions_uncut.Fill(event.CTime_ePiCoinTime_ROC1 + 0.25, event.MMpi)
    CTime_ePiCoinTime_vs_beta_pions_uncut.Fill(event.CTime_ePiCoinTime_ROC1 + 0.25, event.P_gtr_beta)
    P_RFTime_vs_MMpi_pions_uncut.Fill(event.P_RF_Dist, event.MMpi)
    CTime_ePiCoinTime_vs_RFTime_pions_uncut.Fill(event.CTime_ePiCoinTime_ROC1 + 0.25, event.P_RF_Dist)
    P_HGC_xy_npe_pions_uncut.Fill(event.P_hgcer_yAtCer,event.P_hgcer_xAtCer,event.P_hgcer_npeSum)
    P_Aero_xy_npe_pions_uncut.Fill(event.P_aero_yAtAero,event.P_aero_xAtAero,event.P_aero_npeSum)
    P_NGC_xy_npe_pions_uncut.Fill(event.P_ngcer_yAtCer,event.P_ngcer_xAtCer,event.P_ngcer_npeSum)
    P_RFTime_vs_gtr_dp_pions_uncut.Fill(event.P_RF_Dist, event.P_gtr_dp)

# Fill Accpt Cut Hitograms
for event in Cut_Pion_Events_Accpt_Data_tree:
    H_gtr_beta_pions_data_accpt_cut_all.Fill(event.H_gtr_beta)
    H_cal_etottracknorm_pions_data_accpt_cut_all.Fill(event.H_cal_etottracknorm)
    H_cer_npeSum_pions_data_accpt_cut_all.Fill(event.H_cer_npeSum)
    H_RFTime_Dist_pions_data_accpt_cut_all.Fill(event.H_RF_Dist)
    P_gtr_beta_pions_data_accpt_cut_all.Fill(event.P_gtr_beta)
    P_gtr_dp_pions_data_accpt_cut_all.Fill(event.P_gtr_dp)
    P_cal_etottracknorm_pions_data_accpt_cut_all.Fill(event.P_cal_etottracknorm)
    P_hgcer_npeSum_pions_data_accpt_cut_all.Fill(event.P_hgcer_npeSum)
    P_hgcer_xAtCer_pions_data_accpt_cut_all.Fill(event.P_hgcer_xAtCer)
    P_hgcer_yAtCer_pions_data_accpt_cut_all.Fill(event.P_hgcer_yAtCer)
    P_ngcer_npeSum_pions_data_accpt_cut_all.Fill(event.P_ngcer_npeSum)
    P_ngcer_xAtCer_pions_data_accpt_cut_all.Fill(event.P_ngcer_xAtCer)
    P_ngcer_yAtCer_pions_data_accpt_cut_all.Fill(event.P_ngcer_yAtCer)
    P_aero_npeSum_pions_data_accpt_cut_all.Fill(event.P_aero_npeSum)
    P_aero_xAtAero_pions_data_accpt_cut_all.Fill(event.P_aero_xAtAero)
    P_aero_yAtAero_pions_data_accpt_cut_all.Fill(event.P_aero_yAtAero)
    P_kin_MMpi_pions_data_accpt_cut_all.Fill(event.MMpi)
    P_RFTime_Dist_pions_data_accpt_cut_all.Fill(event.P_RF_Dist)
    CTime_ePiCoinTime_ROC1_pions_data_accpt_cut_all.Fill(event.CTime_ePiCoinTime_ROC1 + 0.25)
    H_cal_etottracknorm_vs_cer_npeSum_pions_accpt_cut_all.Fill(event.H_cal_etottracknorm, event.H_cer_npeSum)
    P_hgcer_vs_aero_npe_pions_accpt_cut_all.Fill(event.P_hgcer_npeSum, event.P_aero_npeSum)
    P_ngcer_vs_hgcer_npe_pions_accpt_cut_all.Fill(event.P_ngcer_npeSum, event.P_hgcer_npeSum)
    P_ngcer_vs_aero_npe_pions_accpt_cut_all.Fill(event.P_ngcer_npeSum, event.P_aero_npeSum)
    P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_accpt_cut_all.Fill(event.P_hgcer_yAtCer, event.P_hgcer_xAtCer)
    P_aero_yAtAero_vs_aero_xAtAero_pions_accpt_cut_all.Fill(event.P_aero_yAtAero, event.P_aero_xAtAero)
    P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_accpt_cut_all.Fill(event.P_ngcer_yAtCer, event.P_ngcer_xAtCer)
    P_cal_etottracknorm_vs_ngcer_npe_pions_accpt_cut_all.Fill(event.P_cal_etottracknorm, event.P_ngcer_npeSum)
    CTime_ePiCoinTime_vs_MMpi_pions_accpt_cut_all.Fill(event.CTime_ePiCoinTime_ROC1 + 0.25, event.MMpi)
    CTime_ePiCoinTime_vs_beta_pions_accpt_cut_all.Fill(event.CTime_ePiCoinTime_ROC1 + 0.25, event.P_gtr_beta)
    P_RFTime_vs_MMpi_pions_accpt_cut_all.Fill(event.P_RF_Dist, event.MMpi)
    CTime_ePiCoinTime_vs_RFTime_pions_accpt_cut_all.Fill(event.CTime_ePiCoinTime_ROC1 + 0.25, event.P_RF_Dist)
    P_HGC_xy_npe_pions_accpt_cut_all.Fill(event.P_hgcer_yAtCer,event.P_hgcer_xAtCer,event.P_hgcer_npeSum)
    P_Aero_xy_npe_pions_accpt_cut_all.Fill(event.P_aero_yAtAero,event.P_aero_xAtAero,event.P_aero_npeSum)
    P_NGC_xy_npe_pions_accpt_cut_all.Fill(event.P_ngcer_yAtCer,event.P_ngcer_xAtCer,event.P_ngcer_npeSum)
    P_RFTime_vs_gtr_dp_pions_accpt_cut_all.Fill(event.P_RF_Dist, event.P_gtr_dp)

# Fill Accpt + Prompt Cut Hitograms
for event in Cut_Pion_Events_Prompt_Data_tree:
    HMS_PID_Cut = (event.H_cer_npeSum > H_cer_npeSum_cut_value) & (event.H_cal_etottracknorm > H_cal_etottracknorm_cut_value)
#    SHMS_PID_Cut = (event.P_aero_npeSum > P_aero_npeSum_cut_value)
    SHMS_PID_Cut = (event.P_RF_Dist > P_RF_Dist_low_cut_value) & (event.P_RF_Dist < P_RF_Dist_high_cut_value) & (event.P_aero_npeSum > P_aero_npeSum_cut_value)
    if (HMS_PID_Cut & SHMS_PID_Cut):
        H_gtr_beta_pions_data_prompt_cut_all.Fill(event.H_gtr_beta)
        H_cal_etottracknorm_pions_data_prompt_cut_all.Fill(event.H_cal_etottracknorm)
        H_cer_npeSum_pions_data_prompt_cut_all.Fill(event.H_cer_npeSum)
        H_RFTime_Dist_pions_data_prompt_cut_all.Fill(event.H_RF_Dist)
        P_gtr_beta_pions_data_prompt_cut_all.Fill(event.P_gtr_beta)
        P_gtr_dp_pions_data_prompt_cut_all.Fill(event.P_gtr_dp)
        P_cal_etottracknorm_pions_data_prompt_cut_all.Fill(event.P_cal_etottracknorm)
        P_hgcer_npeSum_pions_data_prompt_cut_all.Fill(event.P_hgcer_npeSum)
        P_hgcer_xAtCer_pions_data_prompt_cut_all.Fill(event.P_hgcer_xAtCer)
        P_hgcer_yAtCer_pions_data_prompt_cut_all.Fill(event.P_hgcer_yAtCer)
        P_ngcer_npeSum_pions_data_prompt_cut_all.Fill(event.P_ngcer_npeSum)
        P_ngcer_xAtCer_pions_data_prompt_cut_all.Fill(event.P_ngcer_xAtCer)
        P_ngcer_yAtCer_pions_data_prompt_cut_all.Fill(event.P_ngcer_yAtCer)
        P_aero_npeSum_pions_data_prompt_cut_all.Fill(event.P_aero_npeSum)
        P_aero_xAtAero_pions_data_prompt_cut_all.Fill(event.P_aero_xAtAero)
        P_aero_yAtAero_pions_data_prompt_cut_all.Fill(event.P_aero_yAtAero)
        P_kin_MMpi_pions_data_prompt_cut_all.Fill(event.MMpi)
        P_RFTime_Dist_pions_data_prompt_cut_all.Fill(event.P_RF_Dist)
        CTime_ePiCoinTime_ROC1_pions_data_prompt_cut_all.Fill(event.CTime_ePiCoinTime_ROC1 + 0.25)
        H_cal_etottracknorm_vs_cer_npeSum_pions_prompt_cut_all.Fill(event.H_cal_etottracknorm, event.H_cer_npeSum)
        P_hgcer_vs_aero_npe_pions_prompt_cut_all.Fill(event.P_hgcer_npeSum, event.P_aero_npeSum)
        P_ngcer_vs_hgcer_npe_pions_prompt_cut_all.Fill(event.P_ngcer_npeSum, event.P_hgcer_npeSum)
        P_ngcer_vs_aero_npe_pions_prompt_cut_all.Fill(event.P_ngcer_npeSum, event.P_aero_npeSum)
        P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_prompt_cut_all.Fill(event.P_hgcer_yAtCer, event.P_hgcer_xAtCer)
        P_aero_yAtAero_vs_aero_xAtAero_pions_prompt_cut_all.Fill(event.P_aero_yAtAero, event.P_aero_xAtAero)
        P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_prompt_cut_all.Fill(event.P_ngcer_yAtCer, event.P_ngcer_xAtCer)
        P_cal_etottracknorm_vs_ngcer_npe_pions_prompt_cut_all.Fill(event.P_cal_etottracknorm, event.P_ngcer_npeSum)
        CTime_ePiCoinTime_vs_MMpi_pions_prompt_cut_all.Fill(event.CTime_ePiCoinTime_ROC1 + 0.25, event.MMpi)
        CTime_ePiCoinTime_vs_beta_pions_prompt_cut_all.Fill(event.CTime_ePiCoinTime_ROC1 + 0.25, event.P_gtr_beta)
        P_RFTime_vs_MMpi_pions_prompt_cut_all.Fill(event.P_RF_Dist, event.MMpi)
        CTime_ePiCoinTime_vs_RFTime_pions_prompt_cut_all.Fill(event.CTime_ePiCoinTime_ROC1 + 0.25, event.P_RF_Dist)
        P_HGC_xy_npe_pions_prompt_cut_all.Fill(event.P_hgcer_yAtCer,event.P_hgcer_xAtCer,event.P_hgcer_npeSum)
        P_Aero_xy_npe_pions_prompt_cut_all.Fill(event.P_aero_yAtAero,event.P_aero_xAtAero,event.P_aero_npeSum)
        P_NGC_xy_npe_pions_prompt_cut_all.Fill(event.P_ngcer_yAtCer,event.P_ngcer_xAtCer,event.P_ngcer_npeSum)
        P_RFTime_vs_gtr_dp_pions_prompt_cut_all.Fill(event.P_RF_Dist, event.P_gtr_dp)


# Fill Accpt + Random Cut Hitograms
for event in Cut_Pion_Events_Random_Data_tree:
    HMS_PID_Cut = (event.H_cer_npeSum > H_cer_npeSum_cut_value) & (event.H_cal_etottracknorm > H_cal_etottracknorm_cut_value)
#    SHMS_PID_Cut = (event.P_aero_npeSum > P_aero_npeSum_cut_value)
    SHMS_PID_Cut = (event.P_RF_Dist > P_RF_Dist_low_cut_value) & (event.P_RF_Dist < P_RF_Dist_high_cut_value) & (event.P_aero_npeSum > P_aero_npeSum_cut_value)
    if (HMS_PID_Cut & SHMS_PID_Cut):
        H_gtr_beta_pions_data_random_cut_all.Fill(event.H_gtr_beta)
        H_cal_etottracknorm_pions_data_random_cut_all.Fill(event.H_cal_etottracknorm)
        H_cer_npeSum_pions_data_random_cut_all.Fill(event.H_cer_npeSum)
        H_RFTime_Dist_pions_data_random_cut_all.Fill(event.H_RF_Dist)
        P_gtr_beta_pions_data_random_cut_all.Fill(event.P_gtr_beta)
        P_gtr_dp_pions_data_random_cut_all.Fill(event.P_gtr_dp)
        P_cal_etottracknorm_pions_data_random_cut_all.Fill(event.P_cal_etottracknorm)
        P_hgcer_npeSum_pions_data_random_cut_all.Fill(event.P_hgcer_npeSum)
        P_hgcer_xAtCer_pions_data_random_cut_all.Fill(event.P_hgcer_xAtCer)
        P_hgcer_yAtCer_pions_data_random_cut_all.Fill(event.P_hgcer_yAtCer)
        P_ngcer_npeSum_pions_data_random_cut_all.Fill(event.P_ngcer_npeSum)
        P_ngcer_xAtCer_pions_data_random_cut_all.Fill(event.P_ngcer_xAtCer)
        P_ngcer_yAtCer_pions_data_random_cut_all.Fill(event.P_ngcer_yAtCer)
        P_aero_npeSum_pions_data_random_cut_all.Fill(event.P_aero_npeSum)
        P_aero_xAtAero_pions_data_random_cut_all.Fill(event.P_aero_xAtAero)
        P_aero_yAtAero_pions_data_random_cut_all.Fill(event.P_aero_yAtAero)
        P_kin_MMpi_pions_data_random_cut_all.Fill(event.MMpi)
        P_RFTime_Dist_pions_data_random_cut_all.Fill(event.P_RF_Dist)
        CTime_ePiCoinTime_ROC1_pions_data_random_cut_all.Fill(event.CTime_ePiCoinTime_ROC1 + 0.25)
        CTime_ePiCoinTime_ROC1_pions_data_random_unsub_cut_all.Fill(event.CTime_ePiCoinTime_ROC1 + 0.25)
        H_cal_etottracknorm_vs_cer_npeSum_pions_random_cut_all.Fill(event.H_cal_etottracknorm, event.H_cer_npeSum)
        P_hgcer_vs_aero_npe_pions_random_cut_all.Fill(event.P_hgcer_npeSum, event.P_aero_npeSum)
        P_ngcer_vs_hgcer_npe_pions_random_cut_all.Fill(event.P_ngcer_npeSum, event.P_hgcer_npeSum)
        P_ngcer_vs_aero_npe_pions_random_cut_all.Fill(event.P_ngcer_npeSum, event.P_aero_npeSum)
        P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_random_cut_all.Fill(event.P_hgcer_yAtCer, event.P_hgcer_xAtCer)
        P_aero_yAtAero_vs_aero_xAtAero_pions_random_cut_all.Fill(event.P_aero_yAtAero, event.P_aero_xAtAero)
        P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_random_cut_all.Fill(event.P_ngcer_yAtCer, event.P_ngcer_xAtCer)
        P_cal_etottracknorm_vs_ngcer_npe_pions_random_cut_all.Fill(event.P_cal_etottracknorm, event.P_ngcer_npeSum)
        CTime_ePiCoinTime_vs_MMpi_pions_random_cut_all.Fill(event.CTime_ePiCoinTime_ROC1 + 0.25, event.MMpi)
        CTime_ePiCoinTime_vs_beta_pions_random_cut_all.Fill(event.CTime_ePiCoinTime_ROC1 + 0.25, event.P_gtr_beta)
        P_RFTime_vs_MMpi_pions_random_cut_all.Fill(event.P_RF_Dist, event.MMpi)
        CTime_ePiCoinTime_vs_RFTime_pions_random_cut_all.Fill(event.CTime_ePiCoinTime_ROC1 + 0.25, event.P_RF_Dist)
        P_HGC_xy_npe_pions_random_cut_all.Fill(event.P_hgcer_yAtCer,event.P_hgcer_xAtCer,event.P_hgcer_npeSum)
        P_Aero_xy_npe_pions_random_cut_all.Fill(event.P_aero_yAtAero,event.P_aero_xAtAero,event.P_aero_npeSum)
        P_NGC_xy_npe_pions_random_cut_all.Fill(event.P_ngcer_yAtCer,event.P_ngcer_xAtCer,event.P_ngcer_npeSum)
        P_RFTime_vs_gtr_dp_pions_random_cut_all.Fill(event.P_RF_Dist, event.P_gtr_dp)

# Histogram for PID Study
for event in Cut_Pion_Events_Prompt_Data_tree:
    SHMS_PID_Cut = (event.P_aero_npeSum > P_aero_npeSum_cut_value)
    HMS_PID_Cut = (event.H_cer_npeSum > H_cer_npeSum_cut_value) & (event.H_cal_etottracknorm > H_cal_etottracknorm_cut_value)
    if (HMS_PID_Cut & SHMS_PID_Cut):
       P_kin_MMpi_pions_data_prompt_aerocut_all.Fill(event.MMpi)

for event in Cut_Pion_Events_Random_Data_tree:
    SHMS_PID_Cut = (event.P_aero_npeSum > P_aero_npeSum_cut_value)
    HMS_PID_Cut = (event.H_cer_npeSum > H_cer_npeSum_cut_value) & (event.H_cal_etottracknorm > H_cal_etottracknorm_cut_value)
    if (HMS_PID_Cut & SHMS_PID_Cut):
       P_kin_MMpi_pions_data_random_aerocut_all.Fill(event.MMpi)

for event in Cut_Pion_Events_Prompt_Data_tree:
    SHMS_PID_Cut = (event.P_RF_Dist > P_RF_Dist_low_cut_value) & (event.P_RF_Dist < P_RF_Dist_high_cut_value)
    HMS_PID_Cut = (event.H_cer_npeSum > H_cer_npeSum_cut_value) & (event.H_cal_etottracknorm > H_cal_etottracknorm_cut_value)
    if (HMS_PID_Cut & SHMS_PID_Cut):
       P_kin_MMpi_pions_data_prompt_RFcut_all.Fill(event.MMpi)

for event in Cut_Pion_Events_Random_Data_tree:
    SHMS_PID_Cut = (event.P_RF_Dist > P_RF_Dist_low_cut_value) & (event.P_RF_Dist < P_RF_Dist_high_cut_value)
    HMS_PID_Cut = (event.H_cer_npeSum > H_cer_npeSum_cut_value) & (event.H_cal_etottracknorm > H_cal_etottracknorm_cut_value)
    if (HMS_PID_Cut & SHMS_PID_Cut):
       P_kin_MMpi_pions_data_random_RFcut_all.Fill(event.MMpi)

for event in Cut_Pion_Events_Accpt_Data_tree:
    HMS_PID_Cut = (event.H_cer_npeSum > H_cer_npeSum_cut_value) & (event.H_cal_etottracknorm > H_cal_etottracknorm_cut_value)
#    SHMS_PID_Cut = (event.P_aero_npeSum > P_aero_npeSum_cut_value)
    SHMS_PID_Cut = (event.P_RF_Dist > P_RF_Dist_low_cut_value) & (event.P_RF_Dist < P_RF_Dist_high_cut_value) & (event.P_aero_npeSum > P_aero_npeSum_cut_value)
    if (HMS_PID_Cut & SHMS_PID_Cut):
       CTime_ePiCoinTime_ROC1_pions_data_pid_cut_all.Fill(event.CTime_ePiCoinTime_ROC1 + 0.25)

print("Histograms filled")

#################################################################################################################################################

# Random subtraction from missing mass
#for event in Cut_Pion_Events_Random_tree:
#    P_kin_MMp_pions_cut_random_scaled.Fill(event.MMp)
#    P_kin_MMp_pions_cut_random_scaled.Scale(1.0/nWindows)
#P_kin_MMp_pions_cut_random_sub.Add(P_kin_MMp_pions_cut_prompt, P_kin_MMp_pions_cut_random_scaled, 1, -1)

H_gtr_beta_pions_data_random_cut_all.Scale(1.0/nWindows)
H_cal_etottracknorm_pions_data_random_cut_all.Scale(1.0/nWindows)
H_cer_npeSum_pions_data_random_cut_all.Scale(1.0/nWindows)
H_RFTime_Dist_pions_data_random_cut_all.Scale(1.0/nWindows)
P_gtr_beta_pions_data_random_cut_all.Scale(1.0/nWindows)
P_gtr_dp_pions_data_random_cut_all.Scale(1.0/nWindows)
P_cal_etottracknorm_pions_data_random_cut_all.Scale(1.0/nWindows)
P_hgcer_npeSum_pions_data_random_cut_all.Scale(1.0/nWindows)
P_hgcer_xAtCer_pions_data_random_cut_all.Scale(1.0/nWindows)
P_hgcer_yAtCer_pions_data_random_cut_all.Scale(1.0/nWindows)
P_ngcer_npeSum_pions_data_random_cut_all.Scale(1.0/nWindows)
P_ngcer_xAtCer_pions_data_random_cut_all.Scale(1.0/nWindows)
P_ngcer_yAtCer_pions_data_random_cut_all.Scale(1.0/nWindows)
P_aero_npeSum_pions_data_random_cut_all.Scale(1.0/nWindows)
P_aero_xAtAero_pions_data_random_cut_all.Scale(1.0/nWindows)
P_aero_yAtAero_pions_data_random_cut_all.Scale(1.0/nWindows)
P_RFTime_Dist_pions_data_random_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_ROC1_pions_data_random_cut_all.Scale(1.0/nWindows)
H_cal_etottracknorm_vs_cer_npeSum_pions_random_cut_all.Scale(1.0/nWindows)
P_hgcer_vs_aero_npe_pions_random_cut_all.Scale(1.0/nWindows)
P_ngcer_vs_hgcer_npe_pions_random_cut_all.Scale(1.0/nWindows)
P_ngcer_vs_aero_npe_pions_random_cut_all.Scale(1.0/nWindows)
P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_random_cut_all.Scale(1.0/nWindows)
P_aero_yAtAero_vs_aero_xAtAero_pions_random_cut_all.Scale(1.0/nWindows)
P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_random_cut_all.Scale(1.0/nWindows)
P_cal_etottracknorm_vs_ngcer_npe_pions_random_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_vs_MMpi_pions_random_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_vs_beta_pions_random_cut_all.Scale(1.0/nWindows)
P_RFTime_vs_MMpi_pions_random_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_vs_RFTime_pions_random_cut_all.Scale(1.0/nWindows)
P_HGC_xy_npe_pions_random_cut_all.Scale(1.0/nWindows)
P_Aero_xy_npe_pions_random_cut_all.Scale(1.0/nWindows)
P_NGC_xy_npe_pions_random_cut_all.Scale(1.0/nWindows)
P_RFTime_vs_gtr_dp_pions_random_cut_all.Scale(1.0/nWindows)
P_kin_MMpi_pions_data_random_cut_all.Scale(1.0/nWindows)
# PID Test
P_kin_MMpi_pions_data_random_aerocut_all.Scale(1.0/nWindows)
P_kin_MMpi_pions_data_random_RFcut_all.Scale(1.0/nWindows)


H_gtr_beta_pions_data_cut_all.Add(H_gtr_beta_pions_data_prompt_cut_all, H_gtr_beta_pions_data_random_cut_all, 1, -1)
H_cal_etottracknorm_pions_data_cut_all.Add(H_cal_etottracknorm_pions_data_prompt_cut_all, H_cal_etottracknorm_pions_data_random_cut_all, 1, -1)
H_cer_npeSum_pions_data_cut_all.Add(H_cer_npeSum_pions_data_prompt_cut_all, H_cer_npeSum_pions_data_random_cut_all, 1, -1)
H_RFTime_Dist_pions_data_cut_all.Add(H_RFTime_Dist_pions_data_prompt_cut_all, H_RFTime_Dist_pions_data_random_cut_all, 1, -1)
P_gtr_beta_pions_data_cut_all.Add(P_gtr_beta_pions_data_prompt_cut_all, P_gtr_beta_pions_data_random_cut_all, 1, -1)
P_gtr_dp_pions_data_cut_all.Add(P_gtr_dp_pions_data_prompt_cut_all, P_gtr_dp_pions_data_random_cut_all, 1, -1)
P_cal_etottracknorm_pions_data_cut_all.Add(P_cal_etottracknorm_pions_data_prompt_cut_all, P_cal_etottracknorm_pions_data_random_cut_all, 1, -1)
P_hgcer_npeSum_pions_data_cut_all.Add(P_hgcer_npeSum_pions_data_prompt_cut_all, P_hgcer_npeSum_pions_data_random_cut_all, 1, -1)
P_hgcer_xAtCer_pions_data_cut_all.Add(P_hgcer_xAtCer_pions_data_prompt_cut_all, P_hgcer_xAtCer_pions_data_random_cut_all, 1, -1)
P_hgcer_yAtCer_pions_data_cut_all.Add(P_hgcer_yAtCer_pions_data_prompt_cut_all, P_hgcer_yAtCer_pions_data_random_cut_all, 1, -1)
P_ngcer_npeSum_pions_data_cut_all.Add(P_ngcer_npeSum_pions_data_prompt_cut_all, P_ngcer_npeSum_pions_data_random_cut_all, 1, -1)
P_ngcer_xAtCer_pions_data_cut_all.Add(P_ngcer_xAtCer_pions_data_prompt_cut_all, P_ngcer_xAtCer_pions_data_random_cut_all, 1, -1)
P_ngcer_yAtCer_pions_data_cut_all.Add(P_ngcer_yAtCer_pions_data_prompt_cut_all, P_ngcer_yAtCer_pions_data_random_cut_all, 1, -1)
P_aero_npeSum_pions_data_cut_all.Add(P_aero_npeSum_pions_data_prompt_cut_all, P_aero_npeSum_pions_data_random_cut_all, 1, -1)
P_aero_xAtAero_pions_data_cut_all.Add(P_aero_xAtAero_pions_data_prompt_cut_all, P_aero_xAtAero_pions_data_random_cut_all, 1, -1)
P_aero_yAtAero_pions_data_cut_all.Add(P_aero_yAtAero_pions_data_prompt_cut_all, P_aero_yAtAero_pions_data_random_cut_all, 1, -1)
P_RFTime_Dist_pions_data_cut_all.Add(P_RFTime_Dist_pions_data_prompt_cut_all, P_RFTime_Dist_pions_data_random_cut_all, 1, -1)
CTime_ePiCoinTime_ROC1_pions_data_cut_all.Add(CTime_ePiCoinTime_ROC1_pions_data_prompt_cut_all, CTime_ePiCoinTime_ROC1_pions_data_random_cut_all, 1, -1)
H_cal_etottracknorm_vs_cer_npeSum_pions_cut_all.Add(H_cal_etottracknorm_vs_cer_npeSum_pions_prompt_cut_all, H_cal_etottracknorm_vs_cer_npeSum_pions_random_cut_all, 1, -1)
P_hgcer_vs_aero_npe_pions_cut_all.Add(P_hgcer_vs_aero_npe_pions_prompt_cut_all, P_hgcer_vs_aero_npe_pions_random_cut_all, 1, -1)
P_ngcer_vs_hgcer_npe_pions_cut_all.Add(P_ngcer_vs_hgcer_npe_pions_prompt_cut_all, P_ngcer_vs_hgcer_npe_pions_random_cut_all, 1, -1)
P_ngcer_vs_aero_npe_pions_cut_all.Add(P_ngcer_vs_aero_npe_pions_prompt_cut_all, P_ngcer_vs_aero_npe_pions_random_cut_all, 1, -1)
P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_cut_all.Add(P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_prompt_cut_all, P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_random_cut_all, 1, -1)
P_aero_yAtAero_vs_aero_xAtAero_pions_cut_all.Add(P_aero_yAtAero_vs_aero_xAtAero_pions_prompt_cut_all, P_aero_yAtAero_vs_aero_xAtAero_pions_random_cut_all, 1, -1)
P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_cut_all.Add(P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_prompt_cut_all, P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_random_cut_all, 1, -1)
P_cal_etottracknorm_vs_ngcer_npe_pions_cut_all.Add(P_cal_etottracknorm_vs_ngcer_npe_pions_prompt_cut_all, P_cal_etottracknorm_vs_ngcer_npe_pions_random_cut_all, 1, -1)
CTime_ePiCoinTime_vs_MMpi_pions_cut_all.Add(CTime_ePiCoinTime_vs_MMpi_pions_prompt_cut_all, CTime_ePiCoinTime_vs_MMpi_pions_random_cut_all, 1, -1)
CTime_ePiCoinTime_vs_beta_pions_cut_all.Add(CTime_ePiCoinTime_vs_beta_pions_prompt_cut_all, CTime_ePiCoinTime_vs_beta_pions_random_cut_all, 1, -1)
P_RFTime_vs_MMpi_pions_cut_all.Add(P_RFTime_vs_MMpi_pions_prompt_cut_all, P_RFTime_vs_MMpi_pions_random_cut_all, 1, -1)
CTime_ePiCoinTime_vs_RFTime_pions_cut_all.Add(CTime_ePiCoinTime_vs_RFTime_pions_prompt_cut_all, CTime_ePiCoinTime_vs_RFTime_pions_random_cut_all, 1, -1)
P_HGC_xy_npe_pions_cut_all.Add(P_HGC_xy_npe_pions_prompt_cut_all, P_HGC_xy_npe_pions_random_cut_all, 1, -1)
P_Aero_xy_npe_pions_cut_all.Add(P_Aero_xy_npe_pions_prompt_cut_all, P_Aero_xy_npe_pions_random_cut_all, 1, -1)
P_NGC_xy_npe_pions_cut_all.Add(P_NGC_xy_npe_pions_prompt_cut_all, P_NGC_xy_npe_pions_random_cut_all, 1, -1)
P_RFTime_vs_gtr_dp_pions_cut_all.Add(P_RFTime_vs_gtr_dp_pions_prompt_cut_all, P_RFTime_vs_gtr_dp_pions_random_cut_all, 1, -1)
P_kin_MMpi_pions_data_cut_all.Add(P_kin_MMpi_pions_data_prompt_cut_all, P_kin_MMpi_pions_data_random_cut_all, 1, -1)
# PID Test
P_kin_MMpi_pions_data_aerocut_all.Add(P_kin_MMpi_pions_data_prompt_aerocut_all, P_kin_MMpi_pions_data_random_aerocut_all, 1, -1)
P_kin_MMpi_pions_data_RFcut_all.Add(P_kin_MMpi_pions_data_prompt_RFcut_all, P_kin_MMpi_pions_data_random_RFcut_all, 1, -1)

############################################################################################################################################

# HGC/NGC/Aero XY Projection vs npe for pions.
HGC_proj_yx_pions_uncut = ROOT.TProfile2D(P_HGC_xy_npe_pions_uncut.Project3DProfile("yx"))
NGC_proj_yx_pions_uncut = ROOT.TProfile2D(P_NGC_xy_npe_pions_uncut.Project3DProfile("yx"))
Aero_proj_yx_pions_uncut = ROOT.TProfile2D(P_Aero_xy_npe_pions_uncut.Project3DProfile("yx"))
HGC_proj_yx_pions_accpt_cut_all = ROOT.TProfile2D(P_HGC_xy_npe_pions_accpt_cut_all.Project3DProfile("yx"))
NGC_proj_yx_pions_accpt_cut_all = ROOT.TProfile2D(P_NGC_xy_npe_pions_accpt_cut_all.Project3DProfile("yx"))
Aero_proj_yx_pions_accpt_cut_all = ROOT.TProfile2D(P_Aero_xy_npe_pions_accpt_cut_all.Project3DProfile("yx"))
HGC_proj_yx_pions_prompt_cut_all = ROOT.TProfile2D(P_HGC_xy_npe_pions_prompt_cut_all.Project3DProfile("yx"))
NGC_proj_yx_pions_prompt_cut_all = ROOT.TProfile2D(P_NGC_xy_npe_pions_prompt_cut_all.Project3DProfile("yx"))
Aero_proj_yx_pions_prompt_cut_all = ROOT.TProfile2D(P_Aero_xy_npe_pions_prompt_cut_all.Project3DProfile("yx"))
HGC_proj_yx_pions_random_cut_all = ROOT.TProfile2D(P_HGC_xy_npe_pions_random_cut_all.Project3DProfile("yx"))
NGC_proj_yx_pions_random_cut_all = ROOT.TProfile2D(P_NGC_xy_npe_pions_random_cut_all.Project3DProfile("yx"))
Aero_proj_yx_pions_random_cut_all = ROOT.TProfile2D(P_Aero_xy_npe_pions_random_cut_all.Project3DProfile("yx"))
HGC_proj_yx_pions_cut_all = ROOT.TProfile2D(P_HGC_xy_npe_pions_cut_all.Project3DProfile("yx"))
NGC_proj_yx_pions_cut_all = ROOT.TProfile2D(P_NGC_xy_npe_pions_cut_all.Project3DProfile("yx"))
Aero_proj_yx_pions_cut_all = ROOT.TProfile2D(P_Aero_xy_npe_pions_cut_all.Project3DProfile("yx"))

############################################################################################################################################

# Removes stat box
ROOT.gStyle.SetOptStat(0)

# Saving histograms in PDF
c1_pid1 = TCanvas("c1_pid1", "MM and Detector Distributions", 100, 0, 1400, 1400)
c1_pid1.Divide(2,3)
c1_pid1.cd(1)
# Format the text string with the "p" format
#text_str = 'Beam Energy = {}, Q^2 = {}, W = {}, SHMS_theta = {}'.format(BEAM_ENERGY, Q2, W, ptheta)
text_str = '{}, {}, {}, {}'.format(BEAM_ENERGY, Q2, W, ptheta)
c1_pid1_text_lines = [
    ROOT.TText(0.5, 0.9, "Pion Physics Production Setting"),
    ROOT.TText(0.5, 0.8, text_str),
    ROOT.TText(0.5, 0.6, "PID Cuts"),
    ROOT.TText(0.5, 0.5, 'H_cer_npeSum > {}'.format(H_cer_npeSum_cut_value)),
    ROOT.TText(0.5, 0.4, 'H_cal_etottracknorm > {}'.format(H_cal_etottracknorm_cut_value)),
    ROOT.TText(0.5, 0.3, 'P_aero_npeSum > {}'.format(P_aero_npeSum_cut_value)),
    ROOT.TText(0.5, 0.2, '{} < P_RF_Dist < {}'.format(P_RF_Dist_low_cut_value, P_RF_Dist_high_cut_value)),
#    ROOT.TText(0.5, 0.1, 'P_hgcer_npeSum > {}'.format(P_hgcer_npeSum_cut_value)),

]
for c1_pid1_text in c1_pid1_text_lines:
    c1_pid1_text.SetTextSize(0.07)
    c1_pid1_text.SetTextAlign(22)
#    c1_pid1_text.SetTextColor(ROOT.kGreen + 4)
#    if c1_pid1_text.GetTitle() == "Red = SIMC":
#       c1_pid1_text.SetTextColor(ROOT.kRed)  # Setting text color to red
#    if c1_pid1_text.GetTitle() == "Blue = DATA":
#       c1_pid1_text.SetTextColor(ROOT.kBlue)  # Setting text color to red
    c1_pid1_text.Draw()
c1_pid1.cd(2)
P_kin_MMpi_pions_data_uncut.SetLineColor(1)
P_kin_MMpi_pions_data_uncut.Draw("hist")
y_min1 = 0  # Set your y-axis minimum range
y_max1 = 15000  # Set your y-axis maximum range
bin_events_h1 = 0.0
bin_events_l1 = 0.0
integral_bkregion1 = 0.0  # Initialize the manual integral
shaded_peak1 = P_kin_MMpi_pions_data_uncut.Clone()
shaded_peak1.SetFillColor(2)  # Color for shaded peak
shaded_peak1.SetFillStyle(3244)  # Fill style for shaded peak
shaded_peak1.GetXaxis().SetRangeUser(minbin, maxbin)  # Set the x-axis range
shaded_peak1.Draw("same hist")  # Draw the shaded peak
shaded_bkregion1 = P_kin_MMpi_pions_data_uncut.Clone()
shaded_bkregion1.SetFillColor(8)  # Color for shaded background
shaded_bkregion1.SetFillStyle(3244)  # Fill style for shaded background
shaded_bkregion1.GetXaxis().SetRangeUser(minbin, maxbin)  # Set the x-axis range
shaded_bkregion1.GetYaxis().SetRangeUser(y_min1, y_max1)  # Set the y-axis range
shaded_bkregion1.Draw("same hist")  # Draw the shaded background
BinLow1 = shaded_peak1.GetXaxis().FindBin(minbin)  # Find bin for minbin
BinHigh1 = shaded_peak1.GetXaxis().FindBin(maxbin)  # Find bin for maxbin
integral_peak1 = shaded_peak1.Integral(BinLow1, BinHigh1)  # Integral for shaded peak
#integral_bkregion1 = shaded_bkregion1.Integral(BinLow, BinHigh)  # Integral for shaded background region
# Loop through the bins and sum the contents if within the y range
for bin_num1 in range(BinLow1, BinHigh1 + 1):
    bin_content1 = shaded_bkregion1.GetBinContent(bin_num1)
    # Check if the bin content is within the specified y range
    if bin_content1 > y_max1:
        bin_events_h1 += y_max1
    if y_min1 <= bin_content1 <= y_max1:
        bin_events_l1 += bin_content1
integral_bkregion1 = bin_events_h1 + bin_events_l1
BinIntegral_pions1 = integral_peak1 - integral_bkregion1
NeutronEvt_pions1 = TPaveText(0.55,0.65,0.90,0.85,"NDC")
NeutronEvt_pions1.SetLineColor(2)
NeutronEvt_pions1.AddText("e #pi n Events (Peak+Bkgd): %i" %(integral_peak1))
NeutronEvt_pions1.AddText("e #pi n Events (Bkgd): %i" %(integral_bkregion1))
NeutronEvt_pions1.AddText("e #pi n Events (Peak): %i" %(BinIntegral_pions1))
NeutronEvt_pions1.Draw()
line1 = ROOT.TLine(minbin, y_max1, maxbin, y_max1)
line1.SetLineColor(4)
line1.SetLineWidth(3)
line1.Draw("same")
'''
shadedpeak_pions1 = P_kin_MMpi_pions_data_uncut.Clone()
shadedpeak_pions1.SetFillColor(2)
shadedpeak_pions1.SetFillStyle(3244)
shadedpeak_pions1.GetXaxis().SetRangeUser(minbin, maxbin)
shadedpeak_pions1.Draw("samehist")
NeutronEvt_pions1 = TPaveText(0.58934,0.675,0.95,0.75,"NDC")
BinLow_pions1 = P_kin_MMpi_pions_data_uncut.GetXaxis().FindBin(minbin)
BinHigh_pions1 = P_kin_MMpi_pions_data_uncut.GetXaxis().FindBin(maxbin)
BinIntegral_pions1 = int(P_kin_MMpi_pions_data_uncut.Integral(BinLow_pions1, BinHigh_pions1))
NeutronEvt_pions1.SetLineColor(2)
NeutronEvt_pions1.AddText("e #pi n Events: %i" %(BinIntegral_pions1))
NeutronEvt_pions1.Draw()
#legend1_pions = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
#legend1_pions.AddEntry("P_kin_MMpi_pions_data_uncut", "without cuts", "l")
#legend1_pions.AddEntry("P_kin_MMpi_pions_data_accpt_cut_all", "with cuts (acpt)", "l")
#legend1_pions.AddEntry("P_kin_MMpi_pions_data_cut_all", "with cuts (acpt/CT/PID)", "l")
#legend1_pions.Draw("same")
'''
c1_pid1.cd(3)
P_kin_MMpi_pions_data_accpt_cut_all.SetLineColor(2)
P_kin_MMpi_pions_data_accpt_cut_all.Draw("hist")
y_min2 = 0  # Set your y-axis minimum range
y_max2 = 9500  # Set your y-axis maximum range
bin_events_h2 = 0.0
bin_events_l2 = 0.0
integral_bkregion2 = 0.0  # Initialize the manual integral
shaded_peak2 = P_kin_MMpi_pions_data_accpt_cut_all.Clone()
shaded_peak2.SetFillColor(2)  # Color for shaded peak
shaded_peak2.SetFillStyle(3244)  # Fill style for shaded peak
shaded_peak2.GetXaxis().SetRangeUser(minbin, maxbin)  # Set the x-axis range
shaded_peak2.Draw("same hist")  # Draw the shaded peak
shaded_bkregion2 = P_kin_MMpi_pions_data_accpt_cut_all.Clone()
shaded_bkregion2.SetFillColor(8)  # Color for shaded background
shaded_bkregion2.SetFillStyle(3244)  # Fill style for shaded background
shaded_bkregion2.GetXaxis().SetRangeUser(minbin, maxbin)  # Set the x-axis range
shaded_bkregion2.GetYaxis().SetRangeUser(y_min2, y_max2)  # Set the y-axis range
shaded_bkregion2.Draw("same hist")  # Draw the shaded background
BinLow2 = shaded_peak2.GetXaxis().FindBin(minbin)  # Find bin for minbin
BinHigh2 = shaded_peak2.GetXaxis().FindBin(maxbin)  # Find bin for maxbin
integral_peak2 = shaded_peak2.Integral(BinLow2, BinHigh2)  # Integral for shaded peak
for bin_num2 in range(BinLow2, BinHigh2 + 1):
    bin_content2 = shaded_bkregion2.GetBinContent(bin_num2)
#    print("Bin content of shaded BK region")
#    print(bin_content2)
    # Check if the bin content is within the specified y range
    if bin_content2 > y_max2:
        bin_events_h2 += y_max2
#        print("Bin content in bin_events_h")
#        print(bin_events_h, "\n")
    if y_min2 <= bin_content2 <= y_max2:
        bin_events_l2 += bin_content2
#        print("Bin content in bin_events_l")
#        print(bin_events_l, "\n")
integral_bkregion2 = bin_events_h2 + bin_events_l2
BinIntegral_pions2 = integral_peak2 - integral_bkregion2
NeutronEvt_pions2 = TPaveText(0.55,0.65,0.90,0.85,"NDC")
NeutronEvt_pions2.SetLineColor(2)
NeutronEvt_pions2.AddText("e #pi n Events (Peak+Bkgd): %i" %(integral_peak2))
NeutronEvt_pions2.AddText("e #pi n Events (Bkgd): %i" %(integral_bkregion2))
NeutronEvt_pions2.AddText("e #pi n Events (Peak): %i" %(BinIntegral_pions2))
NeutronEvt_pions2.Draw()
line2 = ROOT.TLine(minbin, y_max2, maxbin, y_max2)
line2.SetLineColor(4)
line2.SetLineWidth(3)
line2.Draw("same")
#print("integrals")
#print(integral_peak2)
#print(integral_bkregion2)
#print(BinIntegral_pions2)
'''
shadedpeak_pions2 = P_kin_MMpi_pions_data_accpt_cut_all.Clone()
shadedpeak_pions2.SetFillColor(1)
shadedpeak_pions2.SetFillStyle(3244)
shadedpeak_pions2.GetXaxis().SetRangeUser(minbin, maxbin)
shadedpeak_pions2.Draw("samehist")
NeutronEvt_pions2 = TPaveText(0.58934,0.675,0.95,0.75,"NDC")
BinLow_pions2 = P_kin_MMpi_pions_data_accpt_cut_all.GetXaxis().FindBin(minbin)
BinHigh_pions2 = P_kin_MMpi_pions_data_accpt_cut_all.GetXaxis().FindBin(maxbin)
BinIntegral_pions2 = int(P_kin_MMpi_pions_data_accpt_cut_all.Integral(BinLow_pions2, BinHigh_pions2))
NeutronEvt_pions2.SetLineColor(2)
NeutronEvt_pions2.AddText("e #pi n Events: %i" %(BinIntegral_pions2))
NeutronEvt_pions2.Draw()
'''
c1_pid1.cd(4)
P_kin_MMpi_pions_data_aerocut_all.SetLineColor(4)
P_kin_MMpi_pions_data_aerocut_all.Draw("hist")
y_min3 = 0  # Set your y-axis minimum range
y_max3 = 400  # Set your y-axis maximum range
bin_events_h3 = 0.0
bin_events_l3 = 0.0
integral_bkregion3 = 0.0  # Initialize the manual integral
shaded_peak3 = P_kin_MMpi_pions_data_aerocut_all.Clone()
shaded_peak3.SetFillColor(2)  # Color for shaded peak
shaded_peak3.SetFillStyle(3244)  # Fill style for shaded peak
shaded_peak3.GetXaxis().SetRangeUser(minbin, maxbin)  # Set the x-axis range
shaded_peak3.Draw("same hist")  # Draw the shaded peak
shaded_bkregion3 = P_kin_MMpi_pions_data_aerocut_all.Clone()
shaded_bkregion3.SetFillColor(8)  # Color for shaded background
shaded_bkregion3.SetFillStyle(3244)  # Fill style for shaded background
shaded_bkregion3.GetXaxis().SetRangeUser(minbin, maxbin)  # Set the x-axis range
shaded_bkregion3.GetYaxis().SetRangeUser(y_min3, y_max3)  # Set the y-axis range
shaded_bkregion3.Draw("same hist")  # Draw the shaded background
BinLow3 = shaded_peak3.GetXaxis().FindBin(minbin)  # Find bin for minbin
BinHigh3 = shaded_peak3.GetXaxis().FindBin(maxbin)  # Find bin for maxbin
integral_peak3 = shaded_peak3.Integral(BinLow3, BinHigh3)  # Integral for shaded peak
#integral_bkregion1 = shaded_bkregion1.Integral(BinLow, BinHigh)  # Integral for shaded background region
# Loop through the bins and sum the contents if within the y range
for bin_num3 in range(BinLow3, BinHigh3 + 1):
    bin_content3 = shaded_bkregion3.GetBinContent(bin_num3)
    # Check if the bin content is within the specified y range
    if bin_content3 > y_max3:
        bin_events_h3 += y_max3
    if y_min3 <= bin_content3 <= y_max3:
        bin_events_l3 += bin_content3
integral_bkregion3 = bin_events_h3 + bin_events_l3
BinIntegral_pions3 = integral_peak3 - integral_bkregion3
NeutronEvt_pions3 = TPaveText(0.55,0.65,0.90,0.85,"NDC")
NeutronEvt_pions3.SetLineColor(2)
NeutronEvt_pions3.AddText("e #pi n Events (Peak+Bkgd): %i" %(integral_peak3))
NeutronEvt_pions3.AddText("e #pi n Events (Bkgd): %i" %(integral_bkregion3))
NeutronEvt_pions3.AddText("e #pi n Events (Peak): %i" %(BinIntegral_pions3))
NeutronEvt_pions3.Draw()
line3 = ROOT.TLine(minbin, y_max3, maxbin, y_max3)
line3.SetLineColor(1)
line3.SetLineWidth(3)
line3.Draw("same")
'''
shadedpeak_pions3 = P_kin_MMpi_pions_data_aerocut_all.Clone()
shadedpeak_pions3.SetFillColor(2)
shadedpeak_pions3.SetFillStyle(3244)
shadedpeak_pions3.GetXaxis().SetRangeUser(minbin, maxbin)
shadedpeak_pions3.Draw("samehist")
NeutronEvt_pions3 = TPaveText(0.58934,0.675,0.95,0.75,"NDC")
BinLow_pions3 = P_kin_MMpi_pions_data_aerocut_all.GetXaxis().FindBin(minbin)
BinHigh_pions3 = P_kin_MMpi_pions_data_aerocut_all.GetXaxis().FindBin(maxbin)
BinIntegral_pions3 = int(P_kin_MMpi_pions_data_aerocut_all.Integral(BinLow_pions3, BinHigh_pions3))
NeutronEvt_pions3.SetLineColor(2)
NeutronEvt_pions3.AddText("e #pi n Events: %i" %(BinIntegral_pions3))
NeutronEvt_pions3.Draw()
'''
c1_pid1.cd(5)
P_kin_MMpi_pions_data_RFcut_all.SetLineColor(4)
P_kin_MMpi_pions_data_RFcut_all.Draw("hist")
y_min4 = 0  # Set your y-axis minimum range
y_max4 = 400  # Set your y-axis maximum range
bin_events_h4 = 0.0
bin_events_l4 = 0.0
integral_bkregion4 = 0.0  # Initialize the manual integral
shaded_peak4 = P_kin_MMpi_pions_data_RFcut_all.Clone()
shaded_peak4.SetFillColor(2)  # Color for shaded peak
shaded_peak4.SetFillStyle(3244)  # Fill style for shaded peak
shaded_peak4.GetXaxis().SetRangeUser(minbin, maxbin)  # Set the x-axis range
shaded_peak4.Draw("same hist")  # Draw the shaded peak
shaded_bkregion4 = P_kin_MMpi_pions_data_RFcut_all.Clone()
shaded_bkregion4.SetFillColor(8)  # Color for shaded background
shaded_bkregion4.SetFillStyle(3244)  # Fill style for shaded background
shaded_bkregion4.GetXaxis().SetRangeUser(minbin, maxbin)  # Set the x-axis range
shaded_bkregion4.GetYaxis().SetRangeUser(y_min4, y_max4)  # Set the y-axis range
shaded_bkregion4.Draw("same hist")  # Draw the shaded background
BinLow4 = shaded_peak4.GetXaxis().FindBin(minbin)  # Find bin for minbin
BinHigh4 = shaded_peak4.GetXaxis().FindBin(maxbin)  # Find bin for maxbin
integral_peak4 = shaded_peak4.Integral(BinLow4, BinHigh4)  # Integral for shaded peak
#integral_bkregion1 = shaded_bkregion1.Integral(BinLow, BinHigh)  # Integral for shaded background region
# Loop through the bins and sum the contents if within the y range
for bin_num4 in range(BinLow4, BinHigh4 + 1):
    bin_content4 = shaded_bkregion4.GetBinContent(bin_num4)
    # Check if the bin content is within the specified y range
    if bin_content4 > y_max4:
        bin_events_h4 += y_max4
    if y_min4 <= bin_content4 <= y_max4:
        bin_events_l4 += bin_content4
integral_bkregion4 = bin_events_h4 + bin_events_l4
BinIntegral_pions4 = integral_peak4 - integral_bkregion4
NeutronEvt_pions4 = TPaveText(0.55,0.65,0.90,0.85,"NDC")
NeutronEvt_pions4.SetLineColor(2)
NeutronEvt_pions4.AddText("e #pi n Events (Peak+Bkgd): %i" %(integral_peak4))
NeutronEvt_pions4.AddText("e #pi n Events (Bkgd): %i" %(integral_bkregion4))
NeutronEvt_pions4.AddText("e #pi n Events (Peak): %i" %(BinIntegral_pions4))
NeutronEvt_pions4.Draw()
line4 = ROOT.TLine(minbin, y_max4, maxbin, y_max4)
line4.SetLineColor(1)
line4.SetLineWidth(3)
line4.Draw("same")
'''
shadedpeak_pions4 = P_kin_MMpi_pions_data_aerocut_all.Clone()
shadedpeak_pions4.SetFillColor(2)
shadedpeak_pions4.SetFillStyle(3244)
shadedpeak_pions4.GetXaxis().SetRangeUser(minbin, maxbin)
shadedpeak_pions4.Draw("samehist")
NeutronEvt_pions4 = TPaveText(0.58934,0.675,0.95,0.75,"NDC")
BinLow_pions4 = P_kin_MMpi_pions_data_RFcut_all.GetXaxis().FindBin(minbin)
BinHigh_pions4 = P_kin_MMpi_pions_data_RFcut_all.GetXaxis().FindBin(maxbin)
BinIntegral_pions4 = int(P_kin_MMpi_pions_data_RFcut_all.Integral(BinLow_pions4, BinHigh_pions4))
NeutronEvt_pions4.SetLineColor(2)
NeutronEvt_pions4.AddText("e #pi n Events: %i" %(BinIntegral_pions4))
NeutronEvt_pions4.Draw()
'''
c1_pid1.cd(6)
P_kin_MMpi_pions_data_cut_all.SetLineColor(4)
P_kin_MMpi_pions_data_cut_all.Draw("hist")
y_min5 = 0  # Set your y-axis minimum range
y_max5 = 400  # Set your y-axis maximum range
bin_events_h5 = 0.0
bin_events_l5 = 0.0
integral_bkregion5 = 0.0  # Initialize the manual integral
shaded_peak5 = P_kin_MMpi_pions_data_cut_all.Clone()
shaded_peak5.SetFillColor(2)  # Color for shaded peak
shaded_peak5.SetFillStyle(3244)  # Fill style for shaded peak
shaded_peak5.GetXaxis().SetRangeUser(minbin, maxbin)  # Set the x-axis range
shaded_peak5.Draw("same hist")  # Draw the shaded peak
shaded_bkregion5 = P_kin_MMpi_pions_data_cut_all.Clone()
shaded_bkregion5.SetFillColor(8)  # Color for shaded background
shaded_bkregion5.SetFillStyle(3244)  # Fill style for shaded background
shaded_bkregion5.GetXaxis().SetRangeUser(minbin, maxbin)  # Set the x-axis range
shaded_bkregion5.GetYaxis().SetRangeUser(y_min5, y_max5)  # Set the y-axis range
shaded_bkregion5.Draw("same hist")  # Draw the shaded background
BinLow5 = shaded_peak5.GetXaxis().FindBin(minbin)  # Find bin for minbin
BinHigh5 = shaded_peak5.GetXaxis().FindBin(maxbin)  # Find bin for maxbin
integral_peak5 = shaded_peak5.Integral(BinLow5, BinHigh5)  # Integral for shaded peak
#integral_bkregion1 = shaded_bkregion1.Integral(BinLow, BinHigh)  # Integral for shaded background region
# Loop through the bins and sum the contents if within the y range
for bin_num5 in range(BinLow5, BinHigh5 + 1):
    bin_content5 = shaded_bkregion5.GetBinContent(bin_num5)
    # Check if the bin content is within the specified y range
    if bin_content5 > y_max5:
        bin_events_h5 += y_max5
    if y_min5 <= bin_content5 <= y_max5:
        bin_events_l5 += bin_content5
integral_bkregion5 = bin_events_h5 + bin_events_l5
BinIntegral_pions5 = integral_peak5 - integral_bkregion5
NeutronEvt_pions5 = TPaveText(0.55,0.65,0.90,0.85,"NDC")
NeutronEvt_pions5.SetLineColor(2)
NeutronEvt_pions5.AddText("e #pi n Events (Peak+Bkgd): %i" %(integral_peak5))
NeutronEvt_pions5.AddText("e #pi n Events (Bkgd): %i" %(integral_bkregion5))
NeutronEvt_pions5.AddText("e #pi n Events (Peak): %i" %(BinIntegral_pions5))
NeutronEvt_pions5.Draw()
line5 = ROOT.TLine(minbin, y_max5, maxbin, y_max5)
line5.SetLineColor(1)
line5.SetLineWidth(3)
line5.Draw("same")
'''
shadedpeak_pions5 = P_kin_MMpi_pions_data_cut_all.Clone()
shadedpeak_pions5.SetFillColor(2)
shadedpeak_pions5.SetFillStyle(3244)
shadedpeak_pions5.GetXaxis().SetRangeUser(minbin, maxbin)
shadedpeak_pions5.Draw("samehist")
NeutronEvt_pions5 = TPaveText(0.58934,0.675,0.95,0.75,"NDC")
BinLow_pions5 = P_kin_MMpi_pions_data_cut_all.GetXaxis().FindBin(minbin)
BinHigh_pions5 = P_kin_MMpi_pions_data_cut_all.GetXaxis().FindBin(maxbin)
BinIntegral_pions5 = int(P_kin_MMpi_pions_data_cut_all.Integral(BinLow_pions5, BinHigh_pions5))
NeutronEvt_pions5.SetLineColor(2)
NeutronEvt_pions5.AddText("e #pi n Events: %i" %(BinIntegral_pions5))
NeutronEvt_pions5.Draw()
'''
c1_pid1.Print(Pion_Analysis_Distributions + '(')

c1_pid2 = TCanvas("c1_pid2", "2D Detector Distributions", 100, 0, 1400, 1400)
c1_pid2.Divide(2,3)
c1_pid2.cd(1)
gPad.SetLogy()
P_ngcer_npeSum_pions_data_uncut.SetLineColor(1)
P_ngcer_npeSum_pions_data_uncut.Draw()
c1_pid2.cd(2)
gPad.SetLogy()
P_ngcer_npeSum_pions_data_accpt_cut_all.SetLineColor(2)
P_ngcer_npeSum_pions_data_accpt_cut_all.Draw()
P_ngcer_npeSum_pions_data_cut_all.SetLineColor(4)
P_ngcer_npeSum_pions_data_cut_all.Draw("same")
legend3_pions = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
legend3_pions.AddEntry("P_ngcer_npeSum_pions_data_accpt_cut_all", "with cuts (acpt)", "l")
legend3_pions.AddEntry("P_ngcer_npeSum_pions_data_cut_all", "with cuts (acpt/CT/PID)", "l")
legend3_pions.Draw("same")
c1_pid2.cd(3)
gPad.SetLogy()
P_hgcer_npeSum_pions_data_uncut.SetLineColor(1)
P_hgcer_npeSum_pions_data_uncut.Draw()
c1_pid2.cd(4)
gPad.SetLogy()
P_hgcer_npeSum_pions_data_accpt_cut_all.SetLineColor(2)
P_hgcer_npeSum_pions_data_accpt_cut_all.Draw()
P_hgcer_npeSum_pions_data_cut_all.SetLineColor(4)
P_hgcer_npeSum_pions_data_cut_all.Draw("same")
legend4_pions = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
legend4_pions.AddEntry("P_hgcer_npeSum_pions_data_accpt_cut_all", "with cuts (acpt)", "l")
legend4_pions.AddEntry("P_hgcer_npeSum_pions_data_cut_all", "with cuts (acpt/CT/PID)", "l")
legend4_pions.Draw("same")
c1_pid2.cd(5)
gPad.SetLogy()
P_aero_npeSum_pions_data_uncut.SetLineColor(1)
P_aero_npeSum_pions_data_uncut.Draw()
c1_pid2.cd(6)
gPad.SetLogy()
P_aero_npeSum_pions_data_accpt_cut_all.SetLineColor(2)
P_aero_npeSum_pions_data_accpt_cut_all.Draw()
P_aero_npeSum_pions_data_cut_all.SetLineColor(4)
P_aero_npeSum_pions_data_cut_all.Draw("same")
legend5_pions = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
legend5_pions.AddEntry("P_aero_npeSum_pions_data_accpt_cut_all", "with cuts (acpt)", "l")
legend5_pions.AddEntry("P_aero_npeSum_pions_data_cut_all", "with cuts (acpt/CT/PID)", "l")
legend5_pions.Draw("same")
c1_pid2.Print(Pion_Analysis_Distributions)

c1_pid3 = TCanvas("c1_pid3", "2D Detector Distributions", 100, 0, 1400, 1400)
c1_pid3.Divide(2,3)
c1_pid3.cd(1)
gPad.SetLogz()
P_hgcer_vs_aero_npe_pions_uncut.Draw("COLZ")
c1_pid3.cd(2)
gPad.SetLogz()
P_hgcer_vs_aero_npe_pions_cut_all.Draw("COLZ")
c1_pid3.cd(3)
gPad.SetLogz()
P_ngcer_vs_hgcer_npe_pions_uncut.Draw("COLZ")
c1_pid3.cd(4)
gPad.SetLogz()
P_ngcer_vs_hgcer_npe_pions_cut_all.Draw("COLZ")
c1_pid3.cd(5)
gPad.SetLogz()
P_ngcer_vs_aero_npe_pions_uncut.Draw("COLZ")
c1_pid3.cd(6)
gPad.SetLogz()
P_ngcer_vs_aero_npe_pions_cut_all.Draw("COLZ")
c1_pid3.Print(Pion_Analysis_Distributions)

c1_pid4 = TCanvas("c1_pid4", "CT Distributions", 100, 0, 1400, 1400)
c1_pid4.Divide(2,3)
c1_pid4.cd(1)
CTime_ePiCoinTime_ROC1_pions_data_uncut.SetLineColor(1)
CTime_ePiCoinTime_ROC1_pions_data_uncut.Draw()
c1_pid4.cd(2)
#CTime_ePiCoinTime_ROC1_pions_data_accpt_cut_all.SetLineColor(4)
#CTime_ePiCoinTime_ROC1_pions_data_accpt_cut_all.Draw()
CTime_ePiCoinTime_ROC1_pions_data_pid_cut_all.SetLineColor(4)
CTime_ePiCoinTime_ROC1_pions_data_pid_cut_all.Draw()
# Find the x position of maxima
max_bin_CT = CTime_ePiCoinTime_ROC1_pions_data_pid_cut_all.GetMaximumBin()
x_max_CT = CTime_ePiCoinTime_ROC1_pions_data_pid_cut_all.GetBinCenter(max_bin_CT)
print(f"X position of CT maximum: {x_max_CT}")
CTime_ePiCoinTime_ROC1_pions_data_prompt_cut_all.SetLineColor(6)
CTime_ePiCoinTime_ROC1_pions_data_prompt_cut_all.Draw("same")
CTime_ePiCoinTime_ROC1_pions_data_random_unsub_cut_all.SetLineColor(8)
CTime_ePiCoinTime_ROC1_pions_data_random_unsub_cut_all.Draw("same")
legend6_pions = ROOT.TLegend(0.1, 0.815, 0.48, 0.9)
legend6_pions.AddEntry("ePiCoinTime_pions_pidcut", "CT with acpt+PID cuts", "l")
legend6_pions.AddEntry("ePiCoinTime_pions_cut_prompt", "CT_prompt with cuts (acpt/PID)", "l")
legend6_pions.AddEntry("ePiCoinTime_pions_cut_randm", "CT_randoms with cuts (acpt/PID)", "l")
legend6_pions.Draw("same")
c1_pid4.cd(3)
gPad.SetLogz()
CTime_ePiCoinTime_vs_MMpi_pions_uncut.Draw("COLZ")
LowerPrompt1_pions = TLine(PromptWindow[0],gPad.GetUymin(),PromptWindow[0],2)
LowerPrompt1_pions.SetLineColor(2)
LowerPrompt1_pions.SetLineWidth(2)
LowerPrompt1_pions.Draw("same")
UpperPrompt1_pions = TLine(PromptWindow[1],gPad.GetUymin(),PromptWindow[1],2)
UpperPrompt1_pions.SetLineColor(2)
UpperPrompt1_pions.SetLineWidth(2)
UpperPrompt1_pions.Draw("same")
LowerRandomL1_pions = TLine(RandomWindows[0],gPad.GetUymin(),RandomWindows[0],2)
LowerRandomL1_pions.SetLineColor(8)
LowerRandomL1_pions.SetLineWidth(2)
LowerRandomL1_pions.Draw("same")
UpperRandomL1_pions = TLine(RandomWindows[1],gPad.GetUymin(),RandomWindows[1],2)
UpperRandomL1_pions.SetLineColor(8)
UpperRandomL1_pions.SetLineWidth(2)
UpperRandomL1_pions.Draw("same")
LowerRandomR1_pions = TLine(RandomWindows[2],gPad.GetUymin(),RandomWindows[2],2)
LowerRandomR1_pions.SetLineColor(8)
LowerRandomR1_pions.SetLineWidth(2)
LowerRandomR1_pions.Draw("same")
UpperRandomR1_pions = TLine(RandomWindows[3],gPad.GetUymin(),RandomWindows[3],2)
UpperRandomR1_pions.SetLineColor(8)
UpperRandomR1_pions.SetLineWidth(2)
UpperRandomR1_pions.Draw("same")
c1_pid4.cd(4)
gPad.SetLogz()
CTime_ePiCoinTime_vs_MMpi_pions_cut_all.GetYaxis().SetRangeUser(0.4, 1.4)
CTime_ePiCoinTime_vs_MMpi_pions_cut_all.GetXaxis().SetRangeUser(-5, 5)
CTime_ePiCoinTime_vs_MMpi_pions_cut_all.Draw("COLZ")
c1_pid4.cd(5)
P_RFTime_Dist_pions_data_accpt_cut_all.SetLineColor(2)
P_RFTime_Dist_pions_data_accpt_cut_all.Draw()
P_RFTime_Dist_pions_data_cut_all.SetLineColor(4)
P_RFTime_Dist_pions_data_cut_all.Draw("same")
# Find the x position of maxima
max_bin_RF = P_RFTime_Dist_pions_data_cut_all.GetMaximumBin()
x_max_RF = P_RFTime_Dist_pions_data_cut_all.GetBinCenter(max_bin_RF)
print(f"X position of RF maximum: {x_max_RF}")
legend7_pions = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
legend7_pions.AddEntry("P_RFTime_Dist_pions_data_accpt_cut_all", "with cuts (acpt)", "l")
legend7_pions.AddEntry("P_RFTime_Dist_pions_data_cut_all", "with cuts (acpt/CT/PID)", "l")
legend7_pions.Draw("same")
c1_pid4.cd(6)
gPad.SetLogz()
#CTime_ePiCoinTime_vs_RFTime_pions_cut_all.GetXaxis().SetRangeUser(-5, 5)
#CTime_ePiCoinTime_vs_RFTime_pions_cut_all.Draw("COLZ")
P_RFTime_vs_gtr_dp_pions_cut_all.GetYaxis().SetRangeUser(-15, 25)
P_RFTime_vs_gtr_dp_pions_cut_all.Draw("COLZ")
c1_pid4.Print(Pion_Analysis_Distributions)

c1_pid5 = TCanvas("c1_pid5", "RF Distributions", 100, 0, 1400, 1400)
c1_pid5.Divide(2,3)
c1_pid5.cd(1)
gPad.SetLogz()
P_RFTime_vs_MMpi_pions_uncut.Draw("COLZ")
c1_pid5.cd(2)
gPad.SetLogz()
P_RFTime_vs_MMpi_pions_cut_all.GetYaxis().SetRangeUser(0.4, 1.4)
P_RFTime_vs_MMpi_pions_cut_all.Draw("COLZ")
c1_pid5.cd(3)
P_cal_etottracknorm_pions_data_uncut.SetLineColor(2)
P_cal_etottracknorm_pions_data_uncut.Draw()
P_cal_etottracknorm_pions_data_cut_all.SetLineColor(4)
P_cal_etottracknorm_pions_data_cut_all.Draw("same")
legend8_pions = ROOT.TLegend(0.1, 0.815, 0.48, 0.9)
legend8_pions.AddEntry("P_cal_etottracknorm_pions_data_uncut", "without cuts", "l")
legend8_pions.AddEntry("P_cal_etottracknorm_pions_data_cut_all", "with cuts (acpt/CT/PID)", "l")
legend8_pions.Draw("same")
c1_pid5.cd(4)
gPad.SetLogy()
H_cal_etottracknorm_pions_data_uncut.SetLineColor(2)
H_cal_etottracknorm_pions_data_uncut.Draw()
H_cal_etottracknorm_pions_data_cut_all.SetLineColor(4)
H_cal_etottracknorm_pions_data_cut_all.Draw("same")
legend9_pions = ROOT.TLegend(0.1, 0.815, 0.48, 0.9)
legend9_pions.AddEntry("H_cal_etottracknorm_pions_data_uncut", "without cuts", "l")
legend9_pions.AddEntry("H_cal_etottracknorm_pions_data_cut_all", "with cuts (acpt/CT/PID)", "l")
legend9_pions.Draw("same")
c1_pid5.cd(5)
gPad.SetLogy()
H_cer_npeSum_pions_data_uncut.SetLineColor(2)
H_cer_npeSum_pions_data_uncut.Draw()
H_cer_npeSum_pions_data_cut_all.SetLineColor(4)
H_cer_npeSum_pions_data_cut_all.Draw("same")
legend11_pions = ROOT.TLegend(0.1, 0.815, 0.48, 0.9)
legend11_pions.AddEntry("H_cer_npeSum_pions_data_uncut", "without cuts", "l")
legend11_pions.AddEntry("H_cer_npeSum_pions_data_cut_all", "with cuts (acpt/CT/PID)", "l")
legend11_pions.Draw("same")
c1_pid5.cd(6)
gPad.SetLogz()
H_cal_etottracknorm_vs_cer_npeSum_pions_random_cut_all.Draw("COLZ")
c1_pid5.Print(Pion_Analysis_Distributions)

c1_pid6 = TCanvas("c1_pid6", "Detector XY Distributions", 100, 0, 1400, 1400)
c1_pid6.Divide(2,3)
c1_pid6.cd(1)
gPad.SetLogz()
HGC_proj_yx_pions_uncut.Draw("COLZ")
c1_pid6.cd(2)
gPad.SetLogz()
HGC_proj_yx_pions_cut_all.Draw("COLZ")
c1_pid6.cd(3)
gPad.SetLogz()
NGC_proj_yx_pions_uncut.Draw("COLZ")
c1_pid6.cd(4)
gPad.SetLogz()
NGC_proj_yx_pions_cut_all.Draw("COLZ")
c1_pid6.cd(5)
gPad.SetLogz()
Aero_proj_yx_pions_uncut.Draw("COLZ")
c1_pid6.cd(6)
gPad.SetLogz()
Aero_proj_yx_pions_cut_all.Draw("COLZ")
c1_pid6.Print(Pion_Analysis_Distributions)

c1_pid7 = TCanvas("c1_pid7", "MMPi Distributions", 100, 0, 1400, 1400)
c1_pid7.Divide(1,1)
c1_pid7.cd(1)
P_kin_MMpi_pions_data_uncut.SetLineColor(1)
P_kin_MMpi_pions_data_uncut.Draw("hist")
P_kin_MMpi_pions_data_accpt_cut_all.SetLineColor(4)
P_kin_MMpi_pions_data_accpt_cut_all.Draw("same hist")
P_kin_MMpi_pions_data_aerocut_all.SetLineColor(6)
P_kin_MMpi_pions_data_aerocut_all.Draw("same hist")
P_kin_MMpi_pions_data_RFcut_all.SetLineColor(8)
P_kin_MMpi_pions_data_RFcut_all.Draw("same hist")
P_kin_MMpi_pions_data_cut_all.SetLineColor(2)
P_kin_MMpi_pions_data_cut_all.Draw("same hist")
#P_kin_MMpi_pions_data_prompt_cut_all.Draw("same hist")
#P_kin_MMpi_pions_data_random_cut_all.Draw("same hist")
legend7_pions = ROOT.TLegend(0.55,0.65,0.90,0.85)
legend7_pions.AddEntry("P_kin_MMpi_pions_data_uncut", "MMpi without cuts", "l")
legend7_pions.AddEntry("P_kin_MMpi_pions_data_accpt_cut_all", "MMpi with acpt cuts only", "l")
legend7_pions.AddEntry("P_kin_MMpi_pions_data_aerocut_all", "MMpi with acpt+Aero+HMSPID cuts", "l")
legend7_pions.AddEntry("P_kin_MMpi_pions_data_RFcut_all", "MMpi with acpt+RF+HMSPID cuts", "l")
legend7_pions.AddEntry("P_kin_MMpi_pions_data_cut_all", "MMpi with acpt+Aero+RF+HMSPID cuts", "l")
legend7_pions.Draw("same")

c1_pid7.Print(Pion_Analysis_Distributions + ')')


#############################################################################################################################################

# Making directories in output file
outHistFile = ROOT.TFile.Open("%s/%s_%s_%s_%s_%s_ProdCoin_PID_Output_Data.root" % (OUTPATH, BEAM_ENERGY, Q2, W, ptheta, MaxEvent) , "RECREATE")
d_Uncut_Pion_Events_Data = outHistFile.mkdir("Uncut_Pion_Events_Data")
d_Cut_Pion_Events_Accpt_Data = outHistFile.mkdir("Cut_Pion_Events_Accpt_Data")
d_Cut_Pion_Events_Prompt_Data = outHistFile.mkdir("Cut_Pion_Events_Prompt_Data")
d_Cut_Pion_Events_Random_Data = outHistFile.mkdir("Cut_Pion_Events_Random_Data")
d_Cut_Pion_Events_All_Data = outHistFile.mkdir("Cut_Pion_Events_All_Data")

# Writing Histograms for pions
d_Uncut_Pion_Events_Data.cd()
H_gtr_beta_pions_data_uncut.Write()
H_cal_etottracknorm_pions_data_uncut.Write()
H_cer_npeSum_pions_data_uncut.Write()
H_RFTime_Dist_pions_data_uncut.Write()
P_gtr_beta_pions_data_uncut.Write()
P_gtr_dp_pions_data_uncut.Write()
P_cal_etottracknorm_pions_data_uncut.Write()
P_hgcer_npeSum_pions_data_uncut.Write()
P_hgcer_xAtCer_pions_data_uncut.Write()
P_hgcer_yAtCer_pions_data_uncut.Write()
P_ngcer_npeSum_pions_data_uncut.Write()
P_ngcer_xAtCer_pions_data_uncut.Write()
P_ngcer_yAtCer_pions_data_uncut.Write()
P_aero_npeSum_pions_data_uncut.Write()
P_aero_xAtAero_pions_data_uncut.Write()
P_aero_yAtAero_pions_data_uncut.Write()
P_kin_MMpi_pions_data_uncut.Write()
P_RFTime_Dist_pions_data_uncut.Write()
CTime_ePiCoinTime_ROC1_pions_data_uncut.Write()
H_cal_etottracknorm_vs_cer_npeSum_pions_uncut.Write()
P_hgcer_vs_aero_npe_pions_uncut.Write()
P_ngcer_vs_hgcer_npe_pions_uncut.Write()
P_ngcer_vs_aero_npe_pions_uncut.Write()
P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_uncut.Write()
P_aero_yAtAero_vs_aero_xAtAero_pions_uncut.Write()
P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_uncut.Write()
P_cal_etottracknorm_vs_ngcer_npe_pions_uncut.Write()
CTime_ePiCoinTime_vs_MMpi_pions_uncut.Write()
CTime_ePiCoinTime_vs_beta_pions_uncut.Write()
P_RFTime_vs_MMpi_pions_uncut.Write()
CTime_ePiCoinTime_vs_RFTime_pions_uncut.Write()
P_HGC_xy_npe_pions_uncut.Write()
P_Aero_xy_npe_pions_uncut.Write()
P_NGC_xy_npe_pions_uncut.Write()
P_RFTime_vs_gtr_dp_pions_uncut.Write()

d_Cut_Pion_Events_Accpt_Data.cd()
H_gtr_beta_pions_data_accpt_cut_all.Write()
H_cal_etottracknorm_pions_data_accpt_cut_all.Write()
H_cer_npeSum_pions_data_accpt_cut_all.Write()
H_RFTime_Dist_pions_data_accpt_cut_all.Write()
P_gtr_beta_pions_data_accpt_cut_all.Write()
P_gtr_dp_pions_data_accpt_cut_all.Write()
P_cal_etottracknorm_pions_data_accpt_cut_all.Write()
P_hgcer_npeSum_pions_data_accpt_cut_all.Write()
P_hgcer_xAtCer_pions_data_accpt_cut_all.Write()
P_hgcer_yAtCer_pions_data_accpt_cut_all.Write()
P_ngcer_npeSum_pions_data_accpt_cut_all.Write()
P_ngcer_xAtCer_pions_data_accpt_cut_all.Write()
P_ngcer_yAtCer_pions_data_accpt_cut_all.Write()
P_aero_npeSum_pions_data_accpt_cut_all.Write()
P_aero_xAtAero_pions_data_accpt_cut_all.Write()
P_aero_yAtAero_pions_data_accpt_cut_all.Write()
P_kin_MMpi_pions_data_accpt_cut_all.Write()
P_RFTime_Dist_pions_data_accpt_cut_all.Write()
CTime_ePiCoinTime_ROC1_pions_data_accpt_cut_all.Write()
H_cal_etottracknorm_vs_cer_npeSum_pions_accpt_cut_all.Write()
P_hgcer_vs_aero_npe_pions_accpt_cut_all.Write()
P_ngcer_vs_hgcer_npe_pions_accpt_cut_all.Write()
P_ngcer_vs_aero_npe_pions_accpt_cut_all.Write()
P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_accpt_cut_all.Write()
P_aero_yAtAero_vs_aero_xAtAero_pions_accpt_cut_all.Write()
P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_accpt_cut_all.Write()
P_cal_etottracknorm_vs_ngcer_npe_pions_accpt_cut_all.Write()
CTime_ePiCoinTime_vs_MMpi_pions_accpt_cut_all.Write()
CTime_ePiCoinTime_vs_beta_pions_accpt_cut_all.Write()
P_RFTime_vs_MMpi_pions_accpt_cut_all.Write()
CTime_ePiCoinTime_vs_RFTime_pions_accpt_cut_all.Write()
P_HGC_xy_npe_pions_accpt_cut_all.Write()
P_Aero_xy_npe_pions_accpt_cut_all.Write()
P_NGC_xy_npe_pions_accpt_cut_all.Write()
P_RFTime_vs_gtr_dp_pions_accpt_cut_all.Write()

d_Cut_Pion_Events_Prompt_Data.cd()
H_gtr_beta_pions_data_prompt_cut_all.Write()
H_cal_etottracknorm_pions_data_prompt_cut_all.Write()
H_cer_npeSum_pions_data_prompt_cut_all.Write()
H_RFTime_Dist_pions_data_prompt_cut_all.Write()
P_gtr_beta_pions_data_prompt_cut_all.Write()
P_gtr_dp_pions_data_prompt_cut_all.Write()
P_cal_etottracknorm_pions_data_prompt_cut_all.Write()
P_hgcer_npeSum_pions_data_prompt_cut_all.Write()
P_hgcer_xAtCer_pions_data_prompt_cut_all.Write()
P_hgcer_yAtCer_pions_data_prompt_cut_all.Write()
P_ngcer_npeSum_pions_data_prompt_cut_all.Write()
P_ngcer_xAtCer_pions_data_prompt_cut_all.Write()
P_ngcer_yAtCer_pions_data_prompt_cut_all.Write()
P_aero_npeSum_pions_data_prompt_cut_all.Write()
P_aero_xAtAero_pions_data_prompt_cut_all.Write()
P_aero_yAtAero_pions_data_prompt_cut_all.Write()
P_kin_MMpi_pions_data_prompt_cut_all.Write()
P_RFTime_Dist_pions_data_prompt_cut_all.Write()
CTime_ePiCoinTime_ROC1_pions_data_prompt_cut_all.Write()
H_cal_etottracknorm_vs_cer_npeSum_pions_prompt_cut_all.Write()
P_hgcer_vs_aero_npe_pions_prompt_cut_all.Write()
P_ngcer_vs_hgcer_npe_pions_prompt_cut_all.Write()
P_ngcer_vs_aero_npe_pions_prompt_cut_all.Write()
P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_prompt_cut_all.Write()
P_aero_yAtAero_vs_aero_xAtAero_pions_prompt_cut_all.Write()
P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_prompt_cut_all.Write()
P_cal_etottracknorm_vs_ngcer_npe_pions_prompt_cut_all.Write()
CTime_ePiCoinTime_vs_MMpi_pions_prompt_cut_all.Write()
CTime_ePiCoinTime_vs_beta_pions_prompt_cut_all.Write()
P_RFTime_vs_MMpi_pions_prompt_cut_all.Write()
CTime_ePiCoinTime_vs_RFTime_pions_prompt_cut_all.Write()
P_HGC_xy_npe_pions_prompt_cut_all.Write()
P_Aero_xy_npe_pions_prompt_cut_all.Write()
P_NGC_xy_npe_pions_prompt_cut_all.Write()
P_RFTime_vs_gtr_dp_pions_prompt_cut_all.Write()

d_Cut_Pion_Events_Random_Data.cd()
H_gtr_beta_pions_data_random_cut_all.Write()
H_cal_etottracknorm_pions_data_random_cut_all.Write()
H_cer_npeSum_pions_data_random_cut_all.Write()
H_RFTime_Dist_pions_data_random_cut_all.Write()
P_gtr_beta_pions_data_random_cut_all.Write()
P_gtr_dp_pions_data_random_cut_all.Write()
P_cal_etottracknorm_pions_data_random_cut_all.Write()
P_hgcer_npeSum_pions_data_random_cut_all.Write()
P_hgcer_xAtCer_pions_data_random_cut_all.Write()
P_hgcer_yAtCer_pions_data_random_cut_all.Write()
P_ngcer_npeSum_pions_data_random_cut_all.Write()
P_ngcer_xAtCer_pions_data_random_cut_all.Write()
P_ngcer_yAtCer_pions_data_random_cut_all.Write()
P_aero_npeSum_pions_data_random_cut_all.Write()
P_aero_xAtAero_pions_data_random_cut_all.Write()
P_aero_yAtAero_pions_data_random_cut_all.Write()
P_kin_MMpi_pions_data_random_cut_all.Write()
P_RFTime_Dist_pions_data_random_cut_all.Write()
CTime_ePiCoinTime_ROC1_pions_data_random_cut_all.Write()
H_cal_etottracknorm_vs_cer_npeSum_pions_random_cut_all.Write()
P_hgcer_vs_aero_npe_pions_random_cut_all.Write()
P_ngcer_vs_hgcer_npe_pions_random_cut_all.Write()
P_ngcer_vs_aero_npe_pions_random_cut_all.Write()
P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_random_cut_all.Write()
P_aero_yAtAero_vs_aero_xAtAero_pions_random_cut_all.Write()
P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_random_cut_all.Write()
P_cal_etottracknorm_vs_ngcer_npe_pions_random_cut_all.Write()
CTime_ePiCoinTime_vs_MMpi_pions_random_cut_all.Write()
CTime_ePiCoinTime_vs_beta_pions_random_cut_all.Write()
P_RFTime_vs_MMpi_pions_random_cut_all.Write()
CTime_ePiCoinTime_vs_RFTime_pions_random_cut_all.Write()
P_HGC_xy_npe_pions_random_cut_all.Write()
P_Aero_xy_npe_pions_random_cut_all.Write()
P_NGC_xy_npe_pions_random_cut_all.Write()
P_RFTime_vs_gtr_dp_pions_random_cut_all.Write()

d_Cut_Pion_Events_All_Data.cd()
H_gtr_beta_pions_data_cut_all.Write()
H_cal_etottracknorm_pions_data_cut_all.Write()
H_cer_npeSum_pions_data_cut_all.Write()
H_RFTime_Dist_pions_data_cut_all.Write()
P_gtr_beta_pions_data_cut_all.Write()
P_gtr_dp_pions_data_cut_all.Write()
P_cal_etottracknorm_pions_data_cut_all.Write()
P_hgcer_npeSum_pions_data_cut_all.Write()
P_hgcer_xAtCer_pions_data_cut_all.Write()
P_hgcer_yAtCer_pions_data_cut_all.Write()
P_ngcer_npeSum_pions_data_cut_all.Write()
P_ngcer_xAtCer_pions_data_cut_all.Write()
P_ngcer_yAtCer_pions_data_cut_all.Write()
P_aero_npeSum_pions_data_cut_all.Write()
P_aero_xAtAero_pions_data_cut_all.Write()
P_aero_yAtAero_pions_data_cut_all.Write()
P_kin_MMpi_pions_data_cut_all.Write()
P_RFTime_Dist_pions_data_cut_all.Write()
CTime_ePiCoinTime_ROC1_pions_data_cut_all.Write()
H_cal_etottracknorm_vs_cer_npeSum_pions_cut_all.Write()
P_hgcer_vs_aero_npe_pions_cut_all.Write()
P_ngcer_vs_hgcer_npe_pions_cut_all.Write()
P_ngcer_vs_aero_npe_pions_cut_all.Write()
P_hgcer_yAtCer_vs_hgcer_xAtCer_pions_cut_all.Write()
P_aero_yAtAero_vs_aero_xAtAero_pions_cut_all.Write()
P_ngcer_yAtCer_vs_ngcer_xAtCer_pions_cut_all.Write()
P_cal_etottracknorm_vs_ngcer_npe_pions_cut_all.Write()
CTime_ePiCoinTime_vs_MMpi_pions_cut_all.Write()
CTime_ePiCoinTime_vs_beta_pions_cut_all.Write()
P_RFTime_vs_MMpi_pions_cut_all.Write()
CTime_ePiCoinTime_vs_RFTime_pions_cut_all.Write()
P_HGC_xy_npe_pions_cut_all.Write()
P_Aero_xy_npe_pions_cut_all.Write()
P_NGC_xy_npe_pions_cut_all.Write()
P_RFTime_vs_gtr_dp_pions_cut_all.Write()
P_kin_MMpi_pions_data_aerocut_all.Write()
P_kin_MMpi_pions_data_RFcut_all.Write()
CTime_ePiCoinTime_ROC1_pions_data_pid_cut_all.Write()

##################################################################################################################################################################################3

infile_DATA.Close() 
outHistFile.Close()

print ("Processing Complete")
