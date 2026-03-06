#! /usr/bin/python
#
# Description:
# ================================================================
# Time-stamp: "2025-03-13 01:29:19 junaid"
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
from ROOT import TCanvas, TList, TPaveLabel, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TLegend, TGaxis, TLine, TMath, TLatex, TPaveText, TArc, TGraph, TGraphPolar, TText, TString
from ROOT import kBlack, kCyan, kRed, kGreen, kMagenta, kBlue
from functools import reduce
import math as ma
import ctypes

##################################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=1:
    print("!!!!! ERROR !!!!!\n Expected 1 arguments\n Usage is with - PHY_SETTING MaxEvents Suffix RunList\n!!!!! ERROR !!!!!")
    sys.exit(1)
# Extract the first three words from PHY_SETTING for the CSV file name
##################################################################################################################################################

# Input params - run number and max number of events
PHY_SETTING = sys.argv[1]
# Extract the first three words from PHY_SETTING for the CSV file name
setting_name = "_".join(PHY_SETTING.split("_")[:4])

# Constructing inout file names
DATA_ROOTFILE_SUFFIX = "%s_-1_ProdCoin_Analysed_Data" % (setting_name)
DUMMY_ROOTFILE_SUFFIX = "%s_-1_ProdCoin_Analysed_Dummy_Data" % (setting_name)

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
from ltsep import Root

lt=Root(os.path.realpath(__file__), "Plot_ProdCoin")

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

# Input file location and variables taking
rootFile_DATA = "%s/%s.root" % (OUTPATH, DATA_ROOTFILE_SUFFIX)
rootFile_DUMMY = "%s/%s.root" % (OUTPATH, DUMMY_ROOTFILE_SUFFIX)

###################################################################################################################################################
##########################################################################################################################################################################################################

# Read stuff from the main event tree
infile_DATA = ROOT.TFile.Open(rootFile_DATA, "READ")
infile_DUMMY = ROOT.TFile.Open(rootFile_DUMMY, "READ")

# Grab the trees
Cut_Pion_Events_Accpt_Data_tree = infile_DATA.Get("Cut_Pion_Events_Accpt")
Cut_Pion_Events_Accpt_Dummy_tree = infile_DUMMY.Get("Cut_Pion_Events_Accpt")
Cut_Pion_Events_Prompt_Data_tree = infile_DATA.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Prompt_Dummy_tree = infile_DUMMY.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Random_Data_tree = infile_DATA.Get("Cut_Pion_Events_Random")
Cut_Pion_Events_Random_Dummy_tree = infile_DUMMY.Get("Cut_Pion_Events_Random")

###################################################################################################################################################

# Grab runlist
runlist_file = "%s/UTIL_BATCH/InputRunLists/PionLT_2021_2022/%s_center" % (REPLAYPATH, setting_name)

# Function to read run numbers from a given run list file
def read_run_numbers(run_list_file):
    try:
        with open(run_list_file, 'r') as f:
            return [int(line.strip()) for line in f if line.strip().isdigit()]
    except FileNotFoundError:
        print(f"!!!!! ERROR !!!!!\nRun list file not found: {run_list_file}\n!!!!! ERROR !!!!!")
        sys.exit(2)

# Read run numbers from the specified files
runlist = read_run_numbers(runlist_file)
# assign first run number
run_number = runlist[0]

# Section for grabbing Prompt/Random selection parameters from PARAM file
PARAMPATH = UTILPATH+"/DB/PARAM/PionLT"
TimingCutFile = "%s/PionLT_Timing_Parameters.csv" % PARAMPATH # This should match the param file actually being used!

TimingCutf = open(TimingCutFile)
try:
    TimingCutFile
except NameError:
    print("!!!!! ERRROR !!!!!\n One (or more) of the cut files not found!\n!!!!! ERRORR !!!!!")
    sys.exit(2)
print("\n Reading timing cuts from %s to select random windows \n" % TimingCutFile)

PromptWindow = [0, 0]
RandomWindows = [0, 0, 0, 0]
linenum = 0 # Count line number we're on
TempPar = -1 # To check later
for line in TimingCutf: # Read all lines in the cut file
    linenum += 1 # Add one to line number at start of loop
    if(linenum > 1): # Skip first line
        line = line.partition('#')[0] # Treat anything after a # as a comment and ignore it
        line = line.rstrip()
        array_list = line.split(",") # Convert line into an array_list, anything after a comma is a new entry

        # Validate entries and check if run_number falls in the range [start, end]
        try:
            start = int(array_list[0].strip())
            end = int(array_list[1].strip())
        except (IndexError, ValueError):
            continue

        if run_number in range(start, end + 1):
            TempPar += 2  # If run number is in range, set to non -1 value
            BunchSpacing = float(array_list[2])
            CoinOffset = float(array_list[3])  # Coin offset value
            nSkip = float(array_list[4])  # Number of random windows skipped
            nWindows = float(array_list[5])  # Total number of random windows
            PromptPeak = float(array_list[6])  # Pion CT prompt peak position
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

print("Using the following timing windows (in ns):")
print("Prompt Window: [%.2f, %.2f]" % (PromptWindow[0], PromptWindow[1]))
print("Number of Random Windows: %.2f" % nWindows)
print("Random Windows: [%.2f, %.2f] and [%.2f, %.2f]" % (RandomWindows[0], RandomWindows[1], RandomWindows[2], RandomWindows[3]))
TimingCutf.close()

CoinTime_Prompt_Cut = lambda event: (PromptWindow[0] <= event.CTime_ePiCoinTime_ROC2 <= PromptWindow[1])
CoinTime_Random_Cut = lambda event: ((RandomWindows[0] <= event.CTime_ePiCoinTime_ROC2 <= RandomWindows[1]) and (RandomWindows[2] <= event.CTime_ePiCoinTime_ROC2 <= RandomWindows[3]))

#############################################################################################################################################################################################################

# HMS Detector Cuts
HMS_Calo_cut = 0.7
HMS_Cer_cut = 1.4
HMS_Detector_Cut = lambda event: (event.H_cal_etottracknorm > HMS_Calo_cut and event.H_cer_npeSum > HMS_Cer_cut)

# RF cut
RF_Cut_lowvalue = 1.2
RF_Cut_highvalue = 3.4
RF_Cut = lambda event: (RF_Cut_lowvalue <= event.P_RF_Dist <= RF_Cut_highvalue)

# SHMS Detector Cuts
SHMS_Aero_cut = 1.5
Aero_Detector_Cut = lambda event: (event.P_aero_npeSum > SHMS_Aero_cut)   

#############################################################################################################################################################################################################

# Defining Histogram for PID study
MM_nbins = 280
MM_min = 0.0
MM_max = 1.4
CT_nbins = 1000
CT_min = -50.0
CT_max = 50.0
Aero_nbins = 500
Aero_min = 0.0
Aero_max = 50.0

P_kin_MMpi_pions_data_prompt_noRF_noAero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_prompt_noRF_noAero_cut_all", "MIssing Mass data (NoRF NoAero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_data_prompt_noRF_noAero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_data_prompt_noRF_noAero_cut_all", "Electron-Pion CTime (NoRF NoAero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_data_prompt_noRF_noAero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_prompt_noRF_noAero_cut_all", "SHMS aero npeSum (NoRF NoAero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_data_random_noRF_noAero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_random_noRF_noAero_cut_all", "MIssing Mass data (NoRF NoAero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_data_random_noRF_noAero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_data_random_noRF_noAero_cut_all", "Electron-Pion CTime (NoRF NoAero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_data_random_noRF_noAero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_random_noRF_noAero_cut_all", "SHMS aero npeSum (NoRF NoAero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_data_randomsub_noRF_noAero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_randomsub_noRF_noAero_cut_all", "MIssing Mass data (NoRF NoAero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_data_randomsub_noRF_noAero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_data_randomsub_noRF_noAero_cut_all", "Electron-Pion CTime (NoRF NoAero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_data_randomsub_noRF_noAero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_randomsub_noRF_noAero_cut_all", "SHMS aero npeSum (NoRF NoAero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_dummy_prompt_noRF_noAero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_dummy_prompt_noRF_noAero_cut_all", "MIssing Mass data (NoRF NoAero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_dummy_prompt_noRF_noAero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_dummy_prompt_noRF_noAero_cut_all", "Electron-Pion CTime (NoRF NoAero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_dummy_prompt_noRF_noAero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_dummy_prompt_noRF_noAero_cut_all", "SHMS aero npeSum (NoRF NoAero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_dummy_random_noRF_noAero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_dummy_random_noRF_noAero_cut_all", "MIssing Mass data (NoRF NoAero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_dummy_random_noRF_noAero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_dummy_random_noRF_noAero_cut_all", "Electron-Pion CTime (NoRF NoAero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_dummy_random_noRF_noAero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_dummy_random_noRF_noAero_cut_all", "SHMS aero npeSum (NoRF NoAero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_dummy_randomsub_noRF_noAero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_dummy_randomsub_noRF_noAero_cut_all", "MIssing Mass data (NoRF NoAero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_dummy_randomsub_noRF_noAero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_dummy_randomsub_noRF_noAero_cut_all", "Electron-Pion CTime (NoRF NoAero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_dummy_randomsub_noRF_noAero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_dummy_randomsub_noRF_noAero_cut_all", "SHMS aero npeSum (NoRF NoAero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_data_dummysub_noRF_noAero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_dummysub_noRF_noAero_cut_all", "MIssing Mass data (NoRF NoAero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_data_dummysub_noRF_noAero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_data_dummysub_noRF_noAero_cut_all", "Electron-Pion CTime (NoRF NoAero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_data_dummysub_noRF_noAero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_dummysub_noRF_noAero_cut_all", "SHMS aero npeSum (NoRF NoAero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)

P_kin_MMpi_pions_data_prompt_noRF_Aero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_prompt_noRF_Aero_cut_all", "MIssing Mass data (NoRF Aero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_data_prompt_noRF_Aero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_data_prompt_noRF_Aero_cut_all", "Electron-Pion CTime (NoRF Aero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_data_prompt_noRF_Aero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_prompt_noRF_Aero_cut_all", "SHMS aero npeSum (NoRF Aero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_data_random_noRF_Aero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_random_noRF_Aero_cut_all", "MIssing Mass data (NoRF Aero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_data_random_noRF_Aero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_data_random_noRF_Aero_cut_all", "Electron-Pion CTime (NoRF Aero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_data_random_noRF_Aero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_random_noRF_Aero_cut_all", "SHMS aero npeSum (NoRF Aero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_data_randomsub_noRF_Aero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_randomsub_noRF_Aero_cut_all", "MIssing Mass data (NoRF Aero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_data_randomsub_noRF_Aero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_data_randomsub_noRF_Aero_cut_all", "Electron-Pion CTime (NoRF Aero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_data_randomsub_noRF_Aero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_randomsub_noRF_Aero_cut_all", "SHMS aero npeSum (NoRF Aero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_dummy_prompt_noRF_Aero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_dummy_prompt_noRF_Aero_cut_all", "MIssing Mass data (NoRF Aero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_dummy_prompt_noRF_Aero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_dummy_prompt_noRF_Aero_cut_all", "Electron-Pion CTime (NoRF Aero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_dummy_prompt_noRF_Aero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_dummy_prompt_noRF_Aero_cut_all", "SHMS aero npeSum (NoRF Aero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_dummy_random_noRF_Aero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_dummy_random_noRF_Aero_cut_all", "MIssing Mass data (NoRF Aero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_dummy_random_noRF_Aero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_dummy_random_noRF_Aero_cut_all", "Electron-Pion CTime (NoRF Aero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_dummy_random_noRF_Aero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_dummy_random_noRF_Aero_cut_all", "SHMS aero npeSum (NoRF Aero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_dummy_randomsub_noRF_Aero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_dummy_randomsub_noRF_Aero_cut_all", "MIssing Mass data (NoRF Aero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_dummy_randomsub_noRF_Aero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_dummy_randomsub_noRF_Aero_cut_all", "Electron-Pion CTime (NoRF Aero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_dummy_randomsub_noRF_Aero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_dummy_randomsub_noRF_Aero_cut_all", "SHMS aero npeSum (NoRF Aero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_data_dummysub_noRF_Aero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_dummysub_noRF_Aero_cut_all", "MIssing Mass data (NoRF Aero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_data_dummysub_noRF_Aero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_data_dummysub_noRF_Aero_cut_all", "Electron-Pion CTime (NoRF Aero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_data_dummysub_noRF_Aero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_dummysub_noRF_Aero_cut_all", "SHMS aero npeSum (NoRF Aero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)

P_kin_MMpi_pions_data_prompt_RF_noAero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_prompt_RF_noAero_cut_all", "MIssing Mass data (RF NoAero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_data_prompt_RF_noAero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_data_prompt_RF_noAero_cut_all", "Electron-Pion CTime (RF NoAero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_data_prompt_RF_noAero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_prompt_RF_noAero_cut_all", "SHMS aero npeSum (RF NoAero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_data_random_RF_noAero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_random_RF_noAero_cut_all", "MIssing Mass data (RF NoAero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_data_random_RF_noAero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_data_random_RF_noAero_cut_all", "Electron-Pion CTime (RF NoAero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_data_random_RF_noAero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_random_RF_noAero_cut_all", "SHMS aero npeSum (RF NoAero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_data_randomsub_RF_noAero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_randomsub_RF_noAero_cut_all", "MIssing Mass data (RF NoAero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_data_randomsub_RF_noAero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_data_randomsub_RF_noAero_cut_all", "Electron-Pion CTime (RF NoAero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_data_randomsub_RF_noAero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_randomsub_RF_noAero_cut_all", "SHMS aero npeSum (RF NoAero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_dummy_prompt_RF_noAero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_dummy_prompt_RF_noAero_cut_all", "MIssing Mass data (RF NoAero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_dummy_prompt_RF_noAero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_dummy_prompt_RF_noAero_cut_all", "Electron-Pion CTime (RF NoAero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_dummy_prompt_RF_noAero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_dummy_prompt_RF_noAero_cut_all", "SHMS aero npeSum (RF NoAero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_dummy_random_RF_noAero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_dummy_random_RF_noAero_cut_all", "MIssing Mass data (RF NoAero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_dummy_random_RF_noAero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_dummy_random_RF_noAero_cut_all", "Electron-Pion CTime (RF NoAero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_dummy_random_RF_noAero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_dummy_random_RF_noAero_cut_all", "SHMS aero npeSum (RF NoAero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_dummy_randomsub_RF_noAero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_dummy_randomsub_RF_noAero_cut_all", "MIssing Mass data (RF NoAero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_dummy_randomsub_RF_noAero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_dummy_randomsub_RF_noAero_cut_all", "Electron-Pion CTime (RF NoAero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_dummy_randomsub_RF_noAero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_dummy_randomsub_RF_noAero_cut_all", "SHMS aero npeSum (RF NoAero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_data_dummysub_RF_noAero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_dummysub_RF_noAero_cut_all", "MIssing Mass data (RF NoAero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_data_dummysub_RF_noAero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_data_dummysub_RF_noAero_cut_all", "Electron-Pion CTime (RF NoAero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_data_dummysub_RF_noAero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_dummysub_RF_noAero_cut_all", "SHMS aero npeSum (RF NoAero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)

P_kin_MMpi_pions_data_prompt_RF_Aero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_prompt_RF_Aero_cut_all", "MIssing Mass data (RF Aero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_data_prompt_RF_Aero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_data_prompt_RF_Aero_cut_all", "Electron-Pion CTime (RF Aero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_data_prompt_RF_Aero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_prompt_RF_Aero_cut_all", "SHMS aero npeSum (RF Aero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_data_random_RF_Aero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_random_RF_Aero_cut_all", "MIssing Mass data (RF Aero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_data_random_RF_Aero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_data_random_RF_Aero_cut_all", "Electron-Pion CTime (RF Aero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_data_random_RF_Aero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_random_RF_Aero_cut_all", "SHMS aero npeSum (RF Aero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_data_randomsub_RF_Aero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_randomsub_RF_Aero_cut_all", "MIssing Mass data (RF Aero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_data_randomsub_RF_Aero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_data_randomsub_RF_Aero_cut_all", "Electron-Pion CTime (RF Aero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_data_randomsub_RF_Aero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_randomsub_RF_Aero_cut_all", "SHMS aero npeSum (RF Aero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_dummy_prompt_RF_Aero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_dummy_prompt_RF_Aero_cut_all", "MIssing Mass data (RF Aero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_dummy_prompt_RF_Aero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_dummy_prompt_RF_Aero_cut_all", "Electron-Pion CTime (RF Aero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_dummy_prompt_RF_Aero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_dummy_prompt_RF_Aero_cut_all", "SHMS aero npeSum (RF Aero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_dummy_random_RF_Aero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_dummy_random_RF_Aero_cut_all", "MIssing Mass data (RF Aero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_dummy_random_RF_Aero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_dummy_random_RF_Aero_cut_all", "Electron-Pion CTime (RF Aero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_dummy_random_RF_Aero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_dummy_random_RF_Aero_cut_all", "SHMS aero npeSum (RF Aero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_dummy_randomsub_RF_Aero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_dummy_randomsub_RF_Aero_cut_all", "MIssing Mass data (RF Aero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_dummy_randomsub_RF_Aero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_dummy_randomsub_RF_Aero_cut_all", "Electron-Pion CTime (RF Aero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_dummy_randomsub_RF_Aero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_dummy_randomsub_RF_Aero_cut_all", "SHMS aero npeSum (RF Aero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_pions_data_dummysub_RF_Aero_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_dummysub_RF_Aero_cut_all", "MIssing Mass data (RF Aero); MM_{\pi}; Counts", MM_nbins, MM_min, MM_max)
CTime_ePiCoinTime_ROC2_pions_data_dummysub_RF_Aero_cut_all = ROOT.TH1D("CTime_ePiCoinTime_ROC2_pions_data_dummysub_RF_Aero_cut_all", "Electron-Pion CTime (RF Aero); e pi Coin_Time; Counts", CT_nbins, CT_min, CT_max)
P_aero_npeSum_pions_data_dummysub_RF_Aero_cut_all = ROOT.TH1D("P_aero_npeSum_pions_data_dummysub_RF_Aero_cut_all", "SHMS aero npeSum (RF Aero); SHMS_aero_npeSum; Counts", Aero_nbins, Aero_min, Aero_max)

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

#2D histograms for 
P_kin_MMpi_vs_aero_npeSum_pions_data_prompt_noRF_noAero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_data_prompt_noRF_noAero_cut_all", "Missing Mass vs Aero npeSum Data Prompt (NoRF NoAero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_data_random_noRF_noAero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_data_random_noRF_noAero_cut_all", "Missing Mass vs Aero npeSum Data Random (NoRF NoAero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_prompt_noRF_noAero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_prompt_noRF_noAero_cut_all", "e-pi Coin Time vs Aero npeSum Data Prompt (NoRF NoAero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_random_noRF_noAero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_random_noRF_noAero_cut_all", "e-pi Coin Time vs Aero npeSum Data Random (NoRF NoAero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_dummy_prompt_noRF_noAero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_dummy_prompt_noRF_noAero_cut_all", "Missing Mass vs Aero npeSum Dummy Prompt (NoRF NoAero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_dummy_random_noRF_noAero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_dummy_random_noRF_noAero_cut_all", "Missing Mass vs Aero npeSum Dummy Random (NoRF NoAero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_prompt_noRF_noAero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_prompt_noRF_noAero_cut_all", "e-pi Coin Time vs Aero npeSum Dummy Prompt (NoRF NoAero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_random_noRF_noAero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_random_noRF_noAero_cut_all", "e-pi Coin Time vs Aero npeSum Dummy Random (NoRF NoAero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_data_randomsub_noRF_noAero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_data_randomsub_noRF_noAero_cut_all", "Missing Mass vs Aero npeSum Data Random Subtracted (NoRF NoAero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_randomsub_noRF_noAero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_randomsub_noRF_noAero_cut_all", "e-pi Coin Time vs Aero npeSum Data Random Subtracted (NoRF NoAero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_dummy_randomsub_noRF_noAero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_dummy_randomsub_noRF_noAero_cut_all", "Missing Mass vs Aero npeSum Dummy Random Subtracted (NoRF NoAero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_randomsub_noRF_noAero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_randomsub_noRF_noAero_cut_all", "e-pi Coin Time vs Aero npeSum Dummy Random Subtracted (NoRF NoAero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_data_dummysub_noRF_noAero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_data_dummysub_noRF_noAero_cut_all", "Missing Mass vs Aero npeSum Data Dummy Subtracted (NoRF NoAero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_dummysub_noRF_noAero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_dummysub_noRF_noAero_cut_all", "e-pi Coin Time vs Aero npeSum Data Dummy Subtracted (NoRF NoAero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)

P_kin_MMpi_vs_aero_npeSum_pions_data_prompt_noRF_Aero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_data_prompt_noRF_Aero_cut_all", "Missing Mass vs Aero npeSum Data Prompt (NoRF Aero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_data_random_noRF_Aero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_data_random_noRF_Aero_cut_all", "Missing Mass vs Aero npeSum Data Random (NoRF Aero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_prompt_noRF_Aero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_prompt_noRF_Aero_cut_all", "e-pi Coin Time vs Aero npeSum Data Prompt (NoRF Aero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_random_noRF_Aero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_random_noRF_Aero_cut_all", "e-pi Coin Time vs Aero npeSum Data Random (NoRF Aero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_dummy_prompt_noRF_Aero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_dummy_prompt_noRF_Aero_cut_all", "Missing Mass vs Aero npeSum Dummy Prompt (NoRF Aero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_dummy_random_noRF_Aero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_dummy_random_noRF_Aero_cut_all", "Missing Mass vs Aero npeSum Dummy Random (NoRF Aero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_prompt_noRF_Aero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_prompt_noRF_Aero_cut_all", "e-pi Coin Time vs Aero npeSum Dummy Prompt (NoRF Aero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_random_noRF_Aero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_random_noRF_Aero_cut_all", "e-pi Coin Time vs Aero npeSum Dummy Random (NoRF Aero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_data_randomsub_noRF_Aero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_data_randomsub_noRF_Aero_cut_all", "Missing Mass vs Aero npeSum Data Random Subtracted (NoRF Aero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_randomsub_noRF_Aero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_randomsub_noRF_Aero_cut_all", "e-pi Coin Time vs Aero npeSum Data Random Subtracted (NoRF Aero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_dummy_randomsub_noRF_Aero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_dummy_randomsub_noRF_Aero_cut_all", "Missing Mass vs Aero npeSum Dummy Random Subtracted (NoRF Aero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_randomsub_noRF_Aero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_randomsub_noRF_Aero_cut_all", "e-pi Coin Time vs Aero npeSum Dummy Random Subtracted (NoRF Aero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_data_dummysub_noRF_Aero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_data_dummysub_noRF_Aero_cut_all", "Missing Mass vs Aero npeSum Data Dummy Subtracted (NoRF Aero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_dummysub_noRF_Aero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_dummysub_noRF_Aero_cut_all", "e-pi Coin Time vs Aero npeSum Data Dummy Subtracted (NoRF Aero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)

P_kin_MMpi_vs_aero_npeSum_pions_data_prompt_RF_noAero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_data_prompt_RF_noAero_cut_all", "Missing Mass vs Aero npeSum Data Prompt (RF NoAero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_data_random_RF_noAero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_data_random_RF_noAero_cut_all", "Missing Mass vs Aero npeSum Data Random (RF NoAero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_prompt_RF_noAero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_prompt_RF_noAero_cut_all", "e-pi Coin Time vs Aero npeSum Data Prompt (RF NoAero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_random_RF_noAero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_random_RF_noAero_cut_all", "e-pi Coin Time vs Aero npeSum Data Random (RF NoAero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_dummy_prompt_RF_noAero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_dummy_prompt_RF_noAero_cut_all", "Missing Mass vs Aero npeSum Dummy Prompt (RF NoAero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_dummy_random_RF_noAero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_dummy_random_RF_noAero_cut_all", "Missing Mass vs Aero npeSum Dummy Random (RF NoAero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_prompt_RF_noAero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_prompt_RF_noAero_cut_all", "e-pi Coin Time vs Aero npeSum Dummy Prompt (RF NoAero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_random_RF_noAero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_random_RF_noAero_cut_all", "e-pi Coin Time vs Aero npeSum Dummy Random (RF NoAero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_data_randomsub_RF_noAero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_data_randomsub_RF_noAero_cut_all", "Missing Mass vs Aero npeSum Data Random Subtracted (RF NoAero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_randomsub_RF_noAero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_randomsub_RF_noAero_cut_all", "e-pi Coin Time vs Aero npeSum Data Random Subtracted (RF NoAero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_dummy_randomsub_RF_noAero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_dummy_randomsub_RF_noAero_cut_all", "Missing Mass vs Aero npeSum Dummy Random Subtracted (RF NoAero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_randomsub_RF_noAero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_randomsub_RF_noAero_cut_all", "e-pi Coin Time vs Aero npeSum Dummy Random Subtracted (RF NoAero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_data_dummysub_RF_noAero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_data_dummysub_RF_noAero_cut_all", "Missing Mass vs Aero npeSum Data Dummy Subtracted (RF NoAero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_dummysub_RF_noAero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_dummysub_RF_noAero_cut_all", "e-pi Coin Time vs Aero npeSum Data Dummy Subtracted (RF NoAero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)

P_kin_MMpi_vs_aero_npeSum_pions_data_prompt_RF_Aero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_data_prompt_RF_Aero_cut_all", "Missing Mass vs Aero npeSum Data Prompt (RF Aero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_data_random_RF_Aero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_data_random_RF_Aero_cut_all", "Missing Mass vs Aero npeSum Data Random (RF Aero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_prompt_RF_Aero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_prompt_RF_Aero_cut_all", "e-pi Coin Time vs Aero npeSum Data Prompt (RF Aero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_random_RF_Aero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_random_RF_Aero_cut_all", "e-pi Coin Time vs Aero npeSum Data Random (RF Aero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_dummy_prompt_RF_Aero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_dummy_prompt_RF_Aero_cut_all", "Missing Mass vs Aero npeSum Dummy Prompt (RF Aero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_dummy_random_RF_Aero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_dummy_random_RF_Aero_cut_all", "Missing Mass vs Aero npeSum Dummy Random (RF Aero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_prompt_RF_Aero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_prompt_RF_Aero_cut_all", "e-pi Coin Time vs Aero npeSum Dummy Prompt (RF Aero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_random_RF_Aero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_random_RF_Aero_cut_all", "e-pi Coin Time vs Aero npeSum Dummy Random (RF Aero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_data_randomsub_RF_Aero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_data_randomsub_RF_Aero_cut_all", "Missing Mass vs Aero npeSum Data Random Subtracted (RF Aero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_randomsub_RF_Aero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero  npeSum_pions_data_randomsub_RF_Aero_cut_all", "e-pi Coin Time vs Aero npeSum Data Random Subtracted (RF Aero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_dummy_randomsub_RF_Aero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_dummy_randomsub_RF_Aero_cut_all", "Missing Mass vs Aero npeSum Dummy Random Subtracted (RF Aero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_randomsub_RF_Aero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_randomsub_RF_Aero_cut_all", "e-pi Coin Time vs Aero npeSum Dummy Random Subtracted (RF Aero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)
P_kin_MMpi_vs_aero_npeSum_pions_data_dummysub_RF_Aero_cut_all = ROOT.TH2D("P_kin_MMpi_vs_aero_npeSum_pions_data_dummysub_RF_Aero_cut_all", "Missing Mass vs Aero npeSum Data Dummy Subtracted (RF Aero); MM_{\pi}; SHMS_aero_npeSum", MM_nbins, MM_min, MM_max, Aero_nbins, Aero_min, Aero_max)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_dummysub_RF_Aero_cut_all = ROOT.TH2D("CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_dummysub_RF_Aero_cut_all", "e-pi Coin Time vs Aero npeSum Data Dummy Subtracted (RF Aero); e pi Coin_Time; SHMS_aero_npeSum", CT_nbins, CT_min, CT_max, Aero_nbins, Aero_min, Aero_max)

##########################################################################################################################################################################################################

#Fill histograms for Cut All Data
for event in Cut_Pion_Events_Prompt_Data_tree:
    P_kin_MMpi_pions_data_prompt_noRF_noAero_cut_all.Fill(event.MMpi)
    CTime_ePiCoinTime_ROC2_pions_data_prompt_noRF_noAero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
    P_aero_npeSum_pions_data_prompt_noRF_noAero_cut_all.Fill(event.P_aero_npeSum)
    P_kin_MMpi_vs_aero_npeSum_pions_data_prompt_noRF_noAero_cut_all.Fill(event.MMpi, event.P_aero_npeSum)
    CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_prompt_noRF_noAero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2, event.P_aero_npeSum)
for event in Cut_Pion_Events_Random_Data_tree:
    P_kin_MMpi_pions_data_random_noRF_noAero_cut_all.Fill(event.MMpi)
    CTime_ePiCoinTime_ROC2_pions_data_random_noRF_noAero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
    P_aero_npeSum_pions_data_random_noRF_noAero_cut_all.Fill(event.P_aero_npeSum)
    P_kin_MMpi_vs_aero_npeSum_pions_data_random_noRF_noAero_cut_all.Fill(event.MMpi, event.P_aero_npeSum)
    CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_random_noRF_noAero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2, event.P_aero_npeSum)
for event in Cut_Pion_Events_Prompt_Dummy_tree:
    P_kin_MMpi_pions_dummy_prompt_noRF_noAero_cut_all.Fill(event.MMpi)
    CTime_ePiCoinTime_ROC2_pions_dummy_prompt_noRF_noAero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
    P_aero_npeSum_pions_dummy_prompt_noRF_noAero_cut_all.Fill(event.P_aero_npeSum)
    P_kin_MMpi_vs_aero_npeSum_pions_dummy_prompt_noRF_noAero_cut_all.Fill(event.MMpi, event.P_aero_npeSum)
    CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_prompt_noRF_noAero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2, event.P_aero_npeSum)
for event in Cut_Pion_Events_Random_Dummy_tree:
    P_kin_MMpi_pions_dummy_random_noRF_noAero_cut_all.Fill(event.MMpi)
    CTime_ePiCoinTime_ROC2_pions_dummy_random_noRF_noAero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
    P_aero_npeSum_pions_dummy_random_noRF_noAero_cut_all.Fill(event.P_aero_npeSum)
    P_kin_MMpi_vs_aero_npeSum_pions_dummy_random_noRF_noAero_cut_all.Fill(event.MMpi, event.P_aero_npeSum)
    CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_random_noRF_noAero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2, event.P_aero_npeSum)

for event in Cut_Pion_Events_Prompt_Data_tree:
    if Aero_Detector_Cut(event):
        P_kin_MMpi_pions_data_prompt_noRF_Aero_cut_all.Fill(event.MMpi)
        CTime_ePiCoinTime_ROC2_pions_data_prompt_noRF_Aero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
        P_aero_npeSum_pions_data_prompt_noRF_Aero_cut_all.Fill(event.P_aero_npeSum)
        P_kin_MMpi_vs_aero_npeSum_pions_data_prompt_noRF_Aero_cut_all.Fill(event.MMpi, event.P_aero_npeSum)
        CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_prompt_noRF_Aero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2, event.P_aero_npeSum)
for event in Cut_Pion_Events_Random_Data_tree:
    if Aero_Detector_Cut(event):
        P_kin_MMpi_pions_data_random_noRF_Aero_cut_all.Fill(event.MMpi)
        CTime_ePiCoinTime_ROC2_pions_data_random_noRF_Aero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
        P_aero_npeSum_pions_data_random_noRF_Aero_cut_all.Fill(event.P_aero_npeSum)
        P_kin_MMpi_vs_aero_npeSum_pions_data_random_noRF_Aero_cut_all.Fill(event.MMpi, event.P_aero_npeSum)
        CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_random_noRF_Aero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2, event.P_aero_npeSum)
for event in Cut_Pion_Events_Prompt_Dummy_tree:
    if Aero_Detector_Cut(event):
        P_kin_MMpi_pions_dummy_prompt_noRF_Aero_cut_all.Fill(event.MMpi)
        CTime_ePiCoinTime_ROC2_pions_dummy_prompt_noRF_Aero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
        P_aero_npeSum_pions_dummy_prompt_noRF_Aero_cut_all.Fill(event.P_aero_npeSum)
        P_kin_MMpi_vs_aero_npeSum_pions_dummy_prompt_noRF_Aero_cut_all.Fill(event.MMpi, event.P_aero_npeSum)
        CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_prompt_noRF_Aero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2, event.P_aero_npeSum)
for event in Cut_Pion_Events_Random_Dummy_tree:
    if Aero_Detector_Cut(event):
        P_kin_MMpi_pions_dummy_random_noRF_Aero_cut_all.Fill(event.MMpi)
        CTime_ePiCoinTime_ROC2_pions_dummy_random_noRF_Aero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
        P_aero_npeSum_pions_dummy_random_noRF_Aero_cut_all.Fill(event.P_aero_npeSum)
        P_kin_MMpi_vs_aero_npeSum_pions_dummy_random_noRF_Aero_cut_all.Fill(event.MMpi, event.P_aero_npeSum)
        CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_random_noRF_Aero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2, event.P_aero_npeSum)

for event in Cut_Pion_Events_Prompt_Data_tree:
    if RF_Cut(event):
        P_kin_MMpi_pions_data_prompt_RF_noAero_cut_all.Fill(event.MMpi)
        CTime_ePiCoinTime_ROC2_pions_data_prompt_RF_noAero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
        P_aero_npeSum_pions_data_prompt_RF_noAero_cut_all.Fill(event.P_aero_npeSum)
        P_kin_MMpi_vs_aero_npeSum_pions_data_prompt_RF_noAero_cut_all.Fill(event.MMpi, event.P_aero_npeSum)
        CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_prompt_RF_noAero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2, event.P_aero_npeSum)
for event in Cut_Pion_Events_Random_Data_tree:
    if RF_Cut(event):
        P_kin_MMpi_pions_data_random_RF_noAero_cut_all.Fill(event.MMpi)
        CTime_ePiCoinTime_ROC2_pions_data_random_RF_noAero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
        P_aero_npeSum_pions_data_random_RF_noAero_cut_all.Fill(event.P_aero_npeSum)
        P_kin_MMpi_vs_aero_npeSum_pions_data_random_RF_noAero_cut_all.Fill(event.MMpi, event.P_aero_npeSum)
        CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_random_RF_noAero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2, event.P_aero_npeSum)
for event in Cut_Pion_Events_Prompt_Dummy_tree:
    if RF_Cut(event):
        P_kin_MMpi_pions_dummy_prompt_RF_noAero_cut_all.Fill(event.MMpi)
        CTime_ePiCoinTime_ROC2_pions_dummy_prompt_RF_noAero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
        P_aero_npeSum_pions_dummy_prompt_RF_noAero_cut_all.Fill(event.P_aero_npeSum)
        P_kin_MMpi_vs_aero_npeSum_pions_dummy_prompt_RF_noAero_cut_all.Fill(event.MMpi, event.P_aero_npeSum)
        CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_prompt_RF_noAero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2, event.P_aero_npeSum)
for event in Cut_Pion_Events_Random_Dummy_tree:
    if RF_Cut(event):
        P_kin_MMpi_pions_dummy_random_RF_noAero_cut_all.Fill(event.MMpi)
        CTime_ePiCoinTime_ROC2_pions_dummy_random_RF_noAero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
        P_aero_npeSum_pions_dummy_random_RF_noAero_cut_all.Fill(event.P_aero_npeSum)
        P_kin_MMpi_vs_aero_npeSum_pions_dummy_random_RF_noAero_cut_all.Fill(event.MMpi, event.P_aero_npeSum)
        CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_random_RF_noAero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2, event.P_aero_npeSum)

for event in Cut_Pion_Events_Prompt_Data_tree:
    if Aero_Detector_Cut(event) and RF_Cut(event):
        P_kin_MMpi_pions_data_prompt_RF_Aero_cut_all.Fill(event.MMpi)
        CTime_ePiCoinTime_ROC2_pions_data_prompt_RF_Aero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
        P_aero_npeSum_pions_data_prompt_RF_Aero_cut_all.Fill(event.P_aero_npeSum)
        P_kin_MMpi_vs_aero_npeSum_pions_data_prompt_RF_Aero_cut_all.Fill(event.MMpi, event.P_aero_npeSum)
        CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_prompt_RF_Aero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2, event.P_aero_npeSum)
for event in Cut_Pion_Events_Random_Data_tree:
    if Aero_Detector_Cut(event) and RF_Cut(event):
        P_kin_MMpi_pions_data_random_RF_Aero_cut_all.Fill(event.MMpi)
        CTime_ePiCoinTime_ROC2_pions_data_random_RF_Aero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
        P_aero_npeSum_pions_data_random_RF_Aero_cut_all.Fill(event.P_aero_npeSum)
        P_kin_MMpi_vs_aero_npeSum_pions_data_random_RF_Aero_cut_all.Fill(event.MMpi, event.P_aero_npeSum)
        CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_random_RF_Aero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2, event.P_aero_npeSum)
for event in Cut_Pion_Events_Prompt_Dummy_tree:
    if Aero_Detector_Cut(event) and RF_Cut(event):
        P_kin_MMpi_pions_dummy_prompt_RF_Aero_cut_all.Fill(event.MMpi)
        CTime_ePiCoinTime_ROC2_pions_dummy_prompt_RF_Aero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
        P_aero_npeSum_pions_dummy_prompt_RF_Aero_cut_all.Fill(event.P_aero_npeSum)
        P_kin_MMpi_vs_aero_npeSum_pions_dummy_prompt_RF_Aero_cut_all.Fill(event.MMpi, event.P_aero_npeSum)
        CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_prompt_RF_Aero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2, event.P_aero_npeSum)
for event in Cut_Pion_Events_Random_Dummy_tree:
    if Aero_Detector_Cut(event) and RF_Cut(event):
        P_kin_MMpi_pions_dummy_random_RF_Aero_cut_all.Fill(event.MMpi)
        CTime_ePiCoinTime_ROC2_pions_dummy_random_RF_Aero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2)
        P_aero_npeSum_pions_dummy_random_RF_Aero_cut_all.Fill(event.P_aero_npeSum)
        P_kin_MMpi_vs_aero_npeSum_pions_dummy_random_RF_Aero_cut_all.Fill(event.MMpi, event.P_aero_npeSum)
        CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_random_RF_Aero_cut_all.Fill(event.CTime_ePiCoinTime_ROC2, event.P_aero_npeSum)

print("####################################")
print("###### Histogram filling done ######")
print("####################################\n")

#################################################################################################################################################

# Random subtraction from Histograms
Clone_CTime_ePiCoinTime_ROC2_pions_data_random_noRF_noAero_cut_all = CTime_ePiCoinTime_ROC2_pions_data_random_noRF_noAero_cut_all.Clone("Clone_CTime_ePiCoinTime_ROC2_pions_data_random_noRF_noAero_cut_all")
Clone_CTime_ePiCoinTime_ROC2_pions_data_random_noRF_Aero_cut_all = CTime_ePiCoinTime_ROC2_pions_data_random_noRF_Aero_cut_all.Clone("Clone_CTime_ePiCoinTime_ROC2_pions_data_random_noRF_Aero_cut_all")
Clone_CTime_ePiCoinTime_ROC2_pions_data_random_RF_noAero_cut_all = CTime_ePiCoinTime_ROC2_pions_data_random_RF_noAero_cut_all.Clone("Clone_CTime_ePiCoinTime_ROC2_pions_data_random_RF_noAero_cut_all")
Clone_CTime_ePiCoinTime_ROC2_pions_data_random_RF_Aero_cut_all = CTime_ePiCoinTime_ROC2_pions_data_random_RF_Aero_cut_all.Clone("Clone_CTime_ePiCoinTime_ROC2_pions_data_random_RF_Aero_cut_all")

P_kin_MMpi_pions_data_random_noRF_noAero_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_ROC2_pions_data_random_noRF_noAero_cut_all.Scale(1.0/nWindows)
P_aero_npeSum_pions_data_random_noRF_noAero_cut_all.Scale(1.0/nWindows)
P_kin_MMpi_pions_data_random_noRF_Aero_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_ROC2_pions_data_random_noRF_Aero_cut_all.Scale(1.0/nWindows)
P_aero_npeSum_pions_data_random_noRF_Aero_cut_all.Scale(1.0/nWindows)
P_kin_MMpi_pions_data_random_RF_noAero_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_ROC2_pions_data_random_RF_noAero_cut_all.Scale(1.0/nWindows)
P_aero_npeSum_pions_data_random_RF_noAero_cut_all.Scale(1.0/nWindows)
P_kin_MMpi_pions_data_random_RF_Aero_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_ROC2_pions_data_random_RF_Aero_cut_all.Scale(1.0/nWindows)
P_aero_npeSum_pions_data_random_RF_Aero_cut_all.Scale(1.0/nWindows)
P_kin_MMpi_vs_aero_npeSum_pions_data_random_noRF_noAero_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_random_noRF_noAero_cut_all.Scale(1.0/nWindows)
P_kin_MMpi_vs_aero_npeSum_pions_data_random_noRF_Aero_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_random_noRF_Aero_cut_all.Scale(1.0/nWindows)
P_kin_MMpi_vs_aero_npeSum_pions_data_random_RF_noAero_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_random_RF_noAero_cut_all.Scale(1.0/nWindows)
P_kin_MMpi_vs_aero_npeSum_pions_data_random_RF_Aero_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_random_RF_Aero_cut_all.Scale(1.0/nWindows)

P_kin_MMpi_pions_dummy_random_noRF_noAero_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_ROC2_pions_dummy_random_noRF_noAero_cut_all.Scale(1.0/nWindows)
P_aero_npeSum_pions_dummy_random_noRF_noAero_cut_all.Scale(1.0/nWindows)
P_kin_MMpi_pions_dummy_random_noRF_Aero_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_ROC2_pions_dummy_random_noRF_Aero_cut_all.Scale(1.0/nWindows)
P_aero_npeSum_pions_dummy_random_noRF_Aero_cut_all.Scale(1.0/nWindows)
P_kin_MMpi_pions_dummy_random_RF_noAero_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_ROC2_pions_dummy_random_RF_noAero_cut_all.Scale(1.0/nWindows)
P_aero_npeSum_pions_dummy_random_RF_noAero_cut_all.Scale(1.0/nWindows)
P_kin_MMpi_pions_dummy_random_RF_Aero_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_ROC2_pions_dummy_random_RF_Aero_cut_all.Scale(1.0/nWindows)
P_aero_npeSum_pions_dummy_random_RF_Aero_cut_all.Scale(1.0/nWindows)
P_kin_MMpi_vs_aero_npeSum_pions_dummy_random_noRF_noAero_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_random_noRF_noAero_cut_all.Scale(1.0/nWindows)
P_kin_MMpi_vs_aero_npeSum_pions_dummy_random_noRF_Aero_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_random_noRF_Aero_cut_all.Scale(1.0/nWindows)
P_kin_MMpi_vs_aero_npeSum_pions_dummy_random_RF_noAero_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_random_RF_noAero_cut_all.Scale(1.0/nWindows)
P_kin_MMpi_vs_aero_npeSum_pions_dummy_random_RF_Aero_cut_all.Scale(1.0/nWindows)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_random_RF_Aero_cut_all.Scale(1.0/nWindows)

P_kin_MMpi_pions_data_randomsub_noRF_noAero_cut_all.Add(P_kin_MMpi_pions_data_prompt_noRF_noAero_cut_all, P_kin_MMpi_pions_data_random_noRF_noAero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_pions_data_randomsub_noRF_noAero_cut_all.Add(CTime_ePiCoinTime_ROC2_pions_data_prompt_noRF_noAero_cut_all, CTime_ePiCoinTime_ROC2_pions_data_random_noRF_noAero_cut_all, 1, -1)
P_aero_npeSum_pions_data_randomsub_noRF_noAero_cut_all.Add(P_aero_npeSum_pions_data_prompt_noRF_noAero_cut_all, P_aero_npeSum_pions_data_random_noRF_noAero_cut_all, 1, -1)
P_kin_MMpi_pions_data_randomsub_noRF_Aero_cut_all.Add(P_kin_MMpi_pions_data_prompt_noRF_Aero_cut_all, P_kin_MMpi_pions_data_random_noRF_Aero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_pions_data_randomsub_noRF_Aero_cut_all.Add(CTime_ePiCoinTime_ROC2_pions_data_prompt_noRF_Aero_cut_all, CTime_ePiCoinTime_ROC2_pions_data_random_noRF_Aero_cut_all, 1, -1)
P_aero_npeSum_pions_data_randomsub_noRF_Aero_cut_all.Add(P_aero_npeSum_pions_data_prompt_noRF_Aero_cut_all, P_aero_npeSum_pions_data_random_noRF_Aero_cut_all, 1, -1)
P_kin_MMpi_pions_data_randomsub_RF_noAero_cut_all.Add(P_kin_MMpi_pions_data_prompt_RF_noAero_cut_all, P_kin_MMpi_pions_data_random_RF_noAero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_pions_data_randomsub_RF_noAero_cut_all.Add(CTime_ePiCoinTime_ROC2_pions_data_prompt_RF_noAero_cut_all, CTime_ePiCoinTime_ROC2_pions_data_random_RF_noAero_cut_all, 1, -1)
P_aero_npeSum_pions_data_randomsub_RF_noAero_cut_all.Add(P_aero_npeSum_pions_data_prompt_RF_noAero_cut_all, P_aero_npeSum_pions_data_random_RF_noAero_cut_all, 1, -1)
P_kin_MMpi_pions_data_randomsub_RF_Aero_cut_all.Add(P_kin_MMpi_pions_data_prompt_RF_Aero_cut_all, P_kin_MMpi_pions_data_random_RF_Aero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_pions_data_randomsub_RF_Aero_cut_all.Add(CTime_ePiCoinTime_ROC2_pions_data_prompt_RF_Aero_cut_all, CTime_ePiCoinTime_ROC2_pions_data_random_RF_Aero_cut_all, 1, -1)
P_aero_npeSum_pions_data_randomsub_RF_Aero_cut_all.Add(P_aero_npeSum_pions_data_prompt_RF_Aero_cut_all, P_aero_npeSum_pions_data_random_RF_Aero_cut_all, 1, -1)
P_kin_MMpi_vs_aero_npeSum_pions_data_randomsub_noRF_noAero_cut_all.Add(P_kin_MMpi_vs_aero_npeSum_pions_data_prompt_noRF_noAero_cut_all, P_kin_MMpi_vs_aero_npeSum_pions_data_random_noRF_noAero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_randomsub_noRF_noAero_cut_all.Add(CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_prompt_noRF_noAero_cut_all, CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_random_noRF_noAero_cut_all, 1, -1)
P_kin_MMpi_vs_aero_npeSum_pions_data_randomsub_noRF_Aero_cut_all.Add(P_kin_MMpi_vs_aero_npeSum_pions_data_prompt_noRF_Aero_cut_all, P_kin_MMpi_vs_aero_npeSum_pions_data_random_noRF_Aero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_randomsub_noRF_Aero_cut_all.Add(CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_prompt_noRF_Aero_cut_all, CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_random_noRF_Aero_cut_all, 1, -1)
P_kin_MMpi_vs_aero_npeSum_pions_data_randomsub_RF_noAero_cut_all.Add(P_kin_MMpi_vs_aero_npeSum_pions_data_prompt_RF_noAero_cut_all, P_kin_MMpi_vs_aero_npeSum_pions_data_random_RF_noAero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_randomsub_RF_noAero_cut_all.Add(CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_prompt_RF_noAero_cut_all, CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_random_RF_noAero_cut_all, 1, -1)
P_kin_MMpi_vs_aero_npeSum_pions_data_randomsub_RF_Aero_cut_all.Add(P_kin_MMpi_vs_aero_npeSum_pions_data_prompt_RF_Aero_cut_all, P_kin_MMpi_vs_aero_npeSum_pions_data_random_RF_Aero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_randomsub_RF_Aero_cut_all.Add(CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_prompt_RF_Aero_cut_all, CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_random_RF_Aero_cut_all, 1, -1)

P_kin_MMpi_pions_dummy_randomsub_noRF_noAero_cut_all.Add(P_kin_MMpi_pions_dummy_prompt_noRF_noAero_cut_all, P_kin_MMpi_pions_dummy_random_noRF_noAero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_pions_dummy_randomsub_noRF_noAero_cut_all.Add(CTime_ePiCoinTime_ROC2_pions_dummy_prompt_noRF_noAero_cut_all, CTime_ePiCoinTime_ROC2_pions_dummy_random_noRF_noAero_cut_all, 1, -1)
P_aero_npeSum_pions_dummy_randomsub_noRF_noAero_cut_all.Add(P_aero_npeSum_pions_dummy_prompt_noRF_noAero_cut_all, P_aero_npeSum_pions_dummy_random_noRF_noAero_cut_all, 1, -1)
P_kin_MMpi_pions_dummy_randomsub_noRF_Aero_cut_all.Add(P_kin_MMpi_pions_dummy_prompt_noRF_Aero_cut_all, P_kin_MMpi_pions_dummy_random_noRF_Aero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_pions_dummy_randomsub_noRF_Aero_cut_all.Add(CTime_ePiCoinTime_ROC2_pions_dummy_prompt_noRF_Aero_cut_all, CTime_ePiCoinTime_ROC2_pions_dummy_random_noRF_Aero_cut_all, 1, -1)
P_aero_npeSum_pions_dummy_randomsub_noRF_Aero_cut_all.Add(P_aero_npeSum_pions_dummy_prompt_noRF_Aero_cut_all, P_aero_npeSum_pions_dummy_random_noRF_Aero_cut_all, 1, -1)
P_kin_MMpi_pions_dummy_randomsub_RF_noAero_cut_all.Add(P_kin_MMpi_pions_dummy_prompt_RF_noAero_cut_all, P_kin_MMpi_pions_dummy_random_RF_noAero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_pions_dummy_randomsub_RF_noAero_cut_all.Add(CTime_ePiCoinTime_ROC2_pions_dummy_prompt_RF_noAero_cut_all, CTime_ePiCoinTime_ROC2_pions_dummy_random_RF_noAero_cut_all, 1, -1)
P_aero_npeSum_pions_dummy_randomsub_RF_noAero_cut_all.Add(P_aero_npeSum_pions_dummy_prompt_RF_noAero_cut_all, P_aero_npeSum_pions_dummy_random_RF_noAero_cut_all, 1, -1)
P_kin_MMpi_pions_dummy_randomsub_RF_Aero_cut_all.Add(P_kin_MMpi_pions_dummy_prompt_RF_Aero_cut_all, P_kin_MMpi_pions_dummy_random_RF_Aero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_pions_dummy_randomsub_RF_Aero_cut_all.Add(CTime_ePiCoinTime_ROC2_pions_dummy_prompt_RF_Aero_cut_all, CTime_ePiCoinTime_ROC2_pions_dummy_random_RF_Aero_cut_all, 1, -1)
P_aero_npeSum_pions_dummy_randomsub_RF_Aero_cut_all.Add(P_aero_npeSum_pions_dummy_prompt_RF_Aero_cut_all, P_aero_npeSum_pions_dummy_random_RF_Aero_cut_all, 1, -1)
P_kin_MMpi_vs_aero_npeSum_pions_dummy_randomsub_noRF_noAero_cut_all.Add(P_kin_MMpi_vs_aero_npeSum_pions_dummy_prompt_noRF_noAero_cut_all, P_kin_MMpi_vs_aero_npeSum_pions_dummy_random_noRF_noAero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_randomsub_noRF_noAero_cut_all.Add(CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_prompt_noRF_noAero_cut_all, CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_random_noRF_noAero_cut_all, 1, -1)
P_kin_MMpi_vs_aero_npeSum_pions_dummy_randomsub_noRF_Aero_cut_all.Add(P_kin_MMpi_vs_aero_npeSum_pions_dummy_prompt_noRF_Aero_cut_all, P_kin_MMpi_vs_aero_npeSum_pions_dummy_random_noRF_Aero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_randomsub_noRF_Aero_cut_all.Add(CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_prompt_noRF_Aero_cut_all, CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_random_noRF_Aero_cut_all, 1, -1)
P_kin_MMpi_vs_aero_npeSum_pions_dummy_randomsub_RF_noAero_cut_all.Add(P_kin_MMpi_vs_aero_npeSum_pions_dummy_prompt_RF_noAero_cut_all, P_kin_MMpi_vs_aero_npeSum_pions_dummy_random_RF_noAero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_randomsub_RF_noAero_cut_all.Add(CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_prompt_RF_noAero_cut_all, CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_random_RF_noAero_cut_all, 1, -1)
P_kin_MMpi_vs_aero_npeSum_pions_dummy_randomsub_RF_Aero_cut_all.Add(P_kin_MMpi_vs_aero_npeSum_pions_dummy_prompt_RF_Aero_cut_all, P_kin_MMpi_vs_aero_npeSum_pions_dummy_random_RF_Aero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_randomsub_RF_Aero_cut_all.Add(CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_prompt_RF_Aero_cut_all, CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_random_RF_Aero_cut_all, 1, -1)

print("######################################")
print("###### Random substraction done ######")
print("######################################\n")

############################################################################################################################################

# Dummy Subtraction
P_kin_MMpi_pions_data_dummysub_noRF_noAero_cut_all.Add(P_kin_MMpi_pions_data_randomsub_noRF_noAero_cut_all, P_kin_MMpi_pions_dummy_randomsub_noRF_noAero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_pions_data_dummysub_noRF_noAero_cut_all.Add(CTime_ePiCoinTime_ROC2_pions_data_randomsub_noRF_noAero_cut_all, CTime_ePiCoinTime_ROC2_pions_dummy_randomsub_noRF_noAero_cut_all, 1, -1)
P_aero_npeSum_pions_data_dummysub_noRF_noAero_cut_all.Add(P_aero_npeSum_pions_data_randomsub_noRF_noAero_cut_all, P_aero_npeSum_pions_dummy_randomsub_noRF_noAero_cut_all, 1, -1)
P_kin_MMpi_vs_aero_npeSum_pions_data_dummysub_noRF_noAero_cut_all.Add(P_kin_MMpi_vs_aero_npeSum_pions_data_randomsub_noRF_noAero_cut_all, P_kin_MMpi_vs_aero_npeSum_pions_dummy_randomsub_noRF_noAero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_dummysub_noRF_noAero_cut_all.Add(CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_randomsub_noRF_noAero_cut_all, CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_randomsub_noRF_noAero_cut_all, 1, -1)

P_kin_MMpi_pions_data_dummysub_noRF_Aero_cut_all.Add(P_kin_MMpi_pions_data_randomsub_noRF_Aero_cut_all, P_kin_MMpi_pions_dummy_randomsub_noRF_Aero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_pions_data_dummysub_noRF_Aero_cut_all.Add(CTime_ePiCoinTime_ROC2_pions_data_randomsub_noRF_Aero_cut_all, CTime_ePiCoinTime_ROC2_pions_dummy_randomsub_noRF_Aero_cut_all, 1, -1)
P_aero_npeSum_pions_data_dummysub_noRF_Aero_cut_all.Add(P_aero_npeSum_pions_data_randomsub_noRF_Aero_cut_all, P_aero_npeSum_pions_dummy_randomsub_noRF_Aero_cut_all, 1, -1)
P_kin_MMpi_vs_aero_npeSum_pions_data_dummysub_noRF_Aero_cut_all.Add(P_kin_MMpi_vs_aero_npeSum_pions_data_randomsub_noRF_Aero_cut_all, P_kin_MMpi_vs_aero_npeSum_pions_dummy_randomsub_noRF_Aero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_dummysub_noRF_Aero_cut_all.Add(CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_randomsub_noRF_Aero_cut_all, CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_randomsub_noRF_Aero_cut_all, 1, -1)

P_kin_MMpi_pions_data_dummysub_RF_noAero_cut_all.Add(P_kin_MMpi_pions_data_randomsub_RF_noAero_cut_all, P_kin_MMpi_pions_dummy_randomsub_RF_noAero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_pions_data_dummysub_RF_noAero_cut_all.Add(CTime_ePiCoinTime_ROC2_pions_data_randomsub_RF_noAero_cut_all, CTime_ePiCoinTime_ROC2_pions_dummy_randomsub_RF_noAero_cut_all, 1, -1)
P_aero_npeSum_pions_data_dummysub_RF_noAero_cut_all.Add(P_aero_npeSum_pions_data_randomsub_RF_noAero_cut_all, P_aero_npeSum_pions_dummy_randomsub_RF_noAero_cut_all, 1, -1)
P_kin_MMpi_vs_aero_npeSum_pions_data_dummysub_RF_noAero_cut_all.Add(P_kin_MMpi_vs_aero_npeSum_pions_data_randomsub_RF_noAero_cut_all, P_kin_MMpi_vs_aero_npeSum_pions_dummy_randomsub_RF_noAero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_dummysub_RF_noAero_cut_all.Add(CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_randomsub_RF_noAero_cut_all, CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_randomsub_RF_noAero_cut_all, 1, -1)

P_kin_MMpi_pions_data_dummysub_RF_Aero_cut_all.Add(P_kin_MMpi_pions_data_randomsub_RF_Aero_cut_all, P_kin_MMpi_pions_dummy_randomsub_RF_Aero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_pions_data_dummysub_RF_Aero_cut_all.Add(CTime_ePiCoinTime_ROC2_pions_data_randomsub_RF_Aero_cut_all, CTime_ePiCoinTime_ROC2_pions_dummy_randomsub_RF_Aero_cut_all, 1, -1)
P_aero_npeSum_pions_data_dummysub_RF_Aero_cut_all.Add(P_aero_npeSum_pions_data_randomsub_RF_Aero_cut_all, P_aero_npeSum_pions_dummy_randomsub_RF_Aero_cut_all, 1, -1)
P_kin_MMpi_vs_aero_npeSum_pions_data_dummysub_RF_Aero_cut_all.Add(P_kin_MMpi_vs_aero_npeSum_pions_data_randomsub_RF_Aero_cut_all, P_kin_MMpi_vs_aero_npeSum_pions_dummy_randomsub_RF_Aero_cut_all, 1, -1)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_dummysub_RF_Aero_cut_all.Add(CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_randomsub_RF_Aero_cut_all, CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_dummy_randomsub_RF_Aero_cut_all, 1, -1)

print("####################################")
print("###### Dummy subtraction done ######")
print("####################################\n")

#############################################################################################################################################################

# RF efficiency calculation
#Aero_Ndid = Aero_Eff_Ndid_dummysub_data_cut_all.Integral()
#Aero_Nshould = Aero_Eff_Nshould_dummysub_data_cut_all.Integral()

Aero_Ndid = P_aero_npeSum_pions_data_dummysub_RF_Aero_cut_all.Integral()
Aero_Nshould = P_aero_npeSum_pions_data_dummysub_RF_noAero_cut_all.Integral()

Aero_Eff = Aero_Ndid / Aero_Nshould

Aero_EFF_Err = ma.sqrt(((Aero_Nshould * Aero_Ndid) - (Aero_Ndid**2)) / (Aero_Nshould**3))

print("="*50)
print("###### RF Efficiency Calculation #########")
print("Cut applied on Aerogel detector: P_aero_npeSum > %.2f" % SHMS_Aero_cut)
print("="*50)
print("RF Ndid: %.2f" % Aero_Ndid)
print("RF Nshould: %.2f" % Aero_Nshould)
print("RF Efficiency: %.5f +/- %.5f" % (Aero_Eff, Aero_EFF_Err))
print("="*50)

# Writing output to CSV file
csv_Aero_Eff = os.path.join(UTILPATH, "efficiencies", "%s_PionLT_coin_prod_Aero_efficiency_data.csv" % setting_name)
try:
    with open(csv_Aero_Eff, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["physics setting", "Aero_eff", "Aero_eff_error"])
        writer.writerow([PHY_SETTING, "{:.5f}".format(Aero_Eff), "{:.5f}".format(Aero_EFF_Err)])
    print("Wrote RF efficiency to %s" % csv_Aero_Eff)
except Exception as e:
    print("ERROR writing CSV file: %s" % e)

##################################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Save PID hostograms to output PDF
PID_output_pdf_file = os.path.join(UTILPATH, "scripts/efficiency/OUTPUTS/plots", "%s_PionLT_coin_prod_SHMS_PID.pdf" % setting_name)

c1_PID1 = TCanvas("c1_PID1", "PID Efficiency", 100, 0, 1200, 1800)
c1_PID1.Divide(2,3)
c1_PID1.cd(1)
P_kin_MMpi_pions_data_dummysub_noRF_noAero_cut_all.SetLineColor(kBlue)
P_kin_MMpi_pions_data_dummysub_noRF_noAero_cut_all.SetLineWidth(2)
P_kin_MMpi_pions_data_dummysub_noRF_noAero_cut_all.Draw("Hist")
c1_PID1.cd(2)
P_kin_MMpi_pions_data_dummysub_noRF_Aero_cut_all.SetLineColor(kRed)
P_kin_MMpi_pions_data_dummysub_noRF_Aero_cut_all.SetLineWidth(2)
P_kin_MMpi_pions_data_dummysub_noRF_Aero_cut_all.Draw("Hist")
c1_PID1.cd(3)
CTime_ePiCoinTime_ROC2_pions_data_dummysub_noRF_noAero_cut_all.SetLineColor(kBlue)
CTime_ePiCoinTime_ROC2_pions_data_dummysub_noRF_noAero_cut_all.SetLineWidth(2)
CTime_ePiCoinTime_ROC2_pions_data_dummysub_noRF_noAero_cut_all.Draw("Hist")
Clone_CTime_ePiCoinTime_ROC2_pions_data_random_noRF_noAero_cut_all.SetLineColor(kBlue)
Clone_CTime_ePiCoinTime_ROC2_pions_data_random_noRF_noAero_cut_all.SetLineWidth(2)
Clone_CTime_ePiCoinTime_ROC2_pions_data_random_noRF_noAero_cut_all.Draw("Hist Same")
c1_PID1.cd(4)
CTime_ePiCoinTime_ROC2_pions_data_dummysub_noRF_Aero_cut_all.SetLineColor(kRed)
CTime_ePiCoinTime_ROC2_pions_data_dummysub_noRF_Aero_cut_all.SetLineWidth(2)
CTime_ePiCoinTime_ROC2_pions_data_dummysub_noRF_Aero_cut_all.Draw("Hist")
Clone_CTime_ePiCoinTime_ROC2_pions_data_random_noRF_Aero_cut_all.SetLineColor(kRed)
Clone_CTime_ePiCoinTime_ROC2_pions_data_random_noRF_Aero_cut_all.SetLineWidth(2)
Clone_CTime_ePiCoinTime_ROC2_pions_data_random_noRF_Aero_cut_all.Draw("Hist Same")
c1_PID1.cd(5)
P_aero_npeSum_pions_data_dummysub_noRF_noAero_cut_all.SetLineColor(kBlue)
P_aero_npeSum_pions_data_dummysub_noRF_noAero_cut_all.SetLineWidth(2)
P_aero_npeSum_pions_data_dummysub_noRF_noAero_cut_all.Draw()
c1_PID1.cd(6)
P_aero_npeSum_pions_data_dummysub_noRF_Aero_cut_all.SetLineColor(kRed)
P_aero_npeSum_pions_data_dummysub_noRF_Aero_cut_all.SetLineWidth(2)
P_aero_npeSum_pions_data_dummysub_noRF_Aero_cut_all.Draw()
c1_PID1.Print(PID_output_pdf_file + "(")

c1_PID3 = TCanvas("c1_PID3", "PID Efficiency", 100, 0, 1200, 1200)
c1_PID3.Divide(2,2)
c1_PID3.cd(1)
P_kin_MMpi_vs_aero_npeSum_pions_data_dummysub_noRF_noAero_cut_all.Draw("COLZ")
c1_PID3.cd(2)
P_kin_MMpi_vs_aero_npeSum_pions_data_dummysub_noRF_Aero_cut_all.Draw("COLZ")
c1_PID3.cd(3)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_dummysub_noRF_noAero_cut_all.Draw("COLZ")
c1_PID3.cd(4)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_dummysub_noRF_Aero_cut_all.Draw("COLZ")
c1_PID3.Print(PID_output_pdf_file)

c1_PID2 = TCanvas("c1_PID2", "PID Efficiency", 100, 0, 1200, 1800)
c1_PID2.Divide(2,3)
c1_PID2.cd(1)
P_kin_MMpi_pions_data_dummysub_RF_noAero_cut_all.SetLineColor(kBlue)
P_kin_MMpi_pions_data_dummysub_RF_noAero_cut_all.SetLineWidth(2)
P_kin_MMpi_pions_data_dummysub_RF_noAero_cut_all.Draw("Hist")
c1_PID2.cd(2)
P_kin_MMpi_pions_data_dummysub_RF_Aero_cut_all.SetLineColor(kRed)
P_kin_MMpi_pions_data_dummysub_RF_Aero_cut_all.SetLineWidth(2)
P_kin_MMpi_pions_data_dummysub_RF_Aero_cut_all.Draw("Hist")
c1_PID2.cd(3)
CTime_ePiCoinTime_ROC2_pions_data_dummysub_RF_noAero_cut_all.SetLineColor(kBlue)
CTime_ePiCoinTime_ROC2_pions_data_dummysub_RF_noAero_cut_all.SetLineWidth(2)
CTime_ePiCoinTime_ROC2_pions_data_dummysub_RF_noAero_cut_all.Draw("Hist")
Clone_CTime_ePiCoinTime_ROC2_pions_data_random_RF_noAero_cut_all.SetLineColor(kBlue)
Clone_CTime_ePiCoinTime_ROC2_pions_data_random_RF_noAero_cut_all.SetLineWidth(2)
Clone_CTime_ePiCoinTime_ROC2_pions_data_random_RF_noAero_cut_all.Draw("Hist Same")
c1_PID2.cd(4)
CTime_ePiCoinTime_ROC2_pions_data_dummysub_RF_Aero_cut_all.SetLineColor(kRed)
CTime_ePiCoinTime_ROC2_pions_data_dummysub_RF_Aero_cut_all.SetLineWidth(2)
CTime_ePiCoinTime_ROC2_pions_data_dummysub_RF_Aero_cut_all.Draw("Hist")
Clone_CTime_ePiCoinTime_ROC2_pions_data_random_RF_Aero_cut_all.SetLineColor(kRed)
Clone_CTime_ePiCoinTime_ROC2_pions_data_random_RF_Aero_cut_all.SetLineWidth(2)
Clone_CTime_ePiCoinTime_ROC2_pions_data_random_RF_Aero_cut_all.Draw("Hist Same")
c1_PID2.cd(5)
P_aero_npeSum_pions_data_dummysub_RF_noAero_cut_all.SetLineColor(kBlue)
P_aero_npeSum_pions_data_dummysub_RF_noAero_cut_all.SetLineWidth(2)
P_aero_npeSum_pions_data_dummysub_RF_noAero_cut_all.Draw()
c1_PID2.cd(6)
P_aero_npeSum_pions_data_dummysub_RF_Aero_cut_all.SetLineColor(kRed)
P_aero_npeSum_pions_data_dummysub_RF_Aero_cut_all.SetLineWidth(2)
P_aero_npeSum_pions_data_dummysub_RF_Aero_cut_all.Draw()
c1_PID2.Print(PID_output_pdf_file)

c1_PID4 = TCanvas("c1_PID4", "PID Efficiency", 100, 0, 1200, 1200)
c1_PID4.Divide(2,2)
c1_PID4.cd(1)
P_kin_MMpi_vs_aero_npeSum_pions_data_dummysub_RF_noAero_cut_all.Draw("COLZ")
c1_PID4.cd(2)
P_kin_MMpi_vs_aero_npeSum_pions_data_dummysub_RF_noAero_cut_all.Draw("COLZ")
c1_PID4.cd(3)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_dummysub_RF_Aero_cut_all.Draw("COLZ")
c1_PID4.cd(4)
CTime_ePiCoinTime_ROC2_vs_aero_npeSum_pions_data_dummysub_RF_Aero_cut_all.Draw("COLZ")
c1_PID4.Print(PID_output_pdf_file + ")")

#############################################################################################################################################################

# Close input files
infile_DATA.Close()
infile_DUMMY.Close()

print ("Processing Complete")