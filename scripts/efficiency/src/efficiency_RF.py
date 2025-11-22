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
setting_name = "_".join(PHY_SETTING.split("_")[:3])

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
#Pion_Analysis_Distributions = "%s_ProdCoin_Pion_Analysis_scatter_Distributions.pdf" % (PHY_SETTING)

# Input file location and variables taking
rootFile_DATA = "%s/%s.root" % (OUTPATH, DATA_ROOTFILE_SUFFIX)
rootFile_DUMMY = "%s/%s.root" % (OUTPATH, DUMMY_ROOTFILE_SUFFIX)

###################################################################################################################################################

nbins = 200
RF_min = 0.0
RF_max = 4.0

# Defining Histograms for Pions
RF_Eff_Ndid_data_prompt_cut_all = ROOT.TH1D("RF_Eff_Ndid_data_prompt_cut_all", "RF Efficiency Ndid Distribution; RF Time (ns); Counts", nbins, RF_min, RF_max)
RF_Eff_Ndid_data_random_cut_all = ROOT.TH1D("RF_Eff_Ndid_data_random_cut_all", "RF Efficiency Ndid Distribution; RF Time (ns); Counts", nbins, RF_min, RF_max) 
RF_Eff_Nshould_data_prompt_cut_all = ROOT.TH1D("RF_Eff_Nshould_data_prompt_cut_all", "RF Efficiency Nshould Distribution; RF Time (ns); Counts", nbins, RF_min, RF_max)
RF_Eff_Nshould_data_random_cut_all = ROOT.TH1D("RF_Eff_Nshould_data_random_cut_all", "RF Efficiency Nshould Distribution; RF Time (ns); Counts", nbins, RF_min, RF_max)

RF_Eff_Ndid_dummy_prompt_cut_all = ROOT.TH1D("RF_Eff_Ndid_dummy_prompt_cut_all", "RF Efficiency Ndid Distribution; RF Time (ns); Counts", nbins, RF_min, RF_max)
RF_Eff_Ndid_dummy_random_cut_all = ROOT.TH1D("RF_Eff_Ndid_dummy_random_cut_all", "RF Efficiency Ndid Distribution; RF Time (ns); Counts", nbins, RF_min, RF_max) 
RF_Eff_Nshould_dummy_prompt_cut_all = ROOT.TH1D("RF_Eff_Nshould_dummy_prompt_cut_all", "RF Efficiency Nshould Distribution; RF Time (ns); Counts", nbins, RF_min, RF_max)
RF_Eff_Nshould_dummy_random_cut_all = ROOT.TH1D("RF_Eff_Nshould_dummy_random_cut_all", "RF Efficiency Nshould Distribution; RF Time (ns); Counts", nbins, RF_min, RF_max)

Rf_Eff_Ndid_randomsub_data_cut_all = ROOT.TH1D("Rf_Eff_Ndid_randomsub_data_cut_all", "RF Efficiency Ndid Distribution; RF Time (ns); Counts", nbins, RF_min, RF_max)
Rf_Eff_Nshould_randomsub_data_cut_all = ROOT.TH1D("Rf_Eff_Nshould_randomsub_data_cut_all", "RF Efficiency Nshould Distribution; RF Time (ns); Counts", nbins, RF_min, RF_max)
Rf_Eff_Ndid_randomsub_dummy_cut_all = ROOT.TH1D("Rf_Eff_Ndid_randomsub_dummy_cut_all", "RF Efficiency Ndid Distribution; RF Time (ns); Counts", nbins, RF_min, RF_max)
Rf_Eff_Nshould_randomsub_dummy_cut_all = ROOT.TH1D("Rf_Eff_Nshould_randomsub_dummy_cut_all", "RF Efficiency Nshould Distribution; RF Time (ns); Counts", nbins, RF_min, RF_max)

Rf_Eff_Ndid_dummysub_data_cut_all = ROOT.TH1D("Rf_Eff_Ndid_dummysub_data_cut_all", "RF Efficiency Ndid Distribution; RF Time (ns); Counts", nbins, RF_min, RF_max)
Rf_Eff_Nshould_dummysub_data_cut_all = ROOT.TH1D("Rf_Eff_Nshould_dummysub_data_cut_all", "RF Efficiency Nshould Distribution; RF Time (ns); Counts", nbins, RF_min, RF_max)

##########################################################################################################################################################################################################

# Read stuff from the main event tree
infile_DATA = ROOT.TFile.Open(rootFile_DATA, "READ")
infile_DUMMY = ROOT.TFile.Open(rootFile_DUMMY, "READ")

# Grab the trees
Cut_Pion_Events_Data_tree = infile_DATA.Get("Cut_Pion_Events_Accpt")
Cut_Pion_Events_Dummy_tree = infile_DUMMY.Get("Cut_Pion_Events_Accpt")

###################################################################################################################################################

# Grab runlist
runlist_file = "%s/UTIL_BATCH/InputRunLists/PionLT_2021_2022/%s_loweps_center" % (REPLAYPATH, setting_name)

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

# RF cut
RF_Cut_lowvalue = 1.2
RF_Cut_highvalue = 3.4
RF_Cut = lambda event: (RF_Cut_lowvalue <= event.P_RF_Dist <= RF_Cut_highvalue)

# HMS Detector Cuts
HMS_Calo_cut = 0.7
HMS_Cer_cut = 1.5
HMS_Detector_Cut = lambda event: (event.H_cal_etottracknorm > HMS_Calo_cut and event.H_cer_npeSum > HMS_Cer_cut)

# SHMS Detector Cuts
SHMS_HGC_Cut = 1.5
SHMS_Aero_cut = 2.5
SHMS_Detector_Cut = lambda event: (event.P_hgcer_npeSum > SHMS_HGC_Cut and event.P_aero_npeSum > SHMS_Aero_cut)   

#Fill histograms for Cut All Data
for event in Cut_Pion_Events_Data_tree:
    if HMS_Detector_Cut(event) and SHMS_Detector_Cut(event) and CoinTime_Prompt_Cut(event) and RF_Cut(event):
        RF_Eff_Ndid_data_prompt_cut_all.Fill(event.P_RF_Dist)
for event in Cut_Pion_Events_Data_tree:
    if HMS_Detector_Cut(event) and SHMS_Detector_Cut(event) and CoinTime_Random_Cut(event) and RF_Cut(event):
        RF_Eff_Ndid_data_random_cut_all.Fill(event.P_RF_Dist)

for event in Cut_Pion_Events_Data_tree:
    if HMS_Detector_Cut(event) and SHMS_Detector_Cut(event) and CoinTime_Prompt_Cut(event):
        RF_Eff_Nshould_data_prompt_cut_all.Fill(event.P_RF_Dist)
for event in Cut_Pion_Events_Data_tree:
    if HMS_Detector_Cut(event) and SHMS_Detector_Cut(event) and CoinTime_Random_Cut(event):
        RF_Eff_Nshould_data_random_cut_all.Fill(event.P_RF_Dist)

for event in Cut_Pion_Events_Dummy_tree:
    if HMS_Detector_Cut(event) and SHMS_Detector_Cut(event) and CoinTime_Prompt_Cut(event) and RF_Cut(event):
        RF_Eff_Ndid_dummy_prompt_cut_all.Fill(event.P_RF_Dist)
for event in Cut_Pion_Events_Dummy_tree:
    if HMS_Detector_Cut(event) and SHMS_Detector_Cut(event) and CoinTime_Random_Cut(event) and RF_Cut(event):
        RF_Eff_Ndid_dummy_random_cut_all.Fill(event.P_RF_Dist)

for event in Cut_Pion_Events_Dummy_tree:
    if HMS_Detector_Cut(event) and SHMS_Detector_Cut(event) and CoinTime_Prompt_Cut(event):
        RF_Eff_Nshould_dummy_prompt_cut_all.Fill(event.P_RF_Dist)
for event in Cut_Pion_Events_Dummy_tree:
    if HMS_Detector_Cut(event) and SHMS_Detector_Cut(event) and CoinTime_Random_Cut(event):
        RF_Eff_Nshould_dummy_random_cut_all.Fill(event.P_RF_Dist)

print("####################################")
print("###### Histogram filling done ######")
print("####################################\n")

#################################################################################################################################################

# Random subtraction from Histograms
RF_Eff_Ndid_data_random_cut_all.Scale(1.0/nWindows)
RF_Eff_Nshould_data_random_cut_all.Scale(1.0/nWindows)
RF_Eff_Ndid_dummy_random_cut_all.Scale(1.0/nWindows)
RF_Eff_Nshould_dummy_random_cut_all.Scale(1.0/nWindows)

Rf_Eff_Ndid_randomsub_data_cut_all.Add(RF_Eff_Ndid_data_prompt_cut_all, RF_Eff_Ndid_data_random_cut_all, 1, -1)
Rf_Eff_Nshould_randomsub_data_cut_all.Add(RF_Eff_Nshould_data_prompt_cut_all, RF_Eff_Nshould_data_random_cut_all, 1, -1)
Rf_Eff_Ndid_randomsub_dummy_cut_all.Add(RF_Eff_Ndid_dummy_prompt_cut_all, RF_Eff_Ndid_dummy_random_cut_all, 1, -1)
Rf_Eff_Nshould_randomsub_dummy_cut_all.Add(RF_Eff_Nshould_dummy_prompt_cut_all, RF_Eff_Nshould_dummy_random_cut_all, 1, -1)

print("######################################")
print("###### Random substraction done ######")
print("######################################\n")

############################################################################################################################################

# Dummy Subtraction
Rf_Eff_Ndid_dummysub_data_cut_all.Add(Rf_Eff_Ndid_randomsub_data_cut_all, Rf_Eff_Ndid_randomsub_dummy_cut_all, 1, -1)
Rf_Eff_Nshould_dummysub_data_cut_all.Add(Rf_Eff_Nshould_randomsub_data_cut_all, Rf_Eff_Nshould_randomsub_dummy_cut_all, 1, -1)

print("####################################")
print("###### Dummy subtraction done ######")
print("####################################\n")

#############################################################################################################################################################
# RF efficiency calculation
RF_Ndid = Rf_Eff_Ndid_dummysub_data_cut_all.Integral()
RF_Nshould = Rf_Eff_Nshould_dummysub_data_cut_all.Integral()

RF_Eff = RF_Ndid / RF_Nshould

RF_EFF_Err = ma.sqrt(((RF_Nshould * RF_Ndid) - (RF_Ndid**2)) / (RF_Nshould**3))

print("="*50)
print("###### RF Efficiency Calculation #########")
print("="*50)
print("RF Ndid: %.2f" % RF_Ndid)
print("RF Nshould: %.2f" % RF_Nshould)
print("RF Efficiency: %.5f +/- %.5f" % (RF_Eff, RF_EFF_Err))
print("="*50)

# Writing output to CSV file
csv_RF_Eff = os.path.join(UTILPATH, "efficiencies", "%s_PionLT_coin_prod_RF_efficiency_data.csv" % setting_name)
try:
    with open(csv_RF_Eff, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["physics setting", "RF_eff", "RF_eff_error"])
        writer.writerow([PHY_SETTING, "{:.5f}".format(RF_Eff), "{:.5f}".format(RF_EFF_Err)])
    print("Wrote RF efficiency to %s" % csv_RF_Eff)
except Exception as e:
    print("ERROR writing CSV file: %s" % e)
    
#############################################################################################################################################################

##################################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################
# Save histograms to output PDF
output_pdf_file = os.path.join(UTILPATH, "scripts/efficiency/OUTPUTS/plots", "%s_PionLT_coin_prod_RF_efficiency.pdf" % setting_name)

c1_RF = TCanvas("c1_RF", "RF Efficiency", 100, 0, 800, 1200)
c1_RF.Divide(1,2)
c1_RF.cd(1)
Rf_Eff_Ndid_dummysub_data_cut_all.SetLineColor(kBlue)
Rf_Eff_Ndid_dummysub_data_cut_all.SetLineWidth(2)
Rf_Eff_Ndid_dummysub_data_cut_all.Draw()
c1_RF.cd(2)
Rf_Eff_Nshould_dummysub_data_cut_all.SetLineColor(kRed)
Rf_Eff_Nshould_dummysub_data_cut_all.SetLineWidth(2)
Rf_Eff_Nshould_dummysub_data_cut_all.Draw()
c1_RF.Print(output_pdf_file)

#############################################################################################################################################################

# Close input files
infile_DATA.Close()
infile_DUMMY.Close()

print ("Processing Complete")