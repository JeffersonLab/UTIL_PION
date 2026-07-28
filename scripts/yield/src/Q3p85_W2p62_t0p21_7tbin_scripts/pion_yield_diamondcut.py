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
from ROOT import TCanvas, TList, TPaveLabel, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TLegend, TGaxis, TLine, TMath, TLatex, TPaveText, TArc, TGraphPolar, TText, TString
from ROOT import kBlack, kCyan, kRed, kGreen, kMagenta, kBlue
from functools import reduce
import math as ma
import csv

##################################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=8:
    print("!!!!! ERROR !!!!!\n Expected 7 arguments\n Usage is with - PHY_SETTING MaxEvents Suffix RunList CVSFile\n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################

# Input params - run number and max number of events
PHY_SETTING = sys.argv[1]
MaxEvent = sys.argv[2]
DATA_Suffix = sys.argv[3]
DUMMY_Suffix = sys.argv[4]
SIMC_Suffix = sys.argv[5]
DATA_RUN_LIST = sys.argv[6]
DUMMY_RUN_LIST = sys.argv[7]
CSV_FILE = sys.argv[8]

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
RUNLISTPATH = "/u/group/c-pionlt/USERS/%s/hallc_replay_lt/UTIL_BATCH/InputRunLists/PionLT_2021_2022" % (USER)
MMCUT_CSV   = "/u/group/c-pionlt/USERS/%s/hallc_replay_lt/UTIL_PION/LTSep_CSVs/mm_offset_cut_csv" % (USER)
EFF_CSV     = "/u/group/c-pionlt/USERS/%s/hallc_replay_lt/UTIL_PION/efficiencies" % (USER)

# Extract the first three words from PHY_SETTING for the CSV file name
setting_name = "_".join(PHY_SETTING.split("_")[:3])
physet_dir_name = "%s_std" % (setting_name)
SIMCPATH = "/volatile/hallc/c-pionlt/%s/OUTPUT/Analysis/SIMC/%s/" % (USER, physet_dir_name)

#################################################################################################################################################

# Output PDF File Name
print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))
Pion_Analysis_Distributions = "%s/%s_%s_ProdCoin_Pion_Analysis_Diamondcut_Distributions.pdf" % (OUTPATH, PHY_SETTING, MaxEvent)

# Input file location and variables taking
rootFile_DATA = "%s/%s_%s_%s.root" % (OUTPATH, PHY_SETTING, MaxEvent, DATA_Suffix)
rootFile_DUMMY = "%s/%s_%s_%s.root" % (OUTPATH, PHY_SETTING, MaxEvent, DUMMY_Suffix)
rootFile_SIMC = "%s/%s.root" % (SIMCPATH, SIMC_Suffix)
data_run_list = "%s/%s" % (RUNLISTPATH, DATA_RUN_LIST)
dummy_run_list = "%s/%s" % (RUNLISTPATH, DUMMY_RUN_LIST)
csv_file = "%s/%s.csv" % (EFF_CSV, CSV_FILE)
mmcut_csv_file = "%s/%s/%s_mm_offsets_cuts_parameters.csv" % (MMCUT_CSV, physet_dir_name, setting_name)

###################################################################################################################################################

# Define the cuts
# SIMC Cuts for Pions Selection
HMS_Acceptance = lambda event: (event.hsdelta >= -8.0) & (event.hsdelta <= 8.0) & (event.hsxpfp >= -0.08) & (event.hsxpfp <= 0.08) & (event.hsypfp >= -0.045) & (event.hsypfp <= 0.045)
SHMS_Acceptance = lambda event: (event.ssdelta >= -10.0) & (event.ssdelta <= 20.0) & (event.ssxpfp >= -0.06) & (event.ssxpfp <= 0.06) & (event.ssypfp >= -0.04) & (event.ssypfp <= 0.04)
#SHMS_Aero_Cut = lambda event: (event.paero_x_det > -55.0) & (event.paero_x_det < 55.0) & (event.paero_y_det > -50) & (event.paero_y_det < 50) # Aerogel tray n = 1.030
SHMS_Aero_Cut = lambda event: (event.paero_x_det > -45.0) & (event.paero_x_det < 45.0) & (event.paero_y_det > -30) & (event.paero_y_det < 30) # Aerogel tray n = 1.011

# Read the MMpi cut values from the CSV file
try:
    with open(mmcut_csv_file, mode='r', newline='') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            if row[csv_reader.fieldnames[0]].strip() == PHY_SETTING:  # Match first column (Physics_Setting)
                MM_Offset = float(row["MM_Offset"].strip())
                MM_Cut_lowvalue = float(row["MM_Cut_low"].strip())
                MM_Cut_highvalue = float(row["MM_Cut_high"].strip())
                break
        else:
            raise ValueError(f"No matching Physics_Setting '{PHY_SETTING}' found in {mmcut_csv_file}.")

except (FileNotFoundError, ValueError, KeyError) as e:
    print(f"Error: {e}")
    sys.exit(1)

# Print the assigned values
#print(f"MMpi_Offset = {MMpi_Offset:.6f}")
#print(f"MMpi_Cut_lowvalue = {MMpi_Cut_lowvalue}")
#print(f"MMpi_Cut_highvalue = {MMpi_Cut_highvalue}")
DATA_MMpi_Cut = lambda event: (MM_Cut_lowvalue <= (event.MMpi + MM_Offset) <= MM_Cut_highvalue)
SIMC_MMpi_Cut = lambda event: (MM_Cut_lowvalue <= event.missmass <= MM_Cut_highvalue)

###############################################################################################################################################

# Section for grabing Prompt/Random selection parameters from PARAM filePARAMPATH = UTILPATH+"/DB/PARAM"
# Function to read run numbers from a given run list file
def read_run_numbers(run_list_file):
    try:
        with open(run_list_file, 'r') as f:
            return [int(line.strip()) for line in f if line.strip().isdigit()]
    except FileNotFoundError:
        print(f"!!!!! ERROR !!!!!\nRun list file not found: {run_list_file}\n!!!!! ERROR !!!!!")
        sys.exit(2)

# Read run numbers from the specified files
data_run_numbers = read_run_numbers(data_run_list)
dummy_run_numbers = read_run_numbers(dummy_run_list)

# Section for grabbing Prompt/Random selection parameters from PARAM file
PARAMPATH = UTILPATH+"/DB/PARAM/PionLT"
TimingCutFile = "%s/PionLT_Timing_Parameters.csv" % PARAMPATH # This should match the param file actually being used!

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
        array_list = line.split(",") # Convert line into an array_list, anything after a comma is a new entry

        # Check if any of the run numbers are within the specified range
        if any(run_num in range(int(array_list[0]), int(array_list[1]) + 1) for run_num in data_run_numbers + dummy_run_numbers):
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

###############################################################################################################################################

print ('\nPhysics Setting = ',PHY_SETTING, '\n')
print("-"*40)

# Grabs simc number of events and normalizaton factor
simc_hist = "%s/%s.hist" % (SIMCPATH, SIMC_Suffix)

f_simc = open(simc_hist)
for line in f_simc:
#    print(line)
    if "Ncontribute" in line:
        val = line.split("=")
        simc_nevents = int(val[1])
    if "normfac" in line:
        val = line.split("=")
        simc_normfactor = float(val[1])
if 'simc_nevents' and 'simc_normfactor' in locals():
    print('\nsimc_nevents = ',simc_nevents,'\nsimc_normfactor = ',simc_normfactor,'\n')
#    print('\n\ndata_charge = {:.4f} +/- {:.4f}'.format(data_charge, data_charge*eff_errProp_data),'\ndummy_charge = {:.4f} +/- {:.4f}'.format(dummy_charge, dummy_charge*eff_errProp_dummy),'\n\n')
else:
    print("ERROR: Invalid simc hist file %s" % simc_hist)
    sys.exit(1)
f_simc.close()
print("-"*40)

# Normalization factor Calculation for SIMC
normfac_simc = (simc_normfactor)/(simc_nevents)

print ("normfac_simc: ", normfac_simc)
print("-"*40)

###################################################################################################################################################
nbins = 200

# Defining Histograms for Pions
# Histograms having Cuts (Acceptance + PID + RF)
W_pions_data_cut_all = ROOT.TH1D("W_pions_data_cut_all", "W Distribution; W; Counts", nbins, 0, 10.0)
Q2_pions_data_cut_all = ROOT.TH1D("Q2_pions_data_cut_all", "Q2 Distribution; Q^{2}; Counts", nbins, 0, 10.0)
Q2_vs_W_pions_data_cut_all = ROOT.TH2D("Q2_vs_W_pions_data_cut_all", "Q^{2} vs W Distribution; Q^{2}; W", nbins, 3.0, 5.0, nbins, 2.0, 3.0)

# Histograms Having Cuts (Acceptance + PID + RF + Prompt Selection) Data
W_pions_data_prompt_cut_all = ROOT.TH1D("W_pions_data_prompt_cut_all", "W Distribution; W; Counts", nbins, 0, 10.0)
Q2_pions_data_prompt_cut_all = ROOT.TH1D("Q2_pions_data_prompt_cut_all", "Q2 Distribution; Q^{2}; Counts", nbins, 0, 10.0)
Q2_vs_W_pions_data_prompt_cut_all = ROOT.TH2D("Q2_vs_W_pions_data_prompt_cut_all", "Q^{2} vs W Distribution; Q^{2}; W", nbins, 3.0, 5.0, nbins, 2.0, 3.0)

# Histograms Having Cuts (Acceptance + PID + RF + Random Selection) Data
W_pions_data_random_cut_all = ROOT.TH1D("W_pions_data_random_cut_all", "W Distribution; W; Counts", nbins, 0, 10.0)
Q2_pions_data_random_cut_all = ROOT.TH1D("Q2_pions_data_random_cut_all", "Q2 Distribution; Q^{2}; Counts", nbins, 0, 10.0)
Q2_vs_W_pions_data_random_cut_all = ROOT.TH2D("Q2_vs_W_pions_data_random_cut_all", "Q^{2} vs W Distribution; Q^{2}; W", nbins, 3.0, 5.0, nbins, 2.0, 3.0)

# Histograms Having Cuts (Acceptance + PID + RF + Prompt Selection) Dummy
W_pions_dummy_prompt_cut_all = ROOT.TH1D("W_pions_dummy_prompt_cut_all", "W Distribution; W; Counts", nbins, 0, 10.0)
Q2_pions_dummy_prompt_cut_all = ROOT.TH1D("Q2_pions_dummy_prompt_cut_all", "Q2 Distribution; Q^{2}; Counts", nbins, 0, 10.0)
Q2_vs_W_pions_dummy_prompt_cut_all = ROOT.TH2D("Q2_vs_W_pions_dummy_prompt_cut_all", "Q^{2} vs W Distribution; Q^{2}; W", nbins, 3.0, 5.0, nbins, 2.0, 3.0)

# Histograms Having Cuts (Acceptance + PID + RF + Random Selection) Dummy
W_pions_dummy_random_cut_all = ROOT.TH1D("W_pions_dummy_random_cut_all", "W Distribution; W; Counts", nbins, 0, 10.0)
Q2_pions_dummy_random_cut_all = ROOT.TH1D("Q2_pions_dummy_random_cut_all", "Q2 Distribution; Q^{2}; Counts", nbins, 0, 10.0)
Q2_vs_W_pions_dummy_random_cut_all = ROOT.TH2D("Q2_vs_W_pions_dummy_random_cut_all", "Q^{2} vs W Distribution; Q^{2}; W", nbins, 3.0, 5.0, nbins, 2.0, 3.0)

# Histograms having Cuts (Acceptance + PID + RF + RandSub) Data
W_pions_randsub_data_cut_all = ROOT.TH1D("W_pions_randsub_data_cut_all", "W Distribution; W; Counts", nbins, 0, 10.0)
Q2_pions_randsub_data_cut_all = ROOT.TH1D("Q2_pions_randsub_data_cut_all", "Q2 Distribution; Q^{2}; Counts", nbins, 0, 10.0)
Q2_vs_W_pions_randsub_data_cut_all = ROOT.TH2D("Q2_vs_W_pions_randsub_data_cut_all", "Q^{2} vs W Distribution; Q^{2}; W", nbins, 3.0, 5.0, nbins, 2.0, 3.0)

# Histograms having Cuts (Acceptance + PID + RF + RandSub) Dummy
W_pions_randsub_dummy_cut_all = ROOT.TH1D("W_pions_randsub_dummy_cut_all", "W Distribution; W; Counts", nbins, 0, 10.0)
Q2_pions_randsub_dummy_cut_all = ROOT.TH1D("Q2_pions_randsub_dummy_cut_all", "Q2 Distribution; Q^{2}; Counts", nbins, 0, 10.0)
Q2_vs_W_pions_randsub_dummy_cut_all = ROOT.TH2D("Q2_vs_W_pions_randsub_dummy_cut_all", "Q^{2} vs W Distribution; Q^{2}; W", nbins, 3.0, 5.0, nbins, 2.0, 3.0)

# Histograms having Cuts (Acceptance + PID + RF + RandSub + DummySub)
W_pions_dummysub_data_cut_all = ROOT.TH1D("W_pions_dummysub_data_cut_all", "W Distribution; W; Counts", nbins, 0, 10.0)
Q2_pions_dummysub_data_cut_all = ROOT.TH1D("Q2_pions_dummysub_data_cut_all", "Q2 Distribution; Q^{2}; Counts", nbins, 0, 10.0)
Q2_vs_W_pions_dummysub_data_cut_all = ROOT.TH2D("Q2_vs_W_pions_dummysub_data_cut_all", "Q^{2} vs W (Accpt+PID+RF+RandDummySub) Distribution; Q^{2}; W", nbins, 3.0, 5.0, nbins, 2.0, 3.0)

# SIMC Histograms with Cuts
Q2_pions_simc_cut_all = ROOT.TH1D("Q2_pions_simc_cut_all", "#Q^2 Distribution; #Q^2; Counts", nbins, 0, 10)
W_pions_simc_cut_all = ROOT.TH1D("W_pions_simc_cut_all", "W Distribution; W; Counts", nbins, 0, 10.0)
Q2_vs_W_pions_simc_cut_all = ROOT.TH2D("Q2_vs_W_pions_simc_cut_all", "Q^{2} vs W SIMC Distribution; Q^{2}; W", nbins, 3.0, 5.0, nbins, 2.0, 3.0)

##########################################################################################################################################################################################################

# Read stuff from the main event tree
infile_DATA = ROOT.TFile.Open(rootFile_DATA, "READ")
infile_DUMMY = ROOT.TFile.Open(rootFile_DUMMY, "READ")
infile_SIMC = ROOT.TFile.Open(rootFile_SIMC, "READ")

#Uncut_Pion_Events_Data_tree = infile_DATA.Get("Uncut_Pion_Events")
#Cut_Pion_Events_Accpt_Data_tree = infile_DATA.Get("Cut_Pion_Events_Accpt")
Cut_Pion_Events_All_Data_tree = infile_DATA.Get("Cut_Pion_Events_All")
Cut_Pion_Events_Prompt_Data_tree = infile_DATA.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Random_Data_tree = infile_DATA.Get("Cut_Pion_Events_Random")
#nEntries_TBRANCH_DATA  = Cut_Pion_Events_Prompt_Data_tree.GetEntries()

#Uncut_Pion_Events_Dummy_tree = infile_DUMMY.Get("Uncut_Pion_Events")
#Cut_Pion_Events_Accpt_Dummy_tree = infile_DUMMY.Get("Cut_Pion_Events_Accpt")
#Cut_Pion_Events_All_Dummy_tree = infile_DUMMY.Get("Cut_Pion_Events_All")
Cut_Pion_Events_Prompt_Dummy_tree = infile_DUMMY.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Random_Dummy_tree = infile_DUMMY.Get("Cut_Pion_Events_Random")
#nEntries_TBRANCH_DUMMY  = Cut_Pion_Events_Prompt_Dummy_tree.GetEntries()

Uncut_Pion_Events_SIMC_tree = infile_SIMC.Get("h10")
#nEntries_TBRANCH_SIMC  = Uncut_Proton_Events_SIMC_tree.GetEntries()

###################################################################################################################################################

#Fill histograms for Cut All Data
for event in Cut_Pion_Events_All_Data_tree:
    if DATA_MMpi_Cut(event):
        W_pions_data_cut_all.Fill(event.W)
        Q2_pions_data_cut_all.Fill(event.Q2)
        Q2_vs_W_pions_data_cut_all.Fill(event.Q2, event.W)
#    ibin += 1

#Fill histograms for Prompt Data
for event in Cut_Pion_Events_Prompt_Data_tree:
    if DATA_MMpi_Cut(event):
        W_pions_data_prompt_cut_all.Fill(event.W)
        Q2_pions_data_prompt_cut_all.Fill(event.Q2)
        Q2_vs_W_pions_data_prompt_cut_all.Fill(event.Q2, event.W)
#    ibin += 1

#Fill histograms for Random Data
for event in Cut_Pion_Events_Random_Data_tree:
    if DATA_MMpi_Cut(event):
        W_pions_data_random_cut_all.Fill(event.W)
        Q2_pions_data_random_cut_all.Fill(event.Q2)
        Q2_vs_W_pions_data_random_cut_all.Fill(event.Q2, event.W)
#    ibin += 1

# Fill histograms from Prompt Dummy
#ibin = 1
for event in Cut_Pion_Events_Prompt_Dummy_tree:
    if DATA_MMpi_Cut(event):
        W_pions_dummy_prompt_cut_all.Fill(event.W)
        Q2_pions_dummy_prompt_cut_all.Fill(event.Q2)
        Q2_vs_W_pions_dummy_prompt_cut_all.Fill(event.Q2, event.W)
#    ibin += 1

# Fill histograms from Random Dummy
#ibin = 1
for event in Cut_Pion_Events_Random_Dummy_tree:
    if DATA_MMpi_Cut(event):
        W_pions_dummy_random_cut_all.Fill(event.W)
        Q2_pions_dummy_random_cut_all.Fill(event.Q2)
        Q2_vs_W_pions_dummy_random_cut_all.Fill(event.Q2, event.W)
#    ibin += 1

# Fill histograms from SIMC ROOT File
for event in Uncut_Pion_Events_SIMC_tree:
    if HMS_Acceptance(event) & SHMS_Acceptance(event) & SHMS_Aero_Cut(event) & SIMC_MMpi_Cut(event):        
            Q2_pions_simc_cut_all.Fill(event.Q2, event.Weight)
            W_pions_simc_cut_all.Fill(event.W, event.Weight)
            Q2_vs_W_pions_simc_cut_all.Fill(event.Q2, event.W, event.Weight)

print("####################################")
print("###### Histogram filling done ######")
print("####################################\n")

#################################################################################################################################################

# Random subtraction from Histograms
W_pions_data_random_cut_all.Scale(1.0/nWindows)
Q2_pions_data_random_cut_all.Scale(1.0/nWindows)
Q2_vs_W_pions_data_random_cut_all.Scale(1.0/nWindows)
W_pions_randsub_data_cut_all.Add(W_pions_data_prompt_cut_all, W_pions_data_random_cut_all, 1, -1)
Q2_pions_randsub_data_cut_all.Add(Q2_pions_data_prompt_cut_all, Q2_pions_data_random_cut_all, 1, -1)
Q2_vs_W_pions_data_cut_all.Add(Q2_vs_W_pions_data_prompt_cut_all, Q2_vs_W_pions_data_random_cut_all, 1, -1)

W_pions_dummy_random_cut_all.Scale(1.0/nWindows)
Q2_pions_dummy_random_cut_all.Scale(1.0/nWindows)
Q2_vs_W_pions_dummy_random_cut_all.Scale(1.0/nWindows)
W_pions_randsub_dummy_cut_all.Add(W_pions_dummy_prompt_cut_all, W_pions_dummy_random_cut_all, 1, -1)
Q2_pions_randsub_dummy_cut_all.Add(Q2_pions_dummy_prompt_cut_all, Q2_pions_dummy_random_cut_all, 1, -1)
Q2_vs_W_pions_randsub_dummy_cut_all.Add(Q2_vs_W_pions_dummy_prompt_cut_all, Q2_vs_W_pions_dummy_random_cut_all, 1, -1)

print("######################################")
print("###### Random substraction done ######")
print("######################################\n")

############################################################################################################################################

# # SIMC Normalization
Q2_pions_simc_cut_all.Scale(normfac_simc)
W_pions_simc_cut_all.Scale(normfac_simc)
Q2_vs_W_pions_simc_cut_all.Scale(normfac_simc)

############################################################################################################################################

# Dummy Subtraction
W_pions_dummysub_data_cut_all.Add(W_pions_randsub_data_cut_all, W_pions_randsub_dummy_cut_all, 1, -1)
Q2_pions_dummysub_data_cut_all.Add(Q2_pions_randsub_data_cut_all, Q2_pions_randsub_dummy_cut_all, 1, -1)
Q2_vs_W_pions_dummysub_data_cut_all.Add(Q2_vs_W_pions_data_cut_all, Q2_vs_W_pions_randsub_dummy_cut_all, 1, -1)

print("####################################")
print("###### Dummy subtraction done ######")
print("####################################\n")

###########################################################################################################################################

# Section for Diamond Cut

# Get the number of bins in X and Y directions
n_bins_x = Q2_vs_W_pions_dummysub_data_cut_all.GetNbinsX()
n_bins_y = Q2_vs_W_pions_dummysub_data_cut_all.GetNbinsY()

# Initialize variables to store the vertex coordinates
vertex1 = [None, None]  # bottom-left
vertex2 = [None, None]  # top-left
vertex3 = [None, None]  # top-right
vertex4 = [None, None]  # bottom-right

# Project the 2D histogram onto the x and y axes
Q2_projection = Q2_vs_W_pions_dummysub_data_cut_all.ProjectionX()
W_projection = Q2_vs_W_pions_dummysub_data_cut_all.ProjectionY()

# Function to find the start and end bins of a projection
def find_start_end_bins(projection):
    start_bin = None
    end_bin = None
    for bin in range(1, projection.GetNbinsX() + 1):
        if projection.GetBinContent(bin) > 0:
            if start_bin is None:
                start_bin = bin
            end_bin = bin
    return start_bin, end_bin

# Find the start and end bins for the x projection (Q2)
Q2_start_bin, Q2_end_bin = find_start_end_bins(Q2_projection)

# Find the start and end bins for the y projection (W)
W_start_bin, W_end_bin = find_start_end_bins(W_projection)

# Calculate the x and y values where the distribution starts and ends
Q2_start = Q2_projection.GetXaxis().GetBinLowEdge(Q2_start_bin)
Q2_end = Q2_projection.GetXaxis().GetBinUpEdge(Q2_end_bin)
W_start = W_projection.GetXaxis().GetBinLowEdge(W_start_bin)
W_end = W_projection.GetXaxis().GetBinUpEdge(W_end_bin)

# Print the results
print(f"Q2 distribution starts at: {Q2_start}, ends at: {Q2_end}")
print(f"W distribution starts at: {W_start}, ends at: {W_end}")

# Hardcoded vertices for the diamond cut
vertex1 = [3.780, 2.600]  # bottom-left
vertex2 = [3.180, 2.775]  # top-left
vertex3 = [3.902, 2.642]  # top-right
vertex4 = [4.520, 2.452]  # bottom-right

# Print vertices
print("Vertices of the populated area:")
print(f"Vertex 1 (bottom-left): {vertex1}")
print(f"Vertex 2 (top-left): {vertex2}")
print(f"Vertex 3 (top-right): {vertex3}")
print(f"Vertex 4 (bottom-right): {vertex4}")

###########################################################################################################################################

# Diamond Cut Test
Q2_vs_W_pions_diamond_cut_all = ROOT.TH2D("Q2_vs_W_pions_diamond_cut_all", "Q^{2} vs W Diamond Cut Distribution; Q^{2}; W", nbins, 3.0, 5.0, nbins, 2.0, 3.0)
# Define the diamond cut
cutg_diamond_simc = ROOT.TCutG("cutg_diamond_simc", 5)
cutg_diamond_simc.SetVarX("Q2")
cutg_diamond_simc.SetVarY("W")
cutg_diamond_simc.SetPoint(0, vertex1[0], vertex1[1])  # bottom-left
cutg_diamond_simc.SetPoint(1, vertex2[0], vertex2[1])  # top-left
cutg_diamond_simc.SetPoint(2, vertex3[0], vertex3[1])  # top-right
cutg_diamond_simc.SetPoint(3, vertex4[0], vertex4[1])  # bottom-right
cutg_diamond_simc.SetPoint(4, vertex1[0], vertex1[1])  # bottom-left again to close the loop
Diamond_Cut = lambda event: (cutg_diamond_simc.IsInside(event.Q2, event.W))

# Fill the histogram with the diamond cut
for event in Cut_Pion_Events_Prompt_Data_tree:
    if DATA_MMpi_Cut(event) & Diamond_Cut(event):
        Q2_vs_W_pions_diamond_cut_all.Fill(event.Q2, event.W)

###########################################################################################################################################

# Removes stat box
ROOT.gStyle.SetOptStat(0)

# Saving histograms in PDF
c1_delta = TCanvas("c1_delta", "Variables Distributions", 100, 0, 1400, 1400)
c1_delta.Divide(2,2)
c1_delta.cd(1)
c1_delta_text_lines = [
    ROOT.TText(0.5, 0.9, "Pion ProdCoin Setting"),
    ROOT.TText(0.5, 0.8, "{}".format(PHY_SETTING)),
    ROOT.TText(0.5, 0.7, "Bottom-left Vertex: ({:.3f}, {:.3f})".format(vertex1[0], vertex1[1])),
    ROOT.TText(0.5, 0.6, "Top-left Vertex: ({:.3f}, {:.3f})".format(vertex2[0], vertex2[1])),
    ROOT.TText(0.5, 0.5, "Top-right Vertex: ({:.3f}, {:.3f})".format(vertex3[0], vertex3[1])),
    ROOT.TText(0.5, 0.4, "Bottom-right Vertex: ({:.3f}, {:.3f})".format(vertex4[0], vertex4[1]))
]
for c1_delta_text in c1_delta_text_lines:
    c1_delta_text.SetTextSize(0.04)
    c1_delta_text.SetTextAlign(22)
    c1_delta_text.SetTextColor(ROOT.kBlack)
    c1_delta_text.Draw()
c1_delta.cd(2)
ROOT.gPad.SetLogz()
Q2_vs_W_pions_diamond_cut_all.GetXaxis().SetRangeUser(3.0, 4.8)
Q2_vs_W_pions_diamond_cut_all.GetYaxis().SetRangeUser(2.4, 2.8)
Q2_vs_W_pions_diamond_cut_all.Draw("colz")
c1_delta.cd(3)
ROOT.gPad.SetLogz()
Q2_vs_W_pions_dummysub_data_cut_all.GetXaxis().SetRangeUser(3.0, 4.8)
Q2_vs_W_pions_dummysub_data_cut_all.GetYaxis().SetRangeUser(2.4, 2.8)
Q2_vs_W_pions_dummysub_data_cut_all.Draw("colz")
# Draw lines to visualize the diamond cuts
line1 = ROOT.TLine(vertex1[0], vertex1[1], vertex2[0], vertex2[1])
line2 = ROOT.TLine(vertex2[0], vertex2[1], vertex3[0], vertex3[1])
line3 = ROOT.TLine(vertex3[0], vertex3[1], vertex4[0], vertex4[1])
line4 = ROOT.TLine(vertex4[0], vertex4[1], vertex1[0], vertex1[1])
# Set line styles
line1.SetLineColor(ROOT.kRed)
line2.SetLineColor(ROOT.kRed)
line3.SetLineColor(ROOT.kRed)
line4.SetLineColor(ROOT.kRed)
line1.SetLineWidth(2)
line2.SetLineWidth(2)
line3.SetLineWidth(2)
line4.SetLineWidth(2)
line1.Draw("same")
line2.Draw("same")
line3.Draw("same")
line4.Draw("same")
c1_delta.cd(4)
ROOT.gPad.SetLogz()
Q2_vs_W_pions_simc_cut_all.GetXaxis().SetRangeUser(3.0, 4.8)
Q2_vs_W_pions_simc_cut_all.GetYaxis().SetRangeUser(2.4, 2.8)
Q2_vs_W_pions_simc_cut_all.Draw("colz")
# Draw lines to visualize the diamond cuts
line1_simc = ROOT.TLine(vertex1[0], vertex1[1], vertex2[0], vertex2[1])
line2_simc = ROOT.TLine(vertex2[0], vertex2[1], vertex3[0], vertex3[1])
line3_simc = ROOT.TLine(vertex3[0], vertex3[1], vertex4[0], vertex4[1])
line4_simc = ROOT.TLine(vertex4[0], vertex4[1], vertex1[0], vertex1[1])
# Set line styles
line1_simc.SetLineColor(ROOT.kRed)
line2_simc.SetLineColor(ROOT.kRed)
line3_simc.SetLineColor(ROOT.kRed)
line4_simc.SetLineColor(ROOT.kRed)
line1_simc.SetLineWidth(2)
line2_simc.SetLineWidth(2)
line3_simc.SetLineWidth(2)
line4_simc.SetLineWidth(2)
line1_simc.Draw("same")
line2_simc.Draw("same")
line3_simc.Draw("same")
line4_simc.Draw("same")
c1_delta.Print(Pion_Analysis_Distributions)

#############################################################################################################################################

# Writing Offsets and Cuts to CSV file
csv_output_path = "%s/LTSep_CSVs/diamond_cut_csv/%s/%s_diamond_cut_parameters.csv" % (UTILPATH, physet_dir_name, setting_name)

# Use the already defined vertices
vertices = {
    "vertex1": vertex1,  # bottom-left
    "vertex2": vertex2,  # top-left
    "vertex3": vertex3,  # top-right
    "vertex4": vertex4,  # bottom-right
}

# Write the vertices to the CSV file
with open(csv_output_path, mode='w', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)
    # Write the header
    csv_writer.writerow(["Vertex", "X", "Y"])
    # Write the vertex data
    for vertex, coordinates in vertices.items():
        csv_writer.writerow([vertex, f"{coordinates[0]:.3f}", f"{coordinates[1]:.3f}"])

print(f"Vertices saved to {csv_output_path}")

#############################################################################################################################################

# Making directories in output file
outHistFile = ROOT.TFile.Open("%s/%s_%s_Diamondcut_Data.root" % (OUTPATH, PHY_SETTING, MaxEvent) , "RECREATE")
d_Cut_Pion_Events_Data_Cut_All = outHistFile.mkdir("Cut_Pion_Events_Data_Cut_All")
d_Cut_Pion_Events_Prompt_Data_Cut_All = outHistFile.mkdir("Cut_Pion_Events_Prompt_Data_Cut_All")
d_Cut_Pion_Events_Random_Data_Cut_All = outHistFile.mkdir("Cut_Pion_Events_Random_Data_Cut_All")
d_Cut_Pion_Events_RandomSub_Data = outHistFile.mkdir("Cut_Pion_Events_RandomSub_Data")
d_Cut_Pion_Events_DummySub_RandomSub_Data = outHistFile.mkdir("Cut_Pion_Events_DummySub_RandomSub_Data")
d_Cut_Pion_Events_Norm_SIMC_Data = outHistFile.mkdir("Cut_Pion_Events_Norm_SIMC_Data")

d_Cut_Pion_Events_Data_Cut_All.cd()
W_pions_data_cut_all.Write()
Q2_pions_data_cut_all.Write()
Q2_vs_W_pions_data_cut_all.Write()

d_Cut_Pion_Events_Prompt_Data_Cut_All.cd()
W_pions_data_prompt_cut_all.Write()
Q2_pions_data_prompt_cut_all.Write()
Q2_vs_W_pions_data_prompt_cut_all.Write()

d_Cut_Pion_Events_Random_Data_Cut_All.cd()
W_pions_data_random_cut_all.Write()
Q2_pions_data_random_cut_all.Write()
Q2_vs_W_pions_data_random_cut_all.Write()

d_Cut_Pion_Events_RandomSub_Data.cd()
W_pions_randsub_data_cut_all.Write()
Q2_pions_randsub_data_cut_all.Write()
Q2_vs_W_pions_randsub_data_cut_all.Write()

d_Cut_Pion_Events_DummySub_RandomSub_Data.cd()
W_pions_dummysub_data_cut_all.Write()
Q2_pions_dummysub_data_cut_all.Write()
Q2_vs_W_pions_dummysub_data_cut_all.Write()
Q2_vs_W_pions_diamond_cut_all.Write()

d_Cut_Pion_Events_Norm_SIMC_Data.cd()
Q2_pions_simc_cut_all.Write()
W_pions_simc_cut_all.Write()
Q2_vs_W_pions_simc_cut_all.Write()

print("####################################")
print("###### Histogram writing done ######")
print("####################################\n")

infile_DATA.Close() 
infile_DUMMY.Close()
infile_SIMC.Close()
outHistFile.Close()

print ("Processing Complete")