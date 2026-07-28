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

##################################################################################################################################################

# Input params - run number and max number of events
PHY_SETTING = sys.argv[1]

# Constructing inout file names
DATA_Suffix_lowepscenter = "{}_loweps_center".format(PHY_SETTING)
DATA_Suffix_lowepsleft = "{}_loweps_left".format(PHY_SETTING)
DATA_Suffix_highepsright = "{}_higheps_right".format(PHY_SETTING)
DATA_Suffix_highepscenter = "{}_higheps_center".format(PHY_SETTING)
DATA_Suffix_highepsleft = "{}_higheps_left".format(PHY_SETTING)

DATA_ROOTFILE_SUFFIX = "-1_ProdCoin_Analysed_Data"
DUMMY_ROOTFILE_SUFFIX = "-1_ProdCoin_Analysed_Dummy_Data"

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
DCUT_CSV    = "/u/group/c-pionlt/USERS/%s/hallc_replay_lt/UTIL_PION/LTSep_CSVs/diamond_cut_csv" % (USER)

#################################################################################################################################################

# Output PDF File Name
print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))
Pion_Analysis_Distributions = "%s_ProdCoin_Pion_Analysis_scatter_Distributions.pdf" % (PHY_SETTING)

# Extract the first three words from PHY_SETTING for the CSV file name
setting_name = "_".join(PHY_SETTING.split("_")[:3])
mmcut_csv_file = "%s/%s_mm_offsets_cuts_parameters.csv" % (MMCUT_CSV, setting_name)
dcut_csv_file = "%s/%s_diamond_cut_parameters.csv" % (DCUT_CSV, setting_name)

# Input file location and variables taking

rootFile_DATA_lowepsleft = "%s/%s_%s.root" % (OUTPATH, DATA_Suffix_lowepsleft, DATA_ROOTFILE_SUFFIX)
rootFile_DATA_lowepscenter = "%s/%s_%s.root" % (OUTPATH, DATA_Suffix_lowepscenter, DATA_ROOTFILE_SUFFIX)
rootFile_DATA_highepsright = "%s/%s_%s.root" % (OUTPATH, DATA_Suffix_highepsright, DATA_ROOTFILE_SUFFIX)
rootFile_DATA_highepscenter = "%s/%s_%s.root" % (OUTPATH, DATA_Suffix_highepscenter, DATA_ROOTFILE_SUFFIX)
rootFile_DATA_highepsleft = "%s/%s_%s.root" % (OUTPATH, DATA_Suffix_highepsleft, DATA_ROOTFILE_SUFFIX)

rootFile_DUMMY_lowepscenter = "%s/%s_%s.root" % (OUTPATH, DATA_Suffix_lowepscenter, DUMMY_ROOTFILE_SUFFIX)
rootFile_DUMMY_lowepsleft = "%s/%s_%s.root" % (OUTPATH, DATA_Suffix_lowepsleft, DUMMY_ROOTFILE_SUFFIX)
rootFile_DUMMY_highepsright = "%s/%s_%s.root" % (OUTPATH, DATA_Suffix_highepsright, DUMMY_ROOTFILE_SUFFIX)
rootFile_DUMMY_highepscenter = "%s/%s_%s.root" % (OUTPATH, DATA_Suffix_highepscenter, DUMMY_ROOTFILE_SUFFIX)
rootFile_DUMMY_highepsleft = "%s/%s_%s.root" % (OUTPATH, DATA_Suffix_highepsleft, DUMMY_ROOTFILE_SUFFIX)

run_list = "%s/%s_Runlist" % (RUNLISTPATH, PHY_SETTING)

###################################################################################################################################################

# Cuts for Pions Selection
# Read the vertices from the CSV file
vertices = {}
try:
    with open(dcut_csv_file, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            vertex_name = row["Vertex"]
            x_value = round(float(row["X"]), 3)  # Round to 3 decimal places
            y_value = round(float(row["Y"]), 3)  # Round to 3 decimal places
            vertices[vertex_name] = [x_value, y_value]
except (FileNotFoundError, ValueError, KeyError) as e:
    print(f"Error reading diamond cut vertices from {dcut_csv_file}: {e}")
    sys.exit(1)

# Assign the vertices
vertex1 = vertices["vertex1"]  # bottom-left
vertex2 = vertices["vertex2"]  # top-left
vertex3 = vertices["vertex3"]  # top-right
vertex4 = vertices["vertex4"]  # bottom-right

# Print the vertex values rounded to 3 decimals
#print("Diamond Cut Vertices (rounded to 3 decimals):")
#print(f"Vertex 1 (bottom-left): [{vertex1[0]:.3f}, {vertex1[1]:.3f}]")
#print(f"Vertex 2 (top-left): [{vertex2[0]:.3f}, {vertex2[1]:.3f}]")
#print(f"Vertex 3 (top-right): [{vertex3[0]:.3f}, {vertex3[1]:.3f}]")
#print(f"Vertex 4 (bottom-right): [{vertex4[0]:.3f}, {vertex4[1]:.3f}]")

# Define the diamond cut
cutg_diamond = ROOT.TCutG("cutg_diamond", 5)
cutg_diamond.SetVarX("Q2")
cutg_diamond.SetVarY("W")
cutg_diamond.SetPoint(0, vertex1[0], vertex1[1])  # bottom-left
cutg_diamond.SetPoint(1, vertex2[0], vertex2[1])  # top-left
cutg_diamond.SetPoint(2, vertex3[0], vertex3[1])  # top-right
cutg_diamond.SetPoint(3, vertex4[0], vertex4[1])  # bottom-right
cutg_diamond.SetPoint(4, vertex1[0], vertex1[1])  # bottom-left again to close the loop

# Read the MMpi cut values from the CSV file
try:
    with open(mmcut_csv_file, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        first_row = next(csv_reader)  # Get the first row
        if first_row:
            MM_Cut_lowvalue = float(first_row["MM_Cut_low"])  # Assign MM_Cut_low
            MM_Cut_highvalue = float(first_row["MM_Cut_high"])  # Assign MM_Cut_high
        else:
            raise ValueError(f"The CSV file {mmcut_csv_file} is empty or missing required columns.")
except (FileNotFoundError, ValueError, KeyError) as e:
    print(f"Error: {e}")
    sys.exit(1)

# Print the assigned values
#print(f"MMpi_Cut_lowvalue = {MMpi_Cut_lowvalue}")
#print(f"MMpi_Cut_highvalue = {MMpi_Cut_highvalue}")

# Read the MMpi offset values from the CSV file
MM_Offset_lowepscenter = None
MM_Offset_lowepsleft = None
MM_Offset_highepsright = None
MM_Offset_highepscenter = None
MM_Offset_highepsleft = None

try:
    with open(mmcut_csv_file, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            if row[csv_reader.fieldnames[0]].strip() == DATA_Suffix_lowepscenter:
                MM_Offset_lowepscenter = float(row["MM_Offset"])
            elif row[csv_reader.fieldnames[0]].strip() == DATA_Suffix_lowepsleft:
                MM_Offset_lowepsleft = float(row["MM_Offset"])
            elif row[csv_reader.fieldnames[0]].strip() == DATA_Suffix_highepsright:
                MM_Offset_highepsright = float(row["MM_Offset"])
            elif row[csv_reader.fieldnames[0]].strip() == DATA_Suffix_highepscenter:
                MM_Offset_highepscenter = float(row["MM_Offset"])
            elif row[csv_reader.fieldnames[0]].strip() == DATA_Suffix_highepsleft:
                MM_Offset_highepsleft = float(row["MM_Offset"])
except (FileNotFoundError, ValueError, KeyError) as e:
    print(f"Error reading MM_Offset values from {mmcut_csv_file}: {e}")
    sys.exit(1)

# Check if all offsets were assigned
if None in [MM_Offset_lowepscenter, MM_Offset_lowepsleft, MM_Offset_highepsright, MM_Offset_highepscenter, MM_Offset_highepsleft]:
    print("Error: One or more MM_Offset values could not be assigned. Please check the CSV file.")
    sys.exit(1)

# Print the assigned values for verification
#print(f"MM_Offset_lowepscenter = {MM_Offset_lowepscenter}")
#print(f"MM_Offset_lowepsleft = {MM_Offset_lowepsleft}")
#print(f"MM_Offset_highepsright = {MM_Offset_highepsright}")
#print(f"MM_Offset_highepscenter = {MM_Offset_highepscenter}")
#print(f"MM_Offset_highepsleft = {MM_Offset_highepsleft}")

DATA_MMpi_Cut_lowepscenter = lambda event: (MM_Cut_lowvalue <= (event.MMpi + MM_Offset_lowepscenter) <= MM_Cut_highvalue)
DATA_MMpi_Cut_lowepsleft = lambda event: (MM_Cut_lowvalue <= (event.MMpi + MM_Offset_lowepsleft) <= MM_Cut_highvalue)
DATA_MMpi_Cut_highepsright = lambda event: (MM_Cut_lowvalue <= (event.MMpi + MM_Offset_highepsright) <= MM_Cut_highvalue)
DATA_MMpi_Cut_highepscenter = lambda event: (MM_Cut_lowvalue <= (event.MMpi + MM_Offset_highepscenter) <= MM_Cut_highvalue)
DATA_MMpi_Cut_highepsleft = lambda event: (MM_Cut_lowvalue <= (event.MMpi + MM_Offset_highepsleft) <= MM_Cut_highvalue)
Diamond_Cut = lambda event: (cutg_diamond.IsInside(event.Q2, event.W))

t_min = 0.10
t_max = 0.50

t_cut = lambda event: (t_min <= -event.MandelT <= t_max)  # MandelT is the t variable in the event tree

###################################################################################################################################################

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
run_numbers = read_run_numbers(run_list)

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
        if any(run_num in range(int(array_list[0]), int(array_list[1]) + 1) for run_num in run_numbers):
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

###################################################################################################################################################
nbins = 500
Q_min = 2.0
Q_max = 5.5
W_min = 2.0
W_max = 3.2
# Defining Histograms for Pions
Q2_vs_W_pions_data_prompt_lowepscenter_cut_all = ROOT.TH2D("Q2_vs_W_pions_data_prompt_lowepscenter_cut_all", "Q^{2} vs W Distribution; Q^{2}; W", nbins, Q_min, Q_max, nbins, W_min, W_max)
Q2_vs_W_pions_data_prompt_highepscenter_cut_all = ROOT.TH2D("Q2_vs_W_pions_data_prompt_highepscenter_cut_all", "Q^{2} vs W Distribution; Q^{2}; W", nbins, Q_min, Q_max, nbins, W_min, W_max)
Q2_vs_W_pions_data_random_lowepscenter_cut_all = ROOT.TH2D("Q2_vs_W_pions_data_random_lowepscenter_cut_all", "Q^{2} vs W Distribution; Q^{2}; W", nbins, Q_min, Q_max, nbins, W_min, W_max)
Q2_vs_W_pions_data_random_highepscenter_cut_all = ROOT.TH2D("Q2_vs_W_pions_data_random_highepscenter_cut_all", "Q^{2} vs W Distribution; Q^{2}; W", nbins, Q_min, Q_max, nbins, W_min, W_max)

Q2_vs_W_pions_dummy_prompt_lowepscenter_cut_all = ROOT.TH2D("Q2_vs_W_pions_dummy_prompt_lowepscenter_cut_all", "Q^{2} vs W Distribution; Q^{2}; W", nbins, Q_min, Q_max, nbins, W_min, W_max)
Q2_vs_W_pions_dummy_prompt_highepscenter_cut_all = ROOT.TH2D("Q2_vs_W_pions_dummy_prompt_highepscenter_cut_all", "Q^{2} vs W Distribution; Q^{2}; W", nbins, Q_min, Q_max, nbins, W_min, W_max)
Q2_vs_W_pions_dummy_random_lowepscenter_cut_all = ROOT.TH2D("Q2_vs_W_pions_dummy_random_lowepscenter_cut_all", "Q^{2} vs W Distribution; Q^{2}; W", nbins, Q_min, Q_max, nbins, W_min, W_max)
Q2_vs_W_pions_dummy_random_highepscenter_cut_all = ROOT.TH2D("Q2_vs_W_pions_dummy_random_highepscenter_cut_all", "Q^{2} vs W Distribution; Q^{2}; W", nbins, Q_min, Q_max, nbins, W_min, W_max)

Q2_vs_W_pions_data_randsub_lowepscenter_cut_all = ROOT.TH2D("Q2_vs_W_pions_data_randsub_lowepscenter_cut_all", "Q^{2} vs W Distribution; Q^{2}; W", nbins, Q_min, Q_max, nbins, W_min, W_max)
Q2_vs_W_pions_data_randsub_highepscenter_cut_all = ROOT.TH2D("Q2_vs_W_pions_data_randsub_highepscenter_cut_all", "Q^{2} vs W Distribution; Q^{2}; W", nbins, Q_min, Q_max, nbins, W_min, W_max)
Q2_vs_W_pions_dummy_randsub_lowepscenter_cut_all = ROOT.TH2D("Q2_vs_W_pions_dummy_randsub_lowepscenter_cut_all", "Q^{2} vs W Distribution; Q^{2}; W", nbins, Q_min, Q_max, nbins, W_min, W_max)
Q2_vs_W_pions_dummy_randsub_highepscenter_cut_all = ROOT.TH2D("Q2_vs_W_pions_dummy_randsub_highepscenter_cut_all", "Q^{2} vs W Distribution; Q^{2}; W", nbins, Q_min, Q_max, nbins, W_min, W_max)
Q2_vs_W_pions_data_dummysub_lowepscenter_cut_all = ROOT.TH2D("Q2_vs_W_pions_data_dummysub_lowepscenter_cut_all", "Q^{2} vs W Distribution; Q^{2}; W", nbins, Q_min, Q_max, nbins, W_min, W_max)
Q2_vs_W_pions_data_dummysub_highepscenter_cut_all = ROOT.TH2D("Q2_vs_W_pions_data_dummysub_highepscenter_cut_all", "Q^{2} vs W Distribution; Q^{2}; W", nbins, Q_min, Q_max, nbins, W_min, W_max)

# Defining variables for t vs phi histograms
nbins_loweps = 350
min_loweps = 0.10
max_loweps = 0.60
nbins_higheps = 350
min_higheps = 0.10
max_higheps = 0.60

phi_vs_t_data_prompt_lowepscenter_cut_all = ROOT.TGraph()
phi_vs_t_data_prompt_lowepsleft_cut_all = ROOT.TGraph()
phi_vs_t_data_prompt_highepsright_cut_all = ROOT.TGraph()
phi_vs_t_data_prompt_highepscenter_cut_all = ROOT.TGraph()
phi_vs_t_data_prompt_highepsleft_cut_all = ROOT.TGraph()

##########################################################################################################################################################################################################

# Read stuff from the main event tree
infile_DATA_lowepscenter = ROOT.TFile.Open(rootFile_DATA_lowepscenter, "READ")
infile_DATA_lowepsleft = ROOT.TFile.Open(rootFile_DATA_lowepsleft, "READ")
infile_DATA_highepsright = ROOT.TFile.Open(rootFile_DATA_highepsright, "READ")
infile_DATA_highepscenter = ROOT.TFile.Open(rootFile_DATA_highepscenter, "READ")
infile_DATA_highepsleft = ROOT.TFile.Open(rootFile_DATA_highepsleft, "READ")

infile_DUMMY_lowepscenter = ROOT.TFile.Open(rootFile_DUMMY_lowepscenter, "READ")
infile_DUMMY_lowepsleft = ROOT.TFile.Open(rootFile_DUMMY_lowepsleft, "READ")
infile_DUMMY_highepsright = ROOT.TFile.Open(rootFile_DUMMY_highepsright, "READ")
infile_DUMMY_highepscenter = ROOT.TFile.Open(rootFile_DUMMY_highepscenter, "READ")
infile_DUMMY_highepsleft = ROOT.TFile.Open(rootFile_DUMMY_highepsleft, "READ")

# Grab the trees
Cut_Pion_Events_Prompt_Data_lowepscenter_tree = infile_DATA_lowepscenter.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Random_Data_lowepscenter_tree = infile_DATA_lowepscenter.Get("Cut_Pion_Events_Random")
Cut_Pion_Events_Prompt_Data_lowepsleft_tree = infile_DATA_lowepsleft.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Random_Data_lowepsleft_tree = infile_DATA_lowepsleft.Get("Cut_Pion_Events_Random")
Cut_Pion_Events_Prompt_Data_highepsright_tree = infile_DATA_highepsright.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Random_Data_highepsright_tree = infile_DATA_highepsright.Get("Cut_Pion_Events_Random")
Cut_Pion_Events_Prompt_Data_highepscenter_tree = infile_DATA_highepscenter.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Random_Data_highepscenter_tree = infile_DATA_highepscenter.Get("Cut_Pion_Events_Random")
Cut_Pion_Events_Prompt_Data_highepsleft_tree = infile_DATA_highepsleft.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Random_Data_highepsleft_tree = infile_DATA_highepsleft.Get("Cut_Pion_Events_Random")

Cut_Pion_Events_Prompt_Dummy_lowepscenter_tree = infile_DUMMY_lowepscenter.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Random_Dummy_lowepscenter_tree = infile_DUMMY_lowepscenter.Get("Cut_Pion_Events_Random")
Cut_Pion_Events_Prompt_Dummy_lowepsleft_tree = infile_DUMMY_lowepsleft.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Random_Dummy_lowepsleft_tree = infile_DUMMY_lowepsleft.Get("Cut_Pion_Events_Random")
Cut_Pion_Events_Prompt_Dummy_highepsright_tree = infile_DUMMY_highepsright.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Random_Dummy_highepsright_tree = infile_DUMMY_highepsright.Get("Cut_Pion_Events_Random")
Cut_Pion_Events_Prompt_Dummy_highepscenter_tree = infile_DUMMY_highepscenter.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Random_Dummy_highepscenter_tree = infile_DUMMY_highepscenter.Get("Cut_Pion_Events_Random")
Cut_Pion_Events_Prompt_Dummy_highepsleft_tree = infile_DUMMY_highepsleft.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Random_Dummy_highepsleft_tree = infile_DUMMY_highepsleft.Get("Cut_Pion_Events_Random")

###################################################################################################################################################

#Fill histograms for Cut All Data
for event in Cut_Pion_Events_Prompt_Data_lowepscenter_tree:
    if DATA_MMpi_Cut_lowepscenter(event):
        Q2_vs_W_pions_data_prompt_lowepscenter_cut_all.Fill(event.Q2, event.W)
for event in Cut_Pion_Events_Prompt_Data_highepscenter_tree:
    if DATA_MMpi_Cut_highepscenter(event):
        Q2_vs_W_pions_data_prompt_highepscenter_cut_all.Fill(event.Q2, event.W)
for event in Cut_Pion_Events_Random_Data_lowepscenter_tree:
    if DATA_MMpi_Cut_lowepscenter(event):
        Q2_vs_W_pions_data_random_lowepscenter_cut_all.Fill(event.Q2, event.W)
for event in Cut_Pion_Events_Random_Data_highepscenter_tree:
    if DATA_MMpi_Cut_highepscenter(event):
        Q2_vs_W_pions_data_random_highepscenter_cut_all.Fill(event.Q2, event.W)
for event in Cut_Pion_Events_Prompt_Dummy_lowepscenter_tree:
    if DATA_MMpi_Cut_lowepscenter(event):
        Q2_vs_W_pions_dummy_prompt_lowepscenter_cut_all.Fill(event.Q2, event.W)
for event in Cut_Pion_Events_Prompt_Dummy_highepscenter_tree:
    if DATA_MMpi_Cut_highepscenter(event):
        Q2_vs_W_pions_dummy_prompt_highepscenter_cut_all.Fill(event.Q2, event.W)
for event in Cut_Pion_Events_Random_Dummy_lowepscenter_tree:
    if DATA_MMpi_Cut_lowepscenter(event):
        Q2_vs_W_pions_dummy_random_lowepscenter_cut_all.Fill(event.Q2, event.W)
for event in Cut_Pion_Events_Random_Dummy_highepscenter_tree:
    if DATA_MMpi_Cut_highepscenter(event):
        Q2_vs_W_pions_dummy_random_highepscenter_cut_all.Fill(event.Q2, event.W)

#-----------------------------------------------------------------------------------------------------------------------

for event in Cut_Pion_Events_Prompt_Data_lowepscenter_tree:
    if DATA_MMpi_Cut_lowepscenter(event) & Diamond_Cut(event) & t_cut(event):
        phi_vs_t_data_prompt_lowepscenter_cut_all.SetPoint(phi_vs_t_data_prompt_lowepscenter_cut_all.GetN(), event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Prompt_Data_lowepsleft_tree:
    if DATA_MMpi_Cut_lowepsleft(event) & Diamond_Cut(event) & t_cut(event):
        phi_vs_t_data_prompt_lowepsleft_cut_all.SetPoint(phi_vs_t_data_prompt_lowepsleft_cut_all.GetN(), event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Prompt_Data_highepsright_tree:
    if DATA_MMpi_Cut_highepsright(event) & Diamond_Cut(event) & t_cut(event):
        phi_vs_t_data_prompt_highepsright_cut_all.SetPoint(phi_vs_t_data_prompt_highepsright_cut_all.GetN(), event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Prompt_Data_highepscenter_tree:
    if DATA_MMpi_Cut_highepscenter(event) & Diamond_Cut(event) & t_cut(event):
        phi_vs_t_data_prompt_highepscenter_cut_all.SetPoint(phi_vs_t_data_prompt_highepscenter_cut_all.GetN(), event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Prompt_Data_highepsleft_tree:
    if DATA_MMpi_Cut_highepsleft(event) & Diamond_Cut(event) & t_cut(event):
        phi_vs_t_data_prompt_highepsleft_cut_all.SetPoint(phi_vs_t_data_prompt_highepsleft_cut_all.GetN(), event.ph_q, -event.MandelT)

print("####################################")
print("###### Histogram filling done ######")
print("####################################\n")

#################################################################################################################################################

# Random subtraction from Histograms
Q2_vs_W_pions_data_random_lowepscenter_cut_all.Scale(1.0/nWindows)
Q2_vs_W_pions_data_random_highepscenter_cut_all.Scale(1.0/nWindows)
Q2_vs_W_pions_data_randsub_lowepscenter_cut_all.Add(Q2_vs_W_pions_data_prompt_lowepscenter_cut_all, Q2_vs_W_pions_data_random_lowepscenter_cut_all, 1, -1)
Q2_vs_W_pions_data_randsub_highepscenter_cut_all.Add(Q2_vs_W_pions_data_prompt_highepscenter_cut_all, Q2_vs_W_pions_data_random_highepscenter_cut_all, 1, -1)

Q2_vs_W_pions_dummy_random_lowepscenter_cut_all.Scale(1.0/nWindows)
Q2_vs_W_pions_dummy_random_highepscenter_cut_all.Scale(1.0/nWindows)
Q2_vs_W_pions_dummy_randsub_lowepscenter_cut_all.Add(Q2_vs_W_pions_dummy_prompt_lowepscenter_cut_all, Q2_vs_W_pions_dummy_random_lowepscenter_cut_all, 1, -1)
Q2_vs_W_pions_dummy_randsub_highepscenter_cut_all.Add(Q2_vs_W_pions_dummy_prompt_highepscenter_cut_all, Q2_vs_W_pions_dummy_random_highepscenter_cut_all, 1, -1)

print("######################################")
print("###### Random substraction done ######")
print("######################################\n")

############################################################################################################################################

# Dummy Subtraction
Q2_vs_W_pions_data_dummysub_lowepscenter_cut_all.Add(Q2_vs_W_pions_data_randsub_lowepscenter_cut_all, Q2_vs_W_pions_dummy_randsub_lowepscenter_cut_all, 1, -1)
Q2_vs_W_pions_data_dummysub_highepscenter_cut_all.Add(Q2_vs_W_pions_data_randsub_highepscenter_cut_all, Q2_vs_W_pions_dummy_randsub_highepscenter_cut_all, 1, -1)

print("####################################")
print("###### Dummy subtraction done ######")
print("####################################\n")

#############################################################################################################################################################

# Plotting Histograms
# Removes stat box
ROOT.gStyle.SetOptStat(0)

# Saving histograms in PDF
c1_delta = TCanvas("c1_delta", "Variables Distributions", 100, 0, 1800, 1600)
c1_delta.Divide(2, 2)
c1_delta.cd(1)
c1_delta_text_lines = [
    ROOT.TText(0.5, 0.9, "Pion ProdCoin Setting"),
    ROOT.TText(0.5, 0.8, "{}".format(PHY_SETTING)),
]
for c1_delta_text in c1_delta_text_lines:
    c1_delta_text.SetTextSize(0.05)
    c1_delta_text.SetTextAlign(22)
    c1_delta_text.SetTextColor(ROOT.kBlack)
    c1_delta_text.Draw()
c1_delta.cd(2)
ROOT.gPad.SetLogz()
Q2_vs_W_pions_data_dummysub_highepscenter_cut_all.Draw("colz")
Q2_vs_W_pions_data_dummysub_lowepscenter_cut_all.SetLineColor(kGreen-2)
Q2_vs_W_pions_data_dummysub_lowepscenter_cut_all.Draw("box,same")
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
# Create colored dummy boxes for the legend
yellow_box = ROOT.TBox(0, 0, 1, 1)
yellow_box.SetFillColor(ROOT.kYellow-6)
yellow_box.SetLineColor(ROOT.kYellow-6)
green_box = ROOT.TBox(0, 0, 1, 1)
green_box.SetFillColor(ROOT.kGreen-2)
green_box.SetLineColor(ROOT.kGreen-2)
legend = ROOT.TLegend(0.60, 0.70, 0.85, 0.88)
legend.AddEntry(yellow_box, "High #epsilon", "f")
legend.AddEntry(green_box, "Low #epsilon", "f")
legend.AddEntry(line1, "cut", "l")
legend.Draw()
c1_delta.cd(3)
rad_min_loweps = 0.0
rad_max_loweps = 1.0
g_phi_vs_t_data_prompt_lowepscenter_cut_all = ROOT.TGraphPolar(phi_vs_t_data_prompt_lowepscenter_cut_all.GetN(), phi_vs_t_data_prompt_lowepscenter_cut_all.GetX(), phi_vs_t_data_prompt_lowepscenter_cut_all.GetY())
g_phi_vs_t_data_prompt_lowepsleft_cut_all = ROOT.TGraphPolar(phi_vs_t_data_prompt_lowepsleft_cut_all.GetN(), phi_vs_t_data_prompt_lowepsleft_cut_all.GetX(), phi_vs_t_data_prompt_lowepsleft_cut_all.GetY())
g_phi_vs_t_data_prompt_lowepscenter_cut_all.SetMarkerColor(3)
g_phi_vs_t_data_prompt_lowepscenter_cut_all.SetMarkerSize(0.5)
g_phi_vs_t_data_prompt_lowepscenter_cut_all.SetMarkerStyle(20)
g_phi_vs_t_data_prompt_lowepscenter_cut_all.SetTitle(" ")
g_phi_vs_t_data_prompt_lowepscenter_cut_all.GetXaxis().SetName("#Phi")
g_phi_vs_t_data_prompt_lowepscenter_cut_all.GetYaxis().SetName("-t")
g_phi_vs_t_data_prompt_lowepscenter_cut_all.GetXaxis().SetNdivisions(510)
g_phi_vs_t_data_prompt_lowepscenter_cut_all.GetYaxis().SetRangeUser(rad_min_loweps, rad_max_loweps)
g_phi_vs_t_data_prompt_lowepscenter_cut_all.SetMinRadial(rad_min_loweps)
g_phi_vs_t_data_prompt_lowepscenter_cut_all.SetMaxRadial(rad_max_loweps)
g_phi_vs_t_data_prompt_lowepscenter_cut_all.Draw("AP")
g_phi_vs_t_data_prompt_lowepsleft_cut_all.SetMarkerColor(4)
g_phi_vs_t_data_prompt_lowepsleft_cut_all.SetMarkerSize(0.5)
g_phi_vs_t_data_prompt_lowepsleft_cut_all.SetMarkerStyle(20)
g_phi_vs_t_data_prompt_lowepsleft_cut_all.SetTitle(" ")
g_phi_vs_t_data_prompt_lowepsleft_cut_all.SetMinRadial(rad_min_loweps)
g_phi_vs_t_data_prompt_lowepsleft_cut_all.SetMaxRadial(rad_max_loweps)
g_phi_vs_t_data_prompt_lowepsleft_cut_all.Draw("PSame")
c1_delta.cd(4)
rad_min_higheps = 0.0
rad_max_higheps = 1.0
g_phi_vs_t_data_prompt_highepsright_cut_all = ROOT.TGraphPolar(phi_vs_t_data_prompt_highepsright_cut_all.GetN(), phi_vs_t_data_prompt_highepsright_cut_all.GetX(), phi_vs_t_data_prompt_highepsright_cut_all.GetY())
g_phi_vs_t_data_prompt_highepscenter_cut_all = ROOT.TGraphPolar(phi_vs_t_data_prompt_highepscenter_cut_all.GetN(), phi_vs_t_data_prompt_highepscenter_cut_all.GetX(), phi_vs_t_data_prompt_highepscenter_cut_all.GetY())
g_phi_vs_t_data_prompt_highepsleft_cut_all = ROOT.TGraphPolar(phi_vs_t_data_prompt_highepsleft_cut_all.GetN(), phi_vs_t_data_prompt_highepsleft_cut_all.GetX(), phi_vs_t_data_prompt_highepsleft_cut_all.GetY())
g_phi_vs_t_data_prompt_highepsright_cut_all.SetMarkerColor(2)
g_phi_vs_t_data_prompt_highepsright_cut_all.SetMarkerSize(0.5)
g_phi_vs_t_data_prompt_highepsright_cut_all.SetMarkerStyle(20)
g_phi_vs_t_data_prompt_highepsright_cut_all.SetTitle(" ")
g_phi_vs_t_data_prompt_highepsright_cut_all.GetXaxis().SetName("#Phi")
g_phi_vs_t_data_prompt_highepsright_cut_all.GetYaxis().SetName("-t")
g_phi_vs_t_data_prompt_highepsright_cut_all.GetXaxis().SetNdivisions(510)
g_phi_vs_t_data_prompt_highepsright_cut_all.GetXaxis().SetRangeUser(rad_min_higheps, rad_max_higheps)
g_phi_vs_t_data_prompt_highepsright_cut_all.SetMinRadial(rad_min_higheps)
g_phi_vs_t_data_prompt_highepsright_cut_all.SetMaxRadial(rad_max_higheps)
g_phi_vs_t_data_prompt_highepsright_cut_all.Draw("AP")
g_phi_vs_t_data_prompt_highepsleft_cut_all.SetMarkerColor(4)
g_phi_vs_t_data_prompt_highepsleft_cut_all.SetMarkerSize(0.5)
g_phi_vs_t_data_prompt_highepsleft_cut_all.SetMarkerStyle(20)
g_phi_vs_t_data_prompt_highepsleft_cut_all.SetTitle(" ")
g_phi_vs_t_data_prompt_highepsleft_cut_all.SetMinRadial(rad_min_higheps)
g_phi_vs_t_data_prompt_highepsleft_cut_all.SetMaxRadial(rad_max_higheps)
g_phi_vs_t_data_prompt_highepsleft_cut_all.Draw("PSame")
g_phi_vs_t_data_prompt_highepscenter_cut_all.SetMarkerColor(3)
g_phi_vs_t_data_prompt_highepscenter_cut_all.SetMarkerSize(0.5)
g_phi_vs_t_data_prompt_highepscenter_cut_all.SetMarkerStyle(20)
g_phi_vs_t_data_prompt_highepscenter_cut_all.SetTitle(" ")
g_phi_vs_t_data_prompt_highepscenter_cut_all.SetMinRadial(rad_min_higheps)
g_phi_vs_t_data_prompt_highepscenter_cut_all.SetMaxRadial(rad_max_higheps)
g_phi_vs_t_data_prompt_highepscenter_cut_all.Draw("PSame")
c1_delta.Print(Pion_Analysis_Distributions)

#======================================================================================================================================================================================================================================================

# Close input files
infile_DATA_lowepscenter.Close()
infile_DATA_lowepsleft.Close()
infile_DATA_highepsright.Close()
infile_DATA_highepscenter.Close()
infile_DATA_highepsleft.Close()
infile_DUMMY_lowepscenter.Close()
infile_DUMMY_lowepsleft.Close()
infile_DUMMY_highepsright.Close()
infile_DUMMY_highepscenter.Close()
infile_DUMMY_highepsleft.Close() 

print ("Processing Complete")