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
import ctypes

##################################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=5:
    print("!!!!! ERROR !!!!!\n Expected 4 arguments\n Usage is with - PHY_SETTING MaxEvents Suffix RunList\n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################

# Input params - run number and max number of events
PHY_SETTING = sys.argv[1]
MaxEvent = sys.argv[2]
DATA_Suffix = sys.argv[3]
DUMMY_Suffix = sys.argv[4]
RUN_LIST = sys.argv[5]

DATA_Suffix_lowepscenter = "{}_loweps_center".format(PHY_SETTING)
DATA_Suffix_lowepsleft = "{}_loweps_left".format(PHY_SETTING)
DATA_Suffix_highepsright = "{}_higheps_right".format(PHY_SETTING)
DATA_Suffix_highepscenter = "{}_higheps_center".format(PHY_SETTING)
DATA_Suffix_highepsleft = "{}_higheps_left".format(PHY_SETTING)

# Extract the first three words from PHY_SETTING for the CSV file name
setting_name = "_".join(PHY_SETTING.split("_")[:3])
physet_dir_name = "%s_std" % (setting_name)

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
SIMC_tRES_CSV = "/u/group/c-pionlt/USERS/%s/hallc_replay_lt/UTIL_PION/LTSep_CSVs/t_resolution_csv" % (USER)

#################################################################################################################################################

# Output PDF File Name
print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))
Pion_Analysis_Distributions = "%s/%s_%s_ProdCoin_Pion_Analysis_tbinning_Distributions.pdf" % (OUTPATH, PHY_SETTING, MaxEvent)

# Extract the first three words from PHY_SETTING for the CSV file name
mmcut_csv_file = "%s/%s/%s_mm_offsets_cuts_parameters.csv" % (MMCUT_CSV, physet_dir_name, setting_name)
dcut_csv_file = "%s/%s/%s_diamond_cut_parameters.csv" % (DCUT_CSV, physet_dir_name, setting_name)
simc_tres_csv_file = "%s/%s/%s_simc_t_resolution_parameters.csv" % (SIMC_tRES_CSV, physet_dir_name, setting_name)

# Input file location and variables taking
rootFile_DATA_lowepscenter = "%s/%s_%s_%s.root" % (OUTPATH, DATA_Suffix_lowepscenter, MaxEvent, DATA_Suffix)
rootFile_DATA_lowepsleft = "%s/%s_%s_%s.root" % (OUTPATH, DATA_Suffix_lowepsleft, MaxEvent, DATA_Suffix)
rootFile_DATA_highepsright = "%s/%s_%s_%s.root" % (OUTPATH, DATA_Suffix_highepsright, MaxEvent, DATA_Suffix)
rootFile_DATA_highepscenter = "%s/%s_%s_%s.root" % (OUTPATH, DATA_Suffix_highepscenter, MaxEvent, DATA_Suffix)
rootFile_DATA_highepsleft = "%s/%s_%s_%s.root" % (OUTPATH, DATA_Suffix_highepsleft, MaxEvent, DATA_Suffix)

rootFile_DUMMY_lowepscenter = "%s/%s_%s_%s.root" % (OUTPATH, DATA_Suffix_lowepscenter, MaxEvent, DUMMY_Suffix)
rootFile_DUMMY_lowepsleft = "%s/%s_%s_%s.root" % (OUTPATH, DATA_Suffix_lowepsleft, MaxEvent, DUMMY_Suffix)
rootFile_DUMMY_highepsright = "%s/%s_%s_%s.root" % (OUTPATH, DATA_Suffix_highepsright, MaxEvent, DUMMY_Suffix)
rootFile_DUMMY_highepscenter = "%s/%s_%s_%s.root" % (OUTPATH, DATA_Suffix_highepscenter, MaxEvent, DUMMY_Suffix)
rootFile_DUMMY_highepsleft = "%s/%s_%s_%s.root" % (OUTPATH, DATA_Suffix_highepsleft, MaxEvent, DUMMY_Suffix)

run_list = "%s/%s" % (RUNLISTPATH, RUN_LIST)

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

# t-resolution values
# Read the t-resolution values from the CSV file
try:
    with open(simc_tres_csv_file, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            if row["Setting"] == "loweps":
                simc_t_resolution_pions_lowe = float(row["Value"])
            elif row["Setting"] == "higheps":
                simc_t_resolution_pions_highe = float(row["Value"])
except (FileNotFoundError, ValueError, KeyError) as e:
    print(f"Error reading t-resolution parameters from {simc_tres_csv_file}: {e}")
    sys.exit(1)

# Print the assigned values for verification
print(f"SIMC t-resolution (loweps): {simc_t_resolution_pions_lowe}")
print(f"SIMC t-resolution (higheps): {simc_t_resolution_pions_highe}")


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

# Defining Histograms for Pions
t_pions_data_prompt_lowepscenter_cut_all = ROOT.TH1F("t_pions_data_prompt_lowepscenter_cut_all", "t_pions_data_prompt_lowepscenter_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_data_prompt_lowepsleft_cut_all = ROOT.TH1F("t_pions_data_prompt_lowepsleft_cut_all", "t_pions_data_prompt_lowepsleft_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_data_prompt_highepsright_cut_all = ROOT.TH1F("t_pions_data_prompt_highepsright_cut_all", "t_pions_data_prompt_highepsright_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_data_prompt_highepscenter_cut_all = ROOT.TH1F("t_pions_data_prompt_highepscenter_cut_all", "t_pions_data_prompt_highepscenter_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_data_prompt_highepsleft_cut_all = ROOT.TH1F("t_pions_data_prompt_highepsleft_cut_all", "t_pions_data_prompt_highepsleft_cut_all; -t; Counts", nbins, 0.0, 1.0)

t_pions_data_random_lowepscenter_cut_all = ROOT.TH1F("t_pions_data_random_lowepscenter_cut_all", "t_pions_data_random_lowepscenter_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_data_random_lowepsleft_cut_all = ROOT.TH1F("t_pions_data_random_lowepsleft_cut_all", "t_pions_data_random_lowepsleft_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_data_random_highepsright_cut_all = ROOT.TH1F("t_pions_data_random_highepsright_cut_all", "t_pions_data_random_highepsright_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_data_random_highepscenter_cut_all = ROOT.TH1F("t_pions_data_random_highepscenter_cut_all", "t_pions_data_random_highepscenter_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_data_random_highepsleft_cut_all = ROOT.TH1F("t_pions_data_random_highepsleft_cut_all", "t_pions_data_random_highepsleft_cut_all; -t; Counts", nbins, 0.0, 1.0)

t_pions_dummy_prompt_lowepscenter_cut_all = ROOT.TH1F("t_pions_dummy_prompt_lowepscenter_cut_all", "t_pions_dummy_prompt_lowepscenter_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_dummy_prompt_lowepsleft_cut_all = ROOT.TH1F("t_pions_dummy_prompt_lowepsleft_cut_all", "t_pions_dummy_prompt_lowepsleft_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_dummy_prompt_highepsright_cut_all = ROOT.TH1F("t_pions_dummy_prompt_highepsright_cut_all", "t_pions_dummy_prompt_highepsright_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_dummy_prompt_highepscenter_cut_all = ROOT.TH1F("t_pions_dummy_prompt_highepscenter_cut_all", "t_pions_dummy_prompt_highepscenter_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_dummy_prompt_highepsleft_cut_all = ROOT.TH1F("t_pions_dummy_prompt_highepsleft_cut_all", "t_pions_dummy_prompt_highepsleft_cut_all; -t; Counts", nbins, 0.0, 1.0)

t_pions_dummy_random_lowepscenter_cut_all = ROOT.TH1F("t_pions_dummy_random_lowepscenter_cut_all", "t_pions_dummy_random_lowepscenter_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_dummy_random_lowepsleft_cut_all = ROOT.TH1F("t_pions_dummy_random_lowepsleft_cut_all", "t_pions_dummy_random_lowepsleft_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_dummy_random_highepsright_cut_all = ROOT.TH1F("t_pions_dummy_random_highepsright_cut_all", "t_pions_dummy_random_highepsright_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_dummy_random_highepscenter_cut_all = ROOT.TH1F("t_pions_dummy_random_highepscenter_cut_all", "t_pions_dummy_random_highepscenter_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_dummy_random_highepsleft_cut_all = ROOT.TH1F("t_pions_dummy_random_highepsleft_cut_all", "t_pions_dummy_random_highepsleft_cut_all; -t; Counts", nbins, 0.0, 1.0)

t_pions_data_randsub_lowepscenter_cut_all = ROOT.TH1F("t_pions_data_randsub_lowepscenter_cut_all", "t_pions_data_randsub_lowepscenter_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_data_randsub_lowepsleft_cut_all = ROOT.TH1F("t_pions_data_randsub_lowepsleft_cut_all", "t_pions_data_randsub_lowepsleft_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_data_randsub_highepsright_cut_all = ROOT.TH1F("t_pions_data_randsub_highepsright_cut_all", "t_pions_data_randsub_highepsright_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_data_randsub_highepscenter_cut_all = ROOT.TH1F("t_pions_data_randsub_highepscenter_cut_all", "t_pions_data_randsub_highepscenter_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_data_randsub_highepsleft_cut_all = ROOT.TH1F("t_pions_data_randsub_highepsleft_cut_all", "t_pions_data_randsub_highepsleft_cut_all; -t; Counts", nbins, 0.0, 1.0)

t_pions_dummy_randsub_lowepscenter_cut_all = ROOT.TH1F("t_pions_dummy_randsub_lowepscenter_cut_all", "t_pions_dummy_randsub_lowepscenter_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_dummy_randsub_lowepsleft_cut_all = ROOT.TH1F("t_pions_dummy_randsub_lowepsleft_cut_all", "t_pions_dummy_randsub_lowepsleft_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_dummy_randsub_highepsright_cut_all = ROOT.TH1F("t_pions_dummy_randsub_highepsright_cut_all", "t_pions_dummy_randsub_highepsright_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_dummy_randsub_highepscenter_cut_all = ROOT.TH1F("t_pions_dummy_randsub_highepscenter_cut_all", "t_pions_dummy_randsub_highepscenter_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_dummy_randsub_highepsleft_cut_all = ROOT.TH1F("t_pions_dummy_randsub_highepsleft_cut_all", "t_pions_dummy_randsub_highepsleft_cut_all; -t; Counts", nbins, 0.0, 1.0)

t_pions_data_dummysub_lowepscenter_cut_all = ROOT.TH1F("t_pions_data_dummysub_lowepscenter_cut_all", "t_pions_data_dummysub_lowepscenter_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_data_dummysub_lowepsleft_cut_all = ROOT.TH1F("t_pions_data_dummysub_lowepsleft_cut_all", "t_pions_data_dummysub_lowepsleft_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_data_dummysub_highepsright_cut_all = ROOT.TH1F("t_pions_data_dummysub_highepsright_cut_all", "t_pions_data_dummysub_highepsright_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_data_dummysub_highepscenter_cut_all = ROOT.TH1F("t_pions_data_dummysub_highepscenter_cut_all", "t_pions_data_dummysub_highepscenter_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_data_dummysub_highepsleft_cut_all = ROOT.TH1F("t_pions_data_dummysub_highepsleft_cut_all", "t_pions_data_dummysub_highepsleft_cut_all; -t; Counts", nbins, 0.0, 1.0)

t_pions_data_dummysub_loweps_cut_all = ROOT.TH1F("t_pions_data_dummysub_loweps_cut_all", "t_pions_data_dummysub_loweps_cut_all; -t; Counts", nbins, 0.0, 1.0)
t_pions_data_dummysub_higheps_cut_all = ROOT.TH1F("t_pions_data_dummysub_higheps_cut_all", "t_pions_data_dummysub_higheps_cut_all; -t; Counts", nbins, 0.0, 1.0)

# Defining variables for t vs phi histograms
nbins_loweps = 350
min_loweps = 0.10
max_loweps = 0.60
nbins_higheps = 350
min_higheps = 0.10
max_higheps = 0.60

# 2D Histograms
phi_vs_t_data_prompt_lowepscenter_cut_all = ROOT.TH2F("phi_vs_t_data_prompt_lowepscenter_cut_all", "phi_vs_t_data_prompt_lowepscenter_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_loweps, min_loweps, max_loweps)
phi_vs_t_data_prompt_lowepsleft_cut_all = ROOT.TH2F("phi_vs_t_data_prompt_lowepsleft_cut_all", "phi_vs_t_data_prompt_lowepsleft_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_loweps, min_loweps, max_loweps)
phi_vs_t_data_prompt_highepsright_cut_all = ROOT.TH2F("phi_vs_t_data_prompt_highepsright_cut_all", "phi_vs_t_data_prompt_highepsright_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)
phi_vs_t_data_prompt_highepscenter_cut_all = ROOT.TH2F("phi_vs_t_data_prompt_highepscenter_cut_all", "phi_vs_t_data_prompt_highepscenter_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)
phi_vs_t_data_prompt_highepsleft_cut_all = ROOT.TH2F("phi_vs_t_data_prompt_highepsleft_cut_all", "phi_vs_t_data_prompt_highepsleft_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)

phi_vs_t_data_random_lowepscenter_cut_all = ROOT.TH2F("phi_vs_t_data_random_lowepscenter_cut_all", "phi_vs_t_data_random_lowepscenter_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_loweps, min_loweps, max_loweps)
phi_vs_t_data_random_lowepsleft_cut_all = ROOT.TH2F("phi_vs_t_data_random_lowepsleft_cut_all", "phi_vs_t_data_random_lowepsleft_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_loweps, min_loweps, max_loweps)
phi_vs_t_data_random_highepsright_cut_all = ROOT.TH2F("phi_vs_t_data_random_highepsright_cut_all", "phi_vs_t_data_random_highepsright_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)
phi_vs_t_data_random_highepscenter_cut_all = ROOT.TH2F("phi_vs_t_data_random_highepscenter_cut_all", "phi_vs_t_data_random_highepscenter_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)
phi_vs_t_data_random_highepsleft_cut_all = ROOT.TH2F("phi_vs_t_data_random_highepsleft_cut_all", "phi_vs_t_data_random_highepsleft_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)

phi_vs_t_dummy_prompt_lowepscenter_cut_all = ROOT.TH2F("phi_vs_t_dummy_prompt_lowepscenter_cut_all", "phi_vs_t_dummy_prompt_lowepscenter_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_loweps, min_loweps, max_loweps)
phi_vs_t_dummy_prompt_lowepsleft_cut_all = ROOT.TH2F("phi_vs_t_dummy_prompt_lowepsleft_cut_all", "phi_vs_t_dummy_prompt_lowepsleft_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_loweps, min_loweps, max_loweps)
phi_vs_t_dummy_prompt_highepsright_cut_all = ROOT.TH2F("phi_vs_t_dummy_prompt_highepsright_cut_all", "phi_vs_t_dummy_prompt_highepsright_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)
phi_vs_t_dummy_prompt_highepscenter_cut_all = ROOT.TH2F("phi_vs_t_dummy_prompt_highepscenter_cut_all", "phi_vs_t_dummy_prompt_highepscenter_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)
phi_vs_t_dummy_prompt_highepsleft_cut_all = ROOT.TH2F("phi_vs_t_dummy_prompt_highepsleft_cut_all", "phi_vs_t_dummy_prompt_highepsleft_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)

phi_vs_t_dummy_random_lowepscenter_cut_all = ROOT.TH2F("phi_vs_t_dummy_random_lowepscenter_cut_all", "phi_vs_t_dummy_random_lowepscenter_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_loweps, min_loweps, max_loweps)
phi_vs_t_dummy_random_lowepsleft_cut_all = ROOT.TH2F("phi_vs_t_dummy_random_lowepsleft_cut_all", "phi_vs_t_dummy_random_lowepsleft_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_loweps, min_loweps, max_loweps)
phi_vs_t_dummy_random_highepsright_cut_all = ROOT.TH2F("phi_vs_t_dummy_random_highepsright_cut_all", "phi_vs_t_dummy_random_highepsright_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)
phi_vs_t_dummy_random_highepscenter_cut_all = ROOT.TH2F("phi_vs_t_dummy_random_highepscenter_cut_all", "phi_vs_t_dummy_random_highepscenter_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)
phi_vs_t_dummy_random_highepsleft_cut_all = ROOT.TH2F("phi_vs_t_dummy_random_highepsleft_cut_all", "phi_vs_t_dummy_random_highepsleft_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)

phi_vs_t_data_randsub_lowepscenter_cut_all = ROOT.TH2F("phi_vs_t_data_randsub_lowepscenter_cut_all", "phi_vs_t_data_randsub_lowepscenter_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_loweps, min_loweps, max_loweps)
phi_vs_t_data_randsub_lowepsleft_cut_all = ROOT.TH2F("phi_vs_t_data_randsub_lowepsleft_cut_all", "phi_vs_t_data_randsub_lowepsleft_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_loweps, min_loweps, max_loweps)
phi_vs_t_data_randsub_highepsright_cut_all = ROOT.TH2F("phi_vs_t_data_randsub_highepsright_cut_all", "phi_vs_t_data_randsub_highepsright_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)
phi_vs_t_data_randsub_highepscenter_cut_all = ROOT.TH2F("phi_vs_t_data_randsub_highepscenter_cut_all", "phi_vs_t_data_randsub_highepscenter_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)
phi_vs_t_data_randsub_highepsleft_cut_all = ROOT.TH2F("phi_vs_t_data_randsub_highepsleft_cut_all", "phi_vs_t_data_randsub_highepsleft_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)

phi_vs_t_dummy_randsub_lowepscenter_cut_all = ROOT.TH2F("phi_vs_t_dummy_randsub_lowepscenter_cut_all", "phi_vs_t_dummy_randsub_lowepscenter_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_loweps, min_loweps, max_loweps)
phi_vs_t_dummy_randsub_lowepsleft_cut_all = ROOT.TH2F("phi_vs_t_dummy_randsub_lowepsleft_cut_all", "phi_vs_t_dummy_randsub_lowepsleft_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_loweps, min_loweps, max_loweps)
phi_vs_t_dummy_randsub_highepsright_cut_all = ROOT.TH2F("phi_vs_t_dummy_randsub_highepsright_cut_all", "phi_vs_t_dummy_randsub_highepsright_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)
phi_vs_t_dummy_randsub_highepscenter_cut_all = ROOT.TH2F("phi_vs_t_dummy_randsub_highepscenter_cut_all", "phi_vs_t_dummy_randsub_highepscenter_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)
phi_vs_t_dummy_randsub_highepsleft_cut_all = ROOT.TH2F("phi_vs_t_dummy_randsub_highepsleft_cut_all", "phi_vs_t_dummy_randsub_highepsleft_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)

phi_vs_t_data_dummysub_lowepscenter_cut_all = ROOT.TH2F("phi_vs_t_data_dummysub_lowepscenter_cut_all", "phi_vs_t_data_dummysub_lowepscenter_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_loweps, min_loweps, max_loweps)
phi_vs_t_data_dummysub_lowepsleft_cut_all = ROOT.TH2F("phi_vs_t_data_dummysub_lowepsleft_cut_all", "phi_vs_t_data_dummysub_lowepsleft_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_loweps, min_loweps, max_loweps)
phi_vs_t_data_dummysub_highepsright_cut_all = ROOT.TH2F("phi_vs_t_data_dummysub_highepsright_cut_all", "phi_vs_t_data_dummysub_highepsright_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)
phi_vs_t_data_dummysub_highepscenter_cut_all = ROOT.TH2F("phi_vs_t_data_dummysub_highepscenter_cut_all", "phi_vs_t_data_dummysub_highepscenter_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)
phi_vs_t_data_dummysub_highepsleft_cut_all = ROOT.TH2F("phi_vs_t_data_dummysub_highepsleft_cut_all", "phi_vs_t_data_dummysub_highepsleft_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)

phi_vs_t_data_dummysub_loweps_cut_all = ROOT.TH2F("phi_vs_t_data_dummysub_loweps_cut_all", "phi_vs_t_data_dummysub_loweps_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_loweps, min_loweps, max_loweps)
phi_vs_t_data_dummysub_higheps_cut_all = ROOT.TH2F("phi_vs_t_data_dummysub_higheps_cut_all", "phi_vs_t_data_dummysub_higheps_cut_all; #phi; -t", 12, -3.14, 3.14, nbins_higheps, min_higheps, max_higheps)

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
    if DATA_MMpi_Cut_lowepscenter(event) & Diamond_Cut(event):
        t_pions_data_prompt_lowepscenter_cut_all.Fill(-event.MandelT)
        phi_vs_t_data_prompt_lowepscenter_cut_all.Fill(event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Prompt_Data_lowepsleft_tree:
    if DATA_MMpi_Cut_lowepsleft(event) & Diamond_Cut(event):
        t_pions_data_prompt_lowepsleft_cut_all.Fill(-event.MandelT)
        phi_vs_t_data_prompt_lowepsleft_cut_all.Fill(event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Prompt_Data_highepsright_tree:
    if DATA_MMpi_Cut_highepsright(event) & Diamond_Cut(event):
        t_pions_data_prompt_highepsright_cut_all.Fill(-event.MandelT)
        phi_vs_t_data_prompt_highepsright_cut_all.Fill(event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Prompt_Data_highepscenter_tree:
    if DATA_MMpi_Cut_highepscenter(event) & Diamond_Cut(event):
        t_pions_data_prompt_highepscenter_cut_all.Fill(-event.MandelT)
        phi_vs_t_data_prompt_highepscenter_cut_all.Fill(event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Prompt_Data_highepsleft_tree:
    if DATA_MMpi_Cut_highepsleft(event) & Diamond_Cut(event):
        t_pions_data_prompt_highepsleft_cut_all.Fill(-event.MandelT)
        phi_vs_t_data_prompt_highepsleft_cut_all.Fill(event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Random_Data_lowepscenter_tree:
    if DATA_MMpi_Cut_lowepscenter(event) & Diamond_Cut(event):
        t_pions_data_random_lowepscenter_cut_all.Fill(-event.MandelT)
        phi_vs_t_data_random_lowepscenter_cut_all.Fill(event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Random_Data_lowepsleft_tree:
    if DATA_MMpi_Cut_lowepsleft(event) & Diamond_Cut(event):
        t_pions_data_random_lowepsleft_cut_all.Fill(-event.MandelT)
        phi_vs_t_data_random_lowepsleft_cut_all.Fill(event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Random_Data_highepsright_tree:
    if DATA_MMpi_Cut_highepsright(event) & Diamond_Cut(event):
        t_pions_data_random_highepsright_cut_all.Fill(-event.MandelT)
        phi_vs_t_data_random_highepsright_cut_all.Fill(event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Random_Data_highepscenter_tree:
    if DATA_MMpi_Cut_highepscenter(event) & Diamond_Cut(event):
        t_pions_data_random_highepscenter_cut_all.Fill(-event.MandelT)
        phi_vs_t_data_random_highepscenter_cut_all.Fill(event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Random_Data_highepsleft_tree:
    if DATA_MMpi_Cut_highepsleft(event) & Diamond_Cut(event):
        t_pions_data_random_highepsleft_cut_all.Fill(-event.MandelT)
        phi_vs_t_data_random_highepsleft_cut_all.Fill(event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Prompt_Dummy_lowepscenter_tree:
    if DATA_MMpi_Cut_lowepscenter(event) & Diamond_Cut(event):
        t_pions_dummy_prompt_lowepscenter_cut_all.Fill(-event.MandelT)
        phi_vs_t_dummy_prompt_lowepscenter_cut_all.Fill(event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Prompt_Dummy_lowepsleft_tree:
    if DATA_MMpi_Cut_lowepsleft(event) & Diamond_Cut(event):
        t_pions_dummy_prompt_lowepsleft_cut_all.Fill(-event.MandelT)
        phi_vs_t_dummy_prompt_lowepsleft_cut_all.Fill(event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Prompt_Dummy_highepsright_tree:
    if DATA_MMpi_Cut_highepsright(event) & Diamond_Cut(event):
        t_pions_dummy_prompt_highepsright_cut_all.Fill(-event.MandelT)
        phi_vs_t_dummy_prompt_highepsright_cut_all.Fill(event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Prompt_Dummy_highepscenter_tree:
    if DATA_MMpi_Cut_highepscenter(event) & Diamond_Cut(event):
        t_pions_dummy_prompt_highepscenter_cut_all.Fill(-event.MandelT)
        phi_vs_t_dummy_prompt_highepscenter_cut_all.Fill(event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Prompt_Dummy_highepsleft_tree:
    if DATA_MMpi_Cut_highepsleft(event) & Diamond_Cut(event):
        t_pions_dummy_prompt_highepsleft_cut_all.Fill(-event.MandelT)
        phi_vs_t_dummy_prompt_highepsleft_cut_all.Fill(event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Random_Dummy_lowepscenter_tree:
    if DATA_MMpi_Cut_lowepscenter(event) & Diamond_Cut(event):
        t_pions_dummy_random_lowepscenter_cut_all.Fill(-event.MandelT)
        phi_vs_t_dummy_random_lowepscenter_cut_all.Fill(event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Random_Dummy_lowepsleft_tree:
    if DATA_MMpi_Cut_lowepsleft(event) & Diamond_Cut(event):
        t_pions_dummy_random_lowepsleft_cut_all.Fill(-event.MandelT)
        phi_vs_t_dummy_random_lowepsleft_cut_all.Fill(event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Random_Dummy_highepsright_tree:
    if DATA_MMpi_Cut_highepsright(event) & Diamond_Cut(event):
        t_pions_dummy_random_highepsright_cut_all.Fill(-event.MandelT)
        phi_vs_t_dummy_random_highepsright_cut_all.Fill(event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Random_Dummy_highepscenter_tree:
    if DATA_MMpi_Cut_highepscenter(event) & Diamond_Cut(event):
        t_pions_dummy_random_highepscenter_cut_all.Fill(-event.MandelT)
        phi_vs_t_dummy_random_highepscenter_cut_all.Fill(event.ph_q, -event.MandelT)

for event in Cut_Pion_Events_Random_Dummy_highepsleft_tree:
    if DATA_MMpi_Cut_highepsleft(event) & Diamond_Cut(event):
        t_pions_dummy_random_highepsleft_cut_all.Fill(-event.MandelT)
        phi_vs_t_dummy_random_highepsleft_cut_all.Fill(event.ph_q, -event.MandelT)

print("####################################")
print("###### Histogram filling done ######")
print("####################################\n")

#################################################################################################################################################

# Random subtraction from Histograms
t_pions_data_random_lowepscenter_cut_all.Scale(1.0/nWindows)
t_pions_data_random_lowepsleft_cut_all.Scale(1.0/nWindows)
t_pions_data_random_highepsright_cut_all.Scale(1.0/nWindows)
t_pions_data_random_highepscenter_cut_all.Scale(1.0/nWindows)
t_pions_data_random_highepsleft_cut_all.Scale(1.0/nWindows)
phi_vs_t_data_random_lowepscenter_cut_all.Scale(1.0/nWindows)
phi_vs_t_data_random_lowepsleft_cut_all.Scale(1.0/nWindows)
phi_vs_t_data_random_highepsright_cut_all.Scale(1.0/nWindows)
phi_vs_t_data_random_highepscenter_cut_all.Scale(1.0/nWindows)
phi_vs_t_data_random_highepsleft_cut_all.Scale(1.0/nWindows)

t_pions_data_randsub_lowepscenter_cut_all.Add(t_pions_data_prompt_lowepscenter_cut_all, t_pions_data_random_lowepscenter_cut_all, 1, -1)
t_pions_data_randsub_lowepsleft_cut_all.Add(t_pions_data_prompt_lowepsleft_cut_all, t_pions_data_random_lowepsleft_cut_all, 1, -1)
t_pions_data_randsub_highepsright_cut_all.Add(t_pions_data_prompt_highepsright_cut_all, t_pions_data_random_highepsright_cut_all, 1, -1)
t_pions_data_randsub_highepscenter_cut_all.Add(t_pions_data_prompt_highepscenter_cut_all, t_pions_data_random_highepscenter_cut_all, 1, -1)
t_pions_data_randsub_highepsleft_cut_all.Add(t_pions_data_prompt_highepsleft_cut_all, t_pions_data_random_highepsleft_cut_all, 1, -1)
phi_vs_t_data_randsub_lowepscenter_cut_all.Add(phi_vs_t_data_prompt_lowepscenter_cut_all, phi_vs_t_data_random_lowepscenter_cut_all, 1, -1)
phi_vs_t_data_randsub_lowepsleft_cut_all.Add(phi_vs_t_data_prompt_lowepsleft_cut_all, phi_vs_t_data_random_lowepsleft_cut_all, 1, -1)
phi_vs_t_data_randsub_highepsright_cut_all.Add(phi_vs_t_data_prompt_highepsright_cut_all, phi_vs_t_data_random_highepsright_cut_all, 1, -1)
phi_vs_t_data_randsub_highepscenter_cut_all.Add(phi_vs_t_data_prompt_highepscenter_cut_all, phi_vs_t_data_random_highepscenter_cut_all, 1, -1)
phi_vs_t_data_randsub_highepsleft_cut_all.Add(phi_vs_t_data_prompt_highepsleft_cut_all, phi_vs_t_data_random_highepsleft_cut_all, 1, -1)

t_pions_dummy_random_lowepscenter_cut_all.Scale(1.0/nWindows)
t_pions_dummy_random_lowepsleft_cut_all.Scale(1.0/nWindows)
t_pions_dummy_random_highepsright_cut_all.Scale(1.0/nWindows)
t_pions_dummy_random_highepscenter_cut_all.Scale(1.0/nWindows)
t_pions_dummy_random_highepsleft_cut_all.Scale(1.0/nWindows)
phi_vs_t_dummy_random_lowepscenter_cut_all.Scale(1.0/nWindows)
phi_vs_t_dummy_random_lowepsleft_cut_all.Scale(1.0/nWindows)
phi_vs_t_dummy_random_highepsright_cut_all.Scale(1.0/nWindows)
phi_vs_t_dummy_random_highepscenter_cut_all.Scale(1.0/nWindows)
phi_vs_t_dummy_random_highepsleft_cut_all.Scale(1.0/nWindows)

t_pions_dummy_randsub_lowepscenter_cut_all.Add(t_pions_dummy_prompt_lowepscenter_cut_all, t_pions_dummy_random_lowepscenter_cut_all, 1, -1)
t_pions_dummy_randsub_lowepsleft_cut_all.Add(t_pions_dummy_prompt_lowepsleft_cut_all, t_pions_dummy_random_lowepsleft_cut_all, 1, -1)
t_pions_dummy_randsub_highepsright_cut_all.Add(t_pions_dummy_prompt_highepsright_cut_all, t_pions_dummy_random_highepsright_cut_all, 1, -1)
t_pions_dummy_randsub_highepscenter_cut_all.Add(t_pions_dummy_prompt_highepscenter_cut_all, t_pions_dummy_random_highepscenter_cut_all, 1, -1)
t_pions_dummy_randsub_highepsleft_cut_all.Add(t_pions_dummy_prompt_highepsleft_cut_all, t_pions_dummy_random_highepsleft_cut_all, 1, -1)
phi_vs_t_dummy_randsub_lowepscenter_cut_all.Add(phi_vs_t_dummy_prompt_lowepscenter_cut_all, phi_vs_t_dummy_random_lowepscenter_cut_all, 1, -1)
phi_vs_t_dummy_randsub_lowepsleft_cut_all.Add(phi_vs_t_dummy_prompt_lowepsleft_cut_all, phi_vs_t_dummy_random_lowepsleft_cut_all, 1, -1)
phi_vs_t_dummy_randsub_highepsright_cut_all.Add(phi_vs_t_dummy_prompt_highepsright_cut_all, phi_vs_t_dummy_random_highepsright_cut_all, 1, -1)
phi_vs_t_dummy_randsub_highepscenter_cut_all.Add(phi_vs_t_dummy_prompt_highepscenter_cut_all, phi_vs_t_dummy_random_highepscenter_cut_all, 1, -1)
phi_vs_t_dummy_randsub_highepsleft_cut_all.Add(phi_vs_t_dummy_prompt_highepsleft_cut_all, phi_vs_t_dummy_random_highepsleft_cut_all, 1, -1)

print("######################################")
print("###### Random substraction done ######")
print("######################################\n")

############################################################################################################################################

# Dummy Subtraction
t_pions_data_dummysub_lowepscenter_cut_all.Add(t_pions_data_randsub_lowepscenter_cut_all, t_pions_dummy_randsub_lowepscenter_cut_all, 1, -1)
t_pions_data_dummysub_lowepsleft_cut_all.Add(t_pions_data_randsub_lowepsleft_cut_all, t_pions_dummy_randsub_lowepsleft_cut_all, 1, -1)
t_pions_data_dummysub_highepsright_cut_all.Add(t_pions_data_randsub_highepsright_cut_all, t_pions_dummy_randsub_highepsright_cut_all, 1, -1)
t_pions_data_dummysub_highepscenter_cut_all.Add(t_pions_data_randsub_highepscenter_cut_all, t_pions_dummy_randsub_highepscenter_cut_all, 1, -1)
t_pions_data_dummysub_highepsleft_cut_all.Add(t_pions_data_randsub_highepsleft_cut_all, t_pions_dummy_randsub_highepsleft_cut_all, 1, -1)

phi_vs_t_data_dummysub_lowepscenter_cut_all.Add(phi_vs_t_data_randsub_lowepscenter_cut_all, phi_vs_t_dummy_randsub_lowepscenter_cut_all, 1, -1)
phi_vs_t_data_dummysub_lowepsleft_cut_all.Add(phi_vs_t_data_randsub_lowepsleft_cut_all, phi_vs_t_dummy_randsub_lowepsleft_cut_all, 1, -1)
phi_vs_t_data_dummysub_highepsright_cut_all.Add(phi_vs_t_data_randsub_highepsright_cut_all, phi_vs_t_dummy_randsub_highepsright_cut_all, 1, -1)
phi_vs_t_data_dummysub_highepscenter_cut_all.Add(phi_vs_t_data_randsub_highepscenter_cut_all, phi_vs_t_dummy_randsub_highepscenter_cut_all, 1, -1)
phi_vs_t_data_dummysub_highepsleft_cut_all.Add(phi_vs_t_data_randsub_highepsleft_cut_all, phi_vs_t_dummy_randsub_highepsleft_cut_all, 1, -1)

print("####################################")
print("###### Dummy subtraction done ######")
print("####################################\n")

###########################################################################################################################################

# Adding Low Eps and high Eps Settings
# Adding Low Eps and High Eps Settings
t_pions_data_dummysub_loweps_cut_all.Add(t_pions_data_dummysub_lowepscenter_cut_all, t_pions_data_dummysub_lowepsleft_cut_all, 1, 1)

t_pions_data_dummysub_higheps_cut_all.Add(t_pions_data_dummysub_highepsright_cut_all, t_pions_data_dummysub_highepscenter_cut_all, 1, 1)
t_pions_data_dummysub_higheps_cut_all.Add(t_pions_data_dummysub_highepsleft_cut_all,)

phi_vs_t_data_dummysub_loweps_cut_all.Add(phi_vs_t_data_dummysub_lowepscenter_cut_all, phi_vs_t_data_dummysub_lowepsleft_cut_all, 1, 1)
#phi_vs_t_data_dummysub_loweps_cut_all.Add(phi_vs_t_data_dummysub_lowepsleft_cut_all)

phi_vs_t_data_dummysub_higheps_cut_all.Add(phi_vs_t_data_dummysub_highepsright_cut_all)
phi_vs_t_data_dummysub_higheps_cut_all.Add(phi_vs_t_data_dummysub_highepscenter_cut_all)
phi_vs_t_data_dummysub_higheps_cut_all.Add(phi_vs_t_data_dummysub_highepsleft_cut_all)

print("####################################")
print("###### Histogram addition done ######")
print("####################################\n")

###########################################################################################################################################

#############################################################################################################################################

# 7 t-bins
#tbin_pions_min = [81,114,127,138,151,167,187]
#tbin_pions_max = [113,126,138,150,166,186,250]

# 5 t-bins
tbin_pions_min = [81,119,137,156,180]
tbin_pions_max = [118,136,155,179,250]

dN_data_pions_lowe = np.array([array.array('d', [0.0])] * len(tbin_pions_min))
dN_data_pions_highe = np.array([array.array('d', [0.0])] * len(tbin_pions_min))

N_data_pions_lowe = [t_pions_data_dummysub_loweps_cut_all.IntegralAndError(tbin_pions_min[i],tbin_pions_max[i],dN_data_pions_lowe[i],"") for i in range(len(tbin_pions_min))]
N_data_pions_highe = [t_pions_data_dummysub_higheps_cut_all.IntegralAndError(tbin_pions_min[i],tbin_pions_max[i],dN_data_pions_highe[i],"") for i in range(len(tbin_pions_min))]

### Writing values into dataframe
df_yields_pions = pd.DataFrame(columns=['t_min','t_max','Yield (loweps)','Error (loweps)','Yield (higheps)','Error (higheps)'])
for i in range(len(tbin_pions_min)):
     if (simc_t_resolution_pions_lowe >= (t_pions_data_dummysub_loweps_cut_all.GetXaxis().GetBinUpEdge(tbin_pions_max[i]) - t_pions_data_dummysub_loweps_cut_all.GetXaxis().GetBinLowEdge(tbin_pions_min[i])) or (simc_t_resolution_pions_highe >= (t_pions_data_dummysub_higheps_cut_all.GetXaxis().GetBinUpEdge(tbin_pions_max[i]) - t_pions_data_dummysub_higheps_cut_all.GetXaxis().GetBinLowEdge(tbin_pions_min[i])))):
        print("SIMC t resolution is smaller than the bin width")
        break
     df_yields_pions.loc[len(df_yields_pions.index)] = [t_pions_data_dummysub_loweps_cut_all.GetXaxis().GetBinLowEdge(tbin_pions_min[i]),t_pions_data_dummysub_loweps_cut_all.GetXaxis().GetBinUpEdge(tbin_pions_max[i]),N_data_pions_lowe[i],dN_data_pions_lowe[i,0],N_data_pions_highe[i],dN_data_pions_highe[i,0]]

df_yields_pions = df_yields_pions.round(3)
print("-"*40)
print(df_yields_pions.head(8))
print("-"*40)

yields_pions_path =  "%s/LTSep_CSVs/t_binning_csv/%s/%s_tbinning_yields_pions.csv" % (UTILPATH, physet_dir_name, PHY_SETTING)
df_yields_pions.to_csv(yields_pions_path, encoding='utf-8', index=False, header=True)

#############################################################################################################################################

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
# Create a TPaveText object for the table
table_pave = ROOT.TPaveText(0.1, 0.1, 0.95, 0.7, "NDC")  # Define the position (Normalized Device Coordinates)
table_pave.SetTextSize(0.035)
table_pave.SetTextAlign(12)  # Left alignment
table_pave.SetFillColor(0)  # Transparent background
table_pave.SetBorderSize(1)  # Border size
# Add the column headers
table_pave.AddText("t_min   t_max    Y(leps)    Err(leps)    Y(heps)    Err(heps)")
# Add the rows from the DataFrame
for index, row in df_yields_pions.iterrows():
    table_pave.AddText(f"{row['t_min']}    {row['t_max']}    {row['Yield (loweps)']}    {row['Error (loweps)']}    {row['Yield (higheps)']}    {row['Error (higheps)']}")
# Draw the table on the canvas
table_pave.Draw()
c1_delta.cd(2)
t_pions_data_dummysub_loweps_cut_all.SetLineColor(ROOT.kRed)
t_pions_data_dummysub_loweps_cut_all.Draw("Hist")
t_pions_data_dummysub_higheps_cut_all.SetLineColor(ROOT.kBlue)
t_pions_data_dummysub_higheps_cut_all.Draw("same, Hist")

# List to store TLine objects
lines = []
# Draw vertical lines for all t_min and t_max values
for t_min, t_max in zip(tbin_pions_min, tbin_pions_max):
    # Get the bin edges for t_min and t_max
    t_min_edge = t_pions_data_dummysub_loweps_cut_all.GetXaxis().GetBinLowEdge(t_min)
    t_max_edge = t_pions_data_dummysub_loweps_cut_all.GetXaxis().GetBinUpEdge(t_max)
#    print(f"Drawing lines for t_min_edge: {t_min_edge}, t_max_edge: {t_max_edge}")
    # Draw the line for t_min
    t_min_line = ROOT.TLine(t_min_edge, 0, t_min_edge, t_pions_data_dummysub_loweps_cut_all.GetMaximum())
    t_min_line.SetLineColor(ROOT.kBlack)
    t_min_line.SetLineStyle(1)  # Dashed line
    t_min_line.SetLineWidth(1)
    t_min_line.Draw("same")
    lines.append(t_min_line)  # Store the line to prevent garbage collection
    # Draw the line for t_max
    t_max_line = ROOT.TLine(t_max_edge, 0, t_max_edge, t_pions_data_dummysub_loweps_cut_all.GetMaximum())
    t_max_line.SetLineColor(ROOT.kBlack)
    t_max_line.SetLineStyle(1)  # Dashed line
    t_min_line.SetLineWidth(1)
    t_max_line.Draw("same")
    lines.append(t_max_line)  # Store the line to prevent garbage collection
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)  # Set legend position
legend.AddEntry(t_pions_data_dummysub_loweps_cut_all, "Low #epsilon Data", "l")  # Add red histogram
legend.AddEntry(t_pions_data_dummysub_higheps_cut_all, "High #epsilon Data", "l")  # Add blue histogram
legend.Draw("same")
c1_delta.cd(3)
# Draw the 2D histogram
#gPad.SetLogz()
# Set the range for the radial axis
tmin_loweps = phi_vs_t_data_dummysub_loweps_cut_all.GetYaxis().GetBinLowEdge(1)
tmax_loweps = phi_vs_t_data_dummysub_loweps_cut_all.GetYaxis().GetBinUpEdge(phi_vs_t_data_dummysub_loweps_cut_all.GetYaxis().GetNbins())
minrangeuser_loweps = tmin_loweps
maxrangeuser_loweps = tmax_loweps
#phi_vs_t_data_dummysub_loweps_cut_all.GetYaxis().SetRangeUser(minrangeuser_loweps, maxrangeuser_loweps)
# Remove the title from the histogram
phi_vs_t_data_dummysub_loweps_cut_all.SetTitle("")
phi_vs_t_data_dummysub_loweps_cut_all.SetStats(0)  # Disable stats box
phi_vs_t_data_dummysub_loweps_cut_all.SetStats(0)  # Disable stats box for the second histogram
#phi_vs_t_data_dummysub_loweps_cut_all.Draw("SURF2 POL")
phi_vs_t_data_dummysub_loweps_cut_all.Draw("SURF2Z1 POL")
# Set the color palette and viewing angles
gStyle.SetPalette(55)
gPad.SetTheta(90)
gPad.SetPhi(180)
# Add title for the plot
tvsphi_title_pion_loweps = TPaveText(0.02, 0.85, 0.29, 0.95, "NDC")  # Adjusted coordinates
tvsphi_title_pion_loweps.AddText("-t vs #phi - loweps")
tvsphi_title_pion_loweps.SetTextSize(0.04)
tvsphi_title_pion_loweps.SetTextAlign(22)
tvsphi_title_pion_loweps.Draw()
# Add labels for phi = 0
ptphizero_pion_loweps = TPaveText(0.923951, 0.513932, 0.993778, 0.574551, "NDC")
ptphizero_pion_loweps.AddText("#phi = 0")
ptphizero_pion_loweps.SetTextSize(0.03)
ptphizero_pion_loweps.SetTextAlign(22)
ptphizero_pion_loweps.SetLineWidth(2)
ptphizero_pion_loweps.Draw()
# Add line and label for phi = π/2
phihalfpi_pion_loweps = TLine(0, 0, 0, 0.6)
phihalfpi_pion_loweps.SetLineColor(ROOT.kBlack)
phihalfpi_pion_loweps.SetLineWidth(2)
phihalfpi_pion_loweps.Draw()
ptphihalfpi_pion_loweps = TPaveText(0.417855, 0.901876, 0.486574, 0.996358, "NDC")
ptphihalfpi_pion_loweps.AddText("#phi = #frac{#pi}{2}")
ptphihalfpi_pion_loweps.SetTextSize(0.03)
ptphihalfpi_pion_loweps.SetTextAlign(22)
ptphihalfpi_pion_loweps.Draw()
# Add line and label for phi = 3π/2
phithreepi_pion_loweps = TLine(0, 0, 0, -0.6)
phithreepi_pion_loweps.SetLineColor(ROOT.kBlack)
phithreepi_pion_loweps.SetLineWidth(2)
phithreepi_pion_loweps.Draw()
ptphithreepi_pion_loweps = TPaveText(0.419517, 0.00514928, 0.487128, 0.0996315, "NDC")
ptphithreepi_pion_loweps.AddText("#phi = #frac{3#pi}{2}")
ptphithreepi_pion_loweps.SetTextSize(0.03)
ptphithreepi_pion_loweps.SetTextAlign(22)
ptphithreepi_pion_loweps.Draw()
# Add line and label for phi = π
phipi_pion_loweps = TLine(0, 0, -0.6, 0)
phipi_pion_loweps.SetLineColor(ROOT.kBlack)
phipi_pion_loweps.SetLineWidth(2)
phipi_pion_loweps.Draw()
ptphipi_pion_loweps = TPaveText(0.0277092, 0.514217, 0.096428, 0.572746, "NDC")
ptphipi_pion_loweps.AddText("#phi = #pi")
ptphipi_pion_loweps.SetTextSize(0.03)
ptphipi_pion_loweps.SetTextAlign(22)
ptphipi_pion_loweps.Draw()
# Draw concentric arcs for t bins
Arc_pion_loweps = TArc()
for k in range(0, 9):
    Arc_pion_loweps.SetFillStyle(0)
#    Arc_pion_loweps.SetFillColor(ROOT.kWhite)  # Set the fill color to white
    Arc_pion_loweps.SetLineColor(ROOT.kBlack)
    Arc_pion_loweps.SetLineWidth(2)
    # To change the arc radius we have to change number 0.825 in the lower line.
#    Arc_pion_loweps.DrawArc(0, 0, 0.95 * (k + 1) / 10, 0., 360., "same") # for 6 arcs
#    Arc_pion_loweps.DrawArc(0, 0, 0.82 * (k + 1) / 10, 0., 360., "same") # for 7 arcs
    Arc_pion_loweps.DrawArc(0, 0, 0.64 * (k + 1) / 10, 0., 360., "same") # for 9 arcs
# Add radial axis
tradius_pion_loweps = TGaxis(0, 0, 0.575, 0, minrangeuser_loweps, maxrangeuser_loweps, 9, "-+N")
tradius_pion_loweps.SetMaxDigits(3)
tradius_pion_loweps.SetLineColor(2)
tradius_pion_loweps.SetLabelColor(2)
tradius_pion_loweps.SetLabelSize(0.03)
tradius_pion_loweps.Draw()
# Add line for phi = 0
phizero_pion_loweps = TLine(0, 0, 0.6, 0)
phizero_pion_loweps.SetLineColor(ROOT.kBlack)
phizero_pion_loweps.SetLineWidth(2)
phizero_pion_loweps.Draw()
c1_delta.cd(4)
# Draw the 2D histogram
#gPad.SetLogz()
# Set the range for the radial axis
tmin_higheps = phi_vs_t_data_dummysub_higheps_cut_all.GetYaxis().GetBinLowEdge(1)
tmax_higheps = phi_vs_t_data_dummysub_higheps_cut_all.GetYaxis().GetBinUpEdge(phi_vs_t_data_dummysub_higheps_cut_all.GetYaxis().GetNbins())
minrangeuser_higheps = tmin_higheps
maxrangeuser_higheps = tmax_higheps
#phi_vs_t_data_dummysub_higheps_cut_all.GetYaxis().SetRangeUser(minrangeuser_higheps, maxrangeuser_higheps)
# Remove the title from the histogram
phi_vs_t_data_dummysub_higheps_cut_all.SetTitle("")
phi_vs_t_data_dummysub_higheps_cut_all.SetStats(0)  # Disable stats box
#phi_vs_t_data_dummysub_higheps_cut_all.Draw("SURF2 POL")
phi_vs_t_data_dummysub_higheps_cut_all.Draw("SURF2Z1 POL")
# Set the color palette and viewing angles
gStyle.SetPalette(55)
gPad.SetTheta(90)
gPad.SetPhi(180)
# Add title for the plot
tvsphi_title_pion_higheps = TPaveText(0.02, 0.85, 0.29, 0.95, "NDC")  # Adjusted coordinates
tvsphi_title_pion_higheps.AddText("-t vs #phi - higheps")
tvsphi_title_pion_higheps.SetTextSize(0.04)
tvsphi_title_pion_higheps.SetTextAlign(22)
tvsphi_title_pion_higheps.Draw()
# Add labels for phi = 0
ptphizero_pion_higheps = TPaveText(0.923951, 0.513932, 0.993778, 0.574551, "NDC")
ptphizero_pion_higheps.AddText("#phi = 0")
ptphizero_pion_higheps.SetTextSize(0.03)
ptphizero_pion_higheps.SetTextAlign(22)
ptphizero_pion_higheps.SetLineWidth(2)
ptphizero_pion_higheps.Draw()
# Add line and label for phi = π/2
phihalfpi_pion_higheps = TLine(0, 0, 0, 0.6)
phihalfpi_pion_higheps.SetLineColor(ROOT.kBlack)
phihalfpi_pion_higheps.SetLineWidth(2)
phihalfpi_pion_higheps.Draw()
ptphihalfpi_pion_higheps = TPaveText(0.417855, 0.901876, 0.486574, 0.996358, "NDC")
ptphihalfpi_pion_higheps.AddText("#phi = #frac{#pi}{2}")
ptphihalfpi_pion_higheps.SetTextSize(0.03)
ptphihalfpi_pion_higheps.SetTextAlign(22)
ptphihalfpi_pion_higheps.Draw()
# Add line and label for phi = 3π/2
phithreepi_pion_higheps = TLine(0, 0, 0, -0.6)
phithreepi_pion_higheps.SetLineColor(ROOT.kBlack)
phithreepi_pion_higheps.SetLineWidth(2)
phithreepi_pion_higheps.Draw()
ptphithreepi_pion_higheps = TPaveText(0.419517, 0.00514928, 0.487128, 0.0996315, "NDC")
ptphithreepi_pion_higheps.AddText("#phi = #frac{3#pi}{2}")
ptphithreepi_pion_higheps.SetTextSize(0.03)
ptphithreepi_pion_higheps.SetTextAlign(22)
ptphithreepi_pion_higheps.Draw()
# Add line and label for phi = π
phipi_pion_higheps = TLine(0, 0, -0.6, 0)
phipi_pion_higheps.SetLineColor(ROOT.kBlack)
phipi_pion_higheps.SetLineWidth(2)
phipi_pion_higheps.Draw()
ptphipi_pion_higheps = TPaveText(0.0277092, 0.514217, 0.096428, 0.572746, "NDC")
ptphipi_pion_higheps.AddText("#phi = #pi")
ptphipi_pion_higheps.SetTextSize(0.03)
ptphipi_pion_higheps.SetTextAlign(22)
ptphipi_pion_higheps.Draw()
# Draw concentric arcs for t bins
Arc_pion_higheps = TArc()
for k in range(0, 9):
    Arc_pion_higheps.SetFillStyle(0)
#    Arc_pion_higheps.SetFillColor(ROOT.kWhite)  # Set the fill color to white
    Arc_pion_higheps.SetLineColor(ROOT.kBlack)
    Arc_pion_higheps.SetLineWidth(2)
    # To change the arc radius we have to change number 0.825 in the higher line.
#    Arc_pion_loweps.DrawArc(0, 0, 0.95 * (k + 1) / 10, 0., 360., "same") # for 6 arcs
#    Arc_pion_loweps.DrawArc(0, 0, 0.82 * (k + 1) / 10, 0., 360., "same") # for 7 arcs
    Arc_pion_loweps.DrawArc(0, 0, 0.64 * (k + 1) / 10, 0., 360., "same") # for 9 arcs
# Add radial axis
tradius_pion_higheps = TGaxis(0, 0, 0.575, 0, minrangeuser_higheps, maxrangeuser_higheps, 9, "-+N")
tradius_pion_higheps.SetMaxDigits(3)
tradius_pion_higheps.SetLineColor(2)
tradius_pion_higheps.SetLabelColor(2)
tradius_pion_higheps.SetLabelSize(0.03)
tradius_pion_higheps.Draw()
# Add line for phi = 0
phizero_pion_higheps = TLine(0, 0, 0.6, 0)
phizero_pion_higheps.SetLineColor(ROOT.kBlack)
phizero_pion_higheps.SetLineWidth(2)
phizero_pion_higheps.Draw()
c1_delta.Print(Pion_Analysis_Distributions)

#############################################################################################################################################

# Making directories in output file
outHistFile = ROOT.TFile.Open("%s/%s_%s_tbinning_Data.root" % (OUTPATH, PHY_SETTING, MaxEvent) , "RECREATE")
d_Cut_Pion_Events_Prompt_Data_Cut_All = outHistFile.mkdir("Cut_Pion_Events_Prompt_Data_Cut_All")
d_Cut_Pion_Events_Random_Data_Cut_All = outHistFile.mkdir("Cut_Pion_Events_Random_Data_Cut_All")
d_Cut_Pion_Events_RandomSub_Data = outHistFile.mkdir("Cut_Pion_Events_RandomSub_Data")
d_Cut_Pion_Events_DummySub_RandomSub_Data = outHistFile.mkdir("Cut_Pion_Events_DummySub_RandomSub_Data")
d_Cut_Pion_Events_SIMC_Data = outHistFile.mkdir("Cut_Pion_Events_SIMC_Data")
d_Cut_Pion_Events_Combined_RandomSub_Data = outHistFile.mkdir("Cut_Pion_Events_Combined_RandomSub_Data")

# Writing Histograms for pions     
d_Cut_Pion_Events_Prompt_Data_Cut_All.cd()
t_pions_data_prompt_lowepscenter_cut_all.Write()
t_pions_data_prompt_lowepsleft_cut_all.Write()
t_pions_data_prompt_highepsright_cut_all.Write()
t_pions_data_prompt_highepscenter_cut_all.Write()
t_pions_data_prompt_highepsleft_cut_all.Write()
t_pions_dummy_prompt_lowepscenter_cut_all.Write()
t_pions_dummy_prompt_lowepsleft_cut_all.Write()
t_pions_dummy_prompt_highepsright_cut_all.Write()
t_pions_dummy_prompt_highepscenter_cut_all.Write()
t_pions_dummy_prompt_highepsleft_cut_all.Write()
phi_vs_t_data_prompt_lowepscenter_cut_all.Write()
phi_vs_t_data_prompt_lowepsleft_cut_all.Write()
phi_vs_t_data_prompt_highepsright_cut_all.Write()
phi_vs_t_data_prompt_highepscenter_cut_all.Write()
phi_vs_t_data_prompt_highepsleft_cut_all.Write()
phi_vs_t_dummy_prompt_lowepscenter_cut_all.Write()
phi_vs_t_dummy_prompt_lowepsleft_cut_all.Write()
phi_vs_t_dummy_prompt_highepsright_cut_all.Write()
phi_vs_t_dummy_prompt_highepscenter_cut_all.Write()
phi_vs_t_dummy_prompt_highepsleft_cut_all.Write()

d_Cut_Pion_Events_Random_Data_Cut_All.cd()
t_pions_data_random_lowepscenter_cut_all.Write()
t_pions_data_random_lowepsleft_cut_all.Write()
t_pions_data_random_highepsright_cut_all.Write()
t_pions_data_random_highepscenter_cut_all.Write()
t_pions_data_random_highepsleft_cut_all.Write()
t_pions_dummy_random_lowepscenter_cut_all.Write()
t_pions_dummy_random_lowepsleft_cut_all.Write()
t_pions_dummy_random_highepsright_cut_all.Write()
t_pions_dummy_random_highepscenter_cut_all.Write()
t_pions_dummy_random_highepsleft_cut_all.Write()
phi_vs_t_data_random_lowepscenter_cut_all.Write()
phi_vs_t_data_random_lowepsleft_cut_all.Write()
phi_vs_t_data_random_highepsright_cut_all.Write()
phi_vs_t_data_random_highepscenter_cut_all.Write()
phi_vs_t_data_random_highepsleft_cut_all.Write()
phi_vs_t_dummy_random_lowepscenter_cut_all.Write()
phi_vs_t_dummy_random_lowepsleft_cut_all.Write()
phi_vs_t_dummy_random_highepsright_cut_all.Write()
phi_vs_t_dummy_random_highepscenter_cut_all.Write()
phi_vs_t_dummy_random_highepsleft_cut_all.Write()

d_Cut_Pion_Events_RandomSub_Data.cd()
t_pions_data_randsub_lowepscenter_cut_all.Write()
t_pions_data_randsub_lowepsleft_cut_all.Write()
t_pions_data_randsub_highepsright_cut_all.Write()
t_pions_data_randsub_highepscenter_cut_all.Write()
t_pions_data_randsub_highepsleft_cut_all.Write()
t_pions_dummy_randsub_lowepscenter_cut_all.Write()
t_pions_dummy_randsub_lowepsleft_cut_all.Write()
t_pions_dummy_randsub_highepsright_cut_all.Write()
t_pions_dummy_randsub_highepscenter_cut_all.Write()
t_pions_dummy_randsub_highepsleft_cut_all.Write()
phi_vs_t_data_randsub_lowepscenter_cut_all.Write()
phi_vs_t_data_randsub_lowepsleft_cut_all.Write()
phi_vs_t_data_randsub_highepsright_cut_all.Write()
phi_vs_t_data_randsub_highepscenter_cut_all.Write()
phi_vs_t_data_randsub_highepsleft_cut_all.Write()
phi_vs_t_dummy_randsub_lowepscenter_cut_all.Write()
phi_vs_t_dummy_randsub_lowepsleft_cut_all.Write()
phi_vs_t_dummy_randsub_highepsright_cut_all.Write()
phi_vs_t_dummy_randsub_highepscenter_cut_all.Write()
phi_vs_t_dummy_randsub_highepsleft_cut_all.Write()

d_Cut_Pion_Events_DummySub_RandomSub_Data.cd()
t_pions_data_dummysub_lowepscenter_cut_all.Write()
t_pions_data_dummysub_lowepsleft_cut_all.Write()
t_pions_data_dummysub_highepsright_cut_all.Write()
t_pions_data_dummysub_highepscenter_cut_all.Write()
t_pions_data_dummysub_highepsleft_cut_all.Write()
phi_vs_t_data_dummysub_lowepscenter_cut_all.Write()
phi_vs_t_data_dummysub_lowepsleft_cut_all.Write()
phi_vs_t_data_dummysub_highepsright_cut_all.Write()
phi_vs_t_data_dummysub_highepscenter_cut_all.Write()
phi_vs_t_data_dummysub_highepsleft_cut_all.Write()

d_Cut_Pion_Events_Combined_RandomSub_Data.cd()
t_pions_data_dummysub_loweps_cut_all.Write()
t_pions_data_dummysub_higheps_cut_all.Write()
phi_vs_t_data_dummysub_loweps_cut_all.Write()
phi_vs_t_data_dummysub_higheps_cut_all.Write()

print("####################################")
print("###### Histogram writing done ######")
print("####################################\n")

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
outHistFile.Close()

print ("Processing Complete")