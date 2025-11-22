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
import uncertainties as u
from ctypes import c_double

##################################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=7:
    print("!!!!! ERROR !!!!!\n Expected 7 arguments\n Usage is with - ROOTfileSuffixs Beam Energy MaxEvents RunList CVSFile\n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################

# Input params - run number and max number of events
PHY_SETTING = sys.argv[1]
MaxEvent = sys.argv[2]
DATA_Suffix = sys.argv[3]
DUMMY_Suffix = sys.argv[4]
DATA_RUN_LIST = sys.argv[5]
DUMMY_RUN_LIST = sys.argv[6]
CSV_FILE = sys.argv[7]

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
EFF_CSV     = "/u/group/c-pionlt/USERS/%s/hallc_replay_lt/UTIL_PION/efficiencies" % (USER)
MMCUT_CSV   = "/u/group/c-pionlt/USERS/%s/hallc_replay_lt/UTIL_PION/LTSep_CSVs/mm_offset_cut_csv" % (USER)
DCUT_CSV    = "/u/group/c-pionlt/USERS/%s/hallc_replay_lt/UTIL_PION/LTSep_CSVs/diamond_cut_csv" % (USER)
TBINCSVPATH = "/u/group/c-pionlt/USERS/%s/hallc_replay_lt/UTIL_PION/LTSep_CSVs/t_binning_csv" % (USER)

# Extract the first three words from PHY_SETTING for the CSV file name
setting_name = "_".join(PHY_SETTING.split("_")[:3])
physet_dir_name = "%s_std" % (setting_name)

#################################################################################################################################################

# Output PDF File Name
print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))
#Pion_Analysis_Distributions = "%s/%s_%s_ProdCoin_Pion_Analysis_DataYield_Distributions.pdf" % (OUTPATH, PHY_SETTING, MaxEvent)

# Input file location and variables taking
rootFile_DATA = "%s/%s_%s_%s.root" % (OUTPATH, PHY_SETTING, MaxEvent, DATA_Suffix)
rootFile_DUMMY = "%s/%s_%s_%s.root" % (OUTPATH, PHY_SETTING, MaxEvent, DUMMY_Suffix)
data_run_list = "%s/%s" % (RUNLISTPATH, DATA_RUN_LIST)
dummy_run_list = "%s/%s" % (RUNLISTPATH, DUMMY_RUN_LIST)
eff_csv_file = "%s/%s.csv" % (EFF_CSV, CSV_FILE)
mmcut_csv_file = "%s/%s/%s_mm_offsets_cuts_parameters.csv" % (MMCUT_CSV, physet_dir_name, setting_name)
dcut_csv_file = "%s/%s/%s_diamond_cut_parameters.csv" % (DCUT_CSV, physet_dir_name, setting_name)
tbin_csv_file  = "%s/%s/%s_tbinning_yields_pions.csv" % (TBINCSVPATH, physet_dir_name, setting_name)

###################################################################################################################################################

print ('\nPhysics Setting = ',PHY_SETTING, '\n')
print("-"*40)

# Cuts for Pions Selection - Change according to your data
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
Diamond_Cut = lambda event: (cutg_diamond.IsInside(event.Q2, event.W))

#------------------------------------------------------------------------------------------------------------------------------------

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

#-----------------------------------------------------------------------------------------------------------------------------------------

# t-binning cut
# Read the CSV file
tbin_df = pd.read_csv(tbin_csv_file)

# Extract the t_min and t_max columns as arrays
t_min_values = tbin_df['t_min'].values
t_max_values = tbin_df['t_max'].values

# Define the cuts using row 1 and row 3 values
tbin_Cut1 = lambda event: (t_min_values[0] <= -event.MandelT <= t_max_values[0])  # Row 1
tbin_Cut2 = lambda event: (t_min_values[1] < -event.MandelT <= t_max_values[1])  # Row 2
tbin_Cut3 = lambda event: (t_min_values[2] < -event.MandelT <= t_max_values[2])  # Row 3
tbin_Cut4 = lambda event: (t_min_values[3] < -event.MandelT <= t_max_values[3])  # Row 4
tbin_Cut5 = lambda event: (t_min_values[4] < -event.MandelT <= t_max_values[4])  # Row 5

# Bundle them into a list
tbin_cuts = [tbin_Cut1, tbin_Cut2, tbin_Cut3, tbin_Cut4, tbin_Cut5]

# Print the t-binning cuts with 3 decimal places
print("\n t-binning cuts:")
print(f"tbin_Cut1: t_min = {t_min_values[0]:.3f}, t_max = {t_max_values[0]:.3f}")
print(f"tbin_Cut2: t_min = {t_min_values[1]:.3f}, t_max = {t_max_values[1]:.3f}")
print(f"tbin_Cut3: t_min = {t_min_values[2]:.3f}, t_max = {t_max_values[2]:.3f}")
print(f"tbin_Cut4: t_min = {t_min_values[3]:.3f}, t_max = {t_max_values[3]:.3f}")
print(f"tbin_Cut5: t_min = {t_min_values[4]:.3f}, t_max = {t_max_values[4]:.3f}")
print("-"*40)

# Define phi bins (15 bins from 0 to 360 degrees, each 24 degrees wide)
phi_bins = [i for i in range(0, 361, 24)]  # 0, 24, 48, ..., 360

# Calculate total number of t-bins and phi bins
total_tbins = len(tbin_cuts)
total_phibins = len(phi_bins) - 1
#print(f"Total t-bins: {total_tbins}, Total phi bins: {total_phibins}")

#########################################################################################################################################

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

###################################################################################################################################################

print("-"*40)
print ('Calculating effective charge for the data and dummy run list')
print("-"*40)

# Read CSV File and calculate total charge in mC for normalization
# Input runlists and csv files
data_run_list_file = (data_run_list)
dummy_run_list_file = (dummy_run_list)
eff_csv_file_name = (eff_csv_file)

# Read run numbers from the run list file
with open(data_run_list_file, 'r') as run_list_file_1:
    data_runs = [line.strip() for line in run_list_file_1 if line.strip()]
with open(dummy_run_list_file, 'r') as run_list_file_2:
    dummy_runs = [line.strip() for line in run_list_file_2 if line.strip()]

# Read CSV file using pandas
df = pd.read_csv(eff_csv_file_name)

# Filter DataFrame to include only rows with run numbers in the run list
filtered_data_df = df[df['Run_Number'].astype(str).str.replace('.0', '', regex=False).isin(data_runs)]
filtered_dummy_df = df[df['Run_Number'].astype(str).str.replace('.0', '', regex=False).isin(dummy_runs)]

# Intializing Variables to calculate the product of efficiency and charge
total_data_effective_charge = 0.0
total_dummy_effective_charge = 0.0
total_data_effective_charge_sum = 0.0
total_dummy_effective_charge_sum = 0.0
total_data_effective_charge_error_sum = 0.0
total_dummy_effective_charge_error_sum = 0.0
total_dummy_effective_charge_cal = 0.0
total_dummy_effective_charge_cal_error = 0.0

# Detector (HMS Cer + HMS Cal) Efficiencies
hms_Cer_detector_efficiency = 0.9981
hms_Cal_detector_efficiency = 0.9981
hms_Cer_detector_efficiency_error = 0.0001
hms_Cal_detector_efficiency_error = 0.0001

if setting_name == "Q3p85_W2p62_t0p21":
    #RF cut Efficiency
    RF_efficiency = 0.998
    RF_efficiency_error = 0.00007
    # Pion Absorption Correction
    pion_absorption_correction = 0.9654
    pion_absorption_correction_error = 0.00005
else:
    print("!!!!!\n Need to provide RF Eff and Absortion Corr for %s\n!!!!!" % (setting_name))
    RF_efficiency = 1.0
    RF_efficiency_error = 0.0
    pion_absorption_correction = 1.0
    pion_absorption_correction_error = 0.0
    print("Assuming RF Eff = 1.0 +/- 0.0 and Pion Absorption Corr = 1.0 +/- 0.0 \n")

row_data = []

# Print charge values for each run
for index, row in filtered_data_df.iterrows():
    data_charge = row['Corrected_Charge'] 
    data_charge_error = row['Corrected_Charge_ERROR']
    data_Boiling_factor = row['Target_BoilingCorr']
    data_Boiling_factor_error = row['Target_BoilingCorr_ERROR']
    data_edtm_livetime_Corr = row['Non_Scaler_EDTM_Live_Time_Corr']
    data_edtm_livetime_Corr_error = row['Non_Scaler_EDTM_Live_Time_Corr_ERROR']
    data_coinblocking_factor = row['CoinBlocking']
    data_coinblocking_factor_error = row['CoinBlocking_ERROR']

    data_hms_tracking_efficiency = row['HMS_Elec_SING_TRACK_EFF'] 
    data_shms_tracking_efficiency = row['SHMS_Pion_SING_TRACK_EFF'] 
    data_hms_tracking_efficiency_error = row['HMS_Elec_SING_TRACK_EFF_ERROR']
    data_shms_tracking_efficiency_error = row['SHMS_Pion_SING_TRACK_EFF_ERROR']

    data_hms_hodo_3_of_4_efficiency = row['HMS_Hodo_3_of_4_EFF']
    data_shms_hodo_3_of_4_efficiency = row['SHMS_Hodo_3_of_4_EFF']
    data_hms_hodo_3_of_4_efficiency_error = row['HMS_Elec_SING_TRACK_EFF_ERROR']
    data_shms_hodo_3_of_4_efficiency_error = row['SHMS_Pion_SING_TRACK_EFF_ERROR']
#    data_hms_hodo_3_of_4_efficiency_error = row['HMS_Hodo_3_of_4_EFF_ERROR']
#    data_shms_hodo_3_of_4_efficiency_error = row['SHMS_Hodo_3_of_4_EFF_ERROR']

    data_shms_aero_detector_efficiency = row['SHMS_Aero_COIN_Pion_Eff']
    data_shms_aero_detector_efficiency_error = row['SHMS_Aero_COIN_Pion_Eff_ERROR']

    data_product = (data_charge * pion_absorption_correction * data_hms_tracking_efficiency * data_shms_tracking_efficiency * RF_efficiency * hms_Cer_detector_efficiency * hms_Cal_detector_efficiency * data_hms_hodo_3_of_4_efficiency * data_shms_hodo_3_of_4_efficiency * data_shms_aero_detector_efficiency * data_edtm_livetime_Corr * data_Boiling_factor * data_coinblocking_factor)
    data_product_error = data_product * (math.sqrt((data_charge_error/data_charge)** 2 +  (pion_absorption_correction_error/pion_absorption_correction)**2 + (data_hms_tracking_efficiency_error/data_hms_tracking_efficiency)** 2 + (data_shms_tracking_efficiency_error/data_shms_tracking_efficiency)** 2 + (RF_efficiency_error/RF_efficiency)**2 + (data_edtm_livetime_Corr_error/data_edtm_livetime_Corr)** 2 + (hms_Cer_detector_efficiency_error/hms_Cer_detector_efficiency)** 2 + (hms_Cal_detector_efficiency_error/hms_Cal_detector_efficiency)** 2 + (data_hms_hodo_3_of_4_efficiency_error/data_hms_hodo_3_of_4_efficiency)** 2 + (data_shms_hodo_3_of_4_efficiency_error/data_shms_hodo_3_of_4_efficiency)**2 + (data_shms_aero_detector_efficiency_error/data_shms_aero_detector_efficiency)**2 + (data_Boiling_factor_error/data_Boiling_factor)** 2 + (data_coinblocking_factor_error/data_coinblocking_factor)**2)) 

    total_data_effective_charge_sum += data_product
    total_data_effective_charge_error_sum += (data_product_error)** 2

    print('Charge for data run: {:<10} Data Charge: {:.3f} Pion_Absorpt_Corr: {:.3f} HMS Tracking Eff: {:.3f} SHMS Tracking Eff: {:.3f} RF Eff: {:.3f} HMS Cer Detector Eff: {:.3f} HMS Cal Detector_Eff: {:.3f} HMS Hodo 3/4 Eff: {:.3f} SHMS Hodo 3/4 Eff: {:.3f} SHMS Aero Detector Eff: {:.3f} EDTM Live Time: {:.3f} Boiling Correction: {:.3f} Coin Blocking: {:.3f} Product: {:.3f}'.format(row["Run_Number"], data_charge, pion_absorption_correction, data_hms_tracking_efficiency, data_shms_tracking_efficiency, RF_efficiency, hms_Cer_detector_efficiency, hms_Cal_detector_efficiency, data_hms_hodo_3_of_4_efficiency, data_shms_hodo_3_of_4_efficiency, data_shms_aero_detector_efficiency, data_edtm_livetime_Corr, data_Boiling_factor, data_coinblocking_factor, data_product))

    row_data.append({'Run_Number':row["Run_Number"], 'charge':data_charge, 'charge_error':data_charge_error, 'pion_absorption_correction':pion_absorption_correction, 'HMS_Tracking_Eff':data_hms_tracking_efficiency, 'HMS_Tracking_Eff_error':data_hms_tracking_efficiency_error, 'SHMS_Tracking_Eff':data_shms_tracking_efficiency, 'SHMS_Tracking_Eff_error':data_shms_tracking_efficiency_error, 'RF_Eff':RF_efficiency, 'RF_Eff_error':RF_efficiency_error, 'HMS_Cer_Detector_Eff':hms_Cer_detector_efficiency, 'HMS_Cer_Detector_Eff_error':hms_Cer_detector_efficiency_error, 'HMS_Cal_Detector_Eff':hms_Cal_detector_efficiency, 'HMS_Cal_Detector_Eff_error':hms_Cal_detector_efficiency_error, 'HMS_Hodo_3_4_Eff':data_hms_hodo_3_of_4_efficiency, 'HMS_Hodo_3_4_Eff_error':data_hms_hodo_3_of_4_efficiency_error, 'SHMS_Hodo_3_4_Eff':data_shms_hodo_3_of_4_efficiency, 'SHMS_Hodo_3_4_Eff_error':data_shms_hodo_3_of_4_efficiency_error, 'SHMS Aerogel Eff':data_shms_aero_detector_efficiency, 'SHMS Aerogel Eff_error':data_shms_aero_detector_efficiency_error, 'EDTM_Live_Time':data_edtm_livetime_Corr, 'EDTM_Live_Time_error':data_edtm_livetime_Corr_error, 'Boiling_factor':data_Boiling_factor, 'Boiling_factor_error':data_Boiling_factor_error, 'Coin_Blocking':data_coinblocking_factor, 'Coin_Blocking_error':data_coinblocking_factor_error, 'effective_charge':data_product, 'effective_charge_error':data_product_error})

print("-"*40)

total_data_effective_charge = total_data_effective_charge_sum
total_data_effective_charge_error = math.sqrt(total_data_effective_charge_error_sum)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

row_dummy = []

for index, row in filtered_dummy_df.iterrows():
    dummy_charge = row['Corrected_Charge'] 
    dummy_charge_error = row['Corrected_Charge_ERROR']
    dummy_edtm_livetime_Corr = row['Non_Scaler_EDTM_Live_Time_Corr']
    dummy_edtm_livetime_Corr_error = row['Non_Scaler_EDTM_Live_Time_Corr_ERROR']
    dummy_coinblocking_factor = row['CoinBlocking']
    dummy_coinblocking_factor_error = row['CoinBlocking_ERROR']

    dummy_hms_tracking_efficiency = row['HMS_Elec_SING_TRACK_EFF'] 
    dummy_shms_tracking_efficiency = row['SHMS_Pion_SING_TRACK_EFF'] 
    dummy_hms_tracking_efficiency_error = row['HMS_Elec_SING_TRACK_EFF_ERROR']
    dummy_shms_tracking_efficiency_error = row['SHMS_Pion_SING_TRACK_EFF_ERROR']

    dummy_hms_hodo_3_of_4_efficiency = row['HMS_Hodo_3_of_4_EFF']
    dummy_shms_hodo_3_of_4_efficiency = row['SHMS_Hodo_3_of_4_EFF']
    dummy_hms_hodo_3_of_4_efficiency_error = row['HMS_Elec_SING_TRACK_EFF_ERROR']
    dummy_shms_hodo_3_of_4_efficiency_error = row['SHMS_Pion_SING_TRACK_EFF_ERROR']
#    dummy_hms_hodo_3_of_4_efficiency_error = row['HMS_Hodo_3_of_4_EFF_ERROR']
#    dummy_shms_hodo_3_of_4_efficiency_error = row['SHMS_Hodo_3_of_4_EFF_ERROR']

    dummy_shms_aero_detector_efficiency = row['SHMS_Aero_COIN_Pion_Eff']
    dummy_shms_aero_detector_efficiency_error = row['SHMS_Aero_COIN_Pion_Eff_ERROR']

    dummy_product = (dummy_charge * dummy_hms_tracking_efficiency * dummy_shms_tracking_efficiency * RF_efficiency * hms_Cer_detector_efficiency * hms_Cal_detector_efficiency * dummy_hms_hodo_3_of_4_efficiency * dummy_shms_hodo_3_of_4_efficiency * dummy_edtm_livetime_Corr * dummy_shms_aero_detector_efficiency * dummy_coinblocking_factor)
    dummy_product_error = dummy_product * (math.sqrt(((dummy_charge_error/dummy_charge) ** 2 + dummy_hms_tracking_efficiency_error/dummy_hms_tracking_efficiency) ** 2 + (dummy_shms_tracking_efficiency_error/dummy_shms_tracking_efficiency) ** 2 + (RF_efficiency_error/RF_efficiency) **2 + (dummy_edtm_livetime_Corr_error/dummy_edtm_livetime_Corr)** 2 + (hms_Cer_detector_efficiency_error/hms_Cer_detector_efficiency) ** 2 + (hms_Cal_detector_efficiency_error/hms_Cal_detector_efficiency) ** 2 + (dummy_hms_hodo_3_of_4_efficiency_error/dummy_hms_hodo_3_of_4_efficiency) ** 2 + (dummy_shms_hodo_3_of_4_efficiency_error/dummy_shms_hodo_3_of_4_efficiency) ** 2 + (dummy_shms_aero_detector_efficiency_error/dummy_shms_aero_detector_efficiency) ** 2 + (dummy_coinblocking_factor_error/dummy_coinblocking_factor)**2))

    total_dummy_effective_charge_sum += dummy_product
    total_dummy_effective_charge_error_sum += (dummy_product_error)** 2

    print('Charge for dummy run: {:<10} Dummy Charge: {:.3f} HMS Tracking Eff: {:.3f} SHMS Tracking Eff: {:.3f} RF Eff: {:.3f} HMS Cer Detector Eff: {:.3f} HMS Cal Detector_Eff: {:.3f} HMS Hodo 3/4 Eff: {:.3f} SHMS Hodo 3/4 Eff: {:.3f} EDTM Live Time: {:.3f} SHMS Aerogel Eff: {:.3f} Coin Blocking: {:.3f} Product: {:.3f}'.format(row["Run_Number"], dummy_charge, dummy_hms_tracking_efficiency, dummy_shms_tracking_efficiency, RF_efficiency, hms_Cer_detector_efficiency, hms_Cal_detector_efficiency, dummy_hms_hodo_3_of_4_efficiency, dummy_shms_hodo_3_of_4_efficiency, dummy_edtm_livetime_Corr, dummy_shms_aero_detector_efficiency, dummy_coinblocking_factor, dummy_product))

    row_dummy.append({'Run_Number':row["Run_Number"], 'charge':dummy_charge, 'charge_error':dummy_charge_error, 'HMS_Tracking_Eff':dummy_hms_tracking_efficiency, 'HMS_Tracking_Eff_error':dummy_hms_tracking_efficiency_error, 'SHMS_Tracking_Eff':dummy_shms_tracking_efficiency, 'SHMS_Tracking_Eff_error':dummy_shms_tracking_efficiency_error, 'RF_Eff':RF_efficiency, 'RF_Eff_error':RF_efficiency_error, 'HMS_Cer_Detector_Eff':hms_Cer_detector_efficiency, 'HMS_Cer_Detector_Eff_error':hms_Cer_detector_efficiency_error, 'HMS_Cal_Detector_Eff':hms_Cal_detector_efficiency, 'HMS_Cal_Detector_Eff_error':hms_Cal_detector_efficiency_error, 'HMS_Hodo_3_4_Eff':dummy_hms_hodo_3_of_4_efficiency, 'HMS_Hodo_3_4_Eff_error':dummy_hms_hodo_3_of_4_efficiency_error, 'SHMS_Hodo_3_4_Eff':dummy_shms_hodo_3_of_4_efficiency, 'SHMS_Hodo_3_4_Eff_error':dummy_shms_hodo_3_of_4_efficiency_error, 'SHMS Aerogel Eff':dummy_shms_aero_detector_efficiency, 'SHMS Aerogel Eff_error':dummy_shms_aero_detector_efficiency_error, 'EDTM_Live_Time':dummy_edtm_livetime_Corr, 'EDTM_Live_Time_error':dummy_edtm_livetime_Corr_error, 'Boiling_factor':'NA', 'Boiling_factor_error':'NA', 'Coin_Blocking':dummy_coinblocking_factor, 'Coin_Blocking_error':dummy_coinblocking_factor_error, 'effective_charge':dummy_product, 'effective_charge_error':dummy_product_error})

print("-"*40)

# Dummy Target Thickness Correction
#dummy_target_corr = 4.8579 # KaonLT
dummy_target_corr = 3.527 # PionLT
dummy_target_corr_error = 0.227 # PionLT

total_dummy_effective_charge_cal = total_dummy_effective_charge_sum
total_dummy_effective_charge_cal_error = math.sqrt(total_dummy_effective_charge_error_sum)

total_dummy_effective_charge = (total_dummy_effective_charge_cal * dummy_target_corr)
total_dummy_effective_charge_error = total_dummy_effective_charge * math.sqrt((total_dummy_effective_charge_cal_error/total_dummy_effective_charge_cal)**2 + (dummy_target_corr_error/dummy_target_corr)** 2)

print("\nTotal effective charge for the data run list: {:.5f} ± {:.5f}".format(total_data_effective_charge, total_data_effective_charge_error))
print("\nTotal effective charge for the dummy run list: {:.5f} ± {:.5f}".format(total_dummy_effective_charge, total_dummy_effective_charge_error))

# Normalization factor Calculation for Data and Dummy
normfac_data = 1.0/(total_data_effective_charge)
normfac_dummy = 1.0/(total_dummy_effective_charge)

print("-"*40)
print ("normfac_data :", normfac_data)
print ("normfac_dummy: ", normfac_dummy)
print("-"*40)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Define the output CSV file name
eff_cal_output_csv_path = "%s/LTSep_CSVs/physics_yields_csv/%s/%s_Physics_Norm_Eff_Charge_Calculation.csv" % (UTILPATH, physet_dir_name, setting_name)

# Prepare the header for the CSV file
header = [
    "Physics_Setting", "Run Type", "Run Number", "Charge", "Charge Error",
    "HMS Tracking Eff", "HMS Tracking Eff Error", "SHMS Tracking Eff", "SHMS Tracking Eff Error",
    "HMS Cer Detector Eff", "HMS Cer Detector Eff Error", "HMS Cal Detector Eff", "HMS Cal Detector Eff Error",
    "HMS Hodo 3/4 Eff", "HMS Hodo 3/4 Eff Error", "SHMS Hodo 3/4 Eff", "SHMS Hodo 3/4 Eff Error",
    "EDTM Live Time", "EDTM Live Time Error", "Boiling Factor", "Boiling Factor Error", "Dummy Target Correction",
    "Dummy Target Correction Error", "Effective Charge", "Effective Charge Error", "Total Effective Charge", "Total Effective Charge Error", "Normfactor"
]

# Combine data and dummy rows
combined_rows = []

# Add data rows
for row in row_data:
    combined_rows.append({
        "Physics_Setting": PHY_SETTING,
        "Run Type": "data",
        "Run Number": row["Run_Number"],
        "Charge": row["charge"],
        "Charge Error": row["charge_error"],
        "HMS Tracking Eff": row["HMS_Tracking_Eff"],
        "HMS Tracking Eff Error": row["HMS_Tracking_Eff_error"],
        "SHMS Tracking Eff": row["SHMS_Tracking_Eff"],
        "SHMS Tracking Eff Error": row["SHMS_Tracking_Eff_error"],
        "HMS Cer Detector Eff": row["HMS_Cer_Detector_Eff"],
        "HMS Cer Detector Eff Error": row["HMS_Cer_Detector_Eff_error"],
        "HMS Cal Detector Eff": row["HMS_Cal_Detector_Eff"],
        "HMS Cal Detector Eff Error": row["HMS_Cal_Detector_Eff_error"],
        "HMS Hodo 3/4 Eff": row["HMS_Hodo_3_4_Eff"],
        "HMS Hodo 3/4 Eff Error": row["HMS_Hodo_3_4_Eff_error"],
        "SHMS Hodo 3/4 Eff": row["SHMS_Hodo_3_4_Eff"],
        "SHMS Hodo 3/4 Eff Error": row["SHMS_Hodo_3_4_Eff_error"],
        "EDTM Live Time": row["EDTM_Live_Time"],
        "EDTM Live Time Error": row["EDTM_Live_Time_error"],
        "Boiling Factor": row["Boiling_factor"],
        "Boiling Factor Error": row["Boiling_factor_error"],
        "Dummy Target Correction": "NA",
        "Dummy Target Correction Error": "NA",
        "Effective Charge": row["effective_charge"],
        "Effective Charge Error": row["effective_charge_error"],
        "Total Effective Charge": total_data_effective_charge,
        "Total Effective Charge Error": total_data_effective_charge_error,
        "Normfactor": normfac_data
    })

# Add dummy rows
for row in row_dummy:
    combined_rows.append({
        "Physics_Setting": PHY_SETTING,
        "Run Type": "dummy",
        "Run Number": row["Run_Number"],
        "Charge": row["charge"],
        "Charge Error": row["charge_error"],
        "HMS Tracking Eff": row["HMS_Tracking_Eff"],
        "HMS Tracking Eff Error": row["HMS_Tracking_Eff_error"],
        "SHMS Tracking Eff": row["SHMS_Tracking_Eff"],
        "SHMS Tracking Eff Error": row["SHMS_Tracking_Eff_error"],
        "HMS Cer Detector Eff": row["HMS_Cer_Detector_Eff"],
        "HMS Cer Detector Eff Error": row["HMS_Cer_Detector_Eff_error"],
        "HMS Cal Detector Eff": row["HMS_Cal_Detector_Eff"],
        "HMS Cal Detector Eff Error": row["HMS_Cal_Detector_Eff_error"],
        "HMS Hodo 3/4 Eff": row["HMS_Hodo_3_4_Eff"],
        "HMS Hodo 3/4 Eff Error": row["HMS_Hodo_3_4_Eff_error"],
        "SHMS Hodo 3/4 Eff": row["SHMS_Hodo_3_4_Eff"],
        "SHMS Hodo 3/4 Eff Error": row["SHMS_Hodo_3_4_Eff_error"],
        "EDTM Live Time": row["EDTM_Live_Time"],
        "EDTM Live Time Error": row["EDTM_Live_Time_error"],
        "Boiling Factor": row["Boiling_factor"],
        "Boiling Factor Error": row["Boiling_factor_error"],
        "Dummy Target Correction": dummy_target_corr,
        "Dummy Target Correction Error": dummy_target_corr_error,
        "Effective Charge": row["effective_charge"],
        "Effective Charge Error": row["effective_charge_error"],
        "Total Effective Charge": total_dummy_effective_charge,
        "Total Effective Charge Error": total_dummy_effective_charge_error,
        "Normfactor": normfac_dummy
    })

# Read the existing CSV file into memory
existing_rows = []
try:
    with open(eff_cal_output_csv_path, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        existing_rows = list(csv_reader)
except FileNotFoundError:
    # If the file doesn't exist, initialize an empty list
    existing_rows = []

# Update or add rows
updated_rows = []
for new_row in combined_rows:
    row_found = False
    for existing_row in existing_rows:
        if (existing_row.get("Physics_Setting") == new_row.get("Physics Setting") and
            existing_row.get("Run Type") == new_row.get("Run Type") and
            existing_row.get("Run Number") == str(new_row.get("Run Number"))):
            # Update existing row
            existing_row.update(new_row)
            row_found = True
            break
    if not row_found:
        # Add new row
        updated_rows.append(new_row)

# Combine updated rows with existing rows
existing_rows.extend(updated_rows)

# Write everything back to the CSV file
with open(eff_cal_output_csv_path, mode='w', newline='') as csv_file:
    writer = csv.DictWriter(csv_file, fieldnames=header)
    writer.writeheader()
    writer.writerows(existing_rows)

print(f"Run factors written to {eff_cal_output_csv_path}")

###################################################################################################################################################

nbins = 200

hist_config = {
    "MMpi":    ("Missing Mass Distribution (Cut_All); MM_{\pi}; Counts", nbins, 0.8, 1.2),
    "Q2":      ("Q^{2} Distribution (Cut_All); Q^{2}; Counts", nbins, 2.0, 5.0),
    "W":       ("W Distribution (Cut_All); W; Counts", nbins, 2.0, 3.0),
    "epsilon": ("Epsilon Distribution (Cut_All); #epsilon; Counts", nbins, 0.0, 1.0),
    "theta":   ("Theta Distribution (Cut_All); #theta; Counts", nbins, 0.0, 0.5),
}

def initialize_histograms(t_min_values, t_max_values, phi_bins, hist_config, category_suffixes):
    histograms_by_tphi_cut = []
    added_t_bin_vars = set()  # Track which (tmin, tmax, suffix, var) already had W/Q2 added

    for i in range(len(t_min_values)):
        tmin = t_min_values[i]
        tmax = t_max_values[i]
        t_label = f"t({tmin:.2f}-{tmax:.2f})".replace('.', 'p')
        for j in range(len(phi_bins) - 1):
            phimin = phi_bins[j]
            phimax = phi_bins[j + 1]
            phi_label = f"phi({phimin}-{phimax})"
            hist_set = {suffix: {} for suffix in category_suffixes}

            for var, (base_title, nbins, xmin, xmax) in hist_config.items():
                for suffix in category_suffixes:
                    if var == "MMpi":
                        # MMpi → t and φ bins
                        hname = f"{var}_{t_label}_{phi_label}_{suffix}_data_cut_all"
                        title = f"{base_title} (t ∈ [{tmin:.2f}, {tmax:.2f}], φ ∈ [{phimin}, {phimax}], Category: {suffix})"
                        hist_set[suffix][var] = ROOT.TH1D(hname, title, nbins, xmin, xmax)
                    elif var in ("Q2", "W", "epsilon", "theta") and (tmin, tmax, suffix, var) not in added_t_bin_vars:
                        # Q2, W → only once per t-bin across full phi
                        hname = f"{var}_{t_label}_{suffix}_data_cut_all"
                        title = f"{base_title} (t ∈ [{tmin:.2f}, {tmax:.2f}], φ ∈ [0, 360], Category: {suffix})"
                        hist_set[suffix][var] = ROOT.TH1D(hname, title, nbins, xmin, xmax)
                        added_t_bin_vars.add((tmin, tmax, suffix, var))

            histograms_by_tphi_cut.append((tmin, tmax, phimin, phimax, hist_set))
    return histograms_by_tphi_cut

categories_data = ["prompt_data", "random_data", "randsub_data", "scaled_randsub_data", "prompt_dummy", "random_dummy", "randsub_dummy", "scaled_randsub_dummy", "dummysub"]

# Initialize histograms for data
data_histograms_by_tphi_cut = initialize_histograms(t_min_values, t_max_values, phi_bins, hist_config, categories_data)

##########################################################################################################################################################################################################

# Read stuff from the main event tree
infile_DATA = ROOT.TFile.Open(rootFile_DATA, "READ")
infile_DUMMY = ROOT.TFile.Open(rootFile_DUMMY, "READ")

#Uncut_Pion_Events_Data_tree = infile_DATA.Get("Uncut_Pion_Events")
#Cut_Pion_Events_Accpt_Data_tree = infile_DATA.Get("Cut_Pion_Events_Accpt")
#Cut_Pion_Events_All_Data_tree = infile_DATA.Get("Cut_Pion_Events_All")
Cut_Pion_Events_Prompt_Data_tree = infile_DATA.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Random_Data_tree = infile_DATA.Get("Cut_Pion_Events_Random")
#nEntries_TBRANCH_DATA  = Cut_Pion_Events_Prompt_Data_tree.GetEntries()

#Uncut_Pion_Events_Dummy_tree = infile_DUMMY.Get("Uncut_Pion_Events")
#Cut_Pion_Events_Accpt_Dummy_tree = infile_DUMMY.Get("Cut_Pion_Events_Accpt")
#Cut_Pion_Events_All_Dummy_tree = infile_DUMMY.Get("Cut_Pion_Events_All")
Cut_Pion_Events_Prompt_Dummy_tree = infile_DUMMY.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Random_Dummy_tree = infile_DUMMY.Get("Cut_Pion_Events_Random")
#nEntries_TBRANCH_DUMMY  = Cut_Pion_Events_Prompt_Dummy_tree.GetEntries()

###################################################################################################################################################

# Fill the histograms with the data
for event in Cut_Pion_Events_Prompt_Data_tree:
    if DATA_MMpi_Cut(event) and Diamond_Cut(event):
        for tbin_cut, (tmin, tmax) in zip(tbin_cuts, zip(t_min_values, t_max_values)):
            if tbin_cut(event):
                # Fill Q² and W (once per t-bin)
                for tmin_hist, tmax_hist, _, _, hist_set in data_histograms_by_tphi_cut:
                    if tmin == tmin_hist and tmax == tmax_hist:
                        if "Q2" in hist_set["prompt_data"]:
                            hist_set["prompt_data"]["Q2"].Fill(event.Q2)
                        if "W" in hist_set["prompt_data"]:
                            hist_set["prompt_data"]["W"].Fill(event.W)
                        if "epsilon" in hist_set["prompt_data"]:
                            hist_set["prompt_data"]["epsilon"].Fill(event.epsilon)
                        if "theta" in hist_set["prompt_data"]:
                            hist_set["prompt_data"]["theta"].Fill(event.th_q) 
                        break  # Only fill Q2/W once for this t-bin
                # Fill MMπ per t–φ bin
                phi_deg = event.ph_q * (180 / math.pi)  # Convert φ to degrees
                if phi_deg < 0:
                    phi_deg += 360  # Ensure φ is in [0, 360]
                for tmin_hist, tmax_hist, phimin_hist, phimax_hist, hist_set in data_histograms_by_tphi_cut:
                    if tmin == tmin_hist and tmax == tmax_hist and phimin_hist <= phi_deg <= phimax_hist:
                        hist_set["prompt_data"]["MMpi"].Fill(event.MMpi + MM_Offset)
                        break  # Only fill in the matching t–φ bin

for event in Cut_Pion_Events_Random_Data_tree:
    if DATA_MMpi_Cut(event) and Diamond_Cut(event):
        for tbin_cut, (tmin, tmax) in zip(tbin_cuts, zip(t_min_values, t_max_values)):
            if tbin_cut(event):
                # Fill Q2 and W once per t-bin
                for tmin_hist, tmax_hist, _, _, hist_set in data_histograms_by_tphi_cut:
                    if tmin == tmin_hist and tmax == tmax_hist:
                        if "Q2" in hist_set["random_data"]:
                            hist_set["random_data"]["Q2"].Fill(event.Q2)
                        if "W" in hist_set["random_data"]:
                            hist_set["random_data"]["W"].Fill(event.W)
                        if "epsilon" in hist_set["random_data"]:
                            hist_set["random_data"]["epsilon"].Fill(event.epsilon)
                        if "theta" in hist_set["random_data"]:
                            hist_set["random_data"]["theta"].Fill(event.th_q)
                        break
                # Fill MMpi for t–φ bin
                phi_deg = event.ph_q * (180 / math.pi)
                if phi_deg < 0:
                    phi_deg += 360  # Ensure φ is in [0, 360]
                for tmin_hist, tmax_hist, phimin_hist, phimax_hist, hist_set in data_histograms_by_tphi_cut:
                    if tmin == tmin_hist and tmax == tmax_hist and phimin_hist <= phi_deg <= phimax_hist:
                        hist_set["random_data"]["MMpi"].Fill(event.MMpi + MM_Offset)
                        break

# --- Prompt Dummy Tree ---
for event in Cut_Pion_Events_Prompt_Dummy_tree:
    if DATA_MMpi_Cut(event) and Diamond_Cut(event):
        for tbin_cut, (tmin, tmax) in zip(tbin_cuts, zip(t_min_values, t_max_values)):
            if tbin_cut(event):
                # Fill Q2 and W once per t-bin
                for tmin_hist, tmax_hist, _, _, hist_set in data_histograms_by_tphi_cut:
                    if tmin == tmin_hist and tmax == tmax_hist:
                        if "Q2" in hist_set["prompt_dummy"]:
                            hist_set["prompt_dummy"]["Q2"].Fill(event.Q2)
                        if "W" in hist_set["prompt_dummy"]:
                            hist_set["prompt_dummy"]["W"].Fill(event.W)
                        if "epsilon" in hist_set["prompt_dummy"]:
                            hist_set["prompt_dummy"]["epsilon"].Fill(event.epsilon)
                        if "theta" in hist_set["prompt_dummy"]:
                            hist_set["prompt_dummy"]["theta"].Fill(event.th_q)
                        break
                # Fill MMpi for t–φ bin
                phi_deg = event.ph_q * (180 / math.pi)
                if phi_deg < 0:
                    phi_deg += 360  # Ensure φ is in [0, 360]
                for tmin_hist, tmax_hist, phimin_hist, phimax_hist, hist_set in data_histograms_by_tphi_cut:
                    if tmin == tmin_hist and tmax == tmax_hist and phimin_hist <= phi_deg <= phimax_hist:
                        hist_set["prompt_dummy"]["MMpi"].Fill(event.MMpi + MM_Offset)
                        break

# --- Random Dummy Tree ---
for event in Cut_Pion_Events_Random_Dummy_tree:
    if DATA_MMpi_Cut(event) and Diamond_Cut(event):
        for tbin_cut, (tmin, tmax) in zip(tbin_cuts, zip(t_min_values, t_max_values)):
            if tbin_cut(event):
                # Fill Q2 and W once per t-bin
                for tmin_hist, tmax_hist, _, _, hist_set in data_histograms_by_tphi_cut:
                    if tmin == tmin_hist and tmax == tmax_hist:
                        if "Q2" in hist_set["random_dummy"]:
                            hist_set["random_dummy"]["Q2"].Fill(event.Q2)
                        if "W" in hist_set["random_dummy"]:
                            hist_set["random_dummy"]["W"].Fill(event.W)
                        if "epsilon" in hist_set["random_dummy"]:
                            hist_set["random_dummy"]["epsilon"].Fill(event.epsilon)
                        if "theta" in hist_set["random_dummy"]:
                            hist_set["random_dummy"]["theta"].Fill(event.th_q)
                        break
                # Fill MMpi for t–φ bin
                phi_deg = event.ph_q * (180 / math.pi)
                if phi_deg < 0:
                    phi_deg += 360  # Ensure φ is in [0, 360]
                for tmin_hist, tmax_hist, phimin_hist, phimax_hist, hist_set in data_histograms_by_tphi_cut:
                    if tmin == tmin_hist and tmax == tmax_hist and phimin_hist <= phi_deg <= phimax_hist:
                        hist_set["random_dummy"]["MMpi"].Fill(event.MMpi + MM_Offset)
                        break

print("####################################")
print("###### Histogram filling done ######")
print("####################################\n")

#################################################################################################################################################

# Scale random histograms according to nWindows
for tmin_hist, tmax_hist, phimin_hist, phimax_hist, hist_set in data_histograms_by_tphi_cut:
    for var in hist_config:
        # Scale random_data
        if "random_data" in hist_set and var in hist_set["random_data"]:
            random_data_hist = hist_set["random_data"][var]
            if random_data_hist:
                random_data_hist.Scale(1 / nWindows)
        # Scale random_dummy
        if "random_dummy" in hist_set and var in hist_set["random_dummy"]:
            random_dummy_hist = hist_set["random_dummy"][var]
            if random_dummy_hist:
                random_dummy_hist.Scale(1 / nWindows)

# Perform random subtraction for each t and phi bin for real data
for tmin_hist, tmax_hist, phimin_hist, phimax_hist, hist_set in data_histograms_by_tphi_cut:
    for var in hist_config:
        if var in hist_set["prompt_data"] and var in hist_set["random_data"]:
            prompt_data_hist = hist_set["prompt_data"][var]
            random_data_hist = hist_set["random_data"][var]
            # Clone and subtract
            randsub_data_hist = prompt_data_hist.Clone()
            randsub_data_hist.Add(random_data_hist, -1)
            # Store result
            hist_set["randsub_data"][var] = randsub_data_hist

# Perform random subtraction for dummy data
for tmin_hist, tmax_hist, phimin_hist, phimax_hist, hist_set in data_histograms_by_tphi_cut:
    for var in hist_config:
        if var in hist_set["prompt_dummy"] and var in hist_set["random_dummy"]:
            prompt_dummy_hist = hist_set["prompt_dummy"][var]
            random_dummy_hist = hist_set["random_dummy"][var]
            # Clone and subtract
            randsub_dummy_hist = prompt_dummy_hist.Clone()
            randsub_dummy_hist.Add(random_dummy_hist, -1)
            # Store result
            hist_set["randsub_dummy"][var] = randsub_dummy_hist

print("######################################")
print("###### Random substraction done ######")
print("######################################\n")

############################################################################################################################################

# Scale randsub_data and randsub_dummy histograms with normalization factors
for tmin_hist, tmax_hist, phimin_hist, phimax_hist, hist_set in data_histograms_by_tphi_cut:
    for var in hist_config:
        # Scale the randsub_data histogram
        if "randsub_data" in hist_set and var in hist_set["randsub_data"]:
            randsub_data_hist = hist_set["randsub_data"][var]
            if randsub_data_hist:
                scaled_randsub_data_hist = randsub_data_hist.Clone()
                scaled_randsub_data_hist.Scale(normfac_data)
                hist_set["scaled_randsub_data"][var] = scaled_randsub_data_hist
        # Scale the randsub_dummy histogram
        if "randsub_dummy" in hist_set and var in hist_set["randsub_dummy"]:
            randsub_dummy_hist = hist_set["randsub_dummy"][var]
            if randsub_dummy_hist:
                scaled_randsub_dummy_hist = randsub_dummy_hist.Clone()
                scaled_randsub_dummy_hist.Scale(normfac_dummy)
                hist_set["scaled_randsub_dummy"][var] = scaled_randsub_dummy_hist

# Dummy Subtraction
for tmin_hist, tmax_hist, phimin_hist, phimax_hist, hist_set in data_histograms_by_tphi_cut:
    for var in hist_config:
        if "scaled_randsub_data" in hist_set and "scaled_randsub_dummy" in hist_set:
            if var in hist_set["scaled_randsub_data"] and var in hist_set["scaled_randsub_dummy"]:
                scaled_randsub_data_hist = hist_set["scaled_randsub_data"][var]
                scaled_randsub_dummy_hist = hist_set["scaled_randsub_dummy"][var]
                if scaled_randsub_data_hist and scaled_randsub_dummy_hist:
                    dummysub_data_hist = scaled_randsub_data_hist.Clone()
                    dummysub_data_hist.Add(scaled_randsub_dummy_hist, -1)
                    hist_set["dummysub"][var] = dummysub_data_hist

                    # Print confirmation of subtraction
                    # Print only t-bin for W and Q²
#                    if var in ["W", "Q2"]:
#                        print(f"[SUBTRACTION DONE] Variable: {var} | t ∈ [{tmin_hist:.3f}, {tmax_hist:.3f}]")
#                    else:  # Print both t-bin and φ-bin for MMpi
#                        print(f"[SUBTRACTION DONE] Variable: {var} | t ∈ [{tmin_hist:.3f}, {tmax_hist:.3f}] | φ ∈ [{phimin_hist}, {phimax_hist}]")

print("####################################")
print("###### Dummy subtraction done ######")
print("####################################\n")

#############################################################################################################################################

print(f"Average Q2 and W values for each t-bin (with uncertainties)")
Beam_Energy_S = filtered_data_df['Beam_Energy'].iloc[0]

print("Mean and Mean Error for W and Q² (per t-bin):")
avg_kin_tbin = {}
for tmin, tmax in zip(t_min_values, t_max_values):
    for tmin_hist, tmax_hist, phimin_hist, _, hist_set in data_histograms_by_tphi_cut:
        if tmin == tmin_hist and tmax == tmax_hist and phimin_hist == 0:
            hist_Q = hist_set["dummysub"]["Q2"]
            hist_W = hist_set["dummysub"]["W"]
            hist_eps = hist_set["dummysub"]["epsilon"]
            hist_theta = hist_set["dummysub"]["theta"]

            avg_Q2 = hist_Q.GetMean()
            err_Q2 = hist_Q.GetMeanError()
            avg_W = hist_W.GetMean()
            err_W = hist_W.GetMeanError()
            avg_eps = hist_eps.GetMean()
            err_eps = hist_eps.GetMeanError()
            avg_theta = hist_theta.GetMean()
            err_theta = hist_theta.GetMeanError()

            tbin_key = f"t({tmin:.2f}-{tmax:.2f})"
            avg_kin_tbin[tbin_key] = {
                "avg_Q2": avg_Q2, "err_Q2": err_Q2,
                "avg_W": avg_W, "err_W": err_W,
                "avg_eps": avg_eps, "err_eps": err_eps,
                "avg_theta": avg_theta, "err_theta": err_theta
            }

            print(f"{tbin_key}: <Q²> = {avg_Q2:.4f} ± {err_Q2:.4f}, <W> = {avg_W:.4f} ± {err_W:.4f}, <eps> = {avg_eps:.4f} ± {err_eps:.4f}, <θ> = {avg_theta:.4f} ± {err_theta:.4f}")
            break

# Define output CSV file path for average kinematics
avg_kinematics_path = "%s/LTSep_CSVs/physics_yields_csv/%s/%s_Physics_Avg_Data_Kinematics.csv" % (UTILPATH, physet_dir_name, setting_name)

# Define the header
avg_kin_header = [
    "Physics_Setting", "Beam_Energy", "total_tbins", "tbin_number", "t_min", "t_max", "t_central",
    "total_phibins", "phi_central", "avg_Q2", "avg_Q2_err", "avg_W", "avg_W_err",
    "avg_eps", "avg_eps_err", "avg_theta", "avg_theta_err"
]
# Read the existing CSV file
rows = []
try:
    with open(avg_kinematics_path, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        rows = list(csv_reader)
except FileNotFoundError:
    rows = []

# Prepare new rows for this Physics_Setting
new_setting_rows = []
for i, (tmin, tmax) in enumerate(zip(t_min_values, t_max_values), start=1):
    tbin_key = f"t({tmin:.2f}-{tmax:.2f})"
    if tbin_key in avg_kin_tbin:
        avg_kin = avg_kin_tbin[tbin_key]
        new_setting_rows.append({
            "Physics_Setting": PHY_SETTING,
            "Beam_Energy": Beam_Energy_S,
            "total_tbins": total_tbins,
            "tbin_number": i,
            "t_min": f"{tmin:.4f}",
            "t_max": f"{tmax:.4f}",
            "t_central": f"{(tmin + tmax) / 2:.4f}",
            "total_phibins": total_phibins,
            "phi_central": f"{(0 + 360) / 2:.4f}",
            "avg_Q2": f"{avg_kin['avg_Q2']:.6f}",
            "avg_Q2_err": f"{avg_kin['err_Q2']:.6f}",
            "avg_W": f"{avg_kin['avg_W']:.6f}",
            "avg_W_err": f"{avg_kin['err_W']:.6f}",
            "avg_eps": f"{avg_kin['avg_eps']:.6f}",
            "avg_eps_err": f"{avg_kin['err_eps']:.6f}",
            "avg_theta": f"{avg_kin['avg_theta']:.6f}",
            "avg_theta_err": f"{avg_kin['err_theta']:.6f}"
        })

setting_exists = any(
    row["Physics_Setting"] == PHY_SETTING and
    float(row["t_min"]) == float(tmin) and
    float(row["t_max"]) == float(tmax)
    for row in rows
)

# Replace or append
if setting_exists:
    rows = [row for row in rows if row["Physics_Setting"] != PHY_SETTING]
    rows.extend(new_setting_rows)
else:
    rows.extend(new_setting_rows)

# Write everything back to CSV
with open(avg_kinematics_path, mode='w', newline='') as csv_file:
    writer = csv.DictWriter(csv_file, fieldnames=avg_kin_header)
    writer.writeheader()
    writer.writerows(rows)

print(f"Average Q² and W values updated for '{PHY_SETTING}' in: {avg_kinematics_path}")

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Physics Yield Calculations
# Dictionary to store integrals and errors for each t-bin and phi-bin
dN_data_MMpi = {}

# Loop over all t-bins and phi-bins
for tmin, tmax, phimin, phimax, hist_set in data_histograms_by_tphi_cut:
    # Create a unique key for the current t-bin and phi-bin
    bin_key = f"t({tmin:.2f}-{tmax:.2f})_phi({phimin}-{phimax})"

    # Initialize the error array
    dN_data = array.array('d', [0.0])

    # Get the dummysub histogram for the current bin
    hist = hist_set["dummysub"]["MMpi"]  # Replace "MMpi" with the variable of interest

    # Calculate the integral and error
    N_data = hist.IntegralAndError(1, nbins, dN_data, "")

    # Store the results in the dictionary
    dN_data_MMpi[bin_key] = {
        "integral": N_data,
        "error": dN_data[0]
    }

    # Print the result for debugging
#    print(f"Physics Yield: {bin_key} = Yield: {N_data:.8f}, Error: {dN_data[0]:.8f}")
#print("-" * 40)

# Dictionary to store integrals and errors for each t-bin and phi-bin
dN_error = {}

# Loop over all t-bins and phi-bins
for tmin, tmax, phimin, phimax, hist_set in data_histograms_by_tphi_cut:
    # Create a unique key for the current t-bin and phi-bin
    bin_key = f"t({tmin:.2f}-{tmax:.2f})_phi({phimin}-{phimax})"

    # Initialize error arrays
    dN_data_error = array.array('d', [0.0])
    dN_dummy_error = array.array('d', [0.0])

    # Get the randsub_data and randsub_dummy histograms for the current bin
    randsub_data_hist = hist_set["randsub_data"]["MMpi"]  # Replace "MMpi" with the variable of interest
    randsub_dummy_hist = hist_set["randsub_dummy"]["MMpi"]  # Replace "MMpi" with the variable of interest

    # Calculate the integrals and errors
    N_data = randsub_data_hist.IntegralAndError(1, nbins, dN_data_error, "")
    N_dummy = randsub_dummy_hist.IntegralAndError(1, nbins, dN_dummy_error, "")
    N_dummysub = N_data - N_dummy

    # Normalize the integrals
    N_data_norm = N_data / total_data_effective_charge
    N_dummy_norm = N_dummy / total_dummy_effective_charge

    # Calculate the normalized errors
    if N_data != 0:
        dN_data_norm = N_data_norm * math.sqrt((dN_data_error[0] / N_data) ** 2 + (total_data_effective_charge_error / total_data_effective_charge) ** 2)
    else:
        dN_data_norm = N_data_norm * math.sqrt((0) ** 2 + (total_data_effective_charge_error / total_data_effective_charge) ** 2)
    if N_dummy != 0:
        dN_dummy_norm = N_dummy_norm * math.sqrt((dN_dummy_error[0] / N_dummy) ** 2 + (total_dummy_effective_charge_error / total_dummy_effective_charge) ** 2)
    else:
        dN_dummy_norm = N_data_norm * math.sqrt((0) ** 2 + (total_dummy_effective_charge_error / total_dummy_effective_charge) ** 2)

    # Ensure dN_data_norm and dN_dummy_norm are not zero
#    if dN_data_norm == 0:
#        dN_data_norm = 1000.0
#    if dN_dummy_norm == 0:
#        dN_dummy_norm = 1000.0

    # Calculate the final yield and its error
    N_data_dummy_sub_norm = N_data_norm - N_dummy_norm
    dN_data_dummy_sub_norm = math.sqrt(dN_data_norm ** 2 + dN_dummy_norm ** 2)

    # Store the results in the dictionary
    dN_error[bin_key] = {
        "Data Counts": N_data,
        "Dummy Counts": N_dummy,
        "N_dummysub": N_dummysub,
        "N_data_norm": N_data_norm,
        "N_dummy_norm": N_dummy_norm,
        "dN_data_norm": dN_data_norm,
        "dN_dummy_norm": dN_dummy_norm,
        "N_data_dummy_sub_norm": N_data_dummy_sub_norm,
        "dN_data_dummy_sub_norm": dN_data_dummy_sub_norm
    }

    # Print the result for debugging
#    print(f"Physics Yield: {bin_key} = Yield: {N_data_dummy_sub_norm:.8f}, Error: {dN_data_dummy_sub_norm:.8f}")
#print("-" * 40)

# Print the results from both loops
print("=" * 40)
print("Physics Yield Results:")
for bin_key in dN_data_MMpi:
    # Use the correct key names from dN_data_MMP
    Y_data = abs(dN_data_MMpi[bin_key]["integral"])
    dY_data_dummy_sub_norm = dN_error[bin_key]["dN_data_dummy_sub_norm"]
    print(f"Physics Yield: {bin_key} = Yield: {Y_data:.8f}, Error: {dY_data_dummy_sub_norm:.8f}")
print("=" * 40)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Define the output CSV file name
yields_pions_path =  "%s/LTSep_CSVs/physics_yields_csv/%s/%s_Physics_Data_Yield.csv" % (UTILPATH, physet_dir_name, setting_name)

# Read the existing CSV file into memory
rows = []
header = ["Physics_Setting", "tbin_number", "t_min", "t_max", "phibin_number", "phi_min", "phi_max", "DummySub_Counts", "physics_yield", "physics_yield_error", "%error/yield"]
try:
    with open(yields_pions_path, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        rows = list(csv_reader)
except FileNotFoundError:
    # If the file doesn't exist, initialize an empty list
    rows = []

# Gather all unique t-bin and phi-bin edges
tbin_edges = sorted(list({(tmin, tmax) for tmin, tmax, _, _, _ in data_histograms_by_tphi_cut}))

# Update or add rows
updated_rows = []
for tbin_number, (tmin, tmax) in enumerate(tbin_edges, start=1):
    # Filter phi bins for this t range, and sort
    phibins = sorted([(phimin, phimax) for tmin2, tmax2, phimin, phimax, _ in data_histograms_by_tphi_cut if tmin2 == tmin and tmax2 == tmax])
    for phibin_number, (phimin, phimax) in enumerate(phibins, start=1):
        bin_key = f"t({tmin:.2f}-{tmax:.2f})_phi({phimin}-{phimax})"
        N_dummysub = dN_error[bin_key]["N_dummysub"]
        Y_data = dN_data_MMpi[bin_key]["integral"]
        dY_data_dummy_sub_norm = dN_error[bin_key]["dN_data_dummy_sub_norm"]

        # See if the row already exists and update it, else append new
        row_found = False
        for row in rows:
            if (
                row["Physics_Setting"] == PHY_SETTING
                and float(row["t_min"]) == tmin
                and float(row["t_max"]) == tmax
                and float(row["phi_min"]) == phimin
                and float(row["phi_max"]) == phimax
            ):
                row.update({
                    "Physics_Setting": PHY_SETTING,
                    "tbin_number": tbin_number,
                    "t_min": tmin,
                    "t_max": tmax,
                    "phibin_number": phibin_number,
                    "phi_min": phimin,
                    "phi_max": phimax,
                    "DummySub_Counts": N_dummysub,
                    "physics_yield": abs(Y_data),
                    "physics_yield_error": dY_data_dummy_sub_norm,
                    "%error/yield": (dY_data_dummy_sub_norm/Y_data)*100 if Y_data != 0 else dY_data_dummy_sub_norm
                })
                row_found = True
                break

        if not row_found:
            updated_rows.append({
                "Physics_Setting": PHY_SETTING,
                "tbin_number": tbin_number,
                "t_min": tmin,
                "t_max": tmax,
                "phibin_number": phibin_number,
                "phi_min": phimin,
                "phi_max": phimax,
                "DummySub_Counts": N_dummysub,
                "physics_yield": abs(Y_data),
                "physics_yield_error": dY_data_dummy_sub_norm,
                "%error/yield": (dY_data_dummy_sub_norm/Y_data)*100 if Y_data != 0 else dY_data_dummy_sub_norm
            })

# Combine updated rows with existing rows
rows.extend(updated_rows)

# Write everything back to the CSV
with open(yields_pions_path, mode='w', newline='') as csv_file:
    writer = csv.DictWriter(csv_file, fieldnames=header)
    writer.writeheader()
    writer.writerows(rows)

print(f"Physics yield results written to {yields_pions_path}")

#############################################################################################################################################

# Making directories in output file
outHistFile = ROOT.TFile.Open("%s/%s_%s_Yield_Data.root" % (OUTPATH, PHY_SETTING, MaxEvent) , "RECREATE")


# Create top-level directories for Q2 and W
q2_dir = outHistFile.mkdir("Cut_Q2_Pion_Events_Norm_DummySub_Data")
w_dir = outHistFile.mkdir("Cut_W_Pion_Events_Norm_DummySub_Data")
epsilon_dir = outHistFile.mkdir("Cut_Epsilon_Pion_Events_Norm_DummySub_Data")
theta_dir = outHistFile.mkdir("Cut_Theta_Pion_Events_Norm_DummySub_Data")

# Write histograms to t-bin-specific directories
for tmin, tmax, phimin, phimax, hist_set in data_histograms_by_tphi_cut:

    # --- Write MMpi (binned in t and phi) ---
    t_dirname_MM = f"Cut_MM_Pion_Events_Norm_DummySub_t{tmin:.2f}-t{tmax:.2f}_Data".replace(".", "p")
    tdir_MM = outHistFile.GetDirectory(t_dirname_MM)
    if not tdir_MM:
        tdir_MM = outHistFile.mkdir(t_dirname_MM)
    tdir_MM.cd()
    var_MM = "MMpi"
    for cat in ["dummysub"]:
        hist = hist_set[cat].get(var_MM)
        if hist:
            hist.SetTitle(f"{var_MM}_t{tmin:.2f}-t{tmax:.2f}_phi{phimin}-phi{phimax}_{cat}".replace(".", "p"))
            hist.SetName(f"{var_MM}_t{tmin:.2f}-t{tmax:.2f}_phi{phimin}-phi{phimax}_{cat}".replace(".", "p"))
            hist.GetYaxis().SetTitle("Counts")
            hist.Write()

    # --- Write Q2 and W once per t-bin (no phi binning) ---
    for var, target_dir in zip(["Q2", "W", "epsilon", "theta"], [q2_dir, w_dir, epsilon_dir, theta_dir]):
        for cat in ["dummysub"]:
            hist = hist_set[cat].get(var)
            if hist:
                target_dir.cd()
                hist.SetTitle(f"{var}_t{tmin:.2f}-t{tmax:.2f}_{cat}".replace(".", "p"))
                hist.SetName(f"{var}_t{tmin:.2f}-t{tmax:.2f}_{cat}".replace(".", "p"))
                hist.GetYaxis().SetTitle("Counts")
                hist.Write()

print("####################################")
print("###### Histogram writing done ######")
print("####################################\n")

infile_DATA.Close() 
infile_DUMMY.Close()
outHistFile.Close()

print ("Processing Complete")