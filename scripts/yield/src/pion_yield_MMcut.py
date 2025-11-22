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

##################################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=8:
    print("!!!!! ERROR !!!!!\n Expected 8 arguments\n Usage is with - PHY_SETTING Beam MaxEvents Suffix RunList CVSFile\n!!!!! ERROR !!!!!")
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
EFF_CSV     = "/u/group/c-pionlt/USERS/%s/hallc_replay_lt/UTIL_PION/efficiencies" % (USER)

# Define paths SIMC
# Extract the first three words from PHY_SETTING for the CSV file name
setting_name = "_".join(PHY_SETTING.split("_")[:3])
physet_dir_name = "%s_std" % (setting_name)
SIMCPATH = "/volatile/hallc/c-pionlt/%s/OUTPUT/Analysis/SIMC/%s/" % (USER, physet_dir_name)

#################################################################################################################################################

# Output PDF File Name
print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))
Pion_Analysis_Distributions = "%s/%s_%s_ProdCoin_Pion_Analysis_MMcut_Distributions.pdf" % (OUTPATH, PHY_SETTING, MaxEvent)

# Input file location and variables taking
rootFile_DATA = "%s/%s_%s_%s.root" % (OUTPATH, PHY_SETTING, MaxEvent, DATA_Suffix)
rootFile_DUMMY = "%s/%s_%s_%s.root" % (OUTPATH, PHY_SETTING, MaxEvent, DUMMY_Suffix)
rootFile_SIMC = "%s/%s.root" % (SIMCPATH, SIMC_Suffix)
data_run_list = "%s/%s" % (RUNLISTPATH, DATA_RUN_LIST)
dummy_run_list = "%s/%s" % (RUNLISTPATH, DUMMY_RUN_LIST)
csv_file = "%s/%s.csv" % (EFF_CSV, CSV_FILE)

###################################################################################################################################################

# SIMC Cuts for Pions Selection
HMS_Acceptance = lambda event: (event.hsdelta >= -8.0) & (event.hsdelta <= 8.0) & (event.hsxpfp >= -0.08) & (event.hsxpfp <= 0.08) & (event.hsypfp >= -0.045) & (event.hsypfp <= 0.045)
SHMS_Acceptance = lambda event: (event.ssdelta >= -10.0) & (event.ssdelta <= 20.0) & (event.ssxpfp >= -0.06) & (event.ssxpfp <= 0.06) & (event.ssypfp >= -0.04) & (event.ssypfp <= 0.04)
#SHMS_Aero_Cut = lambda event: (event.paero_x_det > -55.0) & (event.paero_x_det < 55.0) & (event.paero_y_det > -50) & (event.paero_y_det < 50) # Aerogel tray n = 1.030
SHMS_Aero_Cut = lambda event: (event.paero_x_det > -45.0) & (event.paero_x_det < 45.0) & (event.paero_y_det > -30) & (event.paero_y_det < 30) # Aerogel tray n = 1.011

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

# Read CSV File and calculate total charge in mC for normalization
# Input runlists and csv files
data_run_list_file = (data_run_list)
dummy_run_list_file = (dummy_run_list)
csv_file_name = (csv_file)

print ('\nPhysics Setting = ',PHY_SETTING, '\n')
print("-"*40)

# Read run numbers from the run list file
with open(data_run_list_file, 'r') as run_list_file_1:
    data_runs = [line.strip() for line in run_list_file_1 if line.strip()]
with open(dummy_run_list_file, 'r') as run_list_file_2:
    dummy_runs = [line.strip() for line in run_list_file_2 if line.strip()]

# Read CSV file using pandas
df = pd.read_csv(csv_file_name)

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
nbins_p = 1000

# Defining Histograms for Pions
# Histograms having Cuts (Acceptance)
P_kin_MMpi_pions_data_accpt_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_accpt_cut_all", "MIssing Mass data (dummysub_accpt_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)

# Histograms Having Cuts (Acceptance + PID + RF + Prompt Selection) Data
P_kin_MMpi_pions_data_prompt_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_prompt_cut_all", "MIssing Mass data (prompt_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)

# Histograms Having Cuts (Acceptance + PID + RF + Random Selection) Data
P_kin_MMpi_pions_data_random_cut_all = ROOT.TH1D("P_kin_MMpi_pions_data_random_cut_all", "MIssing Mass data (random_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)

# Histograms Having Cuts (Acceptance + PID + RF + Prompt Selection) Dummy
P_kin_MMpi_pions_dummy_prompt_cut_all = ROOT.TH1D("P_kin_MMpi_pions_dummy_prompt_cut_all", "MIssing Mass dummy (prompt_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)

# Histograms Having Cuts (Acceptance + PID + RF + Random Selection) Dummy
P_kin_MMpi_pions_dummy_random_cut_all = ROOT.TH1D("P_kin_MMpi_pions_dummy_random_cut_all", "MIssing Mass dummy (random_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)

# Histograms having Cuts (Acceptance + PID + RF + RandSub) Data
P_kin_MMpi_pions_randsub_data_cut_all = ROOT.TH1D("P_kin_MMpi_pions_randsub_data_cut_all", "MIssing Mass data (randsub_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)

# Histograms having Cuts (Acceptance + PID + RF + RandSub) Dummy
P_kin_MMpi_pions_randsub_dummy_cut_all = ROOT.TH1D("P_kin_MMpi_pions_randsub_dummy_cut_all", "MIssing Mass dummy (randsub_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)

# Histograms having Cuts (Acceptance + PID + RF + RandSub + DummySub)
P_kin_MMpi_pions_dummysub_data_cut_all = ROOT.TH1D("P_kin_MMpi_pions_dummysub_data_cut_all", "MIssing Mass data (dummysub_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)

# Histograms having Cuts (Acceptance + PID + RF + RandSub + Norm) Data
P_kin_MMpi_pions_normrandsub_data_cut_all = ROOT.TH1D("P_kin_MMpi_pions_normrandsub_data_cut_all", "MIssing Mass data (datasub_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)

# Histograms having Cuts (Acceptance + PID + RF + RandSub + Norm) Dummy
P_kin_MMpi_pions_normrandsub_dummy_cut_all = ROOT.TH1D("P_kin_MMpi_pions_normrandsub_dummy_cut_all", "MIssing Mass data (dummysub_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)

# Histograms having Cuts (Acceptance + PID + RF + Prompt Selection) Norm Dummy Subtraction Data
P_kin_MMpi_pions_normdummysub_data_cut_all = ROOT.TH1D("P_kin_MMpi_pions_normdummysub_data_cut_all", "MIssing Mass data (dummysub_cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)

# SIMC Histograms with Cuts
MMpi_pions_simc_cut_all = ROOT.TH1D("MMpi_pions_simc_cut_all", "MIssing Mass SIMC (cut_all); MM_{pi}; Counts", nbins_p, 0, 2.2)

# Histograms for Error Calculations
P_kin_secondary_pmiss_pions_data_prompt_cut_all_error = ROOT.TH1D("P_kin_secondary_pmiss_pions_data_prompt_cut_all_error", "pmiss Distribution; pmiss; Counts", nbins_p, 0.0, 2.0) 
P_kin_secondary_pmiss_pions_dummy_prompt_cut_all_error = ROOT.TH1D("P_kin_secondary_pmiss_pions_dummy_prompt_cut_all_error", "pmiss Distribution; pmiss; Counts", nbins_p, 0.0, 2.0)
P_kin_secondary_pmiss_pions_data_random_cut_all_error = ROOT.TH1D("P_kin_secondary_pmiss_pions_data_random_cut_all_error", "pmiss Distribution; pmiss; Counts", nbins_p, 0.0, 2.0) 
P_kin_secondary_pmiss_pions_dummy_random_cut_all_error = ROOT.TH1D("P_kin_secondary_pmiss_pions_dummy_random_cut_all_error", "pmiss Distribution; pmiss; Counts", nbins_p, 0.0, 2.0)
P_kin_secondary_pmiss_pions_randsub_data_cut_all_error = ROOT.TH1D("P_kin_secondary_pmiss_pions_randsub_data_cut_all_error", "pmiss Distribution; pmiss; Counts", nbins_p, 0.0, 2.0)
P_kin_secondary_pmiss_pions_randsub_dummy_cut_all_error = ROOT.TH1D("P_kin_secondary_pmiss_pions_randsub_dummy_cut_all_error", "pmiss Distribution; pmiss; Counts", nbins_p, 0.0, 2.0)
pmiss_pions_simc_cut_all_error = ROOT.TH1D("pmiss_pions_simc_cut_all_error", "pmiss Distribution; pmiss; Counts", nbins_p, 0.0, 2.0)
pmiss_pions_simc_cut_all_norm_error = ROOT.TH1D("pmiss_pions_simc_cut_all_norm_error", "pmiss Distribution; pmiss; Counts", nbins_p, 0.0, 2.0)

##########################################################################################################################################################################################################

# Read stuff from the main event tree
infile_DATA = ROOT.TFile.Open(rootFile_DATA, "READ")
infile_DUMMY = ROOT.TFile.Open(rootFile_DUMMY, "READ")
infile_SIMC = ROOT.TFile.Open(rootFile_SIMC, "READ")

#Uncut_Pion_Events_Data_tree = infile_DATA.Get("Uncut_Pion_Events")
Cut_Pion_Events_Accpt_Data_tree = infile_DATA.Get("Cut_Pion_Events_Accpt")
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

Uncut_Pion_Events_SIMC_tree = infile_SIMC.Get("h10")
#nEntries_TBRANCH_SIMC  = Uncut_Proton_Events_SIMC_tree.GetEntries()

###################################################################################################################################################

#Fill histograms for Cut All Data
for event in Cut_Pion_Events_Accpt_Data_tree:
        P_kin_MMpi_pions_data_accpt_cut_all.Fill(event.MMpi)
#    ibin += 1

#Fill histograms for Prompt Data
for event in Cut_Pion_Events_Prompt_Data_tree:
        P_kin_MMpi_pions_data_prompt_cut_all.Fill(event.MMpi)
        P_kin_secondary_pmiss_pions_data_prompt_cut_all_error.Fill(event.pmiss)
#    ibin += 1

#Fill histograms for Random Data
for event in Cut_Pion_Events_Random_Data_tree:
        P_kin_MMpi_pions_data_random_cut_all.Fill(event.MMpi)
        P_kin_secondary_pmiss_pions_data_random_cut_all_error.Fill(event.pmiss)
#    ibin += 1

# Fill histograms from Prompt Dummy
for event in Cut_Pion_Events_Prompt_Dummy_tree:
        P_kin_MMpi_pions_dummy_prompt_cut_all.Fill(event.MMpi)
        P_kin_secondary_pmiss_pions_dummy_prompt_cut_all_error.Fill(event.pmiss)
#    ibin += 1

# Fill histograms from Random Dummy
#ibin = 1
for event in Cut_Pion_Events_Random_Dummy_tree:
        P_kin_MMpi_pions_dummy_random_cut_all.Fill(event.MMpi)
        P_kin_secondary_pmiss_pions_dummy_random_cut_all_error.Fill(event.pmiss)
#    ibin += 1

# Fill histograms from SIMC ROOT File
for event in Uncut_Pion_Events_SIMC_tree:
    if HMS_Acceptance(event) & SHMS_Acceptance(event) & SHMS_Aero_Cut(event):
            MMpi_pions_simc_cut_all.Fill(event.missmass, event.Weight)
            pmiss_pions_simc_cut_all_error.Fill(event.Pm, event.Weight)
            pmiss_pions_simc_cut_all_norm_error.Fill(event.Pm, event.Weight)

print("####################################")
print("###### Histogram filling done ######")
print("####################################\n")

#################################################################################################################################################

# Random subtraction from Histograms
P_kin_MMpi_pions_data_random_cut_all.Scale(1.0/nWindows)
P_kin_secondary_pmiss_pions_data_random_cut_all_error.Scale(1.0/nWindows)
P_kin_MMpi_pions_randsub_data_cut_all.Add(P_kin_MMpi_pions_data_prompt_cut_all, P_kin_MMpi_pions_data_random_cut_all, 1, -1)
P_kin_secondary_pmiss_pions_randsub_data_cut_all_error.Add(P_kin_secondary_pmiss_pions_data_prompt_cut_all_error, P_kin_secondary_pmiss_pions_data_random_cut_all_error, 1, -1)
P_kin_MMpi_pions_normrandsub_data_cut_all.Add(P_kin_MMpi_pions_data_prompt_cut_all, P_kin_MMpi_pions_data_random_cut_all, 1, -1)

P_kin_MMpi_pions_dummy_random_cut_all.Scale(1.0/nWindows)
P_kin_secondary_pmiss_pions_dummy_random_cut_all_error.Scale(1.0/nWindows)
P_kin_MMpi_pions_randsub_dummy_cut_all.Add(P_kin_MMpi_pions_dummy_prompt_cut_all, P_kin_MMpi_pions_dummy_random_cut_all, 1, -1)
P_kin_secondary_pmiss_pions_randsub_dummy_cut_all_error.Add(P_kin_secondary_pmiss_pions_dummy_prompt_cut_all_error, P_kin_secondary_pmiss_pions_dummy_random_cut_all_error, 1, -1)
P_kin_MMpi_pions_normrandsub_dummy_cut_all.Add(P_kin_MMpi_pions_dummy_prompt_cut_all, P_kin_MMpi_pions_dummy_random_cut_all, 1, -1)

print("######################################")
print("###### Random substraction done ######")
print("######################################\n")

############################################################################################################################################

# Data Normalization
P_kin_MMpi_pions_normrandsub_data_cut_all.Scale(normfac_data)

# Dummy Normalization
P_kin_MMpi_pions_normrandsub_dummy_cut_all.Scale(normfac_dummy)

# SIMC Normalization
MMpi_pions_simc_cut_all.Scale(normfac_simc)
pmiss_pions_simc_cut_all_norm_error.Scale(normfac_simc)

############################################################################################################################################

# Dummy Subtraction
P_kin_MMpi_pions_dummysub_data_cut_all.Add(P_kin_MMpi_pions_randsub_data_cut_all, P_kin_MMpi_pions_randsub_dummy_cut_all, 1, -1)
P_kin_MMpi_pions_normdummysub_data_cut_all.Add(P_kin_MMpi_pions_normrandsub_data_cut_all, P_kin_MMpi_pions_normrandsub_dummy_cut_all, 1, -1)

print("####################################")
print("###### Dummy subtraction done ######")
print("####################################\n")

###########################################################################################################################################

# Yield Calculations
dN_data_MMP = array.array('d', [0.0])
dN_simc_MMP = array.array('d', [0.0])

N_data_MMP = P_kin_MMpi_pions_normdummysub_data_cut_all.IntegralAndError(1,nbins_p,dN_data_MMP,"")
N_simc_MMP = MMpi_pions_simc_cut_all.IntegralAndError(1,nbins_p,dN_simc_MMP,"")

# Physics Yield Calculations
physics_yield_MMP = N_data_MMP

# Data/SIMC Ratio Calculations
dataSimcRatio_MMP = N_data_MMP/N_simc_MMP

# Error Propagation
dN_data_error = array.array('d', [0.0])
dN_dummy_error = array.array('d', [0.0])
dN_simc_error = array.array('d', [0.0])
dN_simc_norm_error = array.array('d', [0.0])

N_data = P_kin_secondary_pmiss_pions_randsub_data_cut_all_error.IntegralAndError(1,nbins_p,dN_data_error,"")
N_dummy = P_kin_secondary_pmiss_pions_randsub_dummy_cut_all_error.IntegralAndError(1,nbins_p,dN_dummy_error,"")
N_simc = pmiss_pions_simc_cut_all_error.IntegralAndError(1,nbins_p,dN_simc_error,"")
N_simc_norm = pmiss_pions_simc_cut_all_norm_error.IntegralAndError(1,nbins_p,dN_simc_norm_error,"")

# Data/SIMC Ratio Calculations with errors
N_data_norm = N_data / total_data_effective_charge
N_dummy_norm = N_dummy / total_dummy_effective_charge
#N_simc_norm = N_simc * normfac_simc
N_data_dummy_sub_norm = N_data_norm - N_dummy_norm
dataSimcRatio = (N_data_dummy_sub_norm) / N_simc_norm

dN_data_norm = N_data_norm * ma.sqrt((dN_data_error[0]/N_data)**2  + (total_data_effective_charge_error/total_data_effective_charge)**2)
dN_dummy_norm = N_dummy_norm * ma.sqrt((dN_dummy_error[0]/N_dummy)**2 + (total_dummy_effective_charge_error/total_dummy_effective_charge)**2)
dN_data_dummy_sub_norm = ma.sqrt((dN_data_norm)**2 + (dN_dummy_norm)**2)

dN_simc_norm = N_simc_norm * (dN_simc_error[0]/N_simc)

dataSimcRatio_err = dataSimcRatio * ma.sqrt((dN_data_dummy_sub_norm/N_data_dummy_sub_norm)**2 + (dN_simc_norm/N_simc_norm)**2)

# Physics Yield Calculations with errors
physics_yield = N_data_norm - N_dummy_norm

dN_data_norm_yield = N_data_norm * ma.sqrt((dN_data_error[0]/N_data)**2  + (total_data_effective_charge_error/total_data_effective_charge)**2)
dN_dummy_norm_yield = N_dummy_norm * ma.sqrt((dN_dummy_error[0]/N_dummy)**2 + (total_dummy_effective_charge_error/total_dummy_effective_charge)**2)
dN_data_dummy_sub_norm_yield = ma.sqrt((dN_data_norm_yield)**2 + (dN_dummy_norm_yield)**2)
physics_yield_err = dN_data_dummy_sub_norm_yield

print("="*40)
print("Data/SIMC ratio MMP = {:.3f} +/- {:.3f}".format(dataSimcRatio_MMP, dataSimcRatio_err))
print("="*40)
print("Physics Yield MMP = {:.3f} +/- {:.3f}".format(physics_yield_MMP, physics_yield_err))
print("="*40)

#############################################################################################################################################

# Define a function for fitting a Gaussian with dynamically determined FWHM range
def fit_gaussian(hist, x_min, x_max, dtype):

    print(hist.GetName(),dtype,"-"*25)

    # Find the corresponding bin numbers
    bin_min = hist.GetXaxis().FindBin(x_min)
    bin_max = hist.GetXaxis().FindBin(x_max)

    # Find the maximum value within the specified range
    max_bin = bin_min
    max_value = hist.GetBinContent(max_bin)
    for i in range(bin_min, bin_max):
        if hist.GetBinContent(i) > max_value:
            max_bin = i
            max_value = hist.GetBinContent(i)

    # Print the results
    print("max_bin", max_bin)
    print("max_value", max_value)
    print("bin_center",hist.GetBinCenter(max_bin))

    half_max = max_value*0.75

    # Find left and right bins closest to half-max value
    left_bin = max_bin
    right_bin = max_bin
    while hist.GetBinContent(left_bin) > half_max and left_bin > 1:
        left_bin -= 1
    while hist.GetBinContent(right_bin) > half_max and right_bin < hist.GetNbinsX():
        right_bin += 1

    #min_range = hist.GetBinCenter(max_bin-100)
    #max_range = hist.GetBinCenter(max_bin+100)

    min_range = hist.GetBinCenter(left_bin)
    print("min_range",min_range)
    max_range = hist.GetBinCenter(right_bin)
    print("max_range",max_range)
    print("="*40)

    hist.Fit("gaus", "Q", "", min_range, max_range)
    fit_func = hist.GetFunction('gaus')

    if dtype == "simc":
        fit_func.SetLineColor(kRed)
    if dtype == "data":
        fit_func.SetLineColor(kBlue)
#    if dtype == "dummy":
#        fit_func.SetLineColor(kGreen)

    mean = fit_func.GetParameter(1)
    mean_err = fit_func.GetParError(1)
    print("mean value",mean)
    print("meean error",mean_err)
    print("="*40)
    return [mean, mean_err]

#############################################################################################################################################

# Missing Mass cut check

# Removes stat box
ROOT.gStyle.SetOptStat(0)

# Saving histograms in PDF
c1_delta = TCanvas("c1_delta", "Variables Distributions", 100, 0, 1400, 1800)
c1_delta.Divide(2,3)
c1_delta.cd(2)
tmp_b_mean_MMpi_simc = fit_gaussian(MMpi_pions_simc_cut_all,0.90, 1.0, "simc")
b_mean_MMpi_simc = tmp_b_mean_MMpi_simc[0]
b_mean_err_MMpi_simc = tmp_b_mean_MMpi_simc[1]
tmp_b_mean_MMpi_data = fit_gaussian(P_kin_MMpi_pions_normdummysub_data_cut_all,0.90, 1.0, "data")
b_mean_MMpi_data = tmp_b_mean_MMpi_data[0]
b_mean_err_MMpi_data = tmp_b_mean_MMpi_data[1]
MMpi_pions_simc_cut_all.GetXaxis().SetRangeUser(0.0, 2.2)
MMpi_pions_simc_cut_all.SetLineColor(kRed)
MMpi_pions_simc_cut_all.SetMarkerColor(kRed)
MMpi_pions_simc_cut_all.Draw("E1")
P_kin_MMpi_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(0.0, 2.2)
P_kin_MMpi_pions_normdummysub_data_cut_all.SetLineColor(kBlue)
P_kin_MMpi_pions_normdummysub_data_cut_all.SetMarkerColor(kBlue)
P_kin_MMpi_pions_normdummysub_data_cut_all.Draw("same, E1")
MMpi_fit_func_simc = MMpi_pions_simc_cut_all.GetFunction('gaus')
MMpi_fit_func_simc.SetLineWidth(1)
MMpi_fit_func_data = P_kin_MMpi_pions_normdummysub_data_cut_all.GetFunction('gaus')
MMpi_fit_func_data.SetLineWidth(1)
c1_delta.cd(1)
c1_delta_text_lines = [
    ROOT.TText(0.5, 0.9, "Pion Physics ProdCoin Setting"),
    ROOT.TText(0.5, 0.8, '{} '.format(PHY_SETTING)),
    ROOT.TText(0.5, 0.7, "Red = SIMC"),
    ROOT.TText(0.5, 0.6, "Blue = DATA"),
    ROOT.TText(0.5, 0.4, 'Mean MMpi DATA = {:.4f}'.format(b_mean_MMpi_data)),
    ROOT.TText(0.5, 0.3, 'Mean MMpi SIMC = {:.4f}'.format(b_mean_MMpi_simc)),
    ROOT.TText(0.5, 0.2, 'Diff MMpi (SIMC - DATA) = {:.5f}'.format(b_mean_MMpi_simc - b_mean_MMpi_data))
]
for c1_delta_text in c1_delta_text_lines:
    c1_delta_text.SetTextSize(0.07)
    c1_delta_text.SetTextAlign(22)
    c1_delta_text.SetTextColor(ROOT.kGreen + 4)
    if c1_delta_text.GetTitle() == "Red = SIMC":
       c1_delta_text.SetTextColor(ROOT.kRed)  # Setting text color to red
    if c1_delta_text.GetTitle() == "Blue = DATA":
       c1_delta_text.SetTextColor(ROOT.kBlue)  # Setting text color to red
    c1_delta_text.Draw()
c1_delta.cd(3)
MMpi_pions_simc_cut_all.GetXaxis().SetRangeUser(0.0, 2.2)
MMpi_pions_simc_cut_all.SetLineColor(kRed)
MMpi_pions_simc_cut_all.SetMarkerColor(kRed)
MMpi_pions_simc_cut_all.Draw("hist")
P_kin_MMpi_pions_normdummysub_data_cut_all.GetXaxis().SetRangeUser(0.0, 2.2)
P_kin_MMpi_pions_normdummysub_data_cut_all.SetLineColor(kBlue)
P_kin_MMpi_pions_normdummysub_data_cut_all.SetMarkerColor(kBlue)
P_kin_MMpi_pions_normdummysub_data_cut_all.Draw("same, hist")
c1_delta.cd(4)
MMpi_pions_simc_cut_all.GetXaxis().SetRangeUser(0.8, 1.4)
MMpi_pions_simc_cut_all.SetLineColor(kRed)
MMpi_pions_simc_cut_all.SetMarkerColor(kRed)
MMpi_pions_simc_cut_all.Draw("hist")
# Clone the histogram
P_kin_MMpi_pions_normdummysub_data_cut_all_offset = P_kin_MMpi_pions_normdummysub_data_cut_all.Clone("P_kin_MMpi_pions_normdummysub_data_cut_all_offset")
MM_Offset = b_mean_MMpi_simc - b_mean_MMpi_data
# Apply the offset to the bin centers and content
for bin_idx in range(1, P_kin_MMpi_pions_normdummysub_data_cut_all_offset.GetNbinsX() + 1):
    bin_center = P_kin_MMpi_pions_normdummysub_data_cut_all.GetBinCenter(bin_idx)
    bin_content = P_kin_MMpi_pions_normdummysub_data_cut_all.GetBinContent(bin_idx)
    new_bin_center = bin_center + MM_Offset
    # Find the bin in the offset histogram corresponding to the new center
    new_bin_idx = P_kin_MMpi_pions_normdummysub_data_cut_all_offset.FindBin(new_bin_center)
    P_kin_MMpi_pions_normdummysub_data_cut_all_offset.SetBinContent(new_bin_idx, bin_content)
# Draw the histogram with the applied offset
P_kin_MMpi_pions_normdummysub_data_cut_all_offset.SetLineColor(kBlue)
P_kin_MMpi_pions_normdummysub_data_cut_all_offset.SetMarkerColor(kBlue)
P_kin_MMpi_pions_normdummysub_data_cut_all_offset.Draw("same, hist")
c1_delta.cd(5)
ROOT.gPad.SetLogy()
MMpi_pions_simc_cut_all.GetXaxis().SetRangeUser(0.8, 1.4)
MMpi_pions_simc_cut_all.SetLineColor(kRed)
MMpi_pions_simc_cut_all.SetMarkerColor(kRed)
MMpi_pions_simc_cut_all.Draw("hist")
P_kin_MMpi_pions_normdummysub_data_cut_all_offset.GetXaxis().SetRangeUser(0.8, 1.4)
P_kin_MMpi_pions_normdummysub_data_cut_all_offset.SetLineColor(kBlue)
P_kin_MMpi_pions_normdummysub_data_cut_all_offset.SetMarkerColor(kBlue)
P_kin_MMpi_pions_normdummysub_data_cut_all_offset.Draw("same, hist")
c1_delta.cd(6)
# Clone the histograms
MMpi_pions_simc_cut_all_MMcut = MMpi_pions_simc_cut_all.Clone("MMpi_pions_simc_cut_all_MMcut")
P_kin_MMpi_pions_normdummysub_data_cut_all_MMcut = P_kin_MMpi_pions_normdummysub_data_cut_all_offset.Clone("P_kin_MMpi_pions_normdummysub_data_cut_all_MMcut")
MM_Cut_low = 0.9
MM_Cut_high = 1.06
# Apply the cut to the cloned histograms
for bin_idx in range(1, MMpi_pions_simc_cut_all_MMcut.GetNbinsX() + 1):
    bin_center = MMpi_pions_simc_cut_all_MMcut.GetBinCenter(bin_idx)
    if bin_center <= MM_Cut_low or bin_center >= MM_Cut_high:
        MMpi_pions_simc_cut_all_MMcut.SetBinContent(bin_idx, 0)
for bin_idx in range(1, P_kin_MMpi_pions_normdummysub_data_cut_all_MMcut.GetNbinsX() + 1):
    bin_center = P_kin_MMpi_pions_normdummysub_data_cut_all_MMcut.GetBinCenter(bin_idx)
    if bin_center <= MM_Cut_low or bin_center >= MM_Cut_high:
        P_kin_MMpi_pions_normdummysub_data_cut_all_MMcut.SetBinContent(bin_idx, 0)
# Draw the histograms with the applied cut
ROOT.gPad.SetLogy()
MMpi_pions_simc_cut_all_MMcut.SetLineColor(kRed)
MMpi_pions_simc_cut_all_MMcut.SetMarkerColor(kRed)
MMpi_pions_simc_cut_all_MMcut.Draw("hist")
P_kin_MMpi_pions_normdummysub_data_cut_all_MMcut.SetLineColor(kBlue)
P_kin_MMpi_pions_normdummysub_data_cut_all_MMcut.SetMarkerColor(kBlue)
P_kin_MMpi_pions_normdummysub_data_cut_all_MMcut.Draw("same, hist")
c1_delta.Print(Pion_Analysis_Distributions)

print("="*40)
print("Missing Mass Offset = {:.5f}".format(MM_Offset))
print("="*40)

#############################################################################################################################################

# Writing Offsets and Cuts to CSV file
csv_output_path = "%s/LTSep_CSVs/mm_offset_cut_csv/%s/%s_mm_offsets_cuts_parameters.csv" % (UTILPATH, physet_dir_name, setting_name)

# Data to write
new_row = [PHY_SETTING, f"{MM_Offset:.6f}", MM_Cut_low, MM_Cut_high]

# Check if the file exists
if os.path.exists(csv_output_path):
    # Read the existing file and check for PHY_SETTING
    updated = False
    rows = []
    with open(csv_output_path, mode='r') as csv_file:
        csv_reader = csv.reader(csv_file)
        for row in csv_reader:
            if row and row[0] == PHY_SETTING:
                rows.append(new_row)  # Replace the row with the same PHY_SETTING
                updated = True
            else:
                rows.append(row)
    
    # If PHY_SETTING was not found, append the new row
    if not updated:
        rows.append(new_row)
    
    # Write the updated rows back to the file
    with open(csv_output_path, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerows(rows)
else:
    # If the file does not exist, create it and write the header and new row
    with open(csv_output_path, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(["Physics Setting", "MM_Offset", "MM_Cut_low", "MM_Cut_high"])
        csv_writer.writerow(new_row)

print(f"CSV file updated at: {csv_output_path}")

#############################################################################################################################################

# Making directories in output file
outHistFile = ROOT.TFile.Open("%s/%s_%s_MMcut_Data.root" % (OUTPATH, PHY_SETTING, MaxEvent) , "RECREATE")
d_Cut_Pion_Events_Accpt_Data_Cut_All = outHistFile.mkdir("Cut_Pion_Events_Accpt_Data_Cut_All")
d_Cut_Pion_Events_Prompt_Data_Cut_All = outHistFile.mkdir("Cut_Pion_Events_Prompt_Data_Cut_All")
d_Cut_Pion_Events_Random_Data_Cut_All = outHistFile.mkdir("Cut_Pion_Events_Random_Data_Cut_All")
d_Cut_Pion_Events_RandomSub_Data = outHistFile.mkdir("Cut_Pion_Events_RandomSub_Data")
d_Cut_Pion_Events_DummySub_RandomSub_Data = outHistFile.mkdir("Cut_Pion_Events_DummySub_RandomSub_Data")
d_Cut_Pion_Events_Norm_DummySub_RandomSub_Data = outHistFile.mkdir("Cut_Pion_Events_Norm_DummySub_RandomSub_Data")
d_Cut_Pion_Events_Norm_SIMC_Data = outHistFile.mkdir("Cut_Pion_Events_Norm_SIMC_Data")

# Writing Histograms for pions
d_Cut_Pion_Events_Accpt_Data_Cut_All.cd()
P_kin_MMpi_pions_data_accpt_cut_all.Write()

d_Cut_Pion_Events_Prompt_Data_Cut_All.cd()
P_kin_MMpi_pions_data_prompt_cut_all.Write()

d_Cut_Pion_Events_Random_Data_Cut_All.cd()
P_kin_MMpi_pions_data_random_cut_all.Write()

d_Cut_Pion_Events_RandomSub_Data.cd()
P_kin_MMpi_pions_randsub_data_cut_all.Write()

d_Cut_Pion_Events_DummySub_RandomSub_Data.cd()
P_kin_MMpi_pions_dummysub_data_cut_all.Write()

d_Cut_Pion_Events_Norm_DummySub_RandomSub_Data.cd()
P_kin_MMpi_pions_normdummysub_data_cut_all.Write()

d_Cut_Pion_Events_Norm_SIMC_Data.cd()
MMpi_pions_simc_cut_all.Write()

print("####################################")
print("###### Histogram writing done ######")
print("####################################\n")

infile_DATA.Close() 
infile_DUMMY.Close()
infile_SIMC.Close()
outHistFile.Close()

print ("Processing Complete")