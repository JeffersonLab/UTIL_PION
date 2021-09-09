#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2021-08-31 03:20:46 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
###################################################################################
# Created - 20/July/21, Author - Muhammad Junaid, University of Regina, Canada
###################################################################################
# Python version of the pion plotting script. Now utilises uproot to select event of each type and writes them to a root file.
# Python should allow for easier reading of databases storing diferent variables.
# This version of script is for physics analysis experts
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

# Defining some constants here
minrangeuser = 0 # min range for -t vs phi plot
maxrangeuser = 0.5 # max range for -t vs phi plot
minbin = 0.92 # minbin for selecting neutrons events in missing mass distribution
maxbin = 0.98 # maxbin for selecting neutrons events in missing mass distribution

##################################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=4:
    print("!!!!! ERROR !!!!!\n Expected 4 arguments\n Usage is with - ROOTfilePrefix RunNumber MaxEvents spec \n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################

# Input params - run number and max number of events
runNum = sys.argv[1]
MaxEvent = sys.argv[2]
ROOTPrefix = sys.argv[3]
spec = sys.argv[4]

spec = spec.upper()

USER = subprocess.getstatusoutput("whoami") # Grab user info for file finding
HOST = subprocess.getstatusoutput("hostname")

if ("farm" in HOST[1]):
    REPLAYPATH = "/group/c-pionlt/USERS/%s/hallc_replay_lt" % USER[1]
#    REPLAYPATH = "/group/c-kaonlt/USERS/%s/hallc_replay_lt" % USER[1]

elif ("qcd" in HOST[1]):
    REPLAYPATH = "/group/c-pionlt/USERS/%s/hallc_replay_lt" % USER[1]
#    REPLAYPATH = "/group/c-kaonlt/USERS/%s/hallc_replay_lt" % USER[1]

elif ("phys.uregina" in HOST[1]):
    REPLAYPATH = "/home/%s/work/JLab/hallc_replay_lt" % USER[1]

elif ("cdaq" in HOST[1]):
    REPLAYPATH = "/home/cdaq/hallc-online/hallc_replay_lt"

elif("skynet" in HOST[1]):
    REPLAYPATH = "/home/%s/Work/JLab/hallc_replay_lt" % USER[1]

#################################################################################################################################################

# Add more path setting as needed in a similar manner                                                                                                                                                          
OUTPATH = "%s/UTIL_PION/OUTPUT/Analysis/HeeP" % REPLAYPATH        # Output folder location                                                                                                     
sys.path.insert(0, '%s/UTIL_PION/bin/python/' % REPLAYPATH)
import kaonlt as klt # Import kaonlt module, need the path setting line above prior to importing this                                                                                                         
print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER[1], HOST[1], REPLAYPATH))
Proton_Analysis_Distributions = "%s/%s_%s_sw_heep_%s_Proton_Analysis_Distributions.pdf" % (OUTPATH, runNum, MaxEvent, spec)

#################################################################################################################################################

# Construct the name of the rootfile based upon the info we provided
rootName = "%s/UTIL_PION/OUTPUT/Analysis/HeeP/%s_%s_%s.root" % (REPLAYPATH, runNum, MaxEvent, ROOTPrefix)     # Input file location and variables taking
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
print ("Attempting to process %s" %(rootName))
if os.path.isfile(rootName):
    print ("%s exists, attempting to process" % (rootName))
else:
    print ("%s not found - do you have the correct sym link/folder set up?" % (rootName))
    sys.exit(4)
print("Output path checks out, outputting to %s" % (OUTPATH))

###############################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Read stuff from the main event tree
infile = ROOT.TFile.Open(rootName, "READ")

Uncut_Proton_Events_tree = infile.Get("Uncut_Proton_Events")
Cut_Proton_Events_All_tree = infile.Get("Cut_Proton_Events_All")
Cut_Proton_Events_Prompt_tree = infile.Get("Cut_Proton_Events_Prompt")
Cut_Proton_Events_Random_tree = infile.Get("Cut_Proton_Events_Random")

###################################################################################################################################################

# Defining Histograms for Protons
if spec == "HMS":
    H_gtr_beta_protons_uncut = ROOT.TH1D("H_gtr_beta_protons_uncut", "HMS #beta; HMS_gtr_#beta; Counts", 200, 0.8, 1.2)
    H_gtr_xp_protons_uncut = ROOT.TH1D("H_gtr_xp_protons_uncut", "HMS x'; HMS_gtr_xp; Counts", 200, -0.2, 0.2)
    H_gtr_yp_protons_uncut = ROOT.TH1D("H_gtr_yp_protons_uncut", "HMS y'; HMS_gtr_yp; Counts", 200, -0.2, 0.2)
    H_gtr_dp_protons_uncut = ROOT.TH1D("H_gtr_dp_protons_uncut", "HMS #delta; HMS_gtr_dp; Counts", 200, -15, 15)
    H_gtr_p_protons_uncut = ROOT.TH1D("H_gtr_p_protons_uncut", "HMS p; HMS_gtr_p; Counts", 200, 4, 8)
    H_hod_goodscinhit_protons_uncut = ROOT.TH1D("H_hod_goodscinhit_protons_uncut", "HMS hod goodscinhit; HMS_hod_goodscinhi; Counts", 200, 0.7, 1.3)
    H_hod_goodstarttime_protons_uncut = ROOT.TH1D("H_hod_goodstarttime_protons_uncut", "HMS hod goodstarttime; HMS_hod_goodstarttime; Counts", 200, 0.7, 1.3)
    H_cal_etotnorm_protons_uncut = ROOT.TH1D("H_cal_etotnorm_protons_uncut", "HMS cal etotnorm; HMS_cal_etotnorm; Counts", 200, 0.2, 1.8)
    H_cal_etottracknorm_protons_uncut = ROOT.TH1D("H_cal_etottracknorm_protons_uncut", "HMS cal etottracknorm; HMS_cal_etottracknorm; Counts", 300, 0.2, 1.8)
    H_cer_npeSum_protons_uncut = ROOT.TH1D("H_cer_npeSum_protons_uncut", "HMS cer npeSum; HMS_cer_npeSum; Counts", 200, 0, 50)
    H_kin_MMp_protons_uncut = ROOT.TH1D("H_kin_MMp_protons_uncut", "MIssing Mass; MM_{p}; Counts", 200, 0.5, 1.8)
    H_RFTime_Dist_protons_uncut = ROOT.TH1D("H_RFTime_Dist_protons_uncut", "HMS RFTime; HMS_RFTime; Counts", 200, 0, 4)
    
    H_gtr_beta_protons_cut_all = ROOT.TH1D("H_gtr_beta_protons_cut_all", "HMS #beta; HMS_gtr_#beta; Counts", 200, 0.8, 1.2)
    H_gtr_xp_protons_cut_all = ROOT.TH1D("H_gtr_xp_protons_cut_all", "HMS x'; HMS_gtr_xp; Counts", 200, -0.2, 0.2)
    H_gtr_yp_protons_cut_all = ROOT.TH1D("H_gtr_yp_protons_cut_all", "HMS y'; HMS_gtr_yp; Counts", 200, -0.2, 0.2)
    H_gtr_dp_protons_cut_all = ROOT.TH1D("H_gtr_dp_protons_cut_all", "HMS #delta; HMS_gtr_dp; Counts", 200, -15, 15)
    H_gtr_p_protons_cut_all = ROOT.TH1D("H_gtr_p_protons_cut_all", "HMS p; HMS_gtr_p; Counts", 200, 4, 8)
    H_hod_goodscinhit_protons_cut_all = ROOT.TH1D("H_hod_goodscinhit_protons_cut_all", "HMS hod goodscinhit; HMS_hod_goodscinhit; Counts", 200, 0.7, 1.3)
    H_hod_goodstarttime_protons_cut_all = ROOT.TH1D("H_hod_goodstarttime_protons_cut_all", "HMS hod goodstarttime; HMS_hod_goodstarttime; Counts", 200, 0.7, 1.3)
    H_cal_etotnorm_protons_cut_all = ROOT.TH1D("H_cal_etotnorm_protons_cut_all", "HMS cal etotnorm; HMS_cal_etotnorm; Counts", 200, 0.6, 1.4)
    H_cal_etottracknorm_protons_cut_all = ROOT.TH1D("H_cal_etottracknorm_protons_cut_all", "HMS cal etottracknorm; HMS_cal_etottracknorm; Counts", 200, 0.6, 1.4)
    H_cer_npeSum_protons_cut_all = ROOT.TH1D("H_cer_npeSum_protons_cut_all", "HMS cer npeSum; HMS_cer_npeSum; Counts", 200, 0, 50)
    H_RFTime_Dist_protons_cut_all = ROOT.TH1D("H_RFTime_Dist_protons_cut_all", "HMS RFTime; HMS_RFTime; Counts", 200, 0, 4)
    H_kin_MMp_protons_cut_all = ROOT.TH1D("H_kin_MMpi_protons_cut_all", "Missing Mass; MM_{#pi}; Counts", 200, 0.5, 1.8)

    H_gtr_beta_protons_cut_prompt = ROOT.TH1D("H_gtr_beta_protons_cut_prompt", "SHMS beta; SHMS_#beta; Counts", 200, 0.8, 1.2)
    H_RFTime_Dist_protons_cut_prompt = ROOT.TH1D("H_RFTime_Dist_protons_cut_prompt", "SHMS RFTime; SHMS_RFTime; Counts", 200, 0, 4)
    H_kin_MMp_protons_cut_prompt = ROOT.TH1D("H_kin_MMp_protons_cut_prompt", "Missing Mass; MM_{p}; Counts", 200, 0.5, 1.8)

    H_gtr_beta_protons_cut_random = ROOT.TH1D("H_gtr_beta_protons_cut_random", "SHMS beta; SHMS_#beta; Counts", 200, 0.8, 1.2)
    H_RFTime_Dist_protons_cut_random = ROOT.TH1D("H_RFTime_Dist_protons_cut_random", "SHMS RFTime; SHMS_RFTime; Counts", 200, 0, 4)
    H_kin_MMp_protons_cut_random = ROOT.TH1D("H_kin_MMp_protons_cut_random", "Missing Mass; MM_{p}; Counts", 200, 0.5, 1.8)

    H_kin_MMp_protons_cut_random_scaled = ROOT.TH1D("H_kin_MMp_protons_cut_random_scaled", "Missing Mass; MM_{p}; Counts", 200, 0.5, 1.8)
    H_kin_MMp_protons_cut_random_sub = ROOT.TH1D("H_kin_MMp_protons_cut_random_sub", "Missing Mass Rndm Sub; MM_{p}; Counts", 200, 0.5, 1.8)
    
if spec == "SHMS":
    P_gtr_beta_protons_uncut = ROOT.TH1D("P_gtr_beta_protons_uncut", "SHMS #beta; SHMS_gtr_#beta; Counts", 200, 0.7, 1.3)
    P_gtr_xp_protons_uncut = ROOT.TH1D("P_gtr_xp_protons_uncut", "SHMS x'; SHMS_gtr_xp; Counts", 200, -0.2, 0.2)
    P_gtr_yp_protons_uncut = ROOT.TH1D("P_gtr_yp_protons_uncut", "SHMS y'; SHMS_gtr_yp; Counts", 200, -0.2, 0.2)
    P_gtr_dp_protons_uncut = ROOT.TH1D("P_gtr_dp_protons_uncut", "SHMS delta; SHMS_gtr_dp; Counts", 200, -30, 30)
    P_gtr_p_protons_uncut = ROOT.TH1D("P_gtr_p_protons_uncut", "SHMS p; SHMS_gtr_p; Counts", 200, 4, 8)
    P_hod_goodscinhit_protons_uncut = ROOT.TH1D("P_hod_goodscinhit_protons_uncut", "SHMS hod goodscinhit; SHMS_hod_goodscinhit; Counts", 200, 0.7, 1.3)
    P_hod_goodstarttime_protons_uncut = ROOT.TH1D("P_hod_goodstarttime_protons_uncut", "SHMS hod goodstarttime; SHMS_hod_goodstarttime; Counts", 200, 0.7, 1.3)
    P_cal_etotnorm_protons_uncut = ROOT.TH1D("P_cal_etotnorm_protons_uncut", "SHMS cal etotnorm; SHMS_cal_etotnorm; Counts", 200, 0, 1)
    P_cal_etottracknorm_protons_uncut = ROOT.TH1D("P_cal_etottracknorm_protons_uncut", "SHMS cal etottracknorm; SHMS_cal_etottracknorm; Counts", 200, 0, 1.6)
    P_hgcer_npeSum_protons_uncut = ROOT.TH1D("P_hgcer_npeSum_protons_uncut", "SHMS HGC npeSum; SHMS_hgcer_npeSum; Counts", 200, 0, 50)
    P_hgcer_xAtCer_protons_uncut = ROOT.TH1D("P_hgcer_xAtCer_protons_uncut", "SHMS HGC xAtCer; SHMS_hgcer_xAtCer; Counts", 200, -60, 60)
    P_hgcer_yAtCer_protons_uncut = ROOT.TH1D("P_hgcer_yAtCer_protons_uncut", "SHMS HGC yAtCer; SHMS_hgcer_yAtCer; Counts", 200, -50, 50)
    P_aero_npeSum_protons_uncut = ROOT.TH1D("P_aero_npeSum_protons_uncut", "SHMS aero npeSum; SHMS_aero_npeSum; Counts", 200, 0, 50)
    P_aero_xAtAero_protons_uncut = ROOT.TH1D("P_acero_xAtAero_protons_uncut", "SHMS aero xAtAero; SHMS_aero_xAtAero; Counts", 200, -60, 60)
    P_aero_yAtAero_protons_uncut = ROOT.TH1D("P_aero_yAtAero_protons_uncut", "SHMS aero yAtAero; SHMS_aero_yAtAero; Counts", 200, -50, 50)
    #P_ngcer_npeSum_protons_uncut = ROOT.TH1D("P_ngcer_npeSum_protons_uncut", "SHMS NGC npeSum; SHMS_ngcer_npeSum; Counts", 200, 0, 0.5)
    #P_ngcer_xAtCer_protons_uncut = ROOT.TH1D("P_ngcer_xAtCer_protons_uncut", "SHMS NGC xAtCer; SHMS_ngcer_xAtCer; Counts", 200, -70, 50)
    #P_ngcer_yAtCer_protons_uncut = ROOT.TH1D("P_ngcer_yAtCer_protons_uncut", "SHMS NGC yAtCer; SHMS_ngcer_yAtCer; Counts", 200, -50, 50)
    P_kin_MMp_protons_uncut = ROOT.TH1D("P_kin_MMp_protons_uncut", "MIssing Mass; MM_{p}; Counts", 200, 0.5, 1.8)
    P_RFTime_Dist_protons_uncut = ROOT.TH1D("P_RFTime_Dist_protons_uncut", "SHMS RFTime; SHMS_RFTime; Counts", 200, 0, 4)
    
    P_gtr_beta_protons_cut_all = ROOT.TH1D("P_gtr_beta_protons_cut_all", "SHMS #beta; SHMS_gtr_#beta; Counts", 200, 0.8, 1.2)
    P_gtr_xp_protons_cut_all = ROOT.TH1D("P_gtr_xp_protons_cut_all", "SHMS x'; SHMS_gtr_xp; Counts", 200, -0.2, 0.2)
    P_gtr_yp_protons_cut_all = ROOT.TH1D("P_gtr_yp_protons_cut_all", "SHMS y'; SHMS_gtr_yp; Counts", 200, -0.2, 0.2)
    P_gtr_dp_protons_cut_all = ROOT.TH1D("P_gtr_dp_protons_cut_all", "SHMS #delta; SHMS_gtr_dp; Counts", 200, -15, 15)
    P_gtr_p_protons_cut_all = ROOT.TH1D("P_gtr_p_protons_cut_all", "SHMS p; SHMS_gtr_p; Counts", 200, 4, 8)
    P_hod_goodscinhit_protons_cut_all = ROOT.TH1D("P_hod_goodscinhit_protons_cut_all", "SHMS hod goodscinhit; SHMS_hod_goodscinhit; Counts", 200, 0.7, 1.3)
    P_hod_goodstarttime_protons_cut_all = ROOT.TH1D("P_hod_goodstarttime_protons_cut_all", "SHMS hod goodstarttime; SHMS_hod_goodstarttime; Counts", 200, 0.7, 1.3)
    P_cal_etotnorm_protons_cut_all = ROOT.TH1D("P_cal_etotnorm_protons_cut_all", "SHMS cal etotnorm; SHMS_cal_etotnorm; Counts", 200, 0, 1.2)
    P_cal_etottracknorm_protons_cut_all = ROOT.TH1D("P_cal_etottracknorm_protons_cut_all", "SHMS cal etottracknorm; SHMS_cal_etottracknorm; Counts", 200, 0, 1.6)
    P_hgcer_npeSum_protons_cut_all = ROOT.TH1D("P_hgcer_npeSum_protons_cut_all", "SHMS HGC npeSum; SHMS_hgcer_npeSum; Counts", 200, 0, 50)
    P_hgcer_xAtCer_protons_cut_all = ROOT.TH1D("P_hgcer_xAtCer_protons_cut_all", "SHMS HGC xAtCer; SHMS_hgcer_xAtCer; Counts", 200, -40, 30)
    P_hgcer_yAtCer_protons_cut_all = ROOT.TH1D("P_hgcer_yAtCer_protons_cut_all", "SHMS HGC yAtCer; SHMS_hgcer_yAtCer; Counts", 200, -30, 30)
    P_aero_npeSum_protons_cut_all = ROOT.TH1D("P_aero_npeSum_protons_cut_all", "SHMS aero npeSum; SHMS_aero_npeSum; Counts", 200, 0, 50)
    P_aero_xAtAero_protons_cut_all = ROOT.TH1D("P_acero_xAtAero_protons_cut_all", "SHMS aero xAtAero; SHMS_aero_xAtAero; Counts", 200, -40, 30)
    P_aero_yAtAero_protons_cut_all = ROOT.TH1D("P_aero_yAtAero_protons_cut_all", "SHMS aero yAtAero; SHMS_aero_yAtAero; Counts", 200, -30, 30)
    #P_ngcer_npeSum_protons_cut_all = ROOT.TH1D("P_ngcer_npeSum_protons_cut_all", "SHMS NGC npeSum; SHMS_ngcer_npeSum; Counts", 200, -10, 50)
    #P_ngcer_xAtCer_protons_cut_all = ROOT.TH1D("P_ngcer_xAtCer_protons_cut_all", "SHMS NGC xAtCer; SHMS_ngcer_xAtCer; Counts", 200, -40, 30)
    #P_ngcer_yAtCer_protons_cut_all = ROOT.TH1D("P_ngcer_yAtCer_protons_cut_all", "SHMS NGC yAtCer; SHMS_ngcer_yAtCer; Counts", 200, -30, 30)
    P_RFTime_Dist_protons_cut_all = ROOT.TH1D("P_RFTime_Dist_protons_cut_all", "SHMS RFTime; SHMS_RFTime; Counts", 200, 0, 4)
    P_kin_MMp_protons_cut_all = ROOT.TH1D("P_kin_MMpi_protons_cut_all", "Missing Mass; MM_{#pi}; Counts", 200, 0.5, 1.8)
    
    P_gtr_beta_protons_cut_prompt = ROOT.TH1D("P_gtr_beta_protons_cut_prompt", "SHMS beta; SHMS_#beta; Counts", 200, 0.8, 1.2)
    P_RFTime_Dist_protons_cut_prompt = ROOT.TH1D("P_RFTime_Dist_protons_cut_prompt", "SHMS RFTime; SHMS_RFTime; Counts", 200, 0, 4)
    P_kin_MMp_protons_cut_prompt = ROOT.TH1D("P_kin_MMp_protons_cut_prompt", "Missing Mass; MM_{p}; Counts", 200, 0.5, 1.8)

    P_gtr_beta_protons_cut_random = ROOT.TH1D("P_gtr_beta_protons_cut_random", "SHMS beta; SHMS_#beta; Counts", 200, 0.8, 1.2)
    P_RFTime_Dist_protons_cut_random = ROOT.TH1D("P_RFTime_Dist_protons_cut_random", "SHMS RFTime; SHMS_RFTime; Counts", 200, 0, 4)
    P_kin_MMp_protons_cut_random = ROOT.TH1D("P_kin_MMp_protons_cut_random", "Missing Mass; MM_{p}; Counts", 200, 0.5, 1.8)

    P_kin_MMp_protons_cut_random_scaled = ROOT.TH1D("P_kin_MMp_protons_cut_random_scaled", "Missing Mass; MM_{p}; Counts", 200, 0.5, 1.8)
    P_kin_MMp_protons_cut_random_sub = ROOT.TH1D("P_kin_MMp_protons_cut_random_sub", "Missing Mass Rndm Sub; MM_{p}; Counts", 200, 0.5, 1.8)    
    
    
CTime_epCoinTime_ROC1_protons_uncut = ROOT.TH1D("CTime_epCoinTime_ROC1_protons_uncut", "Electron-Proton CTime; e p Coin_Time; Counts", 200, -60, 60)
CTime_epCoinTime_ROC1_protons_cut_all = ROOT.TH1D("CTime_epCoinTime_ROC1_protons_cut_all", "Electron-Proton CTime; e p Coin_Time; Counts", 200, -50, 50)
CTime_epCoinTime_ROC1_protons_cut_prompt = ROOT.TH1D("CTime_epCoinTime_ROC1_protons_cut_prompt", "Electron-Proton CTime; e p Coin_Time; Counts", 200, -2, 2)
CTime_epCoinTime_ROC1_protons_cut_random = ROOT.TH1D("CTime_epCoinTime_ROC1_protons_cut_random", "Electron-Proton CTime; e p Coin_Time; Counts", 200, -40, 40)

###################################################################################################################################################

# 2D Histograms for protons
if spec == "HMS":
    H_cal_etottracknorm_vs_H_cer_npeSum_protons_uncut = ROOT.TH2D("H_cal_etottracknorm_vs_H_cer_npeSum_protons_uncut","HMS cal etottracknorm vs HMS cer npeSum (no cut); H_cal_etottracknorm; H_cer_npeSum",100, 0, 2, 100, 0, 40)
    CTime_epCoinTime_ROC1_vs_H_kin_MMp_protons_uncut = ROOT.TH2D("CTime_epCoinTime_ROC1_vs_H_kin_MMp_protons_uncut","Electron-Proton CTime vs Missing Mass (no cut); e p Coin_Time; MM_{p}", 200, -40, 40, 200, 0, 2)
    CTime_epCoinTime_ROC1_vs_beta_protons_uncut = ROOT.TH2D("CTime_epCoinTime_ROC1_vs_beta_protons_uncut", "Electron-Proton CTime vs SHMS #beta (no cut); e p Coin_Time; SHMS_#beta", 200, -40, 40, 200, 0, 2)
    H_kin_MMp_vs_H_RFTime_protons_uncut = ROOT.TH2D("H_kin_MMp_vs_H_RFTime_protons_uncut", "Missing Mass vs SHMS RFTime (no cuts); MM_{p}; SHMS_RFTime_Dist", 100, 0, 2, 100, 0, 4)
    
    H_cal_etottracknorm_vs_H_cer_npeSum_protons_cut_all = ROOT.TH2D("H_cal_etottracknorm_vs_H_cer_npeSum_protons_cut_all","HMS cal etottracknorm vs HMS cer npeSum (with cuts); H_cal_etottracknorm; H_cer_npeSum",100, 0.5, 1.5, 100, 0, 40)
    CTime_epCoinTime_ROC1_vs_H_kin_MMp_protons_cut_all = ROOT.TH2D("CTime_epCoinTime_ROC1_vs_H_kin_MMp_protons_cut_all","Electron-Proton CTime vs Missing Mass (with cuts); e p Coin_Time; MM_{p}", 100, -2, 2, 100, 0, 2)
    CTime_epCoinTime_ROC1_vs_beta_protons_cut_all = ROOT.TH2D("CTime_epCoinTime_ROC1_vs_beta_protons_cut_all", "Electron-Proton CTime vs SHMS #beta (with cuts); e p Coin_Time; SHMS_#beta", 100, -2, 2, 100, 0.6, 1.4)
    H_kin_MMp_vs_H_RFTime_protons_cut_all = ROOT.TH2D("H_kin_MMp_vs_H_RFTime_protons_cut_all", "Missing Mass vs SHMS RFTime (with cuts); MM_{p}; SHMS_RFTime_Dist", 100, 0, 2, 100, 0, 4)

    H_kin_MMp_vs_CTime_epCoinTime_ROC1_protons_cut_prompt = ROOT.TH2D("H_kin_MMp_vs_CTime_epCoinTime_ROC1_protons_cut_prompt","Missing Mass vs Electron-Proton CTime; MM_{p}; e p Coin_Time",100, 0, 2, 100, -2, 2)
    
if spec == "SHMS":
    P_hgcer_npeSum_vs_aero_npeSum_protons_uncut = ROOT.TH2D("P_hgcer_npeSum_vs_aero_npeSum_protons_uncut", "SHMS HGC npeSum vs SHMS Aero npeSum (no cut); SHMS_hgcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
    CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_uncut = ROOT.TH2D("CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_uncut","Electron-Proton CTime vs Missing Mass (no cut); e p Coin_Time; MM_{p}", 200, -40, 40, 200, 0, 2)
    P_hgcer_yAtCer_vs_hgcer_xAtCer_protons_uncut = ROOT.TH2D("P_hgcer_yAtCer_vs_hgcer_xAtCer_protons_uncut", "SHMS HGC yAtCer vs SHMS HGC xAtCer (no cut); SHMS_hgcer_yAtCer; SHMS_hgcer_xAtCer", 100, -50, 50, 100, -50, 50)
    P_aero_yAtAero_vs_aero_xAtAero_protons_uncut = ROOT.TH2D("P_aero_yAtAero_vs_aero_xAtAero_protons_uncut", "SHMS aero yAtAero vs SHMS aero xAtAero (no cut); SHMS_aero_yAtAero; SHMS_aero_xAtAero", 100, -50, 50, 100, -50, 50)
    CTime_epCoinTime_ROC1_vs_beta_protons_uncut = ROOT.TH2D("CTime_epCoinTime_ROC1_vs_beta_protons_uncut", "Electron-Proton CTime vs SHMS #beta (no cut); e p Coin_Time; SHMS_#beta", 200, -40, 40, 200, 0, 2)
    P_kin_MMp_vs_P_RFTime_protons_uncut = ROOT.TH2D("P_kin_MMp_vs_P_RFTime_protons_uncut", "Missing Mass vs SHMS RFTime (no cuts); MM_{p}; SHMS_RFTime_Dist", 100, 0, 2, 100, 0, 4)
    #P_cal_etottracknorm_vs_P_ngcer_npeSum_protons_uncut = ROOT.TH2D("P_cal_etottracknorm_vs_P_ngcer_npeSum_protons_uncut", "SHMS cal etottracknorm vs SHMS NGC xAtCer (no cut); SHMS_cal_etottracknorm; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
    #P_ngcer_yAtCer_vs_ngcer_xAtCer_protons_uncut = ROOT.TH2D("P_ngcer_yAtCer_vs_ngcer_xAtCer_protons_uncut", "SHMS NGC yAtCer vs SHMS NGC xAtCer (no cut); SHMS_ngcer_yAtCer; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
    #P_ngcer_npeSum_vs_hgcer_npeSum_protons_uncut = ROOT.TH2D("P_ngcer_npeSum_vs_hgcer_npeSum_protons_uncut", "SHMS NGC npeSum vs SHMS HGC npeSum (no cut); SHMS_ngcer_npeSum; SHMS_hgcer_npeSum", 100, 0, 50, 100, 0, 50)
    #P_ngcer_npeSum_vs_aero_npeSum_protons_uncut = ROOT.TH2D("P_ngcer_npeSum_vs_aero_npeSum_protons_uncut", "SHMS NGC npeSum vs SHMS aero npeSum (no cut); SHMS_ngcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)

    P_hgcer_npeSum_vs_aero_npeSum_protons_cut_all = ROOT.TH2D("P_hgcer_npeSum_vs_aero_npeSum_protons_cut_all", "SHMS HGC npeSum vs SHMS aero npeSum (with cuts); SHMS_hgcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
    CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_cut_all = ROOT.TH2D("CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_cut_all","Electron-Proton CTime vs Missing Mass (with cuts); e p Coin_Time; MM_{p}", 100, -2, 2, 100, 0, 2)
    P_hgcer_yAtCer_vs_hgcer_xAtCer_protons_cut_all = ROOT.TH2D("P_hgcer_yAtCer_vs_hgcer_xAtCer_protons_cut_all", "SHMS HGC yAtCer vs SHMS HGC xAtCer (with cuts); SHMS_hgcer_yAtCer; SHMS_hgcer_xAtCer", 100, -50, 50, 100, -50, 50)
    P_aero_yAtAero_vs_aero_xAtAero_protons_cut_all = ROOT.TH2D("P_aero_yAtAero_vs_aero_xAtAero_protons_cut_all", "SHMS aero yAtAero vs SHMS aero xAtAero (with cuts); SHMS_aero_yAtAero; SHMS_aero_xAtAero", 100, -50, 50, 100, -50, 50)
    CTime_epCoinTime_ROC1_vs_beta_protons_cut_all = ROOT.TH2D("CTime_epCoinTime_ROC1_vs_beta_protons_cut_all", "Electron-Proton CTime vs SHMS #beta (with cuts); e p Coin_Time; SHMS_#beta", 100, -2, 2, 100, 0.6, 1.4)
    P_kin_MMp_vs_P_RFTime_protons_cut_all = ROOT.TH2D("P_kin_MMp_vs_P_RFTime_protons_cut_all", "Missing Mass vs SHMS RFTime (with cuts); MM_{p}; SHMS_RFTime_Dist", 100, 0, 2, 100, 0, 4)
    #P_cal_etottracknorm_vs_P_ngcer_npeSum_protons_cut_all = ROOT.TH2D("P_cal_etottracknorm_vs_P_ngcer_npeSum_protons_cut_all", "P cal etottracknorm vs SHMS NGC xAtCer (with cuts); SHMS_cal_etottracknorm; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
    #P_ngcer_yAtCer_vs_ngcer_xAtCer_protons_cut_all = ROOT.TH2D("P_ngcer_yAtCer_vs_ngcer_xAtCer_protons_cut_all", "SHMS NGC yAtCer vs SHMS NGC xAtCer (with cuts); SHMS_ngcer_yAtCer; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
    #P_ngcer_npeSum_vs_hgcer_npeSum_protons_cut_all = ROOT.TH2D("P_ngcer_npeSum_vs_hgcer_npeSum_protons_cut_all", "SHMS NGC npeSum vs SHMS HGC npeSum (with cuts); SHMS_ngcer_npeSum; SHMS_hgcer_npeSum", 100, 0, 50, 100, 0, 50)
    #P_ngcer_npeSum_vs_aero_npeSum_protons_cut_all = ROOT.TH2D("P_ngcer_npeSum_vs_aero_npeSum_protons_cut_all", "SHMS NGC npeSum vs SHMS aero npeSum (with cuts); SHMS_ngcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)

    P_kin_MMp_vs_CTime_epCoinTime_ROC1_protons_cut_prompt = ROOT.TH2D("P_kin_MMp_vs_CTime_epCoinTime_ROC1_protons_cut_prompt","Missing Mass vs Electron-Proton CTime; MM_{p}; e p Coin_Time",100, 0, 2, 100, -2, 2)

#################################################################################################################################################

# Filling Histograms for Protons
for event in Uncut_Proton_Events_tree:
    if spec == "HMS":
        H_gtr_beta_protons_uncut.Fill(event.H_gtr_beta)
        H_gtr_xp_protons_uncut.Fill(event.H_gtr_xp)
        H_gtr_yp_protons_uncut.Fill(event.H_gtr_yp)
        H_gtr_dp_protons_uncut.Fill(event.H_gtr_dp)
        H_hod_goodscinhit_protons_uncut.Fill(event.H_hod_goodscinhit)
        H_hod_goodstarttime_protons_uncut.Fill(event.H_hod_goodstarttime)
        H_cal_etotnorm_protons_uncut.Fill(event.H_cal_etotnorm)
        H_cal_etottracknorm_protons_uncut.Fill(event.H_cal_etottracknorm)
        H_cer_npeSum_protons_uncut.Fill(event.H_cer_npeSum)
        H_kin_MMp_protons_uncut.Fill(event.MMp)
        H_RFTime_Dist_protons_uncut.Fill(event.H_RF_Dist)
        H_cal_etottracknorm_vs_H_cer_npeSum_protons_uncut.Fill(event.H_cal_etottracknorm, event.H_cer_npeSum)
        CTime_epCoinTime_ROC1_vs_H_kin_MMp_protons_uncut.Fill(event.CTime_epCoinTime_ROC1, event.MMp)
        H_kin_MMp_vs_H_RFTime_protons_uncut.Fill(event.MMp, event.H_RF_Dist)
        CTime_epCoinTime_ROC1_vs_beta_protons_uncut.Fill(event.CTime_epCoinTime_ROC1, event.H_gtr_beta)
        
    if spec == "SHMS":    
        P_gtr_beta_protons_uncut.Fill(event.P_gtr_beta)
        P_gtr_xp_protons_uncut.Fill(event.P_gtr_xp)
        P_gtr_yp_protons_uncut.Fill(event.P_gtr_yp)
        P_gtr_dp_protons_uncut.Fill(event.P_gtr_dp)
        P_gtr_p_protons_uncut.Fill(event.P_gtr_p)
        P_hod_goodscinhit_protons_uncut.Fill(event.P_hod_goodscinhit)
        P_hod_goodstarttime_protons_uncut.Fill(event.P_hod_goodstarttime)
        P_cal_etotnorm_protons_uncut.Fill(event.P_cal_etotnorm)
        P_cal_etottracknorm_protons_uncut.Fill(event.P_cal_etottracknorm)
        P_hgcer_npeSum_protons_uncut.Fill(event.P_hgcer_npeSum)
        P_hgcer_xAtCer_protons_uncut.Fill(event.P_hgcer_xAtCer)
        P_hgcer_yAtCer_protons_uncut.Fill(event.P_hgcer_yAtCer)
        P_aero_npeSum_protons_uncut.Fill(event.P_aero_npeSum)
        P_aero_xAtAero_protons_uncut.Fill(event.P_aero_xAtAero)
        P_aero_yAtAero_protons_uncut.Fill(event.P_aero_yAtAero)
        #    P_ngcer_npeSum_protons_uncut.Fill(event.P_ngcer_npeSum)
        #    P_ngcer_xAtCer_protons_uncut.Fill(event.P_ngcer_xAtCer)
        #    P_ngcer_yAtCer_protons_uncut.Fill(event.P_ngcer_yAtCer)
        P_kin_MMp_protons_uncut.Fill(event.MMp)
        P_RFTime_Dist_protons_uncut.Fill(event.P_RF_Dist)
        P_hgcer_npeSum_vs_aero_npeSum_protons_uncut.Fill(event.P_hgcer_npeSum, event.P_aero_npeSum)
        CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_uncut.Fill(event.CTime_epCoinTime_ROC1, event.MMp)
        P_kin_MMp_vs_P_RFTime_protons_uncut.Fill(event.MMp, event.P_RF_Dist)
        P_hgcer_yAtCer_vs_hgcer_xAtCer_protons_uncut.Fill(event.P_hgcer_yAtCer, event.P_hgcer_xAtCer)
        P_aero_yAtAero_vs_aero_xAtAero_protons_uncut.Fill(event.P_aero_yAtAero, event.P_aero_xAtAero)
        CTime_epCoinTime_ROC1_vs_beta_protons_uncut.Fill(event.CTime_epCoinTime_ROC1, event.P_gtr_beta)
        #    P_cal_etottracknorm_vs_P_ngcer_npeSum_protons_uncut.Fill(event.P_cal_etottracknorm, event.P_ngcer_npeSum)
        #    P_ngcer_yAtCer_vs_ngcer_xAtCer_protons_uncut.Fill(event.P_ngcer_yAtCer, event.P_ngcer_xAtCer)
        #    P_ngcer_npeSum_vs_hgcer_npeSum_protons_uncut.Fill(event.P_ngcer_npeSum, event.P_hgcer_npeSum)
        #    P_ngcer_npeSum_vs_aero_npeSum_protons_uncut.Fill(event.P_ngcer_npeSum, event.P_aero_npeSum)
        
    CTime_epCoinTime_ROC1_protons_uncut.Fill(event.CTime_epCoinTime_ROC1)

for event in Cut_Proton_Events_All_tree:
    if spec == "HMS":
        H_gtr_beta_protons_cut_all.Fill(event.H_gtr_beta)
        H_gtr_xp_protons_cut_all.Fill(event.H_gtr_xp)
        H_gtr_yp_protons_cut_all.Fill(event.H_gtr_yp)
        H_gtr_dp_protons_cut_all.Fill(event.H_gtr_dp)
        H_hod_goodscinhit_protons_cut_all.Fill(event.H_hod_goodscinhit)
        H_hod_goodstarttime_protons_cut_all.Fill(event.H_hod_goodstarttime)
        H_cal_etotnorm_protons_cut_all.Fill(event.H_cal_etotnorm)
        H_cal_etottracknorm_protons_cut_all.Fill(event.H_cal_etottracknorm)
        H_cer_npeSum_protons_cut_all.Fill(event.H_cer_npeSum)
        H_RFTime_Dist_protons_cut_all.Fill(event.H_RF_Dist)
        H_kin_MMp_protons_cut_all.Fill(event.MMp)
        H_cal_etottracknorm_vs_H_cer_npeSum_protons_cut_all.Fill(event.H_cal_etottracknorm, event.H_cer_npeSum)
        H_kin_MMp_vs_H_RFTime_protons_cut_all.Fill(event.MMp, event.H_RF_Dist)
        CTime_epCoinTime_ROC1_vs_beta_protons_cut_all.Fill(event.CTime_epCoinTime_ROC1, event.H_gtr_beta)
        CTime_epCoinTime_ROC1_vs_H_kin_MMp_protons_cut_all.Fill(event.CTime_epCoinTime_ROC1, event.MMp)        
        
    if spec == "SHMS":        
        P_gtr_beta_protons_cut_all.Fill(event.P_gtr_beta)
        P_gtr_xp_protons_cut_all.Fill(event.P_gtr_xp)
        P_gtr_yp_protons_cut_all.Fill(event.P_gtr_yp)
        P_gtr_dp_protons_cut_all.Fill(event.P_gtr_dp)
        P_gtr_p_protons_cut_all.Fill(event.P_gtr_p)
        P_hod_goodscinhit_protons_cut_all.Fill(event.P_hod_goodscinhit)
        P_hod_goodstarttime_protons_cut_all.Fill(event.P_hod_goodstarttime)
        P_cal_etotnorm_protons_cut_all.Fill(event.P_cal_etotnorm)
        P_cal_etottracknorm_protons_cut_all.Fill(event.P_cal_etottracknorm)
        P_hgcer_npeSum_protons_cut_all.Fill(event.P_hgcer_npeSum)
        P_hgcer_xAtCer_protons_cut_all.Fill(event.P_hgcer_xAtCer)
        P_hgcer_yAtCer_protons_cut_all.Fill(event.P_hgcer_yAtCer)
        P_aero_npeSum_protons_cut_all.Fill(event.P_aero_npeSum)
        P_aero_xAtAero_protons_cut_all.Fill(event.P_aero_xAtAero)
        P_aero_yAtAero_protons_cut_all.Fill(event.P_aero_yAtAero)
        #    P_ngcer_npeSum_protons_cut_all.Fill(event.P_ngcer_npeSum)
        #    P_ngcer_xAtCer_protons_cut_all.Fill(event.P_ngcer_xAtCer)
        #    P_ngcer_yAtCer_protons_cut_all.Fill(event.P_ngcer_yAtCer)   
        P_kin_MMp_protons_cut_all.Fill(event.MMp)
        P_RFTime_Dist_protons_cut_all.Fill(event.P_RF_Dist)
        P_hgcer_npeSum_vs_aero_npeSum_protons_cut_all.Fill(event.P_hgcer_npeSum, event.P_aero_npeSum)
        P_hgcer_yAtCer_vs_hgcer_xAtCer_protons_cut_all.Fill(event.P_hgcer_yAtCer, event.P_hgcer_xAtCer)
        P_aero_yAtAero_vs_aero_xAtAero_protons_cut_all.Fill(event.P_aero_yAtAero, event.P_aero_xAtAero)
        P_kin_MMp_vs_P_RFTime_protons_cut_all.Fill(event.MMp, event.P_RF_Dist)
        CTime_epCoinTime_ROC1_vs_beta_protons_cut_all.Fill(event.CTime_epCoinTime_ROC1, event.P_gtr_beta)
        CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_cut_all.Fill(event.CTime_epCoinTime_ROC1, event.MMp)        
        #    P_cal_etottracknorm_vs_P_ngcer_npeSum_protons_cut_all.Fill(event.P_cal_etottracknorm, event.P_ngcer_npeSum)
        #    P_ngcer_yAtCer_vs_ngcer_xAtCer_protons_cut_all.Fill(event.P_ngcer_yAtCer, event.P_ngcer_xAtCer)
        #    P_ngcer_npeSum_vs_hgcer_npeSum_protons_cut_all.Fill(event.P_ngcer_npeSum, event.P_hgcer_npeSum)
        #    P_ngcer_npeSum_vs_aero_npeSum_protons_cut_all.Fill(event.P_ngcer_npeSum, event.P_aero_npeSum) 
        
    CTime_epCoinTime_ROC1_protons_cut_all.Fill(event.CTime_epCoinTime_ROC1)

for event in Cut_Proton_Events_Prompt_tree:
    if spec == "HMS":
        H_gtr_beta_protons_cut_prompt.Fill(event.H_gtr_beta)
        H_RFTime_Dist_protons_cut_prompt.Fill(event.H_RF_Dist)
        H_kin_MMp_protons_cut_prompt.Fill(event.MMp)
        H_kin_MMp_vs_CTime_epCoinTime_ROC1_protons_cut_prompt.Fill(event.MMp, event.CTime_epCoinTime_ROC1)

    if spec == "SHMS":
        P_gtr_beta_protons_cut_prompt.Fill(event.P_gtr_beta)
        P_RFTime_Dist_protons_cut_prompt.Fill(event.P_RF_Dist)
        P_kin_MMp_protons_cut_prompt.Fill(event.MMp)
        P_kin_MMp_vs_CTime_epCoinTime_ROC1_protons_cut_prompt.Fill(event.MMp, event.CTime_epCoinTime_ROC1)
        
    CTime_epCoinTime_ROC1_protons_cut_prompt.Fill(event.CTime_epCoinTime_ROC1)

for event in Cut_Proton_Events_Random_tree:
    if spec == "HMS":
        H_gtr_beta_protons_cut_random.Fill(event.H_gtr_beta)
        H_RFTime_Dist_protons_cut_random.Fill(event.H_RF_Dist)
        H_kin_MMp_protons_cut_random.Fill(event.MMp)

    if spec == "SHMS":
        P_gtr_beta_protons_cut_random.Fill(event.P_gtr_beta)
        P_RFTime_Dist_protons_cut_random.Fill(event.P_RF_Dist)
        P_kin_MMp_protons_cut_random.Fill(event.MMp)

    CTime_epCoinTime_ROC1_protons_cut_random.Fill(event.CTime_epCoinTime_ROC1)
    
#################################################################################################################################################

# Random subtraction from missing mass
for event in Cut_Proton_Events_Random_tree:
    if spec == "HMS":
        H_kin_MMp_protons_cut_random_scaled.Fill(event.MMp)
        H_kin_MMp_protons_cut_random_scaled.Scale(1.0/nWindows)
        H_kin_MMp_protons_cut_random_sub.Add(H_kin_MMp_protons_cut_prompt, H_kin_MMp_protons_cut_random_scaled, 1, -1)

    if spec == "SHMS":
        P_kin_MMp_protons_cut_random_scaled.Fill(event.MMp)
        P_kin_MMp_protons_cut_random_scaled.Scale(1.0/nWindows)
        P_kin_MMp_protons_cut_random_sub.Add(P_kin_MMp_protons_cut_prompt, P_kin_MMp_protons_cut_random_scaled, 1, -1)        
        
############################################################################################################################################

if spec == "HMS":
    # Saving histograms in PDF
    c1_mom = TCanvas("c1_mom", "Momentum Distributions", 100, 0, 1000, 900)
    c1_mom.Divide(1,1)
    # End of polar plotting section
    c1_mom.cd(1)
    H_gtr_p_protons_uncut.Draw("hist")
    H_gtr_p_protons_cut_all.Draw("hist")
    c1_mom.Print(Proton_Analysis_Distributions + '(')

    c1_acpt = TCanvas("c1_H_kin", "Electron-Proton Acceptance Distributions", 100, 0, 1000, 900)
    c1_acpt.Divide(3,1)
    c1_acpt.cd(1)
    gPad.SetLogy()
    H_gtr_xp_protons_uncut.SetLineColor(2)
    H_gtr_xp_protons_uncut.Draw()
    H_gtr_xp_protons_cut_all.SetLineColor(4)
    H_gtr_xp_protons_cut_all.Draw("same")
    c1_acpt.cd(2)
    gPad.SetLogy()
    H_gtr_yp_protons_uncut.SetLineColor(2)
    H_gtr_yp_protons_uncut.Draw()
    H_gtr_yp_protons_cut_all.SetLineColor(4)
    H_gtr_yp_protons_cut_all.Draw("same")
    c1_acpt.cd(3)
    gPad.SetLogy()
    H_gtr_dp_protons_uncut.SetLineColor(2)
    H_gtr_dp_protons_uncut.Draw()
    H_gtr_dp_protons_cut_all.SetLineColor(4)
    H_gtr_dp_protons_cut_all.Draw("same")
    # TLegend (x1, y1, x2, y2) 
    legend2 = ROOT.TLegend(0.115, 0.8, 0.6, 0.9)
    legend2.AddEntry("H_gtr_dp_protons_uncut", "without cuts", "l")
    legend2.AddEntry("H_gtr_dp_protons_cut_all", "with cuts (acpt/RF/PID)", "l")
    legend2.Draw("same")
    c1_acpt.Print(Proton_Analysis_Distributions)
    
    c1_pid = TCanvas("c1_pid", "Electron-Proton CAL Distributions", 100, 0, 1000, 900)
    c1_pid.Divide(2,1)
    c1_pid.cd(1)
    gPad.SetLogy()
    H_cal_etottracknorm_protons_uncut.SetLineColor(2)
    H_cal_etottracknorm_protons_uncut.Draw()
    H_cal_etottracknorm_protons_cut_all.SetLineColor(4)
    H_cal_etottracknorm_protons_cut_all.Draw("same")
    legend7 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
    legend7.AddEntry("H_cal_etottracknorm_protons_uncut", "without cuts", "l")
    legend7.AddEntry("H_cal_etottracknorm_protons_cut_all", "with cuts (acpt/RF/PID)", "l")
    legend7.Draw("same")
    c1_pid.cd(2)
    H_cal_etottracknorm_vs_H_cer_npeSum_protons_cut_all.Draw("COLZ")
    c1_pid.Print(Proton_Analysis_Distributions + ')')
    
if spec == "SHMS":
    # Saving histograms in PDF
    c1_mom = TCanvas("c1_mom", "Momentum Distributions", 100, 0, 1000, 900)
    c1_mom.Divide(1,1)
    # End of polar plotting section
    c1_mom.cd(1)
    P_gtr_p_protons_uncut.Draw("hist")
    P_gtr_p_protons_cut_all.Draw("hist")
    c1_mom.Print(Proton_Analysis_Distributions + '(')

    c1_acpt = TCanvas("c1_H_kin", "Electron-Proton Acceptance Distributions", 100, 0, 1000, 900)
    c1_acpt.Divide(3,1)
    c1_acpt.cd(1)
    gPad.SetLogy()
    P_gtr_xp_protons_uncut.SetLineColor(2)
    P_gtr_xp_protons_uncut.Draw()
    P_gtr_xp_protons_cut_all.SetLineColor(4)
    P_gtr_xp_protons_cut_all.Draw("same")
    c1_acpt.cd(2)
    gPad.SetLogy()
    P_gtr_yp_protons_uncut.SetLineColor(2)
    P_gtr_yp_protons_uncut.Draw()
    P_gtr_yp_protons_cut_all.SetLineColor(4)
    P_gtr_yp_protons_cut_all.Draw("same")
    c1_acpt.cd(3)
    gPad.SetLogy()
    P_gtr_dp_protons_uncut.SetLineColor(2)
    P_gtr_dp_protons_uncut.Draw()
    P_gtr_dp_protons_cut_all.SetLineColor(4)
    P_gtr_dp_protons_cut_all.Draw("same")
    c1_acpt.Print(Proton_Analysis_Distributions)
    
    c1_pid = TCanvas("c1_pid", "Electron-Proton CAL Distributions", 100, 0, 1000, 900)
    c1_pid.Divide(1,1)
    c1_pid.cd(1)
    gPad.SetLogy()
    P_cal_etottracknorm_protons_uncut.SetLineColor(2)
    P_cal_etottracknorm_protons_uncut.Draw()
    P_cal_etottracknorm_protons_cut_all.SetLineColor(4)
    P_cal_etottracknorm_protons_cut_all.Draw("same")
    legend8 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
    legend8.AddEntry("P_cal_etottracknorm_protons_uncut", "without cuts", "l")
    legend8.AddEntry("P_cal_etottracknorm_protons_cut_all", "with cuts (acpt/RF/PID)", "l")
    legend8.Draw("same")
    c1_pid.Print(Proton_Analysis_Distributions)
    
    c2_pid = TCanvas("c2_pid", "Electron-Proton Aero/HGC/NGC PID Distributions", 100, 0, 1000, 900)
    c2_pid.Divide(3,1)
    c2_pid.cd(1)
    gPad.SetLogy()
    P_hgcer_npeSum_protons_uncut.SetLineColor(2)
    P_hgcer_npeSum_protons_uncut.Draw()
    P_hgcer_npeSum_protons_cut_all.SetLineColor(4)
    P_hgcer_npeSum_protons_cut_all.Draw("same")
    legend10 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
    legend10.AddEntry("P_hgcer_npeSum_protons_uncut", "without cuts", "l")
    legend10.AddEntry("P_hgcer_npeSum_protons_cut_all", "with cuts (acpt/RF/PID)", "l")
    legend10.Draw("same")
    c2_pid.cd(2)
    gPad.SetLogy()
    P_aero_npeSum_protons_uncut.SetLineColor(2)
    P_aero_npeSum_protons_uncut.Draw()
    P_aero_npeSum_protons_cut_all.SetLineColor(4)
    P_aero_npeSum_protons_cut_all.Draw("same")
    legend11 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
    legend11.AddEntry("P_aero_npeSum_protons_uncut", "without cuts", "l")
    legend11.AddEntry("P_aero_npeSum_protons_cut_all", "with cuts (acpt/RF/PID)", "l")
    legend11.Draw("same")
    c2_pid.cd(3)
    #P_ngcer_npeSum_protons_uncut.SetLineColor(2)
    #P_ngcer_npeSum_protons_uncut.Draw()
    #P_ngcer_npeSum_protons_cut_all.SetLineColor(4)
    #P_ngcer_npeSum_protons_cut_all.Draw("same") 
    #legend12 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
    #legend12.AddEntry("P_ngcer_npeSum_protons_uncut", "without cuts", "l")
    #legend12.AddEntry("P_ngcer_npeSum_protons_cut_all", "with cuts (acpt/RF/PID)", "l")
    #legend12.Draw("same")
    c2_pid.Print(Proton_Analysis_Distributions)
    
    c3_pid = TCanvas("c3_pid", "Electron-Proton Aero/HGC/NGC PID Distributions", 100, 0, 1000, 900)
    c3_pid.Divide(2,3)
    c3_pid.cd(1)
    gPad.SetLogz()
    P_hgcer_npeSum_vs_aero_npeSum_protons_uncut.Draw("COLZ")
    c3_pid.cd(2)
    P_hgcer_npeSum_vs_aero_npeSum_protons_cut_all.Draw("COLZ")
    c3_pid.cd(3)
    #P_ngcer_npeSum_vs_hgcer_npeSum_protons_uncut.Draw("COLZ")
    c3_pid.cd(4)
    #P_ngcer_npeSum_vs_hgcer_npeSum_protons_cut_all.Draw("COLZ")
    c3_pid.cd(5)
    #P_ngcer_npeSum_vs_aero_npeSum_protons_uncut.Draw("COLZ")
    c3_pid.cd(6)
    #P_ngcer_npeSum_vs_aero_npeSum_protons_cut_all.Draw("COLZ")
    c3_pid.Print(Proton_Analysis_Distributions + ')')

#############################################################################################################################################

# Making directories in output file
outHistFile = ROOT.TFile.Open("%s/%s_%s_%s_Output_Data.root" % (OUTPATH, spec, runNum, MaxEvent), "RECREATE")                                                                                                    
d_Uncut_Proton_Events = outHistFile.mkdir("Uncut_Proton_Events")
d_Cut_Proton_Events_All = outHistFile.mkdir("Cut_Proton_Events_All")
d_Cut_Proton_Events_Prompt = outHistFile.mkdir("Cut_Proton_Events_Prompt")
d_Cut_Proton_Events_Random = outHistFile.mkdir("Cut_Proton_Events_Random")

# Writing Histograms for protons                                                                  
d_Uncut_Proton_Events.cd()
if spec == "HMS":
    H_gtr_beta_protons_uncut.Write()
    H_gtr_xp_protons_uncut.Write()
    H_gtr_yp_protons_uncut.Write()
    H_gtr_dp_protons_uncut.Write()
    H_hod_goodscinhit_protons_uncut.Write()
    H_hod_goodstarttime_protons_uncut.Write()
    H_cal_etotnorm_protons_uncut.Write()
    H_cal_etottracknorm_protons_uncut.Write()
    H_cer_npeSum_protons_uncut.Write()
    H_kin_MMp_protons_uncut.Write()
    H_RFTime_Dist_protons_uncut.Write()
    H_cal_etottracknorm_vs_H_cer_npeSum_protons_uncut.Write()
    CTime_epCoinTime_ROC1_vs_H_kin_MMp_protons_uncut.Write()
    H_kin_MMp_vs_H_RFTime_protons_uncut.Write()

if spec == "SHMS":    
    P_gtr_beta_protons_uncut.Write()
    P_gtr_xp_protons_uncut.Write()
    P_gtr_yp_protons_uncut.Write()
    P_gtr_dp_protons_uncut.Write()
    P_gtr_p_protons_uncut.Write()
    P_hod_goodscinhit_protons_uncut.Write()
    P_hod_goodstarttime_protons_uncut.Write()
    P_cal_etotnorm_protons_uncut.Write()
    P_cal_etottracknorm_protons_uncut.Write()
    P_hgcer_npeSum_protons_uncut.Write()
    P_hgcer_xAtCer_protons_uncut.Write()
    P_hgcer_yAtCer_protons_uncut.Write()
    P_aero_npeSum_protons_uncut.Write()
    P_aero_xAtAero_protons_uncut.Write()
    P_aero_yAtAero_protons_uncut.Write()
    #P_ngcer_npeSum_protons_uncut.Write()
    #P_ngcer_xAtCer_protons_uncut.Write()
    #P_ngcer_yAtCer_protons_uncut.Write() 
    P_kin_MMp_protons_uncut.Write()
    P_RFTime_Dist_protons_uncut.Write()
    P_hgcer_npeSum_vs_aero_npeSum_protons_uncut.Write()
    CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_uncut.Write()
    P_hgcer_yAtCer_vs_hgcer_xAtCer_protons_uncut.Write()
    P_aero_yAtAero_vs_aero_xAtAero_protons_uncut.Write()
    P_kin_MMp_vs_P_RFTime_protons_uncut.Write()
    #P_cal_etottracknorm_vs_P_ngcer_npeSum_protons_uncut.Write()
    #P_ngcer_yAtCer_vs_ngcer_xAtCer_protons_uncut.Write()
    #P_ngcer_npeSum_vs_hgcer_npeSum_protons_uncut.Write()
    #P_ngcer_npeSum_vs_aero_npeSum_protons_uncut.Write()

CTime_epCoinTime_ROC1_protons_uncut.Write()
CTime_epCoinTime_ROC1_vs_beta_protons_uncut.Write()


d_Cut_Proton_Events_All.cd()
if spec == "HMS":
    H_gtr_beta_protons_cut_all.Write()
    H_gtr_xp_protons_cut_all.Write()
    H_gtr_yp_protons_cut_all.Write()
    H_gtr_dp_protons_cut_all.Write()
    H_hod_goodscinhit_protons_cut_all.Write()
    H_hod_goodstarttime_protons_cut_all.Write()
    H_cal_etotnorm_protons_cut_all.Write()
    H_cal_etottracknorm_protons_cut_all.Write()
    H_cer_npeSum_protons_cut_all.Write()
    H_RFTime_Dist_protons_cut_all.Write()
    H_kin_MMp_protons_cut_all.Write()
    H_cal_etottracknorm_vs_H_cer_npeSum_protons_cut_all.Write()
    CTime_epCoinTime_ROC1_vs_H_kin_MMp_protons_cut_all.Write()
    H_kin_MMp_vs_H_RFTime_protons_cut_all.Write()
    H_kin_MMp_protons_cut_random_sub.Write()

if spec == "SHMS":    
    P_gtr_beta_protons_cut_all.Write()
    P_gtr_xp_protons_cut_all.Write()
    P_gtr_yp_protons_cut_all.Write()
    P_gtr_dp_protons_cut_all.Write()
    P_gtr_p_protons_cut_all.Write()
    P_hod_goodscinhit_protons_cut_all.Write()
    P_hod_goodstarttime_protons_cut_all.Write()
    P_cal_etotnorm_protons_cut_all.Write()
    P_cal_etottracknorm_protons_cut_all.Write()
    P_hgcer_npeSum_protons_cut_all.Write()
    P_hgcer_xAtCer_protons_cut_all.Write()
    P_hgcer_yAtCer_protons_cut_all.Write()
    P_aero_npeSum_protons_cut_all.Write()
    P_aero_xAtAero_protons_cut_all.Write()
    P_aero_yAtAero_protons_cut_all.Write()
    #P_ngcer_npeSum_protons_cut_all.Write()
    #P_ngcer_xAtCer_protons_cut_all.Write()
    #P_ngcer_yAtCer_protons_cut_all.Write()
    P_kin_MMp_protons_cut_all.Write()
    P_RFTime_Dist_protons_cut_all.Write()
    P_hgcer_npeSum_vs_aero_npeSum_protons_cut_all.Write()
    CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_cut_all.Write()
    P_hgcer_yAtCer_vs_hgcer_xAtCer_protons_cut_all.Write()
    P_aero_yAtAero_vs_aero_xAtAero_protons_cut_all.Write()
    P_kin_MMp_vs_P_RFTime_protons_cut_all.Write()
    P_kin_MMp_protons_cut_random_sub.Write()
    #P_cal_etottracknorm_vs_P_ngcer_npeSum_protons_cut_all.Write()
    #P_ngcer_yAtCer_vs_ngcer_xAtCer_protons_cut_all.Write()
    #P_ngcer_npeSum_vs_hgcer_npeSum_protons_cut_all.Write()
    #P_ngcer_npeSum_vs_aero_npeSum_protons_cut_all.Write()
    
CTime_epCoinTime_ROC1_protons_cut_all.Write()
CTime_epCoinTime_ROC1_vs_beta_protons_cut_all.Write()


d_Cut_Proton_Events_Prompt.cd()
if spec == "HMS":
    H_gtr_beta_protons_cut_prompt.Write()
    H_RFTime_Dist_protons_cut_prompt.Write()
    CTime_epCoinTime_ROC1_protons_cut_prompt.Write()
    H_kin_MMp_protons_cut_prompt.Write()
    H_kin_MMp_vs_CTime_epCoinTime_ROC1_protons_cut_prompt.Write()
    
if spec == "SHMS":
    P_gtr_beta_protons_cut_prompt.Write()
    P_RFTime_Dist_protons_cut_prompt.Write()
    CTime_epCoinTime_ROC1_protons_cut_prompt.Write()
    P_kin_MMp_protons_cut_prompt.Write()
    P_kin_MMp_vs_CTime_epCoinTime_ROC1_protons_cut_prompt.Write()

    
d_Cut_Proton_Events_Random.cd()
if spec == "HMS":
    H_gtr_beta_protons_cut_random.Write()
    H_RFTime_Dist_protons_cut_random.Write()
    CTime_epCoinTime_ROC1_protons_cut_random.Write()
    H_kin_MMp_protons_cut_random.Write()
    
if spec == "SHMS":
    P_gtr_beta_protons_cut_random.Write()
    P_RFTime_Dist_protons_cut_random.Write()
    CTime_epCoinTime_ROC1_protons_cut_random.Write()
    P_kin_MMp_protons_cut_random.Write() 
    
outHistFile.Close()
infile.Close() 
print ("Processing Complete")
