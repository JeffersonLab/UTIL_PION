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
Analysis_Distributions = "%s/%s_%s_sw_heep_%s_Analysis_Distributions.pdf" % (OUTPATH, runNum, MaxEvent, spec)

#################################################################################################################################################

# Construct the name of the rootfile based upon the info we provided
rootName = "%s/UTIL_PION/OUTPUT/Analysis/HeeP/%s_%s_%s_%s.root" % (REPLAYPATH, spec, runNum, MaxEvent, ROOTPrefix)     # Input file location and variables taking
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

Uncut_Events_tree = infile.Get("Uncut_Events")
Cut_Events_All_tree = infile.Get("Cut_Events_All")

###################################################################################################################################################

 # Defining Histograms for Protons
if spec == "HMS":
    H_gtr_beta_uncut = ROOT.TH1D("H_gtr_beta_uncut", "HMS #beta; HMS_gtr_#beta; Counts", 200, 0.8, 1.2)
    H_gtr_xp_uncut = ROOT.TH1D("H_gtr_xp_uncut", "HMS x'; HMS_gtr_xp; Counts", 200, -0.2, 0.2)
    H_gtr_yp_uncut = ROOT.TH1D("H_gtr_yp_uncut", "HMS y'; HMS_gtr_yp; Counts", 200, -0.2, 0.2)
    H_gtr_dp_uncut = ROOT.TH1D("H_gtr_dp_uncut", "HMS #delta; HMS_gtr_dp; Counts", 200, -15, 15)
    H_gtr_p_uncut = ROOT.TH1D("H_gtr_p_uncut", "HMS p; HMS_gtr_p; Counts", 200, 4, 8)
    H_hod_goodscinhit_uncut = ROOT.TH1D("H_hod_goodscinhit_uncut", "HMS hod goodscinhit; HMS_hod_goodscinhi; Counts", 200, 0.7, 1.3)
    H_hod_goodstarttime_uncut = ROOT.TH1D("H_hod_goodstarttime_uncut", "HMS hod goodstarttime; HMS_hod_goodstarttime; Counts", 200, 0.7, 1.3)
    H_cal_etotnorm_uncut = ROOT.TH1D("H_cal_etotnorm_uncut", "HMS cal etotnorm; HMS_cal_etotnorm; Counts", 200, 0.2, 1.8)
    H_cal_etottracknorm_uncut = ROOT.TH1D("H_cal_etottracknorm_uncut", "HMS cal etottracknorm; HMS_cal_etottracknorm; Counts", 300, 0.2, 1.8)
    H_cer_npeSum_uncut = ROOT.TH1D("H_cer_npeSum_uncut", "HMS cer npeSum; HMS_cer_npeSum; Counts", 200, 0, 50)
    H_RFTime_Dist_uncut = ROOT.TH1D("H_RFTime_Dist_uncut", "HMS RFTime; HMS_RFTime; Counts", 200, 0, 4)
    
    H_gtr_beta_cut_all = ROOT.TH1D("H_gtr_beta_cut_all", "HMS #beta; HMS_gtr_#beta; Counts", 200, 0.8, 1.2)
    H_gtr_xp_cut_all = ROOT.TH1D("H_gtr_xp_cut_all", "HMS x'; HMS_gtr_xp; Counts", 200, -0.2, 0.2)
    H_gtr_yp_cut_all = ROOT.TH1D("H_gtr_yp_cut_all", "HMS y'; HMS_gtr_yp; Counts", 200, -0.2, 0.2)
    H_gtr_dp_cut_all = ROOT.TH1D("H_gtr_dp_cut_all", "HMS #delta; HMS_gtr_dp; Counts", 200, -15, 15)
    H_gtr_p_cut_all = ROOT.TH1D("H_gtr_p_cut_all", "HMS p; HMS_gtr_p; Counts", 200, 4, 8)
    H_hod_goodscinhit_cut_all = ROOT.TH1D("H_hod_goodscinhit_cut_all", "HMS hod goodscinhit; HMS_hod_goodscinhit; Counts", 200, 0.7, 1.3)
    H_hod_goodstarttime_cut_all = ROOT.TH1D("H_hod_goodstarttime_cut_all", "HMS hod goodstarttime; HMS_hod_goodstarttime; Counts", 200, 0.7, 1.3)
    H_cal_etotnorm_cut_all = ROOT.TH1D("H_cal_etotnorm_cut_all", "HMS cal etotnorm; HMS_cal_etotnorm; Counts", 200, 0.6, 1.4)
    H_cal_etottracknorm_cut_all = ROOT.TH1D("H_cal_etottracknorm_cut_all", "HMS cal etottracknorm; HMS_cal_etottracknorm; Counts", 200, 0.6, 1.4)
    H_cer_npeSum_cut_all = ROOT.TH1D("H_cer_npeSum_cut_all", "HMS cer npeSum; HMS_cer_npeSum; Counts", 200, 0, 50)
    H_RFTime_Dist_cut_all = ROOT.TH1D("H_RFTime_Dist_cut_all", "HMS RFTime; HMS_RFTime; Counts", 200, 0, 4)
    
if spec == "SHMS":
    P_gtr_beta_uncut = ROOT.TH1D("P_gtr_beta_uncut", "SHMS #beta; SHMS_gtr_#beta; Counts", 200, 0.7, 1.3)
    P_gtr_xp_uncut = ROOT.TH1D("P_gtr_xp_uncut", "SHMS x'; SHMS_gtr_xp; Counts", 200, -0.2, 0.2)
    P_gtr_yp_uncut = ROOT.TH1D("P_gtr_yp_uncut", "SHMS y'; SHMS_gtr_yp; Counts", 200, -0.2, 0.2)
    P_gtr_dp_uncut = ROOT.TH1D("P_gtr_dp_uncut", "SHMS delta; SHMS_gtr_dp; Counts", 200, -30, 30)
    P_gtr_p_uncut = ROOT.TH1D("P_gtr_p_uncut", "SHMS p; SHMS_gtr_p; Counts", 200, 4, 8)
    P_hod_goodscinhit_uncut = ROOT.TH1D("P_hod_goodscinhit_uncut", "SHMS hod goodscinhit; SHMS_hod_goodscinhit; Counts", 200, 0.7, 1.3)
    P_hod_goodstarttime_uncut = ROOT.TH1D("P_hod_goodstarttime_uncut", "SHMS hod goodstarttime; SHMS_hod_goodstarttime; Counts", 200, 0.7, 1.3)
    P_cal_etotnorm_uncut = ROOT.TH1D("P_cal_etotnorm_uncut", "SHMS cal etotnorm; SHMS_cal_etotnorm; Counts", 200, 0, 1)
    P_cal_etottracknorm_uncut = ROOT.TH1D("P_cal_etottracknorm_uncut", "SHMS cal etottracknorm; SHMS_cal_etottracknorm; Counts", 200, 0, 1.6)
    P_hgcer_npeSum_uncut = ROOT.TH1D("P_hgcer_npeSum_uncut", "SHMS HGC npeSum; SHMS_hgcer_npeSum; Counts", 200, 0, 50)
    P_hgcer_xAtCer_uncut = ROOT.TH1D("P_hgcer_xAtCer_uncut", "SHMS HGC xAtCer; SHMS_hgcer_xAtCer; Counts", 200, -60, 60)
    P_hgcer_yAtCer_uncut = ROOT.TH1D("P_hgcer_yAtCer_uncut", "SHMS HGC yAtCer; SHMS_hgcer_yAtCer; Counts", 200, -50, 50)
    P_aero_npeSum_uncut = ROOT.TH1D("P_aero_npeSum_uncut", "SHMS aero npeSum; SHMS_aero_npeSum; Counts", 200, 0, 50)
    P_aero_xAtAero_uncut = ROOT.TH1D("P_acero_xAtAero_uncut", "SHMS aero xAtAero; SHMS_aero_xAtAero; Counts", 200, -60, 60)
    P_aero_yAtAero_uncut = ROOT.TH1D("P_aero_yAtAero_uncut", "SHMS aero yAtAero; SHMS_aero_yAtAero; Counts", 200, -50, 50)
    P_ngcer_npeSum_uncut = ROOT.TH1D("P_ngcer_npeSum_uncut", "SHMS NGC npeSum; SHMS_ngcer_npeSum; Counts", 200, 0, 50)
    P_ngcer_xAtCer_uncut = ROOT.TH1D("P_ngcer_xAtCer_uncut", "SHMS NGC xAtCer; SHMS_ngcer_xAtCer; Counts", 200, -70, 50)
    P_ngcer_yAtCer_uncut = ROOT.TH1D("P_ngcer_yAtCer_uncut", "SHMS NGC yAtCer; SHMS_ngcer_yAtCer; Counts", 200, -50, 50)
    P_RFTime_Dist_uncut = ROOT.TH1D("P_RFTime_Dist_uncut", "SHMS RFTime; SHMS_RFTime; Counts", 200, 0, 4)
    
    P_gtr_beta_cut_all = ROOT.TH1D("P_gtr_beta_cut_all", "SHMS #beta; SHMS_gtr_#beta; Counts", 200, 0.8, 1.2)
    P_gtr_xp_cut_all = ROOT.TH1D("P_gtr_xp_cut_all", "SHMS x'; SHMS_gtr_xp; Counts", 200, -0.2, 0.2)
    P_gtr_yp_cut_all = ROOT.TH1D("P_gtr_yp_cut_all", "SHMS y'; SHMS_gtr_yp; Counts", 200, -0.2, 0.2)
    P_gtr_dp_cut_all = ROOT.TH1D("P_gtr_dp_cut_all", "SHMS #delta; SHMS_gtr_dp; Counts", 200, -15, 15)
    P_gtr_p_cut_all = ROOT.TH1D("P_gtr_p_cut_all", "SHMS p; SHMS_gtr_p; Counts", 200, 4, 8)
    P_hod_goodscinhit_cut_all = ROOT.TH1D("P_hod_goodscinhit_cut_all", "SHMS hod goodscinhit; SHMS_hod_goodscinhit; Counts", 200, 0.7, 1.3)
    P_hod_goodstarttime_cut_all = ROOT.TH1D("P_hod_goodstarttime_cut_all", "SHMS hod goodstarttime; SHMS_hod_goodstarttime; Counts", 200, 0.7, 1.3)
    P_cal_etotnorm_cut_all = ROOT.TH1D("P_cal_etotnorm_cut_all", "SHMS cal etotnorm; SHMS_cal_etotnorm; Counts", 200, 0, 1.2)
    P_cal_etottracknorm_cut_all = ROOT.TH1D("P_cal_etottracknorm_cut_all", "SHMS cal etottracknorm; SHMS_cal_etottracknorm; Counts", 200, 0, 1.6)
    P_hgcer_npeSum_cut_all = ROOT.TH1D("P_hgcer_npeSum_cut_all", "SHMS HGC npeSum; SHMS_hgcer_npeSum; Counts", 200, 0, 50)
    P_hgcer_xAtCer_cut_all = ROOT.TH1D("P_hgcer_xAtCer_cut_all", "SHMS HGC xAtCer; SHMS_hgcer_xAtCer; Counts", 200, -40, 30)
    P_hgcer_yAtCer_cut_all = ROOT.TH1D("P_hgcer_yAtCer_cut_all", "SHMS HGC yAtCer; SHMS_hgcer_yAtCer; Counts", 200, -30, 30)
    P_aero_npeSum_cut_all = ROOT.TH1D("P_aero_npeSum_cut_all", "SHMS aero npeSum; SHMS_aero_npeSum; Counts", 200, 0, 50)
    P_aero_xAtAero_cut_all = ROOT.TH1D("P_acero_xAtAero_cut_all", "SHMS aero xAtAero; SHMS_aero_xAtAero; Counts", 200, -40, 30)
    P_aero_yAtAero_cut_all = ROOT.TH1D("P_aero_yAtAero_cut_all", "SHMS aero yAtAero; SHMS_aero_yAtAero; Counts", 200, -30, 30)
    P_ngcer_npeSum_cut_all = ROOT.TH1D("P_ngcer_npeSum_cut_all", "SHMS NGC npeSum; SHMS_ngcer_npeSum; Counts", 200, -10, 50)
    P_ngcer_xAtCer_cut_all = ROOT.TH1D("P_ngcer_xAtCer_cut_all", "SHMS NGC xAtCer; SHMS_ngcer_xAtCer; Counts", 200, -40, 30)
    P_ngcer_yAtCer_cut_all = ROOT.TH1D("P_ngcer_yAtCer_cut_all", "SHMS NGC yAtCer; SHMS_ngcer_yAtCer; Counts", 200, -30, 30)
    P_RFTime_Dist_cut_all = ROOT.TH1D("P_RFTime_Dist_cut_all", "SHMS RFTime; SHMS_RFTime; Counts", 200, 0, 4)

###################################################################################################################################################

# 2D Histograms for protons
if spec == "HMS":
    H_cal_etottracknorm_vs_H_cer_npeSum_uncut = ROOT.TH2D("H_cal_etottracknorm_vs_H_cer_npeSum_uncut","HMS cal etottracknorm vs HMS cer npeSum (no cut); H_cal_etottracknorm; H_cer_npeSum",100, 0, 2, 100, 0, 40)
    H_cal_etottracknorm_vs_H_cer_npeSum_cut_all = ROOT.TH2D("H_cal_etottracknorm_vs_H_cer_npeSum_cut_all","HMS cal etottracknorm vs HMS cer npeSum (with cuts); H_cal_etottracknorm; H_cer_npeSum",100, 0.5, 1.5, 100, 0, 40)
    
if spec == "SHMS":
    P_hgcer_npeSum_vs_aero_npeSum_uncut = ROOT.TH2D("P_hgcer_npeSum_vs_aero_npeSum_uncut", "SHMS HGC npeSum vs SHMS Aero npeSum (no cut); SHMS_hgcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
    P_hgcer_yAtCer_vs_hgcer_xAtCer_uncut = ROOT.TH2D("P_hgcer_yAtCer_vs_hgcer_xAtCer_uncut", "SHMS HGC yAtCer vs SHMS HGC xAtCer (no cut); SHMS_hgcer_yAtCer; SHMS_hgcer_xAtCer", 100, -50, 50, 100, -50, 50)
    P_aero_yAtAero_vs_aero_xAtAero_uncut = ROOT.TH2D("P_aero_yAtAero_vs_aero_xAtAero_uncut", "SHMS aero yAtAero vs SHMS aero xAtAero (no cut); SHMS_aero_yAtAero; SHMS_aero_xAtAero", 100, -50, 50, 100, -50, 50)
    P_cal_etottracknorm_vs_P_ngcer_npeSum_uncut = ROOT.TH2D("P_cal_etottracknorm_vs_P_ngcer_npeSum_uncut", "SHMS cal etottracknorm vs SHMS NGC xAtCer (no cut); SHMS_cal_etottracknorm; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
    P_ngcer_yAtCer_vs_ngcer_xAtCer_uncut = ROOT.TH2D("P_ngcer_yAtCer_vs_ngcer_xAtCer_uncut", "SHMS NGC yAtCer vs SHMS NGC xAtCer (no cut); SHMS_ngcer_yAtCer; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
    P_ngcer_npeSum_vs_hgcer_npeSum_uncut = ROOT.TH2D("P_ngcer_npeSum_vs_hgcer_npeSum_uncut", "SHMS NGC npeSum vs SHMS HGC npeSum (no cut); SHMS_ngcer_npeSum; SHMS_hgcer_npeSum", 100, 0, 50, 100, 0, 50)
    P_ngcer_npeSum_vs_aero_npeSum_uncut = ROOT.TH2D("P_ngcer_npeSum_vs_aero_npeSum_uncut", "SHMS NGC npeSum vs SHMS aero npeSum (no cut); SHMS_ngcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)

    P_hgcer_npeSum_vs_aero_npeSum_cut_all = ROOT.TH2D("P_hgcer_npeSum_vs_aero_npeSum_cut_all", "SHMS HGC npeSum vs SHMS aero npeSum (with cuts); SHMS_hgcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
    P_hgcer_yAtCer_vs_hgcer_xAtCer_cut_all = ROOT.TH2D("P_hgcer_yAtCer_vs_hgcer_xAtCer_cut_all", "SHMS HGC yAtCer vs SHMS HGC xAtCer (with cuts); SHMS_hgcer_yAtCer; SHMS_hgcer_xAtCer", 100, -50, 50, 100, -50, 50)
    P_aero_yAtAero_vs_aero_xAtAero_cut_all = ROOT.TH2D("P_aero_yAtAero_vs_aero_xAtAero_cut_all", "SHMS aero yAtAero vs SHMS aero xAtAero (with cuts); SHMS_aero_yAtAero; SHMS_aero_xAtAero", 100, -50, 50, 100, -50, 50)
    P_cal_etottracknorm_vs_P_ngcer_npeSum_cut_all = ROOT.TH2D("P_cal_etottracknorm_vs_P_ngcer_npeSum_cut_all", "P cal etottracknorm vs SHMS NGC xAtCer (with cuts); SHMS_cal_etottracknorm; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
    P_ngcer_yAtCer_vs_ngcer_xAtCer_cut_all = ROOT.TH2D("P_ngcer_yAtCer_vs_ngcer_xAtCer_cut_all", "SHMS NGC yAtCer vs SHMS NGC xAtCer (with cuts); SHMS_ngcer_yAtCer; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
    P_ngcer_npeSum_vs_hgcer_npeSum_cut_all = ROOT.TH2D("P_ngcer_npeSum_vs_hgcer_npeSum_cut_all", "SHMS NGC npeSum vs SHMS HGC npeSum (with cuts); SHMS_ngcer_npeSum; SHMS_hgcer_npeSum", 100, 0, 50, 100, 0, 50)
    P_ngcer_npeSum_vs_aero_npeSum_cut_all = ROOT.TH2D("P_ngcer_npeSum_vs_aero_npeSum_cut_all", "SHMS NGC npeSum vs SHMS aero npeSum (with cuts); SHMS_ngcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)

#################################################################################################################################################

# Filling Histograms for Protons
for event in Uncut_Events_tree:
    if spec == "HMS":
        H_gtr_beta_uncut.Fill(event.H_gtr_beta)
        H_gtr_xp_uncut.Fill(event.H_gtr_xp)
        H_gtr_yp_uncut.Fill(event.H_gtr_yp)
        H_gtr_dp_uncut.Fill(event.H_gtr_dp)
        H_gtr_p_uncut.Fill(event.H_gtr_p)
        H_hod_goodscinhit_uncut.Fill(event.H_hod_goodscinhit)
        H_hod_goodstarttime_uncut.Fill(event.H_hod_goodstarttime)
        H_cal_etotnorm_uncut.Fill(event.H_cal_etotnorm)
        H_cal_etottracknorm_uncut.Fill(event.H_cal_etottracknorm)
        H_cer_npeSum_uncut.Fill(event.H_cer_npeSum)
        H_RFTime_Dist_uncut.Fill(event.H_RF_Dist)
        H_cal_etottracknorm_vs_H_cer_npeSum_uncut.Fill(event.H_cal_etottracknorm, event.H_cer_npeSum)
        
    if spec == "SHMS":    
        P_gtr_beta_uncut.Fill(event.P_gtr_beta)
        P_gtr_xp_uncut.Fill(event.P_gtr_xp)
        P_gtr_yp_uncut.Fill(event.P_gtr_yp)
        P_gtr_dp_uncut.Fill(event.P_gtr_dp)
        P_gtr_p_uncut.Fill(event.P_gtr_p)
        P_hod_goodscinhit_uncut.Fill(event.P_hod_goodscinhit)
        P_hod_goodstarttime_uncut.Fill(event.P_hod_goodstarttime)
        P_cal_etotnorm_uncut.Fill(event.P_cal_etotnorm)
        P_cal_etottracknorm_uncut.Fill(event.P_cal_etottracknorm)
        P_hgcer_npeSum_uncut.Fill(event.P_hgcer_npeSum)
        P_hgcer_xAtCer_uncut.Fill(event.P_hgcer_xAtCer)
        P_hgcer_yAtCer_uncut.Fill(event.P_hgcer_yAtCer)
        P_aero_npeSum_uncut.Fill(event.P_aero_npeSum)
        P_aero_xAtAero_uncut.Fill(event.P_aero_xAtAero)
        P_aero_yAtAero_uncut.Fill(event.P_aero_yAtAero)
        P_ngcer_npeSum_uncut.Fill(event.P_ngcer_npeSum)
        P_ngcer_xAtCer_uncut.Fill(event.P_ngcer_xAtCer)
        P_ngcer_yAtCer_uncut.Fill(event.P_ngcer_yAtCer)
        P_RFTime_Dist_uncut.Fill(event.P_RF_Dist)
        P_hgcer_npeSum_vs_aero_npeSum_uncut.Fill(event.P_hgcer_npeSum, event.P_aero_npeSum)
        P_hgcer_yAtCer_vs_hgcer_xAtCer_uncut.Fill(event.P_hgcer_yAtCer, event.P_hgcer_xAtCer)
        P_aero_yAtAero_vs_aero_xAtAero_uncut.Fill(event.P_aero_yAtAero, event.P_aero_xAtAero)
        P_cal_etottracknorm_vs_P_ngcer_npeSum_uncut.Fill(event.P_cal_etottracknorm, event.P_ngcer_npeSum)
        P_ngcer_yAtCer_vs_ngcer_xAtCer_uncut.Fill(event.P_ngcer_yAtCer, event.P_ngcer_xAtCer)
        P_ngcer_npeSum_vs_hgcer_npeSum_uncut.Fill(event.P_ngcer_npeSum, event.P_hgcer_npeSum)
        P_ngcer_npeSum_vs_aero_npeSum_uncut.Fill(event.P_ngcer_npeSum, event.P_aero_npeSum)
        

for event in Cut_Events_All_tree:
    if spec == "HMS":
        H_gtr_beta_cut_all.Fill(event.H_gtr_beta)
        H_gtr_xp_cut_all.Fill(event.H_gtr_xp)
        H_gtr_yp_cut_all.Fill(event.H_gtr_yp)
        H_gtr_dp_cut_all.Fill(event.H_gtr_dp)
        H_gtr_p_cut_all.Fill(event.H_gtr_p)
        H_hod_goodscinhit_cut_all.Fill(event.H_hod_goodscinhit)
        H_hod_goodstarttime_cut_all.Fill(event.H_hod_goodstarttime)
        H_cal_etotnorm_cut_all.Fill(event.H_cal_etotnorm)
        H_cal_etottracknorm_cut_all.Fill(event.H_cal_etottracknorm)
        H_cer_npeSum_cut_all.Fill(event.H_cer_npeSum)
        H_RFTime_Dist_cut_all.Fill(event.H_RF_Dist)
        H_cal_etottracknorm_vs_H_cer_npeSum_cut_all.Fill(event.H_cal_etottracknorm, event.H_cer_npeSum)
        
    if spec == "SHMS":        
        P_gtr_beta_cut_all.Fill(event.P_gtr_beta)
        P_gtr_xp_cut_all.Fill(event.P_gtr_xp)
        P_gtr_yp_cut_all.Fill(event.P_gtr_yp)
        P_gtr_dp_cut_all.Fill(event.P_gtr_dp)
        P_gtr_p_cut_all.Fill(event.P_gtr_p)
        P_hod_goodscinhit_cut_all.Fill(event.P_hod_goodscinhit)
        P_hod_goodstarttime_cut_all.Fill(event.P_hod_goodstarttime)
        P_cal_etotnorm_cut_all.Fill(event.P_cal_etotnorm)
        P_cal_etottracknorm_cut_all.Fill(event.P_cal_etottracknorm)
        P_hgcer_npeSum_cut_all.Fill(event.P_hgcer_npeSum)
        P_hgcer_xAtCer_cut_all.Fill(event.P_hgcer_xAtCer)
        P_hgcer_yAtCer_cut_all.Fill(event.P_hgcer_yAtCer)
        P_aero_npeSum_cut_all.Fill(event.P_aero_npeSum)
        P_aero_xAtAero_cut_all.Fill(event.P_aero_xAtAero)
        P_aero_yAtAero_cut_all.Fill(event.P_aero_yAtAero)
        P_ngcer_npeSum_cut_all.Fill(event.P_ngcer_npeSum)
        P_ngcer_xAtCer_cut_all.Fill(event.P_ngcer_xAtCer)
        P_ngcer_yAtCer_cut_all.Fill(event.P_ngcer_yAtCer)   
        P_RFTime_Dist_cut_all.Fill(event.P_RF_Dist)
        P_hgcer_npeSum_vs_aero_npeSum_cut_all.Fill(event.P_hgcer_npeSum, event.P_aero_npeSum)
        P_hgcer_yAtCer_vs_hgcer_xAtCer_cut_all.Fill(event.P_hgcer_yAtCer, event.P_hgcer_xAtCer)
        P_aero_yAtAero_vs_aero_xAtAero_cut_all.Fill(event.P_aero_yAtAero, event.P_aero_xAtAero)
        P_cal_etottracknorm_vs_P_ngcer_npeSum_cut_all.Fill(event.P_cal_etottracknorm, event.P_ngcer_npeSum)
        P_ngcer_yAtCer_vs_ngcer_xAtCer_cut_all.Fill(event.P_ngcer_yAtCer, event.P_ngcer_xAtCer)
        P_ngcer_npeSum_vs_hgcer_npeSum_cut_all.Fill(event.P_ngcer_npeSum, event.P_hgcer_npeSum)
        P_ngcer_npeSum_vs_aero_npeSum_cut_all.Fill(event.P_ngcer_npeSum, event.P_aero_npeSum) 
                
############################################################################################################################################

if spec == "HMS":
    # Saving histograms in PDF
    c1_mom = TCanvas("c1_mom", "Momentum Distributions", 100, 0, 1000, 900)
    c1_mom.Divide(1,1)
    c1_mom.cd(1)
    H_gtr_p_uncut.SetLineColor(2)
    H_gtr_p_uncut.Draw("hist")
    H_gtr_p_cut_all.SetLineColor(4)
    H_gtr_p_cut_all.Draw("histsame")
    c1_mom.Print(Analysis_Distributions + '(')

    c1_acpt = TCanvas("c1_H_kin", "Electron-Proton Acceptance Distributions", 100, 0, 1000, 900)
    c1_acpt.Divide(3,1)
    c1_acpt.cd(1)
    gPad.SetLogy()
    H_gtr_xp_uncut.SetLineColor(2)
    H_gtr_xp_uncut.Draw()
    H_gtr_xp_cut_all.SetLineColor(4)
    H_gtr_xp_cut_all.Draw("same")
    c1_acpt.cd(2)
    gPad.SetLogy()
    H_gtr_yp_uncut.SetLineColor(2)
    H_gtr_yp_uncut.Draw()
    H_gtr_yp_cut_all.SetLineColor(4)
    H_gtr_yp_cut_all.Draw("same")
    c1_acpt.cd(3)
    gPad.SetLogy()
    H_gtr_dp_uncut.SetLineColor(2)
    H_gtr_dp_uncut.Draw()
    H_gtr_dp_cut_all.SetLineColor(4)
    H_gtr_dp_cut_all.Draw("same")
    # TLegend (x1, y1, x2, y2) 
    legend2 = ROOT.TLegend(0.115, 0.8, 0.6, 0.9)
    legend2.AddEntry("H_gtr_dp_uncut", "without cuts", "l")
    legend2.AddEntry("H_gtr_dp_cut_all", "with cuts (acpt/RF/PID)", "l")
    legend2.Draw("same")
    c1_acpt.Print(Analysis_Distributions)
    
    c1_pid = TCanvas("c1_pid", "Electron-Proton CAL Distributions", 100, 0, 1000, 900)
    c1_pid.Divide(2,1)
    c1_pid.cd(1)
    gPad.SetLogy()
    H_cal_etottracknorm_uncut.SetLineColor(2)
    H_cal_etottracknorm_uncut.Draw()
    H_cal_etottracknorm_cut_all.SetLineColor(4)
    H_cal_etottracknorm_cut_all.Draw("same")
    legend7 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
    legend7.AddEntry("H_cal_etottracknorm_uncut", "without cuts", "l")
    legend7.AddEntry("H_cal_etottracknorm_cut_all", "with cuts (acpt/RF/PID)", "l")
    legend7.Draw("same")
    c1_pid.cd(2)
    H_cal_etottracknorm_vs_H_cer_npeSum_cut_all.Draw("COLZ")
    c1_pid.Print(Analysis_Distributions + ')')
    
if spec == "SHMS":
    # Saving histograms in PDF
    c1_mom = TCanvas("c1_mom", "Momentum Distributions", 100, 0, 1000, 900)
    c1_mom.Divide(1,1)
    c1_mom.cd(1)
    P_gtr_p_uncut.SetLineColor(2)
    P_gtr_p_uncut.Draw("hist")
    P_gtr_p_cut_all.SetLineColor(4)
    P_gtr_p_cut_all.Draw("histsame")
    c1_mom.Print(Analysis_Distributions + '(')

    c1_acpt = TCanvas("c1_P_kin", "Electron-Proton Acceptance Distributions", 100, 0, 1000, 900)
    c1_acpt.Divide(3,1)
    c1_acpt.cd(1)
    gPad.SetLogy()
    P_gtr_xp_uncut.SetLineColor(2)
    P_gtr_xp_uncut.Draw()
    P_gtr_xp_cut_all.SetLineColor(4)
    P_gtr_xp_cut_all.Draw("same")
    c1_acpt.cd(2)
    gPad.SetLogy()
    P_gtr_yp_uncut.SetLineColor(2)
    P_gtr_yp_uncut.Draw()
    P_gtr_yp_cut_all.SetLineColor(4)
    P_gtr_yp_cut_all.Draw("same")
    c1_acpt.cd(3)
    gPad.SetLogy()
    P_gtr_dp_uncut.SetLineColor(2)
    P_gtr_dp_uncut.Draw()
    P_gtr_dp_cut_all.SetLineColor(4)
    P_gtr_dp_cut_all.Draw("same")
    c1_acpt.Print(Analysis_Distributions)
    
    c1_pid = TCanvas("c1_pid", "Electron-Proton CAL Distributions", 100, 0, 1000, 900)
    c1_pid.Divide(1,1)
    c1_pid.cd(1)
    gPad.SetLogy()
    P_cal_etottracknorm_uncut.SetLineColor(2)
    P_cal_etottracknorm_uncut.Draw()
    P_cal_etottracknorm_cut_all.SetLineColor(4)
    P_cal_etottracknorm_cut_all.Draw("same")
    legend8 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
    legend8.AddEntry("P_cal_etottracknorm_uncut", "without cuts", "l")
    legend8.AddEntry("P_cal_etottracknorm_cut_all", "with cuts (acpt/RF/PID)", "l")
    legend8.Draw("same")
    c1_pid.Print(Analysis_Distributions)
    
    c2_pid = TCanvas("c2_pid", "Electron-Proton Aero/HGC/NGC PID Distributions", 100, 0, 1000, 900)
    c2_pid.Divide(3,1)
    c2_pid.cd(1)
    gPad.SetLogy()
    P_hgcer_npeSum_uncut.SetLineColor(2)
    P_hgcer_npeSum_uncut.Draw()
    P_hgcer_npeSum_cut_all.SetLineColor(4)
    P_hgcer_npeSum_cut_all.Draw("same")
    legend10 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
    legend10.AddEntry("P_hgcer_npeSum_uncut", "without cuts", "l")
    legend10.AddEntry("P_hgcer_npeSum_cut_all", "with cuts (acpt/RF/PID)", "l")
    legend10.Draw("same")
    c2_pid.cd(2)
    gPad.SetLogy()
    P_aero_npeSum_uncut.SetLineColor(2)
    P_aero_npeSum_uncut.Draw()
    P_aero_npeSum_cut_all.SetLineColor(4)
    P_aero_npeSum_cut_all.Draw("same")
    legend11 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
    legend11.AddEntry("P_aero_npeSum_uncut", "without cuts", "l")
    legend11.AddEntry("P_aero_npeSum_cut_all", "with cuts (acpt/RF/PID)", "l")
    legend11.Draw("same")
    c2_pid.cd(3)
    gPad.SetLogy()
    P_ngcer_npeSum_uncut.SetLineColor(2)
    P_ngcer_npeSum_uncut.Draw()
    P_ngcer_npeSum_cut_all.SetLineColor(4)
    P_ngcer_npeSum_cut_all.Draw("same") 
    legend12 = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
    legend12.AddEntry("P_ngcer_npeSum_uncut", "without cuts", "l")
    legend12.AddEntry("P_ngcer_npeSum_cut_all", "with cuts (acpt/RF/PID)", "l")
    legend12.Draw("same")
    c2_pid.Print(Analysis_Distributions)
    
    c3_pid = TCanvas("c3_pid", "Electron-Proton Aero/HGC/NGC PID Distributions", 100, 0, 1000, 900)
    c3_pid.Divide(2,3)
    c3_pid.cd(1)
    gPad.SetLogz()
    P_hgcer_npeSum_vs_aero_npeSum_uncut.Draw("COLZ")
    c3_pid.cd(2)
    P_hgcer_npeSum_vs_aero_npeSum_cut_all.Draw("COLZ")
    c3_pid.cd(3)
    P_ngcer_npeSum_vs_hgcer_npeSum_uncut.Draw("COLZ")
    c3_pid.cd(4)
    P_ngcer_npeSum_vs_hgcer_npeSum_cut_all.Draw("COLZ")
    c3_pid.cd(5)
    P_ngcer_npeSum_vs_aero_npeSum_uncut.Draw("COLZ")
    c3_pid.cd(6)
    P_ngcer_npeSum_vs_aero_npeSum_cut_all.Draw("COLZ")
    c3_pid.Print(Analysis_Distributions + ')')

#############################################################################################################################################

# Making directories in output file
outHistFile = ROOT.TFile.Open("%s/%s_%s_%s_Output_Data.root" % (OUTPATH, spec, runNum, MaxEvent), "RECREATE")                                                                                                    
d_Uncut_Events = outHistFile.mkdir("Uncut_Events")
d_Cut_Events_All = outHistFile.mkdir("Cut_Events_All")
d_Cut_Events_Prompt = outHistFile.mkdir("Cut_Events_Prompt")
d_Cut_Events_Random = outHistFile.mkdir("Cut_Events_Random")

# Writing Histograms for protons                                                                  
d_Uncut_Events.cd()
if spec == "HMS":
    H_gtr_beta_uncut.Write()
    H_gtr_xp_uncut.Write()
    H_gtr_yp_uncut.Write()
    H_gtr_dp_uncut.Write()
    H_gtr_p_uncut.Write()
    H_hod_goodscinhit_uncut.Write()
    H_hod_goodstarttime_uncut.Write()
    H_cal_etotnorm_uncut.Write()
    H_cal_etottracknorm_uncut.Write()
    H_cer_npeSum_uncut.Write()
    H_RFTime_Dist_uncut.Write()
    H_cal_etottracknorm_vs_H_cer_npeSum_uncut.Write()

if spec == "SHMS":    
    P_gtr_beta_uncut.Write()
    P_gtr_xp_uncut.Write()
    P_gtr_yp_uncut.Write()
    P_gtr_dp_uncut.Write()
    P_gtr_p_uncut.Write()
    P_hod_goodscinhit_uncut.Write()
    P_hod_goodstarttime_uncut.Write()
    P_cal_etotnorm_uncut.Write()
    P_cal_etottracknorm_uncut.Write()
    P_hgcer_npeSum_uncut.Write()
    P_hgcer_xAtCer_uncut.Write()
    P_hgcer_yAtCer_uncut.Write()
    P_aero_npeSum_uncut.Write()
    P_aero_xAtAero_uncut.Write()
    P_aero_yAtAero_uncut.Write()
    P_ngcer_npeSum_uncut.Write()
    P_ngcer_xAtCer_uncut.Write()
    P_ngcer_yAtCer_uncut.Write() 
    P_RFTime_Dist_uncut.Write()
    P_hgcer_npeSum_vs_aero_npeSum_uncut.Write()
    P_hgcer_yAtCer_vs_hgcer_xAtCer_uncut.Write()
    P_aero_yAtAero_vs_aero_xAtAero_uncut.Write()
    P_cal_etottracknorm_vs_P_ngcer_npeSum_uncut.Write()
    P_ngcer_yAtCer_vs_ngcer_xAtCer_uncut.Write()
    P_ngcer_npeSum_vs_hgcer_npeSum_uncut.Write()
    P_ngcer_npeSum_vs_aero_npeSum_uncut.Write()

d_Cut_Events_All.cd()
if spec == "HMS":
    H_gtr_beta_cut_all.Write()
    H_gtr_xp_cut_all.Write()
    H_gtr_yp_cut_all.Write()
    H_gtr_dp_cut_all.Write()
    H_gtr_p_cut_all.Write()
    H_hod_goodscinhit_cut_all.Write()
    H_hod_goodstarttime_cut_all.Write()
    H_cal_etotnorm_cut_all.Write()
    H_cal_etottracknorm_cut_all.Write()
    H_cer_npeSum_cut_all.Write()
    H_RFTime_Dist_cut_all.Write()
    H_cal_etottracknorm_vs_H_cer_npeSum_cut_all.Write()

if spec == "SHMS":    
    P_gtr_beta_cut_all.Write()
    P_gtr_xp_cut_all.Write()
    P_gtr_yp_cut_all.Write()
    P_gtr_dp_cut_all.Write()
    P_gtr_p_cut_all.Write()
    P_hod_goodscinhit_cut_all.Write()
    P_hod_goodstarttime_cut_all.Write()
    P_cal_etotnorm_cut_all.Write()
    P_cal_etottracknorm_cut_all.Write()
    P_hgcer_npeSum_cut_all.Write()
    P_hgcer_xAtCer_cut_all.Write()
    P_hgcer_yAtCer_cut_all.Write()
    P_aero_npeSum_cut_all.Write()
    P_aero_xAtAero_cut_all.Write()
    P_aero_yAtAero_cut_all.Write()
    P_ngcer_npeSum_cut_all.Write()
    P_ngcer_xAtCer_cut_all.Write()
    P_ngcer_yAtCer_cut_all.Write()
    P_RFTime_Dist_cut_all.Write()
    P_hgcer_npeSum_vs_aero_npeSum_cut_all.Write()
    P_hgcer_yAtCer_vs_hgcer_xAtCer_cut_all.Write()
    P_aero_yAtAero_vs_aero_xAtAero_cut_all.Write()
    P_cal_etottracknorm_vs_P_ngcer_npeSum_cut_all.Write()
    P_ngcer_yAtCer_vs_ngcer_xAtCer_cut_all.Write()
    P_ngcer_npeSum_vs_hgcer_npeSum_cut_all.Write()
    P_ngcer_npeSum_vs_aero_npeSum_cut_all.Write()
        
outHistFile.Close()
infile.Close() 
print ("Processing Complete")
