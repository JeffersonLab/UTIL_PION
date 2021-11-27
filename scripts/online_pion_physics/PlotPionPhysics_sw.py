#! /usr/bin/python
###########################################################################################################################
# Created - 20/July/21, Author - Muhammad Junaid (mjo147@uregina.ca), University of Regina, Canada (Copyright (c) junaid) #
###########################################################################################################################
# Python version of the pion plotting script. Now utilises uproot to select event of each type and writes them to a root file.
# Python should allow for easier reading of databases storing diferent variables.
# This version of script is for shift workers at JLab
# To run this script, execute: python3 scriptname runnumber

# 05/October/21: Jacob Murphy added in another page of plots with focal plane variables vs beta.
# 11/10/21 - SJDK - There's lots of commented out plots in here, if they aren't being used anymore, can they please just be deleted?
# 16/10/21 - JM - Added in X/Y Calo with projections of ADC hits (SJDK Stop-Gap Measure)

###################################################################################################################################################

# Import relevant packages
import ROOT
import sys, math, os, subprocess
import array
import re # Regexp package - for string manipulation
from ROOT import TCanvas, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TGaxis, TLine, TMath, TPaveText, TArc, TGraphPolar 
from ROOT import kBlack, kBlue, kRed
sys.path.insert(0, 'python/')

##################################################################################################################################################

# Defining some variables here
#minbin = 0.90 # minimum bin for selecting neutrons events in missing mass distribution
#maxbin = 0.98 # maximum bin for selecting neutrons events in missing mass distribution
# SJDK - 25/11/21 - Doubled the range of the MM selection for the LD2 runs at Q2 = 6.0. I expanded this asymmetrically, I raised the upper limit by 0.06 and the lower limit by 0.02
minbin = 0.88 # minimum bin for selecting neutrons events in missing mass distribution
maxbin = 1.04 # maximum bin for selecting neutrons events in missing mass distribution
minrangeuser = 0 # min range for -t vs phi plot
maxrangeuser = 1.5 #  max range for -t vs phi plot     : 2021/11/11 JM updated t bounds - 20/11/2021 - SJDK, adjusted from 1.5 to 0.9 for Q2 = 6.0 middle epsilon setting # NH 2021 11 25 - changed range to 1.5

##################################################################################################################################################

# Check the number of arguments provided to the script
FilenameOverride=False # SJDK 21/09/21 - Added a secret 4th argument so that a full kinematic can be processed, the run number is needed for cuts, but for the kinematic analysis, the filename is not by run number
if len(sys.argv)-1!=3:
    if len(sys.argv)-1!=4:
        print("!!!!! ERROR !!!!!\n Expected 3 arguments\n Usage is with - ROOTfileSuffix RunNumber MaxEvents \n!!!!! ERROR !!!!!")
        sys.exit(1)
    else:
        print ("!!!!! Running with secret 4th argument - FilenameOverride - Taking file name to process as stated EXACTLY in 4th arg !!!!!")
        FilenameOverride=sys.argv[4] # If 4 arguments provided, set the FilenameOverride value to be arg 4

##################################################################################################################################################

# Input params - run number and max number of events
# 23/09/21 - SJDK - Changed ROOTPrefix to ROOTSuffix so that this is actually accurate (and the insturctions printed above make sense)
ROOTSuffix = sys.argv[1]
runNum = sys.argv[2]
MaxEvent = sys.argv[3]

USER = subprocess.getstatusoutput("whoami") # Grab user info for file finding
HOST = subprocess.getstatusoutput("hostname")

if ("farm" in HOST[1]):
    REPLAYPATH = "/group/c-pionlt/online_analysis/hallc_replay_lt"
elif ("qcd" in HOST[1]):
    REPLAYPATH = "/group/c-pionlt/USERS/%s/hallc_replay_lt" % USER[1]
elif ("cdaq" in HOST[1]):
    REPLAYPATH = "/home/cdaq/hallc-online/hallc_replay_lt"
elif ("phys.uregina" in HOST[1]):
    REPLAYPATH = "/home/%s/work/JLab/hallc_replay_lt" % USER[1]
elif("skynet" in HOST[1]):
    REPLAYPATH = "/home/%s/Work/JLab/hallc_replay_lt" % USER[1]

#################################################################################################################################################

# Add more path setting as needed in a similar manner                                                                                                                                                          
OUTPATH = "%s/UTIL_PION/OUTPUT/Analysis/PionLT" % REPLAYPATH        # Output folder location                                                                                                     
sys.path.insert(0, '%s/UTIL_PION/bin/python/' % REPLAYPATH)
import kaonlt as klt # Import kaonlt module, need the path setting line above prior to importing this                                                                                                         
print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER[1], HOST[1], REPLAYPATH))
if (FilenameOverride == False):
    Pion_Analysis_Distributions = "%s/%s_%s_sw_Pion_Analysis_Distributions.pdf" % (OUTPATH, runNum, MaxEvent)
elif (FilenameOverride != False): # If filename override set, format the file name based upon the override file name
    Pion_Analysis_Distributions = "%s/%s_Pion_Analysis_Distributions.pdf" %(OUTPATH, (FilenameOverride.split("_Analysed_Data.root",1)[0]))

if (FilenameOverride != False):
    # Split string on W value, replace p with . and strip the leading Q before converting to a float
    Q2Val = float(((FilenameOverride.split("W")[0]).replace("p",".")).lstrip("Q"))
    # Split string on W, replace p with . and strip leading W, THEN split based on numeric characters and non numeric characters. Select the first entry and convert it to a float with 1 DP precision
    WVal = round(float(re.split('[a-zAZ]',((FilenameOverride.split("W")[1]).replace("p",".")).lstrip("W"))[0]),1)
    Q2min = Q2Val - 2 # Minimum value of Q2 on the Q2 vs W plot
    Q2max = Q2Val + 2 # Maximum value of Q2 on the Q2 vs W plot
    Wmin = WVal - 0.5 # min y-range for Q2vsW plot
    Wmax = WVal + 0.5 # max y-range for Q2vsW plot

#################################################################################################################################################

# Construct the name of the rootfile based upon the info we provided
if (FilenameOverride == False): # Standard running condition, construct file name from run number and max events e.t.c.
    rootName = "%s/UTIL_PION/OUTPUT/Analysis/PionLT/%s_%s_%s.root" % (REPLAYPATH, runNum, MaxEvent, ROOTSuffix)     # Input file location and variables taking
elif (FilenameOverride != False): # Special condition, with 4th arg, use 4th arg as file name
    rootName = "%s/UTIL_PION/OUTPUT/Analysis/PionLT/%s" % (REPLAYPATH, FilenameOverride)
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
if os.path.isfile(rootName):
    print ("%s exists, processing" % (rootName))
else:
    print ("%s not found - do you have the correct sym link/folder set up?" % (rootName))
    sys.exit(4)
print("Output path checks out, outputting to %s" % (OUTPATH))

###############################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Section for grabing Prompt/Random selection parameters from PARAM file
PARAMPATH = "%s/UTIL_PION/DB/PARAM" % REPLAYPATH
print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER[1], HOST[1], REPLAYPATH))
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

###############################################################################################################################################

# Read stuff from the main event tree
infile = ROOT.TFile.Open(rootName, "READ")
Uncut_Pion_Events_tree = infile.Get("Uncut_Pion_Events")
Cut_Pion_Events_noRF_tree = infile.Get("Cut_Pion_Events_noRF")
Cut_Pion_Events_All_tree = infile.Get("Cut_Pion_Events_All")
Cut_Pion_Events_Prompt_tree = infile.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Random_tree = infile.Get("Cut_Pion_Events_Random")
Cut_Pion_Events_Prompt_MM_tree = infile.Get("Cut_Pion_Events_Prompt_MM")

###############################################################################################################################################

# Defining Histograms for Pions
P_RFTime_pions_cut_noRF = ROOT.TH1D("P_RFTime_pions_cut_noRF", "SHMS RFTime; SHMS_RFTime; Counts", 200, 0, 4)

H_xp_pions_uncut = ROOT.TH1D("H_xp_pions_uncut", "HMS x'; HMS_gtr_xp; Counts", 200, -0.2, 0.2)
H_yp_pions_uncut = ROOT.TH1D("H_yp_pions_uncut", "HMS y'; HMS_gtr_yp; Counts", 200, -0.2, 0.2)
H_dp_pions_uncut = ROOT.TH1D("H_dp_pions_uncut", "HMS #delta; HMS_gtr_dp; Counts", 200, -12, 12)
H_cal_etottracknorm_pions_uncut = ROOT.TH1D("H_cal_etottracknorm_pions_uncut", "HMS cal etottracknorm; HMS_cal_etottracknorm; Counts", 300, 0.2, 1.8)
P_xp_pions_uncut = ROOT.TH1D("P_xp_pions_uncut", "SHMS x'; SHMS_gtr_xp; Counts", 200, -0.2, 0.2)
P_yp_pions_uncut = ROOT.TH1D("P_yp_pions_uncut", "SHMS y'; SHMS_gtr_yp; Counts", 200, -0.2, 0.2)
P_dp_pions_uncut = ROOT.TH1D("P_dp_pions_uncut", "SHMS #delta; SHMS_gtr_dp; Counts", 200, -30, 30)
P_cal_etottracknorm_pions_uncut = ROOT.TH1D("P_cal_etottracknorm_pions_uncut", "SHMS cal etottracknorm; SHMS_cal_etottracknorm; Counts", 200, 0, 1.6)
P_cal_fly_numGoodAdcHits_pions_uncut = ROOT.TH1D("P_cal_fly_numGoodAdcHits_pions_uncut", "SHMS Shower Good Occupancy; PMT Number; Number of Good ADC Hits", 224, 0.5, 224.5)
P_hgcer_npe_pions_uncut = ROOT.TH1D("P_hgcer_npe_pions_uncut", "SHMS HGC npeSum; SHMS_hgcer_npeSum; Counts", 200, 0, 50)
P_aero_npe_pions_uncut = ROOT.TH1D("P_aero_npe_pions_uncut", "SHMS aero npeSum; SHMS_aero_npeSum; Counts", 200, 0, 50)
P_ngcer_npe_pions_uncut = ROOT.TH1D("P_ngcer_npe_pions_uncut", "SHMS NGC npeSum; SHMS_ngcer_npeSum; Counts", 200, 0, 50)
P_MMpi_pions_uncut = ROOT.TH1D("P_MMpi_pions_uncut", "MIssing Mass (no cuts); MM_{#pi}; Counts", 200, 0.5, 1.8)
P_RFTime_pions_uncut = ROOT.TH1D("P_RFTime_pions_uncut", "SHMS RFTime; SHMS_RFTime; Counts", 200, 0, 4)
#ePiCoinTime_pions_uncut = ROOT.TH1D("ePiCoinTime_pions_uncut", "Electron-Pion CTime (no cuts); e #pi Coin_Time; Counts", 200, -50, 50)
# SJDK 26/11/21 - Narrowed the range for 2ns running, bin size still the same
ePiCoinTime_pions_uncut = ROOT.TH1D("ePiCoinTime_pions_uncut", "Electron-Pion CTime (no cuts); e #pi Coin_Time; Counts", 100, -25, 25)

H_xp_pions_cut = ROOT.TH1D("H_xp_pions_cut", "HMS x'; HMS_gtr_xp; Counts", 200, -0.2, 0.2)
H_yp_pions_cut = ROOT.TH1D("H_yp_pions_cut", "HMS y'; HMS_gtr_yp; Counts", 200, -0.2, 0.2)
H_dp_pions_cut = ROOT.TH1D("H_dp_pions_cut", "HMS #delta; HMS_gtr_dp; Counts", 200, -12, 12)
H_cal_etottracknorm_pions_cut = ROOT.TH1D("H_cal_etottracknorm_pions_cut", "HMS cal etottracknorm; HMS_cal_etottracknorm; Counts", 200, 0.6, 1.4)
P_xp_pions_cut = ROOT.TH1D("P_xp_pions_cut", "SHMS x'; SHMS_gtr_xp; Counts", 200, -0.2, 0.2)
P_yp_pions_cut = ROOT.TH1D("P_yp_pions_cut", "SHMS y'; SHMS_gtr_yp; Counts", 200, -0.2, 0.2)
P_dp_pions_cut = ROOT.TH1D("P_dp_pions_cut", "SHMS #delta; SHMS_gtr_dp; Counts", 200, -15, 15)
P_cal_etottracknorm_pions_cut = ROOT.TH1D("P_cal_etottracknorm_pions_cut", "SHMS cal etottracknorm; SHMS_cal_etottracknorm; Counts", 200, 0, 1.6)
P_cal_fly_numGoodAdcHits_pions_cut = ROOT.TH1D("P_cal_fly_numGoodAdcHits_pions_cut", "SHMS Shower Good Occupancy; PMT Number; Number of Good ADC Hits", 224, 0.5, 224.5)
P_hgcer_npe_pions_cut = ROOT.TH1D("P_hgcer_npe_pions_cut", "SHMS HGC npeSum; SHMS_hgcer_npeSum; Counts", 200, 0, 50)
P_aero_npe_pions_cut = ROOT.TH1D("P_aero_npe_pions_cut", "SHMS aero npeSum; SHMS_aero_npeSum; Counts", 200, 0, 50)
P_ngcer_npe_pions_cut = ROOT.TH1D("P_ngcer_npe_pions_cut", "SHMS NGC npeSum; SHMS_ngcer_npeSum; Counts", 200, 0, 50)
P_MMpi_pions_cut = ROOT.TH1D("P_MMpi_pions_cut", "Missing Mass (with cuts); MM_{#pi}; Counts", 200, 0.5, 1.8)
P_RFTime_pions_cut = ROOT.TH1D("P_RFTime_pions_cut", "SHMS RFTime; SHMS_RFTime; Counts", 200, 0, 4)
#ePiCoinTime_pions_cut = ROOT.TH1D("ePiCoinTime_pions_cut", "Electron-Pion CTime (with cuts); e #pi Coin_Time; Counts", 200, -50, 50)
# SJDK 26/11/21 - Narrowed the range for 2ns running, bin size still the same
ePiCoinTime_pions_cut = ROOT.TH1D("ePiCoinTime_pions_cut", "Electron-Pion CTime (with cuts); e #pi Coin_Time; Counts", 100, -25, 25)
epsilon_pions_cut = ROOT.TH1D("epsilon_pions_cut", "Epsilon Dist for Prompt Events (Incl MM Cut); epsilon; Counts", 200, 0, 1.0)

ePiCoinTime_pions_cut_prompt = ROOT.TH1D("ePiCoinTime_pions_cut_prompt", "Electron-Pion CTime; e #pi Coin_Time; Counts", 8, -2, 2)
P_MMpi_pions_cut_prompt = ROOT.TH1D("P_MMpi_pions_cut_prompt", "Missing Mass; MM_{#pi}; Counts", 200, 0.5, 1.8)

#ePiCoinTime_pions_cut_randm = ROOT.TH1D("ePiCoinTime_pions_cut_randm", "Electron-Pion CTime; e #pi Coin_Time; Counts", 160, -40, 40)
# SJDK 26/11/21 - Narrowed the range for 2ns running, bin size still the same
ePiCoinTime_pions_cut_randm = ROOT.TH1D("ePiCoinTime_pions_cut_randm", "Electron-Pion CTime; e #pi Coin_Time; Counts", 80, -20, 20)
P_MMpi_pions_cut_randm = ROOT.TH1D("P_MMpi_pions_cut_randm", "Missing Mass; MM_{#pi}; Counts", 200, 0.5, 1.8)

P_MMpi_pions_cut_randm_scaled = ROOT.TH1D("P_MMpi_pions_cut_randm_scaled", "Missing Mass; MM_{#pi}; Counts", 200, 0.5, 1.8)
P_MMpi_pions_cut_randm_sub = ROOT.TH1D("P_MMpi_pions_cut_randm_sub", "Missing Mass Rndm Sub; MM_{#pi}; Counts", 200, 0.5, 1.8)
# SJDK - 26/10/21 - Changed the plot title to be more accurate
phiq_plot = ROOT.TH1D("Phiq", "#phi Dist for Prompt Events (Incl MM Cut); #phi; Counts", 12, -3.14, 3.14) # 2021/10/25 NH - added these at garths request
t_plot = ROOT.TH1D("-t", "-t Dist for Prompt Events (Incl MM Cut); -t; Counts", 48, 0, 1.5) # 2021/11/11 JM changed t range
if (FilenameOverride == False): # Standard running condition, construct file name from run number and max events e.t.c.
    Q2_pions_cut = ROOT.TH1D("Q2_pions_cut", "Q2 Dist for Prompt Events (Incl MM Cut); Q2; Counts", 200, 3, 8) # 2021/11/25 JM - Adjusting Hist bounds for the 8.0 Physics settings 
    W_pions_cut = ROOT.TH1D("W_pions_cut", "W Dist for Prompt Events (Incl MM Cut); W; Counts", 200, 2., 3.)
elif (FilenameOverride != False): # Special case, run with specifc file name, construct histo with ranges based upon filename
    Q2_pions_cut = ROOT.TH1D("Q2_pions_cut", "Q2 Dist for Prompt Events (Incl MM Cut); Q2; Counts", 200, Q2min, Q2max) 
    W_pions_cut = ROOT.TH1D("W_pions_cut", "W Dist for Prompt Events (Incl MM Cut); W; Counts", 200, Wmin, Wmax)  

##############################################################################################################################################

# 2D Histograms for pions
P_hgcer_vs_aero_npe_pions_uncut = ROOT.TH2D("P_hgcer_vs_aero_npe_pions_uncut", "SHMS HGC npeSum vs SHMS Aero npeSum (no cut); SHMS_hgcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
#ePiCoinTime_vs_MMpi_pions_uncut = ROOT.TH2D("ePiCoinTime_vs_MMpi_pions_uncut","Electron-Pion CTime vs Missing Mass (no cut); e #pi Coin_Time; MM_{#pi}", 160, -40, 40, 200, 0, 2)
#ePiCoinTime_vs_beta_pions_uncut = ROOT.TH2D("ePiCoinTime_vs_beta_pions_uncut", "Electron-Pion CTime vs SHMS #beta (no cut); e #pi Coin_Time; SHMS_#beta", 160, -40, 40, 200, 0, 2)
# SJDK 26/11/21 - Narrowed the range for 2ns running, bin size still the same
ePiCoinTime_vs_MMpi_pions_uncut = ROOT.TH2D("ePiCoinTime_vs_MMpi_pions_uncut","Electron-Pion CTime vs Missing Mass (no cut); e #pi Coin_Time; MM_{#pi}", 100, -25, 25, 200, 0, 2)
ePiCoinTime_vs_beta_pions_uncut = ROOT.TH2D("ePiCoinTime_vs_beta_pions_uncut", "Electron-Pion CTime vs SHMS #beta (no cut); e #pi Coin_Time; SHMS_#beta", 100, -25, 25, 200, 0, 2)
P_RFTime_vs_MMpi_pions_uncut = ROOT.TH2D("P_RFTime_vs_MMpi_pions_uncut", "SHMS RFTime vs Missing Mass (no cuts); SHMS_RFTime_Dist; MM_{#pi}", 100, 0, 4, 100, 0, 2)
P_cal_etottracknorm_vs_ngcer_npe_pions_uncut = ROOT.TH2D("P_cal_etottracknorm_vs_ngcer_npe_pions_uncut", "SHMS cal etottracknorm vs SHMS NGC xAtCer (no cut); SHMS_cal_etottracknorm; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
P_cal_fly_numGoodAdcHits_col_row_uncut = ROOT.TH2D("P_cal_fly_numGoodAdcHits_col_row_uncut","SHMS Shower  Occupancy ; Col Number ;  Row Number",14,0,14,16,0,16)
P_ngcer_vs_hgcer_npe_pions_uncut = ROOT.TH2D("P_ngcer_vs_hgcer_npe_pions_uncut", "SHMS NGC npeSum vs SHMS HGC npeSum (no cut); SHMS_ngcer_npeSum; SHMS_hgcer_npeSum", 100, 0, 50, 100, 0, 50)
P_ngcer_vs_aero_npe_pions_uncut = ROOT.TH2D("P_ngcer_vs_aero_npe_pions_uncut", "SHMS NGC npeSum vs SHMS aero npeSum (no cut); SHMS_ngcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
H_dp_vs_beta_pions_uncut = ROOT.TH2D("H_dp_vs_beta_pions_uncut", "HMS #delta vs HMS #beta (no cut); HMS #delta; HMS_#beta", 200, -12, 12, 200, 0, 2)
P_dp_vs_beta_pions_uncut = ROOT.TH2D("P_dp_vs_beta_pions_uncut", "SHMS #delta vs SHMS #beta (no cut); SHMS #delta; HMS_#beta", 200, -12, 12, 200, 0, 2)
H_xfp_vs_beta_pions_uncut = ROOT.TH2D("H_xfp_vs_beta_pions_uncut", "HMS X_{fp} vs HMS #beta (no cut); HMS X_{fp}; HMS_#beta", 160, -40, 40, 200, 0, 2)
P_xfp_vs_beta_pions_uncut = ROOT.TH2D("P_xfp_vs_beta_pions_uncut", "SHMS X_{fp} vs SHMS #beta (no cut); SHMS X_{fp}; SHMS_#beta", 160, -40, 40, 200, 0, 2)
H_yfp_vs_beta_pions_uncut = ROOT.TH2D("H_yfp_vs_beta_pions_uncut", "HMS Y_{fp} vs HMS #beta (no cut); HMS Y_{fp}; HMS_#beta", 100, -25, 25, 200, 0, 2)
P_yfp_vs_beta_pions_uncut = ROOT.TH2D("P_yfp_vs_beta_pions_uncut", "SHMS Y_{fp} vs SHMS #beta (no cut); SHMS Y_{fp}; SHMS_#beta", 160, -40, 40, 200, 0, 2)
H_xpfp_vs_beta_pions_uncut = ROOT.TH2D("H_xpfp_vs_beta_pions_uncut", "HMS X'_{fp} vs HMS #beta (no cut); HMS X'_{fp}; HMS_#beta", 160, -0.4, 0.4, 200, 0, 2)
P_xpfp_vs_beta_pions_uncut = ROOT.TH2D("P_xpfp_vs_beta_pions_uncut", "SHMS X'_{fp} vs SHMS #beta (no cut); SHMS X'_{fp}; SHMS_#beta", 160, -0.4, 0.4, 200, 0, 2)
H_ypfp_vs_beta_pions_uncut = ROOT.TH2D("H_ypfp_vs_beta_pions_uncut", "HMS Y'_{fp} vs HMS #beta (no cut); HMS Y'_{fp}; HMS_#beta", 160, -0.4, 0.4, 200, 0, 2)
P_ypfp_vs_beta_pions_uncut = ROOT.TH2D("P_ypfp_vs_beta_pions_uncut", "SHMS Y'_{fp} vs SHMS #beta (no cut); SHMS Y'_{fp}; SHMS_#beta", 160, -0.4, 0.4, 200, 0, 2)
P_MMpi_vs_beta_pions_uncut = ROOT.TH2D("P_MMpi_vs_beta_pions_uncut", "Missing Mass vs SHMS #beta (no cut); MM_{#pi}; SHMS_#beta", 100, 0, 2, 200, 0, 2)
P_cal_xy_pions_uncut = ROOT.TH2D("P_cal_xy_pions_uncut", "SHMS Calorimeter yCalo vs xCalo (no cuts); cal_yCalo(cm); cal_xCalo(cm)", 14, -62.3, 62.3, 16, -71.2, 71.2)
P_DPexit_xy_pions_uncut = ROOT.TH2D("P_DPexit_xy_pions_uncut", "SHMS Dipole Exit yExit vs xExit (no cuts); yExit(cm); xExit(cm)", 200, -50.0, 50.0, 200, -50.0, 50.0)

H_cal_etottracknorm_vs_cer_npe_pions_cut = ROOT.TH2D("H_cal_etottracknorm_vs_cer_npe_pions_cut","HMS cal etottracknorm vs HMS cer npeSum (with cuts); H_cal_etottracknorm; H_cer_npeSum",100, 0.5, 1.5, 100, 0, 40)
P_hgcer_vs_aero_npe_pions_cut = ROOT.TH2D("P_hgcer_vs_aero_npe_pions_cut", "SHMS HGC npeSum vs SHMS aero npeSum (with cuts); SHMS_hgcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
ePiCoinTime_vs_MMpi_pions_cut = ROOT.TH2D("ePiCoinTime_vs_MMpi_pions_cut","Electron-Pion CTime vs Missing Mass (with PID cuts); e #pi Coin_Time; MM_{#pi}", 100, -2, 2, 100, 0, 2)
ePiCoinTime_vs_beta_pions_cut = ROOT.TH2D("ePiCoinTime_vs_beta_pions_cut", "Electron-Pion CTime vs SHMS #beta (with PID cuts); e #pi Coin_Time; SHMS_#beta", 100, -2, 2, 100, 0.6, 1.4)
P_RFTime_vs_MMpi_pions_cut = ROOT.TH2D("P_RFTime_vs_MMpi_pions_cut", "SHMS RFTime vs Missing Mass (with cuts); SHMS_RFTime_Dist; MM_{#pi}", 100, 0, 4, 100, 0, 2)
P_cal_etottracknorm_vs_ngcer_npe_pions_cut = ROOT.TH2D("P_cal_etottracknorm_vs_ngcer_npe_pions_cut", "P cal etottracknorm vs SHMS NGC xAtCer (with cuts); SHMS_cal_etottracknorm; SHMS_ngcer_xAtCer", 100, -10, 10, 100, -10, 10)
P_cal_fly_numGoodAdcHits_col_row_cut = ROOT.TH2D("P_cal_fly_numGoodAdcHits_col_row_cut","SHMS Shower  Occupancy ; Col Number ;  Row Number",14,0,14,16,0,16)
P_ngcer_vs_hgcer_npe_pions_cut = ROOT.TH2D("P_ngcer_vs_hgcer_npe_pions_cut", "SHMS NGC npeSum vs SHMS HGC npeSum (with cuts); SHMS_ngcer_npeSum; SHMS_hgcer_npeSum", 100, 0, 50, 100, 0, 50)
P_ngcer_vs_aero_npe_pions_cut = ROOT.TH2D("P_ngcer_vs_aero_npe_pions_cut", "SHMS NGC npeSum vs SHMS aero npeSum (with cuts); SHMS_ngcer_npeSum; SHMS_aero_npeSum", 100, 0, 50, 100, 0, 50)
H_dp_vs_beta_pions_cut = ROOT.TH2D("H_dp_vs_beta_pions_cut", "HMS #delta vs HMS #beta (with cut); HMS #delta; HMS_#beta", 200, -12, 12, 200, 0, 2)
P_dp_vs_beta_pions_cut = ROOT.TH2D("P_dp_vs_beta_pions_cut", "SHMS #delta vs SHMS #beta (with cut); SHMS #delta; SHMS_#beta", 200, -30, 30, 200, 0, 2)
H_xfp_vs_beta_pions_cut = ROOT.TH2D("H_xfp_vs_beta_pions_cut", "HMS X_{fp} vs HMS #beta (with cut); HMS X_{fp}; HMS_#beta", 160, -40, 40, 200, 0, 2)
P_xfp_vs_beta_pions_cut = ROOT.TH2D("P_xfp_vs_beta_pions_cut", "SHMS X_{fp} vs SHMS #beta (with cut); SHMS X_{fp}; SHMS_#beta", 160, -40, 40, 200, 0, 2)
H_yfp_vs_beta_pions_cut = ROOT.TH2D("H_yfp_vs_beta_pions_cut", "HMS Y_{fp} vs HMS #beta (with cut); HMS Y_{fp}; HMS_#beta", 100, -25, 25, 200, 0, 2)
P_yfp_vs_beta_pions_cut = ROOT.TH2D("P_yfp_vs_beta_pions_cut", "SHMS Y_{fp} vs SHMS #beta (with cut); SHMS Y_{fp}; SHMS_#beta", 160, -40, 40, 200, 0, 2)
H_xpfp_vs_beta_pions_cut = ROOT.TH2D("H_xpfp_vs_beta_pions_cut", "HMS X'_{fp} vs HMS #beta (with cut); HMS X'_{fp}; HMS_#beta", 160, -40, 40, 200, 0, 2)
P_xpfp_vs_beta_pions_cut = ROOT.TH2D("P_xpfp_vs_beta_pions_cut", "SHMS X'_{fp} vs SHMS #beta (with cut); SHMS X'_{fp}; SHMS_#beta", 160, -40, 40, 200, 0, 2)
H_ypfp_vs_beta_pions_cut = ROOT.TH2D("H_ypfp_vs_beta_pions_cut", "HMS Y'_{fp} vs HMS #beta (with cut); HMS Y'_{fp}; HMS_#beta", 160, -40, 40, 200, 0, 2)
P_ypfp_vs_beta_pions_cut = ROOT.TH2D("P_ypfp_vs_beta_pions_cut", "SHMS Y'_{fp} vs SHMS #beta (with cut); SHMS Y'_{fp}; SHMS_#beta", 160, -40, 40, 200, 0, 2)
P_MMpi_vs_beta_pions_cut = ROOT.TH2D("P_MMpi_vs_beta_pions_cut", "Missing Mass vs SHMS #beta (with cut); MM_{#pi}; SHMS_#beta", 100, 0, 2, 200, 0, 2)
P_cal_xy_pions_cut = ROOT.TH2D("P_cal_xy_pions_cut", "SHMS Calorimeter yCalo vs xCalo (with cuts); cal_yCalo(cm); cal_xCalo(cm)", 14, -62.3, 62.3, 16, -71.2, 71.2)
P_DPexit_xy_pions_cut = ROOT.TH2D("P_DPexit_xy_pions_cut", "SHMS Dipole Exit  yExit vs xExit (with cuts); yExit(cm); xExit(cm)", 200, -50.0, 50.0, 200, -50.0, 50.0)
# SJDK - 26/10/21 - Changed the plot title to be more accurate
if (FilenameOverride == False): # Standard running condition, construct file name from run number and max events e.t.c.
    Q2vsW_pions_cut = ROOT.TH2D("Q2vsW_pions_cut", "Q2 vs W Dist for Prompt Events (Incl MM Cut); Q2; W", 200, 3, 8, 200, 2., 3.) # 2021/11/25 JM Adjusted Q2 bounds
elif (FilenameOverride != False): # Special case, run with specifc file name, construct histo with ranges based upon filename
    Q2vsW_pions_cut = ROOT.TH2D("Q2vsW_pions_cut", "Q2 vs W Dist for Prompt Events (Incl MM Cut); Q2; W", 200, Q2min, Q2max, 200, Wmin, Wmax)
# SJDK - 26/10/21 - Changed the plot title to be more accurate
# NH 2021 11 25 - I changed max value back to 1.5, do not change this or the t marks on the t-phi plots will be wrong!!! Instead ONLY change the range with the constanst on line 31!!
phiqvst_pions_cut = ROOT.TH2D("phiqvst_pions_cut","#phi vs -t Dist for Prompt Events (Incl MM Cut); #phi ;-t", 12, -3.14, 3.14, 48, 0.0, 1.5) #2021 08 12 - NH doubled binning range and # bins : 2021/11/11 JM changed t range - 20/11/21 SJDK - Changed max value from 1.5 to 0.9 
# SJDK 03/11/21 - New phi_q vs theta_q plot that Steve Wood wanted
phiqvsthq_pions_cut = ROOT.TH2D("phiqvsthq_pions_cut","#phi_{q} vs #theta_{q}; #phi_{q}; #theta_{q}", 12, -3.14, 3.14, 20, 0, 0.1)

P_HGC_xy_npe_pions_uncut = ROOT.TH3D("P_HGC_xy_npe_pions_uncut", "SHMS HGC NPE as fn of yAtCer vs SHMS HGC xAtCer (no cuts); HGC_yAtCer(cm); HGC_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_Aero_xy_npe_pions_uncut = ROOT.TH3D("P_Aero_xy_npe_pions_uncut", "SHMS Aerogel NPE as fn of yAtCer vs xAtCer (no cuts); Aero_yAtCer(cm); Aero_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_NGC_xy_npe_pions_uncut = ROOT.TH3D("P_NGC_xy_npe_pions_uncut", "SHMS NGC NPE as fn of yAtCer vs xAtCer (no cuts); NGC_yAtCer(cm); NGC_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_cal_xy_etottracknorm_pions_uncut = ROOT.TH3D("P_cal_xy_etottracknorm_pions_uncut", "SHMS Calorimeter etottracknorm as fn of yCalo vs xCalo (no cuts); cal_yCalo(cm); cal_xCalo(cm); etottracknorm", 14, -62.3, 62.3, 16, -71.2, 71.2, 100, 0.1 , 50)
P_cal_xy_hits_pions_uncut = ROOT.TH3D("P_cal_xy_hits_pions_uncut", "SHMS Calorimeter total ADC hits as fn of yCalo vs xCalo (no cuts); cal_yCalo(cm); cal_xCalo(cm); Adc_Hits", 14, -62.3, 62.3, 16, -71.2, 71.2, 100, 0 , 2000)
P_HGC_xy_npe_pions_cut = ROOT.TH3D("P_HGC_xy_npe_pions_cut", "SHMS HGC NPE as fn of yAtCer vs SHMS HGC xAtCer (with cuts); HGC_yAtCer(cm); HGC_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_Aero_xy_npe_pions_cut = ROOT.TH3D("P_Aero_xy_npe_pions_cut", "SHMS Aerogel NPE as fn of yAtCer vs xAtCer (with cuts); Aero_yAtCer(cm); Aero_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_NGC_xy_npe_pions_cut = ROOT.TH3D("P_NGC_xy_npe_pions_cut", "SHMS NGC NPE as fn of yAtCer vs xAtCer (with cuts); NGC_yAtCer(cm); NGC_xAtCer(cm); NPE", 100, -50, 50, 100, -50, 50, 100, 0.1 , 50)
P_cal_xy_etottracknorm_pions_cut = ROOT.TH3D("P_cal_xy_etottracknorm_pions_cut", "SHMS Calorimeter etottracknorm as fn of yCalo vs xCalo (with cuts); cal_yCalo(cm); cal_xCalo(cm); etottracknorm", 14, -62.3, 62.3, 16, -71.2, 71.2, 100, 0.1 , 50)
P_cal_xy_hits_pions_cut = ROOT.TH3D("P_cal_xy_hits_pions_cut", "SHMS Calorimeter total ADC hits as fn of yCalo vs xCalo (with cuts); cal_yCalo(cm); cal_xCalo(cm); Adc_Hits", 14, -62.3, 62.3, 16, -71.2, 71.2, 100, 0 , 2000)
#################################################################################################################################################

# Filling Histograms for Pions
for event in Cut_Pion_Events_noRF_tree:
    P_RFTime_pions_cut_noRF.Fill(event.P_RF_Dist)

for event in Uncut_Pion_Events_tree:
    H_xp_pions_uncut.Fill(event.H_gtr_xp)
    H_yp_pions_uncut.Fill(event.H_gtr_yp)
    H_dp_pions_uncut.Fill(event.H_gtr_dp)
    H_cal_etottracknorm_pions_uncut.Fill(event.H_cal_etottracknorm)
    P_xp_pions_uncut.Fill(event.P_gtr_xp)
    P_yp_pions_uncut.Fill(event.P_gtr_yp)
    P_dp_pions_uncut.Fill(event.P_gtr_dp)
    P_cal_etottracknorm_pions_uncut.Fill(event.P_cal_etottracknorm)
    P_hgcer_npe_pions_uncut.Fill(event.P_hgcer_npeSum)
    P_aero_npe_pions_uncut.Fill(event.P_aero_npeSum)
    P_ngcer_npe_pions_uncut.Fill(event.P_ngcer_npeSum)
    P_MMpi_pions_uncut.Fill(event.MMpi)
    P_RFTime_pions_uncut.Fill(event.P_RF_Dist)
    ePiCoinTime_pions_uncut.Fill(event.CTime_ePiCoinTime_ROC1)
    P_hgcer_vs_aero_npe_pions_uncut.Fill(event.P_hgcer_npeSum, event.P_aero_npeSum)
    ePiCoinTime_vs_MMpi_pions_uncut.Fill(event.CTime_ePiCoinTime_ROC1, event.MMpi)
    P_RFTime_vs_MMpi_pions_uncut.Fill(event.P_RF_Dist, event.MMpi)
    ePiCoinTime_vs_beta_pions_uncut.Fill(event.CTime_ePiCoinTime_ROC1, event.P_gtr_beta)
    P_cal_etottracknorm_vs_ngcer_npe_pions_uncut.Fill(event.P_cal_etottracknorm, event.P_ngcer_npeSum)
    P_ngcer_vs_hgcer_npe_pions_uncut.Fill(event.P_ngcer_npeSum, event.P_hgcer_npeSum)
    P_ngcer_vs_aero_npe_pions_uncut.Fill(event.P_ngcer_npeSum, event.P_aero_npeSum)
    H_dp_vs_beta_pions_uncut.Fill(event.H_gtr_dp, event.H_gtr_beta) 
    P_dp_vs_beta_pions_uncut.Fill(event.P_gtr_dp, event.P_gtr_beta) 
    H_xfp_vs_beta_pions_uncut.Fill(event.H_dc_xfp, event.H_gtr_beta)
    P_xfp_vs_beta_pions_uncut.Fill(event.P_dc_xfp, event.P_gtr_beta)
    H_xpfp_vs_beta_pions_uncut.Fill(event.H_dc_xpfp, event.H_gtr_beta)
    P_xpfp_vs_beta_pions_uncut.Fill(event.P_dc_xpfp, event.P_gtr_beta)
    H_yfp_vs_beta_pions_uncut.Fill(event.H_dc_yfp, event.H_gtr_beta)
    P_yfp_vs_beta_pions_uncut.Fill(event.P_dc_yfp, event.P_gtr_beta)
    H_ypfp_vs_beta_pions_uncut.Fill(event.H_dc_ypfp, event.H_gtr_beta)
    P_ypfp_vs_beta_pions_uncut.Fill(event.P_dc_ypfp, event.P_gtr_beta)
    P_MMpi_vs_beta_pions_uncut.Fill(event.MMpi, event.P_gtr_beta)
    P_cal_xy_pions_uncut.Fill(event.yCalo, event.xCalo)
    P_DPexit_xy_pions_uncut.Fill(event.yExit, event.xExit)
    P_HGC_xy_npe_pions_uncut.Fill(event.P_hgcer_yAtCer,event.P_hgcer_xAtCer,event.P_hgcer_npeSum)
    P_Aero_xy_npe_pions_uncut.Fill(event.P_aero_yAtAero,event.P_aero_xAtAero,event.P_aero_npeSum)    
    P_NGC_xy_npe_pions_uncut.Fill(event.P_ngcer_yAtCer,event.P_ngcer_xAtCer,event.P_ngcer_npeSum)
    P_cal_xy_etottracknorm_pions_uncut.Fill(event.yCalo,event.xCalo,event.P_cal_etottracknorm)
    P_cal_xy_hits_pions_uncut.Fill(event.yCalo,event.xCalo,event.Cal_Adc_Hits)

# SJDK - 26/10/21 - Moved filling of kinematic quantity distributions to a different loop (one over tree with MM cut applied)
for event in Cut_Pion_Events_All_tree:
    H_xp_pions_cut.Fill(event.H_gtr_xp)
    H_yp_pions_cut.Fill(event.H_gtr_yp)
    H_dp_pions_cut.Fill(event.H_gtr_dp)
    H_cal_etottracknorm_pions_cut.Fill(event.H_cal_etottracknorm)
    P_xp_pions_cut.Fill(event.P_gtr_xp)
    P_yp_pions_cut.Fill(event.P_gtr_yp)
    P_dp_pions_cut.Fill(event.P_gtr_dp)
    P_cal_etottracknorm_pions_cut.Fill(event.P_cal_etottracknorm)
    P_hgcer_npe_pions_cut.Fill(event.P_hgcer_npeSum)
    P_aero_npe_pions_cut.Fill(event.P_aero_npeSum)
    P_ngcer_npe_pions_cut.Fill(event.P_ngcer_npeSum)
    P_MMpi_pions_cut.Fill(event.MMpi)
    P_RFTime_pions_cut.Fill(event.P_RF_Dist)
    ePiCoinTime_pions_cut.Fill(event.CTime_ePiCoinTime_ROC1)
    H_cal_etottracknorm_vs_cer_npe_pions_cut.Fill(event.H_cal_etottracknorm, event.H_cer_npeSum)    
    P_hgcer_vs_aero_npe_pions_cut.Fill(event.P_hgcer_npeSum, event.P_aero_npeSum)
    ePiCoinTime_vs_beta_pions_cut.Fill(event.CTime_ePiCoinTime_ROC1, event.P_gtr_beta)
    ePiCoinTime_vs_MMpi_pions_cut.Fill(event.CTime_ePiCoinTime_ROC1, event.MMpi)
    P_RFTime_vs_MMpi_pions_cut.Fill(event.P_RF_Dist, event.MMpi)    
    P_cal_etottracknorm_vs_ngcer_npe_pions_cut.Fill(event.P_cal_etottracknorm, event.P_ngcer_npeSum)
    P_ngcer_vs_hgcer_npe_pions_cut.Fill(event.P_ngcer_npeSum, event.P_hgcer_npeSum)
    P_ngcer_vs_aero_npe_pions_cut.Fill(event.P_ngcer_npeSum, event.P_aero_npeSum)
    H_dp_vs_beta_pions_cut.Fill(event.H_gtr_dp, event.H_gtr_beta) 
    P_dp_vs_beta_pions_cut.Fill(event.P_gtr_dp, event.P_gtr_beta) 
    H_xfp_vs_beta_pions_cut.Fill(event.H_dc_xfp, event.H_gtr_beta)
    P_xfp_vs_beta_pions_cut.Fill(event.P_dc_xfp, event.P_gtr_beta)
    H_xpfp_vs_beta_pions_cut.Fill(event.H_dc_xpfp, event.H_gtr_beta)
    P_xpfp_vs_beta_pions_cut.Fill(event.P_dc_xpfp, event.P_gtr_beta)
    H_yfp_vs_beta_pions_cut.Fill(event.H_dc_yfp, event.H_gtr_beta)
    P_yfp_vs_beta_pions_cut.Fill(event.P_dc_yfp, event.P_gtr_beta)
    H_ypfp_vs_beta_pions_cut.Fill(event.H_dc_ypfp, event.H_gtr_beta)
    P_ypfp_vs_beta_pions_cut.Fill(event.P_dc_ypfp, event.P_gtr_beta)
    P_MMpi_vs_beta_pions_cut.Fill(event.MMpi, event.P_gtr_beta)
    P_cal_xy_pions_cut.Fill(event.yCalo, event.xCalo)
    P_DPexit_xy_pions_cut.Fill(event.yExit, event.xExit)
    P_HGC_xy_npe_pions_cut.Fill(event.P_hgcer_yAtCer,event.P_hgcer_xAtCer,event.P_hgcer_npeSum)
    P_Aero_xy_npe_pions_cut.Fill(event.P_aero_yAtAero,event.P_aero_xAtAero,event.P_aero_npeSum)
    P_NGC_xy_npe_pions_cut.Fill(event.P_ngcer_yAtCer,event.P_ngcer_xAtCer,event.P_ngcer_npeSum)
    P_cal_xy_etottracknorm_pions_cut.Fill(event.yCalo,event.xCalo,event.P_cal_etottracknorm)
    P_cal_xy_hits_pions_cut.Fill(event.yCalo,event.xCalo,event.Cal_Adc_Hits)

for event in Cut_Pion_Events_Prompt_tree:
    ePiCoinTime_pions_cut_prompt.Fill(event.CTime_ePiCoinTime_ROC1)
    P_MMpi_pions_cut_prompt.Fill(event.MMpi)

for event in Cut_Pion_Events_Random_tree:
    ePiCoinTime_pions_cut_randm.Fill(event.CTime_ePiCoinTime_ROC1)
    P_MMpi_pions_cut_randm.Fill(event.MMpi)

# SJDK 26/10/21 - For kinematic quantities, fill the histograms from the tree with the pion MM cut
for event in Cut_Pion_Events_Prompt_MM_tree:
    phiq_plot.Fill(event.ph_q)
    t_plot.Fill(-event.MandelT) 
    Q2_pions_cut.Fill(event.Q2)
    W_pions_cut.Fill(event.W)
    epsilon_pions_cut.Fill(event.epsilon)
    Q2vsW_pions_cut.Fill(event.Q2, event.W)
    phiqvst_pions_cut.Fill(event.ph_q, -event.MandelT)
    phiqvsthq_pions_cut.Fill(event.ph_q, event.th_q)

print("Histograms filled")

##############################################################################################################################################

# Random subtraction from missing mass and coin_Time
for event in Cut_Pion_Events_Random_tree:
    P_MMpi_pions_cut_randm_scaled.Fill(event.MMpi)
    P_MMpi_pions_cut_randm_scaled.Scale(1.0/nWindows)
P_MMpi_pions_cut_randm_sub.Add(P_MMpi_pions_cut_prompt, P_MMpi_pions_cut_randm_scaled, 1, -1)

###########################################################################################################################################

# HGC/NGC/Aero XY Projection vs npe for pions.
HGC_proj_yx_pions_uncut = ROOT.TProfile2D(P_HGC_xy_npe_pions_uncut.Project3DProfile("yx"))
HGC_proj_yx_pions_cut = ROOT.TProfile2D(P_HGC_xy_npe_pions_cut.Project3DProfile("yx"))
NGC_proj_yx_pions_uncut = ROOT.TProfile2D(P_NGC_xy_npe_pions_uncut.Project3DProfile("yx"))
NGC_proj_yx_pions_cut = ROOT.TProfile2D(P_NGC_xy_npe_pions_cut.Project3DProfile("yx"))
Aero_proj_yx_pions_uncut = ROOT.TProfile2D(P_Aero_xy_npe_pions_uncut.Project3DProfile("yx"))
Aero_proj_yx_pions_cut = ROOT.TProfile2D(P_Aero_xy_npe_pions_cut.Project3DProfile("yx"))


Calo_proj_yx_pions_uncut = ROOT.TProfile2D(P_cal_xy_etottracknorm_pions_uncut.Project3DProfile("yx"))
Calo_proj_yx_pions_cut = ROOT.TProfile2D(P_cal_xy_etottracknorm_pions_cut.Project3DProfile("yx"))
Calo_proj_hits_yx_pions_uncut = ROOT.TProfile2D(P_cal_xy_hits_pions_uncut.Project3DProfile("yx"))
Calo_proj_hits_yx_pions_cut = ROOT.TProfile2D(P_cal_xy_hits_pions_cut.Project3DProfile("yx"))

############################################################################################################################################

# Saving histograms in PDF
c1_pions_kin = TCanvas("c1_pions_kin", "Pions Kinematic Distributions", 100, 0, 1000, 900)
c1_pions_kin.Divide(2,2)
c1_pions_kin.cd(1)
Q2vsW_pions_cut.Draw("COLZ")
c1_pions_kin.cd(2)
epsilon_pions_cut.Draw()
c1_pions_kin.cd(3)
phiqvst_pions_cut.SetStats(0)
phiqvst_pions_cut.GetYaxis().SetRangeUser(minrangeuser,maxrangeuser)
phiqvst_pions_cut.Draw("SURF2 POL")
# Section for polar plotting
gStyle.SetPalette(55)
gPad.SetTheta(90)
gPad.SetPhi(180)
tvsphi_title_pions = TPaveText(0.0277092,0.89779,0.096428,0.991854,"NDC")
tvsphi_title_pions.AddText("-t vs #phi")
tvsphi_title_pions.Draw()
ptphizero_pions = TPaveText(0.923951,0.513932,0.993778,0.574551,"NDC")
ptphizero_pions.AddText("#phi = 0")
ptphizero_pions.Draw()
phihalfpi_pions = TLine(0,0,0,0.6)
phihalfpi_pions.SetLineColor(kBlack)
phihalfpi_pions.SetLineWidth(2)
phihalfpi_pions.Draw()
ptphihalfpi_pions = TPaveText(0.417855,0.901876,0.486574,0.996358,"NDC")
ptphihalfpi_pions.AddText("#phi = #frac{#pi}{2}")
ptphihalfpi_pions.Draw()
phipi_pions = TLine(0,0,-0.6,0)
phipi_pions.SetLineColor(kBlack)
phipi_pions.SetLineWidth(2)
phipi_pions.Draw()
ptphipi_pions = TPaveText(0.0277092,0.514217,0.096428,0.572746,"NDC")
ptphipi_pions.AddText("#phi = #pi")
ptphipi_pions.Draw()
phithreepi_pions = TLine(0,0,0,-0.6)
phithreepi_pions.SetLineColor(kBlack)
phithreepi_pions.SetLineWidth(2)
phithreepi_pions.Draw()
ptphithreepi_pions = TPaveText(0.419517,0.00514928,0.487128,0.0996315,"NDC")
ptphithreepi_pions.AddText("#phi = #frac{3#pi}{2}")
ptphithreepi_pions.Draw()
Arc_pions = TArc()
for k in range(0, 6):
     Arc_pions.SetFillStyle(0)
     Arc_pions.SetLineWidth(2)
     # To change the arc radius we have to change number 0.825 in the lower line.
     Arc_pions.DrawArc(0,0,0.95*(k+1)/(10),0.,360.,"same")
tradius_pions = TGaxis(0,0,0.575,0,minrangeuser,maxrangeuser,6,"-+N") # NH 2021 08 12 - added "N" option which forces there to be 6 divisions, meaning we should be able to change the range of this plot without having to fiddle with the drawing of the arcs!
tradius_pions.SetLineColor(2)
tradius_pions.SetLabelColor(2)
tradius_pions.Draw()
phizero_pions = TLine(0,0,0.6,0)
phizero_pions.SetLineColor(kBlack)
phizero_pions.SetLineWidth(2)
phizero_pions.Draw()
# End of polar plotting section
c1_pions_kin.cd(4)
P_MMpi_pions_cut_randm_sub.Draw("hist")
# Section for Neutron Peak Events Selection
shadedpeak_pions = P_MMpi_pions_cut_randm_sub.Clone()
shadedpeak_pions.SetFillColor(2)
shadedpeak_pions.SetFillStyle(3244)
shadedpeak_pions.GetXaxis().SetRangeUser(minbin, maxbin)
shadedpeak_pions.Draw("samehist")
NeutronEvt_pions = TPaveText(0.58934,0.675,0.95,0.75,"NDC")
BinLow_pions = P_MMpi_pions_cut_randm_sub.GetXaxis().FindBin(minbin)
BinHigh_pions = P_MMpi_pions_cut_randm_sub.GetXaxis().FindBin(maxbin)
BinIntegral_pions = int(P_MMpi_pions_cut_randm_sub.Integral(BinLow_pions, BinHigh_pions))
NeutronEvt_pions.SetLineColor(2)
NeutronEvt_pions.AddText("e #pi n Events: %i" %(BinIntegral_pions))
NeutronEvt_pions.Draw()
# End of Neutron Peak Events Selection Section
c1_pions_kin.Print(Pion_Analysis_Distributions + '(')

c1_pion_kin_pg2 = TCanvas("c1_pions_kin_pg2", "Pion Kinematic Distributions part 2", 100, 0, 1000, 900)
c1_pion_kin_pg2.Divide(2,2)
c1_pion_kin_pg2.cd(1)
phiq_plot.SetMinimum(0) # SJDK 02/11/21 - Added to change the autoscaling of the plot
phiq_plot.Draw("Hist")
c1_pion_kin_pg2.cd(2)
t_plot.Draw("Hist")
c1_pion_kin_pg2.cd(3)
Q2_pions_cut.Draw("Hist")
c1_pion_kin_pg2.cd(4)
W_pions_cut.Draw("Hist")
c1_pion_kin_pg2.Print(Pion_Analysis_Distributions)

c1_pions_acpt = TCanvas("c1_pions_acpt", "Electron-Pion Acceptance Distributions", 100, 0, 1000, 900)
c1_pions_acpt.Divide(3,2)
c1_pions_acpt.cd(1)
gPad.SetLogy()
H_xp_pions_uncut.SetLineColor(2)
H_xp_pions_uncut.Draw()
H_xp_pions_cut.SetLineColor(4)
H_xp_pions_cut.Draw("same")
c1_pions_acpt.cd(2)
gPad.SetLogy()
H_yp_pions_uncut.SetLineColor(2)
H_yp_pions_uncut.Draw()
H_yp_pions_cut.SetLineColor(4)
H_yp_pions_cut.Draw("same")
c1_pions_acpt.cd(3)
gPad.SetLogy()
H_dp_pions_uncut.SetMinimum(0.1*H_dp_pions_cut.GetMinimum()+1) # min of plot should be one order of magnitude below the min bin in cut distribution
H_dp_pions_uncut.SetMaximum(10*H_dp_pions_uncut.GetBinContent(H_dp_pions_uncut.GetMaximumBin())) # Max of plot should be 1 order of magnitude greater than the max bin in uncut distribution
H_dp_pions_uncut.SetLineColor(2)
H_dp_pions_cut.SetLineColor(4)
H_dp_pions_uncut.Draw()
H_dp_pions_cut.Draw("same")
# TLegend (x1, y1, x2, y2) 
legend2_pions = ROOT.TLegend(0.115, 0.8, 0.6, 0.9)
legend2_pions.AddEntry("H_dp_pions_uncut", "without cuts", "l")
legend2_pions.AddEntry("H_dp_pions_cut", "with cuts (acpt/RF/PID)", "l")
legend2_pions.Draw("same")
c1_pions_acpt.cd(4)
gPad.SetLogy()
P_xp_pions_uncut.SetLineColor(2)
P_xp_pions_uncut.Draw()
P_xp_pions_cut.SetLineColor(4)
P_xp_pions_cut.Draw("same")
c1_pions_acpt.cd(5)
gPad.SetLogy()
P_yp_pions_uncut.SetLineColor(2)
P_yp_pions_uncut.Draw()
P_yp_pions_cut.SetLineColor(4)
P_yp_pions_cut.Draw("same")
c1_pions_acpt.cd(6)
gPad.SetLogy()
P_dp_pions_uncut.SetMinimum(0.1*P_dp_pions_cut.GetMinimum()+1) # Implemented same fixed as used above for HMS
P_dp_pions_uncut.SetMaximum(10*P_dp_pions_uncut.GetBinContent(P_dp_pions_uncut.GetMaximumBin()))
P_dp_pions_uncut.SetLineColor(2)
P_dp_pions_uncut.Draw()
P_dp_pions_cut.SetLineColor(4)
P_dp_pions_cut.Draw("same")
c1_pions_acpt.Print(Pion_Analysis_Distributions)

c1_pions_pid = TCanvas("c1_pions_pid", "Electron-Pion CAL Distributions", 100, 0, 1000, 900)
c1_pions_pid.Divide(2,2)
c1_pions_pid.cd(1)
gPad.SetLogy()
H_cal_etottracknorm_pions_uncut.SetLineColor(2)
H_cal_etottracknorm_pions_uncut.Draw()
H_cal_etottracknorm_pions_cut.SetLineColor(4)
H_cal_etottracknorm_pions_cut.Draw("same")
legend7_pions = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
legend7_pions.AddEntry("H_cal_etottracknorm_pions_uncut", "without cuts", "l")
legend7_pions.AddEntry("H_cal_etottracknorm_pions_cut", "with cuts (acpt/RF/PID)", "l")
legend7_pions.Draw("same")
c1_pions_pid.cd(2)
H_cal_etottracknorm_vs_cer_npe_pions_cut.Draw("COLZ")
c1_pions_pid.cd(3)
gPad.SetLogy()
P_cal_etottracknorm_pions_uncut.SetLineColor(2)
P_cal_etottracknorm_pions_uncut.Draw()
P_cal_etottracknorm_pions_cut.SetLineColor(4)
P_cal_etottracknorm_pions_cut.Draw("same")
legend8_pions = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
legend8_pions.AddEntry("P_cal_etottracknorm_pions_uncut", "without cuts", "l")
legend8_pions.AddEntry("P_cal_etottracknorm_pions_cut", "with cuts (acpt/RF/PID)", "l")
legend8_pions.Draw("same")
c1_pions_pid.cd(4)
c1_pions_pid.Print(Pion_Analysis_Distributions)

c1_pions_RF = TCanvas("c1_pions_RF", "Electron-Pion RF Distributions", 100, 0, 1000, 900)
c1_pions_RF.Divide(2,2)
c1_pions_RF.cd(1)
P_RFTime_pions_uncut.SetLineColor(2)
P_RFTime_pions_uncut.Draw()
P_RFTime_pions_cut.SetLineColor(4)
P_RFTime_pions_cut.Draw("same")
legend9_pions = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
legend9_pions.AddEntry("P_RFTime_pions_uncut", "without cuts", "l")
legend9_pions.AddEntry("P_RFTime_pions_cut", "with cuts (acpt/RF/PID)", "l")
legend9_pions.Draw("same")
c1_pions_RF.cd(2)
P_RFTime_pions_cut_noRF.SetLineColor(2)
P_RFTime_pions_cut_noRF.Draw()
P_RFTime_pions_cut.SetLineColor(4)
P_RFTime_pions_cut.Draw("same")
legend_pions = ROOT.TLegend(0.115, 0.835, 0.415, 0.9)
legend_pions.AddEntry("P_RFTime_pions_cut_noRF", "noRF_cuts (acpt/PID)", "l")
legend_pions.AddEntry("P_RFTime_pions_cut", "RF_cuts (acpt/PID)", "l")
legend_pions.Draw("same")
c1_pions_RF.cd(3)
P_RFTime_vs_MMpi_pions_uncut.Draw("COLZ")
c1_pions_RF.cd(4)
P_RFTime_vs_MMpi_pions_cut.Draw("COLZ")
c1_pions_RF.Print(Pion_Analysis_Distributions)

c2_pions_pid = TCanvas("c2_pions_pid", "Electron-Pion Aero/HGC/NGC PID Distributions", 100, 0, 1000, 900)
c2_pions_pid.Divide(2,2)
c2_pions_pid.cd(1)
gPad.SetLogy()
P_hgcer_npe_pions_uncut.SetLineColor(2)
P_hgcer_npe_pions_uncut.Draw()
P_hgcer_npe_pions_cut.SetLineColor(4)
P_hgcer_npe_pions_cut.Draw("same")
legend10_pions = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
legend10_pions.AddEntry("P_hgcer_npe_pions_uncut", "without cuts", "l")
legend10_pions.AddEntry("P_hgcer_npe_pions_cut", "with cuts (acpt/RF/PID)", "l")
legend10_pions.Draw("same")
c2_pions_pid.cd(2)
gPad.SetLogy()
P_aero_npe_pions_uncut.SetLineColor(2)
P_aero_npe_pions_uncut.Draw()
P_aero_npe_pions_cut.SetLineColor(4)
P_aero_npe_pions_cut.Draw("same")
legend11_pions = ROOT.TLegend(0.115, 0.835, 0.43, 0.9)
legend11_pions.AddEntry("P_aero_npe_pions_uncut", "without cuts", "l")
legend11_pions.AddEntry("P_aero_npe_pions_cut", "with cuts (acpt/RF/PID)", "l")
legend11_pions.Draw("same")
c2_pions_pid.cd(3)
gPad.SetLogy()
P_ngcer_npe_pions_uncut.SetLineColor(2)
P_ngcer_npe_pions_uncut.Draw()
P_ngcer_npe_pions_cut.SetLineColor(4)
P_ngcer_npe_pions_cut.Draw("same")
c2_pions_pid.cd(4)
#
c2_pions_pid.Print(Pion_Analysis_Distributions)

c3_pions_pid = TCanvas("c3_pions_pid", "Electron-Pion Aero/HGC/NGC PID Distributions", 100, 0, 1000, 900)
c3_pions_pid.Divide(2,3)
c3_pions_pid.cd(1)
gPad.SetLogz()
P_hgcer_vs_aero_npe_pions_uncut.Draw("COLZ")
c3_pions_pid.cd(2)
gPad.SetLogz()
P_hgcer_vs_aero_npe_pions_cut.Draw("COLZ")
c3_pions_pid.cd(3)
gPad.SetLogz()
P_ngcer_vs_hgcer_npe_pions_uncut.Draw("COLZ")
c3_pions_pid.cd(4)
gPad.SetLogz()
P_ngcer_vs_hgcer_npe_pions_cut.Draw("COLZ")
c3_pions_pid.cd(5)
gPad.SetLogz()
P_ngcer_vs_aero_npe_pions_uncut.Draw("COLZ")
c3_pions_pid.cd(6)
gPad.SetLogz()
P_ngcer_vs_aero_npe_pions_cut.Draw("COLZ")
c3_pions_pid.Print(Pion_Analysis_Distributions)

c1_pions_MM = TCanvas("c1_pions_MM", "Electron-Pion CTime/Missing Mass Distributions", 100, 0, 1000, 900)
c1_pions_MM.Divide(3,2)
c1_pions_MM.cd(1)
ePiCoinTime_pions_uncut.SetLineColor(4)
ePiCoinTime_pions_uncut.Draw()
ePiCoinTime_pions_cut_prompt.SetLineColor(6)
ePiCoinTime_pions_cut_prompt.Draw("same")
ePiCoinTime_pions_cut_randm.SetLineColor(8)
ePiCoinTime_pions_cut_randm.Draw("same")
legend13_pions = ROOT.TLegend(0.1, 0.815, 0.48, 0.9)
legend13_pions.AddEntry("ePiCoinTime_pions_uncut", "CT_without cuts", "l")
legend13_pions.AddEntry("ePiCoinTime_pions_cut_prompt", "CT_prompt with cuts (acpt/RF/PID)", "l")
legend13_pions.AddEntry("ePiCoinTime_pions_cut_randm", "CT_randoms with cuts (acpt/RF/PID)", "l")
legend13_pions.Draw("same")
c1_pions_MM.cd(2)
ePiCoinTime_pions_cut.SetLineColor(4)
ePiCoinTime_pions_cut.Draw()
ePiCoinTime_pions_cut_prompt.SetLineColor(6)
ePiCoinTime_pions_cut_prompt.Draw("same")
ePiCoinTime_pions_cut_randm.SetLineColor(8)
ePiCoinTime_pions_cut_randm.Draw("same")
legend13_pions = ROOT.TLegend(0.1, 0.815, 0.48, 0.9)
legend13_pions.AddEntry("ePiCoinTime_pions_uncut", "CT_without cuts", "l")
legend13_pions.AddEntry("ePiCoinTime_pions_cut_prompt", "CT_prompt with cuts (acpt/RF/PID)", "l")
legend13_pions.AddEntry("ePiCoinTime_pions_cut_randm", "CT_randoms with cuts (acpt/RF/PID)", "l")
legend13_pions.Draw("same")
c1_pions_MM.cd(2)
ePiCoinTime_pions_cut.SetLineColor(4)
ePiCoinTime_pions_cut.Draw()
ePiCoinTime_pions_cut_prompt.SetLineColor(6)
ePiCoinTime_pions_cut_prompt.Draw("same")
ePiCoinTime_pions_cut_randm.SetLineColor(8)
ePiCoinTime_pions_cut_randm.Draw("same")
legend14_pions = ROOT.TLegend(0.1, 0.815, 0.48, 0.9)
legend14_pions.AddEntry("ePiCoinTime_pions_cut", "CT_with cuts (acpt/RF/PID)", "l")
legend14_pions.AddEntry("ePiCoinTime_pions_cut_prompt", "CT_prompt with cuts (acpt/RF/PID)", "l")
legend14_pions.AddEntry("ePiCoinTime_pions_cut_randm", "CT_randoms with cuts (acpt/RF/PID)", "l")
legend14_pions.Draw("same")
c1_pions_MM.cd(3)
P_MMpi_pions_uncut.Draw()
c1_pions_MM.cd(4)
P_MMpi_pions_cut.SetLineColor(4)
P_MMpi_pions_cut.Draw()
P_MMpi_pions_cut_prompt.SetLineColor(6)
P_MMpi_pions_cut_prompt.Draw("same")
P_MMpi_pions_cut_randm.SetLineColor(8)
P_MMpi_pions_cut_randm.Draw("same")
legend15_pions = ROOT.TLegend(0.4, 0.815, 0.78, 0.9)
legend15_pions.AddEntry("P_MMpi_pions_cut", "MM with cuts (acpt/RF/PID)", "l")
legend15_pions.AddEntry("P_MMpi_pions_cut_prompt", "MM_prompt with cuts (acpt/RF/PID)", "l")
legend15_pions.AddEntry("P_MMpi_pions_cut_randm", "MM_randoms with cuts (acpt/RF/PID)", "l")
legend15_pions.Draw("same")
c1_pions_MM.cd(5)
P_MMpi_vs_beta_pions_uncut.Draw("COLZ")
c1_pions_MM.cd(6)
P_MMpi_vs_beta_pions_cut.Draw("COLZ")
c1_pions_MM.Print(Pion_Analysis_Distributions)
c1_pions_CT = TCanvas("c1_pions_CT", "Electron-Pion CTime Distributions", 100, 0, 1000, 900)
c1_pions_CT.Divide(2,2)
c1_pions_CT.cd(1)
ePiCoinTime_vs_beta_pions_uncut.Draw("COLZ")
# TLine (x1, y1, x2, y2)
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
c1_pions_CT.cd(2)
ePiCoinTime_vs_beta_pions_cut.Draw("COLZ")
c1_pions_CT.cd(3)
ePiCoinTime_vs_MMpi_pions_uncut.Draw("COLZ")
LowerPrompt2_pions = TLine(PromptWindow[0],gPad.GetUymin(),PromptWindow[0],2)
LowerPrompt2_pions.SetLineColor(2)
LowerPrompt2_pions.SetLineWidth(2)
LowerPrompt2_pions.Draw("same")
UpperPrompt2_pions = TLine(PromptWindow[1],gPad.GetUymin(),PromptWindow[1],2)
UpperPrompt2_pions.SetLineColor(2)
UpperPrompt2_pions.SetLineWidth(2)
UpperPrompt2_pions.Draw("same")
LowerRandomL2_pions = TLine(RandomWindows[0],gPad.GetUymin(),RandomWindows[0],2)
LowerRandomL2_pions.SetLineColor(8)
LowerRandomL2_pions.SetLineWidth(2)
LowerRandomL2_pions.Draw("same")
UpperRandomL2_pions = TLine(RandomWindows[1],gPad.GetUymin(),RandomWindows[1],2)
UpperRandomL2_pions.SetLineColor(8)
UpperRandomL2_pions.SetLineWidth(2)
UpperRandomL2_pions.Draw("same")
LowerRandomR2_pions = TLine(RandomWindows[2],gPad.GetUymin(),RandomWindows[2],2)
LowerRandomR2_pions.SetLineColor(8)
LowerRandomR2_pions.SetLineWidth(2)
LowerRandomR2_pions.Draw("same")
UpperRandomR2_pions = TLine(RandomWindows[3],gPad.GetUymin(),RandomWindows[3],2)
UpperRandomR2_pions.SetLineColor(8)
UpperRandomR2_pions.SetLineWidth(2)
UpperRandomR2_pions.Draw("same")
c1_pions_CT.cd(4)
ePiCoinTime_vs_MMpi_pions_cut.Draw("COLZ")
c1_pions_CT.Print(Pion_Analysis_Distributions)

c1_pions_delta = TCanvas("c1_pions_delta", "Delta Debugging", 100, 0, 1000, 900)
c1_pions_delta.Divide(2,2)
c1_pions_delta.cd(1)
H_dp_vs_beta_pions_uncut.Draw("COLZ")
c1_pions_delta.cd(2)
H_dp_vs_beta_pions_cut.Draw("COLZ")
c1_pions_delta.cd(3)
P_dp_vs_beta_pions_uncut.Draw("COLZ")
c1_pions_delta.cd(4)
P_dp_vs_beta_pions_cut.Draw("COLZ")
c1_pions_delta.Print(Pion_Analysis_Distributions)

c1_pions_fpH = TCanvas("c1_pions_fpH", "HMS Focal Plane and Beta", 100, 0, 1000, 900)
c1_pions_fpH.Divide(2,2)
c1_pions_fpH.cd(1)
H_xfp_vs_beta_pions_uncut.Draw("COLZ")
c1_pions_fpH.cd(2)
H_xfp_vs_beta_pions_cut.Draw("COLZ")
c1_pions_fpH.cd(3)
H_yfp_vs_beta_pions_uncut.Draw("COLZ")
c1_pions_fpH.cd(4)
H_yfp_vs_beta_pions_cut.Draw("COLZ")
c1_pions_fpH.Print(Pion_Analysis_Distributions)

c1_pions_fpP = TCanvas("c1_pions_fpP", "SHMS Focal Plane and Beta", 100, 0, 1000, 900)
c1_pions_fpP.Divide(2,2)
c1_pions_fpP.cd(1)
P_xfp_vs_beta_pions_uncut.Draw("COLZ")
c1_pions_fpP.cd(2)
P_xfp_vs_beta_pions_cut.Draw("COLZ")
c1_pions_fpP.cd(3)
P_yfp_vs_beta_pions_uncut.Draw("COLZ")
c1_pions_fpP.cd(4)
P_yfp_vs_beta_pions_cut.Draw("COLZ")
c1_pions_fpP.Print(Pion_Analysis_Distributions)

c1_pions_cal_proj = TCanvas("c1_pions_cal_proj", "SHMS Cal XY Projections", 100, 0, 1000, 900)
c1_pions_cal_proj.Divide(2,3)
c1_pions_cal_proj.cd(1)
P_cal_xy_pions_uncut.Draw("COLZ")
c1_pions_cal_proj.cd(2)
P_cal_xy_pions_cut.Draw("COLZ")
c1_pions_cal_proj.cd(3)
Calo_proj_yx_pions_uncut.Draw("COLZ")
c1_pions_cal_proj.cd(4)
Calo_proj_yx_pions_cut.Draw("COLZ")
c1_pions_cal_proj.cd(5)
Calo_proj_hits_yx_pions_uncut.Draw("COLZ")
c1_pions_cal_proj.cd(6)
Calo_proj_hits_yx_pions_cut.Draw("COLZ")
c1_pions_cal_proj.Print(Pion_Analysis_Distributions)


c1_pions_dp_proj = TCanvas("c1_pions_dp_proj", "SHMS Dipole Exit XY Projections", 100, 0, 1000, 900)
c1_pions_dp_proj.Divide(1,2)
c1_pions_dp_proj.cd(1)
P_DPexit_xy_pions_uncut.Draw("COLZ")
c1_pions_dp_proj.cd(2)
P_DPexit_xy_pions_cut.Draw("COLZ")
c1_pions_dp_proj.Print(Pion_Analysis_Distributions)

c1_pions_proj = TCanvas("c1_pions_proj", "HGC/NGC/Aero XY Projection", 100, 0, 1000, 900)
c1_pions_proj.Divide(2,3)
c1_pions_proj.cd(1)
HGC_proj_yx_pions_uncut.Draw("COLZ")
c1_pions_proj.cd(2)
HGC_proj_yx_pions_cut.Draw("COLZ")
c1_pions_proj.cd(3)
NGC_proj_yx_pions_uncut.Draw("COLZ")
c1_pions_proj.cd(4)
NGC_proj_yx_pions_cut.Draw("COLZ")
c1_pions_proj.cd(5)
Aero_proj_yx_pions_uncut.Draw("COLZ")
c1_pions_proj.cd(6)
Aero_proj_yx_pions_cut.Draw("COLZ")
c1_pions_proj.Print(Pion_Analysis_Distributions + ')')

# SJDK - 02/11/21 - Just some stuff I've added in a Steve Wood's request, I'll delete this once I'm done!
# c_PlotsForSteve = TCanvas("c_PlotsForSteve", "Phiq Plot", 100, 0, 2560, 1440)
# phiq_plot.SetMinimum(0) # SJDK 02/11/21 - Added to change the autoscaling of the plot
# phiq_plot.Draw("Hist")
# c_PlotsForSteve.Print("%s/Q2_8p5_Phiq_Plot.pdf" % (OUTPATH))
# c_PlotsForSteve.Print("%s/Q2_8p5_Phiq_Plot.png" % (OUTPATH))
# c_PlotsForSteve2 = TCanvas("c_PlotsForSteve2", "Phiq vs Thetaq Plot", 100, 0, 2560, 1440)
# phiqvsthq_pions_cut.SetStats(0)
# phiqvsthq_pions_cut.Draw("COLZ")
# c_PlotsForSteve2.Print("%s/Q2_8p5_Phiq_Thetaq_Plot.pdf" % (OUTPATH))
# c_PlotsForSteve2.Print("%s/Q2_8p5_Phiq_Thetaq_Plot.png" % (OUTPATH))
# c_PlotsForSteve3 = TCanvas("c_PlotsForSteve3", "Phiq vs Thetaq Plot", 100, 0, 2560, 1440)
# phiqvsthq_pions_cut.SetStats(0)
# phiqvsthq_pions_cut.GetYaxis().SetRangeUser(0,0.1)
# phiqvsthq_pions_cut.Draw("SURF2 POL")
# # Section for polar plotting
# gStyle.SetPalette(55)
# gPad.SetTheta(90)
# gPad.SetPhi(180)
# thetavsphi_title_pions = TPaveText(0.0277092,0.89779,0.096428,0.991854,"NDC")
# thetavsphi_title_pions.AddText("#theta_{q} vs #phi_{q}")
# thetavsphi_title_pions.Draw()
# ptphizero_pions.Draw()
# phihalfpi_pions.Draw()
# ptphihalfpi_pions.Draw()
# phipi_pions.Draw()
# ptphipi_pions.Draw()
# phithreepi_pions.Draw()
# ptphithreepi_pions.Draw()
# Arc_pions_2 = TArc()
# for k in range(0, 5):
#     Arc_pions_2.SetFillStyle(0)
#     Arc_pions_2.SetLineWidth(2)
#     Arc_pions_2.DrawArc(0,0,0.95*(k+1)/(10),0.,360.,"same")
# tradius_pions_2 = TGaxis(0,0,0.575,0,0,0.1,5,"-+N") # NH 2021 08 12 - added "N" option which forces there to be 5 divisions, meaning we should be able to change the range of this plot without having to fiddle with the drawing of the arcs!
# tradius_pions_2.SetLineColor(2)
# tradius_pions_2.SetLabelColor(2)
# tradius_pions_2.Draw()
# phizero_pions.Draw()

# c_PlotsForSteve3.Print("%s/Q2_8p5_Phiq_Thetaq_Plot_Polar.pdf" % (OUTPATH))
# c_PlotsForSteve3.Print("%s/Q2_8p5_Phiq_Thetaq_Plot_Polar.png" % (OUTPATH))

#############################################################################################################################################

# Making directories in output file
if (FilenameOverride == False): # Standard running condition, construct file name from run number and max events e.t.c.
       outHistFile = ROOT.TFile.Open("%s/%s_%s_Output_Data.root" % (OUTPATH, runNum, MaxEvent), "RECREATE")
elif (FilenameOverride != False): # Special condition, with 4th arg, use 4th arg as file name
       outHistFile = ROOT.TFile.Open("%s/%s_%s_Output_Data.root" % (OUTPATH, (FilenameOverride.split("_Analysed_Data.root",1)[0]), MaxEvent), "RECREATE")

d_Uncut_Pion_Events = outHistFile.mkdir("Uncut_Pion_Events")
d_Cut_Pion_Events_All = outHistFile.mkdir("Cut_Pion_Events_All")
d_Cut_Pion_Events_Prompt = outHistFile.mkdir("Cut_Pion_Events_Prompt")
d_Cut_Pion_Events_Random = outHistFile.mkdir("Cut_Pion_Events_Random")

# Writing Histograms in output root file
d_Uncut_Pion_Events.cd()
H_xp_pions_uncut.Write()
H_yp_pions_uncut.Write()
H_dp_pions_uncut.Write()
H_cal_etottracknorm_pions_uncut.Write()
P_xp_pions_uncut.Write()
P_yp_pions_uncut.Write()
P_dp_pions_uncut.Write()
P_cal_etottracknorm_pions_uncut.Write()
P_hgcer_npe_pions_uncut.Write()
P_aero_npe_pions_uncut.Write()
P_ngcer_npe_pions_uncut.Write()
P_MMpi_pions_uncut.Write()
P_RFTime_pions_uncut.Write()
P_RFTime_pions_cut_noRF.Write()
ePiCoinTime_pions_uncut.Write()
P_hgcer_vs_aero_npe_pions_uncut.Write()
ePiCoinTime_vs_MMpi_pions_uncut.Write()
ePiCoinTime_vs_beta_pions_uncut.Write()
P_RFTime_vs_MMpi_pions_uncut.Write()
P_cal_etottracknorm_vs_ngcer_npe_pions_uncut.Write()
P_ngcer_vs_hgcer_npe_pions_uncut.Write()
P_ngcer_vs_aero_npe_pions_uncut.Write()
H_dp_vs_beta_pions_uncut.Write()
P_dp_vs_beta_pions_uncut.Write()
H_xfp_vs_beta_pions_uncut.Write()
H_xpfp_vs_beta_pions_uncut.Write()
P_xfp_vs_beta_pions_uncut.Write()
P_xpfp_vs_beta_pions_uncut.Write()
H_yfp_vs_beta_pions_uncut.Write()
H_ypfp_vs_beta_pions_uncut.Write()
P_yfp_vs_beta_pions_uncut.Write()
P_ypfp_vs_beta_pions_uncut.Write()
P_MMpi_vs_beta_pions_uncut.Write()
P_cal_xy_pions_uncut.Write()
HGC_proj_yx_pions_uncut.Write()
NGC_proj_yx_pions_uncut.Write()
Aero_proj_yx_pions_uncut.Write()
P_HGC_xy_npe_pions_uncut.Write()
P_Aero_xy_npe_pions_uncut.Write()
P_NGC_xy_npe_pions_uncut.Write()
Calo_proj_yx_pions_uncut.Write()
Calo_proj_hits_yx_pions_uncut.Write()
P_DPexit_xy_pions_uncut.Write()

d_Cut_Pion_Events_All.cd()
H_xp_pions_cut.Write()
H_yp_pions_cut.Write()
H_dp_pions_cut.Write()
H_cal_etottracknorm_pions_cut.Write()
P_xp_pions_cut.Write()
P_yp_pions_cut.Write()
P_dp_pions_cut.Write()
P_cal_etottracknorm_pions_cut.Write()
P_hgcer_npe_pions_cut.Write()
P_aero_npe_pions_cut.Write()
P_ngcer_npe_pions_cut.Write()
P_MMpi_pions_cut.Write()
P_RFTime_pions_cut.Write()
ePiCoinTime_pions_cut.Write()
epsilon_pions_cut.Write()
Q2vsW_pions_cut.Write()
phiqvst_pions_cut.Write()
H_cal_etottracknorm_vs_cer_npe_pions_cut.Write()
P_hgcer_vs_aero_npe_pions_cut.Write()
ePiCoinTime_vs_MMpi_pions_cut.Write()
ePiCoinTime_vs_beta_pions_cut.Write()
P_RFTime_vs_MMpi_pions_cut.Write()
P_MMpi_pions_cut_randm_sub.Write()
P_cal_etottracknorm_vs_ngcer_npe_pions_cut.Write()
P_ngcer_vs_hgcer_npe_pions_cut.Write()
P_ngcer_vs_aero_npe_pions_cut.Write()
H_dp_vs_beta_pions_cut.Write()
P_dp_vs_beta_pions_cut.Write()
H_xfp_vs_beta_pions_cut.Write()
H_xpfp_vs_beta_pions_cut.Write()
P_xfp_vs_beta_pions_cut.Write()
P_xpfp_vs_beta_pions_cut.Write()
H_yfp_vs_beta_pions_cut.Write()
H_ypfp_vs_beta_pions_cut.Write()
P_yfp_vs_beta_pions_cut.Write()
P_ypfp_vs_beta_pions_cut.Write()
P_MMpi_vs_beta_pions_cut.Write()
P_cal_xy_pions_cut.Write()
P_HGC_xy_npe_pions_cut.Write()
P_Aero_xy_npe_pions_cut.Write()
P_NGC_xy_npe_pions_cut.Write()
HGC_proj_yx_pions_cut.Write()
NGC_proj_yx_pions_cut.Write()
Aero_proj_yx_pions_cut.Write()
Calo_proj_yx_pions_cut.Write()
Calo_proj_hits_yx_pions_cut.Write()
P_DPexit_xy_pions_cut.Write()
d_Cut_Pion_Events_Prompt.cd()
ePiCoinTime_pions_cut_prompt.Write()
P_MMpi_pions_cut_prompt.Write()

d_Cut_Pion_Events_Random.cd()
ePiCoinTime_pions_cut_randm.Write()
P_MMpi_pions_cut_randm.Write() 

outHistFile.Close()
infile.Close() 
print ("Processing Complete")
print("!!!!!!!!\n %i pi-n events \n!!!!!!!!" % BinIntegral_pions)
