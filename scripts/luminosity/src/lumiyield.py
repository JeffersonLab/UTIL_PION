#! /usr/bin/python
#
# Description: This is where the variables for the yield calculations are formulated.
# Variables calculated: tot_events, h_int_goodscin_evts, p_int_goodscin_evts, SHMSTRIG_cut, HMSTRIG_cut, HMS_track, HMS_track_uncern, SHMS_track, SHMS_track_uncern, accp_edtm
# ================================================================
# Time-stamp: "2021-11-18 05:15:31 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import uproot as up
import numpy as np
import pandas as pd
import scipy
import scipy.integrate as integrate
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, math, os, subprocess

################################################################################################################################################
'''
User Inputs
'''

ROOTPrefix = sys.argv[1]
runNum = sys.argv[2]
MaxEvent=sys.argv[3]

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
import ltsep as lt 

# Add this to all files for more dynamic pathing
USER =  lt.SetPath(os.path.realpath(__file__)).getPath("USER") # Grab user info for file finding
HOST = lt.SetPath(os.path.realpath(__file__)).getPath("HOST")
REPLAYPATH = lt.SetPath(os.path.realpath(__file__)).getPath("REPLAYPATH")
UTILPATH = lt.SetPath(os.path.realpath(__file__)).getPath("UTILPATH")
ANATYPE=lt.SetPath(os.path.realpath(__file__)).getPath("ANATYPE")

################################################################################################################################################
'''
Import scaler script for use in luminosity analysis
'''

# Import scaler table
import scaler

################################################################################################################################################

print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))

# Output for luminosity table
out_f = UTILPATH+"/scripts/luminosity/OUTPUTS/lumi_data.csv"

################################################################################################################################################
'''
Check that root/output paths and files exist for use
'''

# Construct the name of the rootfile based upon the info we provided
OUTPATH = UTILPATH+"/OUTPUT/Analysis/%sLT" % ANATYPE        # Output folder location
rootName = UTILPATH+"/ROOTfiles/Analysis/Lumi/%s_%s_%s.root" % (ROOTPrefix,runNum,MaxEvent)     # Input file location and variables taking
print ("Attempting to process %s" %(rootName))
lt.SetPath(os.path.realpath(__file__)).checkDir(OUTPATH)
lt.SetPath(os.path.realpath(__file__)).checkFile(rootName)
print("Output path checks out, outputting to %s" % (OUTPATH))

################################################################################################################################################
'''
Grab prescale values and tracking efficiencies from report file
'''

# Open report file to grab prescale values and tracking efficiency
report = UTILPATH+"/REPORT_OUTPUT/Analysis/Lumi/%s_%s_%s.report" % (ROOTPrefix,runNum,MaxEvent)
f = open(report)
psList = ['SW_Ps1_factor','SW_Ps2_factor','SW_Ps3_factor','SW_Ps4_factor','SW_Ps5_factor','SW_Ps6_factor']
    
# Prescale input value (psValue) to its actual DAQ understanding (psActual)
psActual = [-1,1,2,3,5,9,17,33,65,129,257,513,1025,2049,4097,8193,16385,32769]
psValue = [-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]

# Search root file for prescale values and tracking efficiency, then save as variables
for line in f:
    data = line.split(':')
    track_data = line.split(':')
    if ('SW_SHMS_Electron_Singles_TRACK_EFF' in track_data[0]):
        SHMS_track_info = track_data[1].split("+-")
    if ('SW_HMS_Electron_Singles_TRACK_EFF' in track_data[0]):
        HMS_track_info = track_data[1].split("+-")
    for i, obj in enumerate(psList) :
        if (psList[i] in data[0]) : 
            if (i == 0) :  
                ps1_tmp = data[1].strip()
            if (i == 1) : 
                ps2_tmp = data[1].strip()
            if (i == 2) :
                ps3_tmp = data[1].strip()
            if (i == 3) :
                ps4_tmp = data[1].strip()
            if (i == 4) :
                ps5_tmp = data[1].strip()
            if (i == 5) :
                ps6_tmp = data[1].strip()
ps1=int(ps1_tmp)
ps2=int(ps2_tmp)
ps3=int(ps3_tmp)
ps4=int(ps4_tmp)
ps5=int(ps5_tmp)
ps6=int(ps6_tmp)
#SHMS_track_eff = float(SHMS_track_info[0]) # Depreciated, define below 
SHMS_track_uncern = float(SHMS_track_info[1])
#HMS_track_eff = float(HMS_track_info[0]) # Depreciated, define below 
HMS_track_uncern = float(HMS_track_info[1])

# Convert the prescale input values to their actual DAQ values
for i,index in enumerate(psActual):
    #psValue
    if (index == ps1) :
        if(index == -1):
            PS1 = 0
        else:
            PS1 = psActual[i]
    if (index == ps2) :
        if(index == -1):
            PS2 = 0
        else:
            PS2 = psActual[i]            
    if (index == ps3) :
        if(index == -1):
            PS3 = 0
        else:
            PS3 = psActual[i]
    if (index == ps4) :
        if(index == -1):
            PS4 = 0
        else:
            PS4 = psActual[i]            
    if (index == ps5) :
        if(index == -1):
            PS5 = 0
        else:
            PS5 = psActual[i]
    if (index == ps6) :
        if(index == -1):
            PS6 = 0
        else:
            PS6 = psActual[i]            
f.close()

################################################################################################################################################
'''
Define prescale variables
'''

print("\nPre-scale values...\nPS1:{0}, PS2:{1}, PS3:{2}, PS4:{3}, PS5:{4}, PS6:{5}\n".format(PS1,PS2,PS3,PS4,PS5,PS6))

# Save only the used prescale triggers to the PS_used list
PS_list = [["PS1",PS1],["PS2",PS2],["PS3",PS3],["PS4",PS4],["PS5",PS5],["PS6",PS6]]
PS_used = []
for val in PS_list:
    if val[1] != 0:
        PS_used.append(val)

# Check if COIN trigger is used by seeing it was saved in the PS_used list
if len(PS_used) > 2:
    PS_names = [PS_used[0][0],PS_used[1][0],PS_used[2][0]]
    SHMS_PS = PS_used[0][1]
    HMS_PS = PS_used[1][1]
    COIN_PS = PS_used[2][1]
else:
    PS_names = [PS_used[0][0],PS_used[1][0]]
    SHMS_PS = PS_used[0][1]
    HMS_PS = PS_used[1][1]

################################################################################################################################################
'''
SCALER TREE, TSP
'''

s_tree = up.open(rootName)["TSP"]

P_BCM4A_scalerCharge = s_tree.array("P.BCM4A.scalerCharge")

################################################################################################################################################
'''
ANALYSIS TREE, T
'''

tree = up.open(rootName)["T"]

if PS_names[1] is "PS3" or PS_names[1] is "PS4":
    W = tree.array("H.kin.primary.W")
    H_cal_etotnorm = tree.array("H.cal.etotnorm")
    H_cer_npeSum = tree.array("H.cer.npeSum")
    H_gtr_dp = tree.array("H.gtr.dp")
    H_tr_tg_th = tree.array("H.gtr.th")
    H_tr_tg_ph = tree.array("H.gtr.ph")
    H_gtr_beta = tree.array("H.gtr.beta")
    H_tr_chi2 = tree.array("H.tr.chi2")
    H_tr_ndof = tree.array("H.tr.ndof")
    H_hod_goodscinhit = tree.array("H.hod.goodscinhit")
    H_hod_betanotrack = tree.array("H.hod.betanotrack")
    H_hod_goodstarttime = tree.array("H.hod.goodstarttime")
    H_dc_ntrack = tree.array("H.dc.ntrack")
    
    H_dc_1x1_nhit = tree.array("H.dc.1x1.nhit")
    H_dc_1u2_nhit = tree.array("H.dc.1u2.nhit")
    H_dc_1u1_nhit = tree.array("H.dc.1u1.nhit")
    H_dc_1v1_nhit = tree.array("H.dc.1v1.nhit")
    H_dc_1x2_nhit = tree.array("H.dc.1x2.nhit")
    H_dc_1v2_nhit = tree.array("H.dc.1v2.nhit")
    H_dc_2x1_nhit = tree.array("H.dc.2x1.nhit")
    H_dc_2u2_nhit = tree.array("H.dc.2u2.nhit")
    H_dc_2u1_nhit = tree.array("H.dc.2u1.nhit")
    H_dc_2v1_nhit = tree.array("H.dc.2v1.nhit")
    H_dc_2x2_nhit = tree.array("H.dc.2x2.nhit")
    H_dc_2v2_nhit = tree.array("H.dc.2v2.nhit")
    
    H_bcm_bcm4a_AvgCurrent = tree.array("H.bcm.bcm4a.AvgCurrent")
    H_cal_etottracknorm = tree.array("H.cal.etottracknorm")

if PS_names[0] is "PS1" or PS_names[0] is "PS2":
    #W = tree.array("P.kin.primary.W")
    P_cal_etotnorm = tree.array("P.cal.etotnorm")
    P_hgcer_npeSum = tree.array("P.hgcer.npeSum")
    P_aero_npeSum = tree.array("P.aero.npeSum")
    P_gtr_dp = tree.array("P.gtr.dp")
    P_gtr_th = tree.array("P.gtr.th")
    P_gtr_ph = tree.array("P.gtr.ph")
    P_gtr_beta = tree.array("P.gtr.beta")
    P_tr_chi2 = tree.array("P.tr.chi2")
    P_tr_ndof = tree.array("P.tr.ndof")
    P_hod_goodscinhit = tree.array("P.hod.goodscinhit")
    P_hod_betanotrack = tree.array("P.hod.betanotrack")
    P_hod_goodstarttime = tree.array("P.hod.goodstarttime")
    P_dc_ntrack = tree.array("P.dc.ntrack")
    if ANATYPE == "Pion":
        P_ngcer_npeSum = tree.array("P.ngcer.npeSum")
    
    P_dc_1x1_nhit = tree.array("P.dc.1x1.nhit")
    P_dc_1u2_nhit = tree.array("P.dc.1u2.nhit")
    P_dc_1u1_nhit = tree.array("P.dc.1u1.nhit")
    P_dc_1v1_nhit = tree.array("P.dc.1v1.nhit")
    P_dc_1x2_nhit = tree.array("P.dc.1x2.nhit")
    P_dc_1v2_nhit = tree.array("P.dc.1v2.nhit")
    P_dc_2x1_nhit = tree.array("P.dc.2x1.nhit")
    P_dc_2u2_nhit = tree.array("P.dc.2u2.nhit")
    P_dc_2u1_nhit = tree.array("P.dc.2u1.nhit")
    P_dc_2v1_nhit = tree.array("P.dc.2v1.nhit")
    P_dc_2x2_nhit = tree.array("P.dc.2x2.nhit")
    P_dc_2v2_nhit = tree.array("P.dc.2v2.nhit")

    P_cal_etottracknorm = tree.array("P.cal.etottracknorm")

if PS_names[0] is "PS1":
    T_coin_pTRIG_SHMS_ROC1_tdcTimeRaw = tree.array("T.coin.pTRIG1_ROC1_tdcTimeRaw")
    T_coin_pTRIG_SHMS_ROC2_tdcTimeRaw = tree.array("T.coin.pTRIG1_ROC2_tdcTimeRaw")
    T_coin_pTRIG_SHMS_ROC1_tdcTime = tree.array("T.coin.pTRIG1_ROC1_tdcTime")
    T_coin_pTRIG_SHMS_ROC2_tdcTime = tree.array("T.coin.pTRIG1_ROC2_tdcTime")

if PS_names[0] is "PS2":
    T_coin_pTRIG_SHMS_ROC1_tdcTimeRaw = tree.array("T.coin.pTRIG2_ROC1_tdcTimeRaw")
    T_coin_pTRIG_SHMS_ROC2_tdcTimeRaw = tree.array("T.coin.pTRIG2_ROC2_tdcTimeRaw")
    T_coin_pTRIG_SHMS_ROC1_tdcTime = tree.array("T.coin.pTRIG2_ROC1_tdcTime")
    T_coin_pTRIG_SHMS_ROC2_tdcTime = tree.array("T.coin.pTRIG2_ROC2_tdcTime")

if PS_names[1] is "PS3":
    T_coin_pTRIG_HMS_ROC1_tdcTimeRaw = tree.array("T.coin.pTRIG3_ROC1_tdcTimeRaw")
    T_coin_pTRIG_HMS_ROC2_tdcTimeRaw = tree.array("T.coin.pTRIG3_ROC2_tdcTimeRaw")
    T_coin_pTRIG_HMS_ROC1_tdcTime = tree.array("T.coin.pTRIG3_ROC1_tdcTime")
    T_coin_pTRIG_HMS_ROC2_tdcTime = tree.array("T.coin.pTRIG3_ROC2_tdcTime")

if PS_names[1] is "PS4":
    T_coin_pTRIG_HMS_ROC1_tdcTimeRaw = tree.array("T.coin.pTRIG4_ROC1_tdcTimeRaw")
    T_coin_pTRIG_HMS_ROC2_tdcTimeRaw = tree.array("T.coin.pTRIG4_ROC2_tdcTimeRaw")
    T_coin_pTRIG_HMS_ROC1_tdcTime = tree.array("T.coin.pTRIG4_ROC1_tdcTime")
    T_coin_pTRIG_HMS_ROC2_tdcTime = tree.array("T.coin.pTRIG4_ROC2_tdcTime")

# Check if COIN trigger is used
if len(PS_used) > 2:
    if PS_names[2] is "PS5":
        T_coin_pTRIG_COIN_ROC1_tdcTimeRaw = tree.array("T.coin.pTRIG5_ROC1_tdcTimeRaw")
        T_coin_pTRIG_COIN_ROC2_tdcTimeRaw = tree.array("T.coin.pTRIG5_ROC2_tdcTimeRaw")
        T_coin_pTRIG_COIN_ROC1_tdcTime = tree.array("T.coin.pTRIG5_ROC1_tdcTime")
        T_coin_pTRIG_COIN_ROC2_tdcTime = tree.array("T.coin.pTRIG5_ROC2_tdcTime")
        
    if PS_names[2] is "PS6":
        T_coin_pTRIG_COIN_ROC1_tdcTimeRaw = tree.array("T.coin.pTRIG6_ROC1_tdcTimeRaw")
        T_coin_pTRIG_COIN_ROC2_tdcTimeRaw = tree.array("T.coin.pTRIG6_ROC2_tdcTimeRaw")
        T_coin_pTRIG_COIN_ROC1_tdcTime = tree.array("T.coin.pTRIG6_ROC1_tdcTime")
        T_coin_pTRIG_COIN_ROC2_tdcTime = tree.array("T.coin.pTRIG6_ROC2_tdcTime")

T_coin_pFADC_TREF_ROC2_adcPed = tree.array("T.coin.pFADC_TREF_ROC2_adcPed")
T_coin_hFADC_TREF_ROC1_adcPed = tree.array("T.coin.hFADC_TREF_ROC1_adcPed")
T_coin_pFADC_TREF_ROC2_adcPulseTimeRaw = tree.array("T.coin.pFADC_TREF_ROC2_adcPulseTimeRaw")
T_coin_hFADC_TREF_ROC1_adcPulseTimeRaw = tree.array("T.coin.hFADC_TREF_ROC1_adcPulseTimeRaw")
T_coin_pEDTM_tdcTimeRaw = tree.array("T.coin.pEDTM_tdcTimeRaw")
T_coin_pEDTM_tdcTime = tree.array("T.coin.pEDTM_tdcTime")
EvtType = tree.array("fEvtHdr.fEvtType")

################################################################################################################################################
'''
Define and set up cuts
'''

fout = UTILPATH+'/DB/CUTS/run_type/lumi.cuts'

if ANATYPE == "Pion":
    cuts = ["h_cal_nt","h_cer_nt","p_cal_nt","p_hgcer_nt","p_aero_nt","p_ngcer_nt","p_ecut_lumi_nt","h_ecut_lumi_nt","c_noedtm","c_edtm","c_ptrigHMS","c_ptrigSHMS","c_curr","h_etrack_lumi_before","h_etrack_lumi_after","p_etrack_lumi_before","p_etrack_lumi_after"]
    # Check if COIN trigger is used
    if len(PS_used) > 2:
        cuts = ["h_cal_nt","h_cer_nt","p_cal_nt","p_hgcer_nt","p_aero_nt","p_ngcer_nt","p_ecut_lumi_nt","h_ecut_lumi_nt","c_noedtm","c_edtm","c_ptrigHMS","c_ptrigSHMS","c_ptrigCOIN","c_curr","h_etrack_lumi_before","h_etrack_lumi_after","p_etrack_lumi_before","p_etrack_lumi_after"]
else:
    cuts = ["h_cal_nt","h_cer_nt","p_cal_nt","p_hgcer_nt","p_aero_nt","p_ecut_lumi_nt","h_ecut_lumi_nt","c_noedtm","c_edtm","c_ptrigHMS","c_ptrigSHMS","c_curr","h_etrack_lumi_before","h_etrack_lumi_after","p_etrack_lumi_before","p_etrack_lumi_after"]
    # Check if COIN trigger is used
    if len(PS_used) > 2:
        cuts = ["h_cal_nt","h_cer_nt","p_cal_nt","p_hgcer_nt","p_aero_nt","p_ecut_lumi_nt","h_ecut_lumi_nt","c_noedtm","c_edtm","c_ptrigHMS","c_ptrigSHMS","c_ptrigCOIN","c_curr","h_etrack_lumi_before","h_etrack_lumi_after","p_etrack_lumi_before","p_etrack_lumi_after"]

cutVals = []
def make_cutDict(cuts,fout,runNum,CURRENT_ENV,DEBUG=False):
    '''
    This method calls several methods in kaonlt package. It is required to create properly formated
    dictionaries. The evaluation must be in the analysis script because the analysis variables (i.e. the
    leaves of interest) are not defined in the kaonlt package. This makes the system more flexible
    overall, but a bit more cumbersome in the analysis script. Perhaps one day a better solution will be
    implimented.
    '''
    print ("\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    print (" Ensure that current values are set for the run you are processing in UTIL_PION/DB/PARAM/Current_Parameters.csv \n")
    print(" You will get errors like - NameError: name 'current' is not defined - if it is not")
    print ("\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")

    # read in cuts file and make dictionary
    importDict = lt.SetCuts(CURRENT_ENV).importDict(cuts,fout,runNum,False)
    for i,cut in enumerate(cuts):
        x = lt.SetCuts(CURRENT_ENV,importDict).booleanDict(cut)
        #######################################################################################
        # Make list of cut strings
        cutVals.append(x)

        # Threshold current
        if cut == "c_curr":
            global thres_curr, report_current
            # e.g. Grabbing threshold current (ie 2.5) from something like this [' {"H_bcm_bcm4a_AvgCurrent" : (abs(H_bcm_bcm4a_AvgCurrent-55) < 2.5)}']
            thres_curr = float(x[0].split(":")[1].split("<")[1].split(")")[0].strip())
            # e.g. Grabbing set current for run (ie 55) from something like this [' {"H_bcm_bcm4a_AvgCurrent" : (abs(H_bcm_bcm4a_AvgCurrent-55) < 2.5)}']
            report_current = float(x[0].split(":")[1].split("<")[0].split(")")[0].split("-")[1].strip())
        #######################################################################################
        print("\n%s" % cut)
        print(x, "\n")
        if i == 0:
            inputDict = {}
        cutDict = lt.SetCuts(CURRENT_ENV,importDict).readDict(cut,inputDict)
        for j,val in enumerate(x):
            cutDict = lt.SetCuts(CURRENT_ENV,importDict).evalDict(cut,eval(x[j]),cutDict)
    return lt.SetCuts(CURRENT_ENV,cutDict)

c = make_cutDict(cuts,fout,runNum,os.path.realpath(__file__),True)

################################################################################################################################################

def pid_cuts():
    '''
    Plots of pid cuts that will be applied to the event selection
    '''

    ###########################
    ######## 1D plots  ########
    ###########################

    f = plt.figure(figsize=(11.69,8.27))
    f.suptitle("Run %s" % runNum)

    ax = f.add_subplot(231)
    ax.hist(H_cal_etotnorm,bins=c.setbin(H_cal_etotnorm,200,0,2.0),label='no cut',histtype='step',alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(H_cal_etotnorm,"h_cal_nt"),bins=c.setbin(H_cal_etotnorm,200,0,2.0),label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('H_cal_etotnorm')
    plt.ylabel('Count')

    ax = f.add_subplot(232)
    ax.hist(H_cer_npeSum,bins=c.setbin(H_cer_npeSum,200,0,60),label='no cut',histtype='step', alpha=0.5,stacked=True, fill=True)
    ax.hist(c.add_cut(H_cer_npeSum,"h_cer_nt"),bins=c.setbin(H_cer_npeSum,200,0,60),label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('H_cer_npeSum')
    plt.ylabel('Count')

    ax = f.add_subplot(233)
    ax.hist(P_cal_etotnorm,bins=c.setbin(P_cal_etotnorm,200,0,4),label='no cut',histtype='step',alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(P_cal_etotnorm,"p_cal_nt"),bins=c.setbin(P_cal_etotnorm,200,0,4),label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('P_cal_etotnorm')
    plt.ylabel('Count')

    ax = f.add_subplot(234)
    ax.hist(P_hgcer_npeSum,bins=c.setbin(P_hgcer_npeSum,200,0,200),label='no cut',histtype='step',alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(P_hgcer_npeSum,"p_hgcer_nt"),bins=c.setbin(P_hgcer_npeSum,200,0,200),label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('P_hgcer_npeSum')
    plt.ylabel('Count')

    ax = f.add_subplot(235)
    ax.hist(P_aero_npeSum,bins=c.setbin(P_aero_npeSum,200,0,400),label='no cut',histtype='step',alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(P_aero_npeSum,"p_aero_nt"),bins=c.setbin(P_aero_npeSum,200,0,400),label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('P_aero_npeSum')
    plt.ylabel('Count')   

    if ANATYPE == "Pion":
        ax = f.add_subplot(236)
        ax.hist(P_ngcer_npeSum,bins=c.setbin(P_ngcer_npeSum,200,0,250),label='no cut',histtype='step',alpha=0.5, stacked=True, fill=True)
        ax.hist(c.add_cut(P_ngcer_npeSum,"p_ngcer_nt"), bins=c.setbin(P_ngcer_npeSum,200,0,250),label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
        plt.yscale('log')
        plt.xlabel('P_ngcer_npeSum')
        plt.ylabel('Count')

    plt.legend(loc="upper right")

    plt.tight_layout(rect=[0,0.03,1,0.95])   
    plt.savefig(UTILPATH+'/scripts/luminosity/OUTPUTS/plots/pid/pid_%s.png' % (runNum))

    ###########################
    ######## 2D plots  ########
    ###########################

    f = plt.figure(figsize=(19.20,8.00))
    f.suptitle("Run %s" % runNum)

    ax = f.add_subplot(241)
    ax.hist2d(H_cal_etotnorm,H_cer_npeSum,bins=[c.setbin(H_cal_etotnorm,400,0,2.0),c.setbin(H_cer_npeSum,400,0,30)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(H_cal_etotnorm,"h_ecut_lumi_nt"), c.add_cut(H_cer_npeSum,"h_ecut_lumi_nt"), bins=[c.setbin(H_cal_etotnorm,400,0,2.0),c.setbin(H_cer_npeSum,400,0,30)],cmin=1,label='cut', alpha=0.5)
    plt.xlabel('H_cal_etotnorm')
    plt.ylabel('H_cer_npeSum')

    ax = f.add_subplot(242)
    ax.hist2d(P_cal_etotnorm,P_hgcer_npeSum,bins=[c.setbin(P_cal_etotnorm,400,0,4),c.setbin(P_hgcer_npeSum,400,0,80)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(P_cal_etotnorm,"p_ecut_lumi_nt"),c.add_cut(P_hgcer_npeSum,"p_ecut_lumi_nt"),bins=[c.setbin(P_cal_etotnorm,400,0,4),c.setbin(P_hgcer_npeSum,400,0,80)],cmin=1,label='cut', alpha=0.5)
    plt.xlabel('P_cal_etotnorm')
    plt.ylabel('P_hgcer_npeSum')

    ax = f.add_subplot(243)
    ax.hist2d(P_cal_etotnorm,P_aero_npeSum,bins=[c.setbin(P_cal_etotnorm,400,0,4),c.setbin(P_aero_npeSum,400,0,100)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(P_cal_etotnorm,"p_ecut_lumi_nt"),c.add_cut(P_aero_npeSum,"p_ecut_lumi_nt"),bins=[c.setbin(P_cal_etotnorm,400,0,4),c.setbin(P_aero_npeSum,400,0,100)],cmin=1,label='cut', alpha=0.5)
    plt.xlabel('P_cal_etotnorm')
    plt.ylabel('P_aero_npeSum')

    if ANATYPE == "Pion":
        ax = f.add_subplot(244)
        ax.hist2d(P_cal_etotnorm,P_ngcer_npeSum,bins=[c.setbin(P_cal_etotnorm,400,0,4),c.setbin(P_ngcer_npeSum,400,0,80)],cmin=1,label='no cut',alpha=0.5)
        ax.hist2d(c.add_cut(P_cal_etotnorm,"p_ecut_lumi_nt"),c.add_cut(P_ngcer_npeSum,"p_ecut_lumi_nt"),bins=[c.setbin(P_cal_etotnorm,400,0,4),c.setbin(P_ngcer_npeSum,400,0,80)],cmin=1,label='cut', alpha=0.5)
        plt.xlabel('P_cal_etotnorm')
        plt.ylabel('P_ngcer_npeSum')

    ax = f.add_subplot(245)
    ax.hist2d(P_aero_npeSum,P_hgcer_npeSum,bins=[c.setbin(P_aero_npeSum,400,0,100),c.setbin(P_hgcer_npeSum,400,0,80)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(P_aero_npeSum,"p_ecut_lumi_nt"),c.add_cut(P_hgcer_npeSum,"p_ecut_lumi_nt"),bins=[c.setbin(P_aero_npeSum,400,0,100),c.setbin(P_hgcer_npeSum,400,0,80)],cmin=1,label='cut', alpha=0.5)
    plt.xlabel('P_aero_npeSum')
    plt.ylabel('P_hgcer_npeSum')

    if ANATYPE == "Pion":
        ax = f.add_subplot(246)
        ax.hist2d(P_ngcer_npeSum,P_hgcer_npeSum,bins=[c.setbin(P_ngcer_npeSum,400,0,80),c.setbin(P_hgcer_npeSum,400,0,80)],cmin=1,label='no cut',alpha=0.5)
        ax.hist2d(c.add_cut(P_ngcer_npeSum,"p_ecut_lumi_nt"),c.add_cut(P_hgcer_npeSum,"p_ecut_lumi_nt"),bins=[c.setbin(P_ngcer_npeSum,400,0,80),c.setbin(P_hgcer_npeSum,400,0,80)],cmin=1,label='cut', alpha=0.5)
        plt.xlabel('P_ngcer_npeSum')
        plt.ylabel('P_hgcer_npeSum')

        ax = f.add_subplot(247)
        ax.hist2d(P_aero_npeSum,P_ngcer_npeSum,bins=[c.setbin(P_aero_npeSum,400,0,100),c.setbin(P_ngcer_npeSum,400,0,80)],cmin=1,label='no cut',alpha=0.5)
        ax.hist2d(c.add_cut(P_aero_npeSum,"p_ecut_lumi_nt"),c.add_cut(P_ngcer_npeSum,"p_ecut_lumi_nt"),bins=[c.setbin(P_aero_npeSum,400,0,100),c.setbin(P_ngcer_npeSum,400,0,80)],cmin=1,label='cut', alpha=0.5)
        plt.xlabel('P_aero_npeSum')
        plt.ylabel('P_ngcer_npeSum')

    ax = f.add_subplot(248)
    plt.axis('off')
    i=0
    plt.text(-0.15,1.00,"HMS cuts...",fontsize=8)
    for cut,val in zip(cuts,cutVals):
        if cut == "h_ecut_lumi_nt":
            for v in val:
                plt.text(-0.15,0.95-i/10," {}".format(v),fontsize=8)
                i+=1
    plt.text(-0.15,0.95-((i)/10+0.05),"SHMS cuts...",fontsize=8)
    for cut,val in zip(cuts,cutVals):
        if cut == "p_ecut_lumi_nt":
            for v in val:
                plt.text(-0.15,0.95-(i+1)/10," {}".format(v),fontsize=8)
                i+=1

    plt.tight_layout(rect=[0,0.03,1,0.95])   
    plt.savefig(UTILPATH+'/scripts/luminosity/OUTPUTS/plots/pid/pid2D_%s.png' % (runNum))

################################################################################################################################################
    
def analysis():
    '''
    Calculate variables for table
    '''

    # Applies edtm window cuts to edtm to get accepted edtm events
    EDTM = c.add_cut(T_coin_pEDTM_tdcTimeRaw,"c_edtm")
    
    # Applies trigger window cuts to trigger to get accepted trigger events
    SHMSTRIG = [x
                for x in c.add_cut(T_coin_pTRIG_SHMS_ROC2_tdcTime,"c_ptrigSHMS")
                if x != 0.0]
    # Applies trigger window cuts to trigger to get accepted trigger events
    HMSTRIG  = [x
                for x in c.add_cut(T_coin_pTRIG_HMS_ROC1_tdcTime,"c_ptrigHMS")
                if x !=0.0]
    # Check if COIN trigger is used
    if len(PS_used) > 2:
        # Applies trigger window cuts to trigger to get accepted trigger events
        COINTRIG  = [x
                     for x in c.add_cut(T_coin_pTRIG_COIN_ROC1_tdcTime,"c_ptrigCOIN")
                     if x !=0.0]
    
    # Cuts event type with current cut, probably pointless?
    EventType = c.add_cut(EvtType,"c_curr")

    # Applies trigger window cuts to trigger to get accepted trigger events
    SHMSTRIG_cut = [trig1
                    for (trig1,evt) in zip(c.add_cut(T_coin_pTRIG_SHMS_ROC2_tdcTime,"c_ptrigSHMS"),EvtType)
                    if evt == 1]
    # Applies trigger window cuts to trigger to get accepted trigger events
    HMSTRIG_cut = [ x
                    for (x, evt) in zip(c.add_cut(T_coin_pTRIG_HMS_ROC1_tdcTime,"c_ptrigHMS"), EvtType)
                    if evt == 2]
    # Check if COIN trigger is used
    if len(PS_used) > 2:
        # Applies trigger window cuts to trigger to get accepted trigger events
        COINTRIG_cut = [ x
                         for (x, evt) in zip(c.add_cut(T_coin_pTRIG_COIN_ROC1_tdcTime,"c_ptrigCOIN"), EvtType)
                         if (evt == 1 or evt == 2)]

    h_et_should = len(c.add_cut(H_cal_etotnorm,"h_etrack_lumi_before"))
    h_et_did = len(c.add_cut(H_cal_etotnorm,"h_etrack_lumi_after"))
    HMS_track_eff = h_et_did/h_et_should

    p_et_should = len(c.add_cut(P_cal_etotnorm,"p_etrack_lumi_before"))
    p_et_did = len(c.add_cut(P_cal_etotnorm,"p_etrack_lumi_after"))
    SHMS_track_eff = p_et_did/p_et_should

    # Applies PID cuts, once integrated this will give the events (no track)
    h_etotnorm = c.add_cut(H_cal_etotnorm,"h_ecut_lumi_nt") 
    p_etotnorm = c.add_cut(P_cal_etotnorm,"p_ecut_lumi_nt")

    # Applies PID cuts, once integrated this will give the events (track)
    h_hadcuts_goodscinhit = c.add_cut(H_hod_goodscinhit,"h_ecut_lumi_nt")
    p_pcuts_goodscinhit = c.add_cut(P_hod_goodscinhit,"p_ecut_lumi_nt")
    
    # Creates a dictionary for the calculated luminosity values 
    track_info = {
        
        "tot_events" : len(EventType),
        "h_int_etotnorm_evts" : (scipy.integrate.simps(h_etotnorm)),
        "p_int_etotnorm_evts" : (scipy.integrate.simps(p_etotnorm)),
        "h_int_goodscin_evts" : scipy.integrate.simps(h_hadcuts_goodscinhit),
        "p_int_goodscin_evts" : scipy.integrate.simps(p_pcuts_goodscinhit),
        "SHMSTRIG_cut" : len(SHMSTRIG_cut),
        "HMSTRIG_cut" : len(HMSTRIG_cut),
        "HMS_track" : HMS_track_eff,
        "HMS_track_uncern" : HMS_track_uncern,
        "SHMS_track" : SHMS_track_eff,
        "SHMS_track_uncern" : SHMS_track_uncern,
        "accp_edtm" : (len(EDTM)),
            
    }
    # Check if COIN trigger is used
    if len(PS_used) > 2:
        # Creates a dictionary for the calculated luminosity values 
        track_info = {
        
            "tot_events" : len(EventType),
            "h_int_etotnorm_evts" : scipy.integrate.simps(h_etotnorm),
            "p_int_etotnorm_evts" : scipy.integrate.simps(p_etotnorm),
            "h_int_goodscin_evts" : scipy.integrate.simps(h_hadcuts_goodscinhit),
            "p_int_goodscin_evts" : scipy.integrate.simps(p_pcuts_goodscinhit),
            "SHMSTRIG_cut" : len(SHMSTRIG_cut),
            "HMSTRIG_cut" : len(HMSTRIG_cut),
            "COINTRIG_cut" : len(COINTRIG_cut),
            "HMS_track" : HMS_track_eff,
            "HMS_track_uncern" : HMS_track_uncern,
            "SHMS_track" : SHMS_track_eff,
            "SHMS_track_uncern" : SHMS_track_uncern,
            "accp_edtm" : (len(EDTM)),
            
        }

    print("\nTerminate","Selection rules have been applied, plotting results")
    print("Total number of events: %.0f" % (track_info["tot_events"]))
    print("Number of EDTM  Events: %.0f" % (track_info["accp_edtm"]))
    print("Number of HMSTRIG Events: %.0f" % (HMS_PS*track_info["HMSTRIG_cut"]))
    print("Number of SHMSTRIG Events: %.0f" % (SHMS_PS*track_info["SHMSTRIG_cut"]))

    print("\nNumber of HMS good events: %.0f +/- %.0f " % ((HMS_PS*track_info["h_int_goodscin_evts"]), math.sqrt(HMS_PS*track_info["h_int_goodscin_evts"])))
    print("Calculated HMS tracking efficiency: %f +/- %f\n" % ((track_info["HMS_track"]), (track_info["HMS_track_uncern"])))

    print("Number of SHMS good events: %.0f +/- %.0f " % ((SHMS_PS*track_info["h_int_goodscin_evts"]), math.sqrt(SHMS_PS*track_info["h_int_goodscin_evts"])))
    print("Calculated SHMS tracking efficiency: %f +/- %f\n" % ((track_info["SHMS_track"]), (track_info["SHMS_track_uncern"])))

    print("============================================================================\n\n")
          
    return track_info

################################################################################################################################################

def main():

    pid_cuts()
    #plt.show()

    # lumi_data = {**scalers , **track_info} # only python 3.5+

    # Import dictionaries
    scalers = scaler.scaler(PS_names, HMS_PS, SHMS_PS, thres_curr, report_current, runNum, MaxEvent, s_tree) 
    track_info = analysis()

    # Merge and sort the two dictionaries of calculations
    data = {}
    for d in (scalers, track_info): 
        data.update(d)
    lumi_data = {i : data[i] for i in sorted(data.keys())}

    # Convert merged dictionary to a pandas dataframe then sort it
    table  = pd.DataFrame([lumi_data], columns=lumi_data.keys())
    table = table.reindex(sorted(table.columns), axis=1)
    
    file_exists = os.path.isfile(out_f)

    # Updates csv file with luminosity calculated values for later analysis (see plot_yield.py)
    if file_exists:
        try:
            out_data = pd.read_csv(out_f)
        except IOError:
            print("Error: %s does not appear to exist." % out_f)
        # Checks if run number is alread in csv and replaces it if it is there
        run_index = out_data.index[out_data["run number"] == int(runNum)].tolist()
        out_data.drop(run_index, inplace=True)
        out_data = out_data.append(table,ignore_index=True)
        print("Output luminosity values\n",out_data)
        out_data.to_csv(out_f, index = False, header=True, mode='w+',)
    else:
        table.to_csv(out_f, index = False, header=True, mode='a',)

if __name__ == '__main__':
    main()
