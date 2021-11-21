#! /usr/bin/python
#
# Description: This is where the variables for the yield calculations are formulated.
# Variables calculated: tot_events, h_int_goodscin_evts, p_int_goodscin_evts, SHMSTRIG_cut, HMSTRIG_cut, HMS_track, HMS_track_uncern, SHMS_track, SHMS_track_uncern, accp_edtm
# ================================================================
# Time-stamp: "2021-11-04 02:47:33 trottar" heinricn
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>, Nathan Heinrich <heinricn@uregina.ca>
#
# This is a version of the lumiyeild.py code, but adapted to work for fADC deadtime studies.
# Richard wrote the origonal code, this version was adapted by Nathan 
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

# Import scaler table
import scaler

################################################################################################################################################

print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))

# Output for luminosity table
out_f = UTILPATH+"/scripts/luminosity/OUTPUTS/fADC_data.csv"

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

# Open report file to grab prescale values and tracking efficiency
report = "%s/UTIL_PION/REPORT_OUTPUT/Analysis/Lumi/%s_%s_%s.report" % (REPLAYPATH,ROOTPrefix,runNum,MaxEvent)
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
SHMS_track_eff = float(SHMS_track_info[0])
SHMS_track_uncern = float(SHMS_track_info[1])
HMS_track_eff = float(HMS_track_info[0])
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


'''
SCALER TREE, TSP
'''

s_tree = up.open(rootName)["TSP"]

P_BCM4A_scalerCharge = s_tree.array("P.BCM4A.scalerCharge")

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
CTime_ePiCoinTime_ROC1 = tree.array("CTime.ePiCoinTime_ROC1")
H_cal_etottracknorm = tree.array("H.cal.etottracknorm")
P_cal_etottracknorm = tree.array("P.cal.etottracknorm")
EvtType = tree.array("fEvtHdr.fEvtType")

#list of all cuts that are used
cuts = ["h_cal", "h_cer", "h_cal_nt", "h_cer_nt", "p_cal", "p_hgcer", "p_aero", "p_cal_nt", "p_hgcer_nt", "p_aero_nt", "p_ngcer_nt", "p_picut_lumi_eff", "p_picut_lumi_nt",  "h_ecut_lumi_eff", "h_ecut_lumi_nt",  "c_noedtm", "c_edtm", "c_ptrigHMS", "c_ptrigSHMS", "c_ptrigCOIN", "coin_pid_only", "c_curr", "coin_pid_notrack", "coin_pid_notrack_rand", "coin_pid_track", "coin_pid_track_rand"]     
fout = REPLAYPATH+'/UTIL_PION/DB/CUTS/run_type/fADCdeadtime.cuts'

#"p_pitrack_lumi_before", "p_pitrack_lumi_after",
#"h_etrack_lumi_before", "h_etrack_lumi_after",

def make_cutDict(cuts,fout,runNum,CURRENT_ENV):
    '''
    This method calls several methods in kaonlt package. It is required to create properly formated
    dictionaries. The evaluation must be in the analysis script because the analysis variables (i.e. the
    leaves of interest) are not defined in the kaonlt package. This makes the system more flexible
    overall, but a bit more cumbersome in the analysis script. Perhaps one day a better solution will be
    implimented.
    '''

    # read in cuts file and make dictionary
    importDict = lt.SetCuts(CURRENT_ENV).importDict(cuts,fout,runNum)
    for i,cut in enumerate(cuts):
        x = lt.SetCuts(CURRENT_ENV,importDict).booleanDict(cut)
        #######################################################################################
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

c = make_cutDict(cuts,fout,runNum,os.path.realpath(__file__))

def pid_cuts():
    '''
    Plots of pid cuts that will be applied to the event selection
    '''

    #### 1D plots ####

    f = plt.figure(figsize=(11.69,8.27))
    ax = f.add_subplot(331)
    ax.hist(H_cal_etotnorm,bins=c.setbin(H_cal_etotnorm,200),label='no cut',histtype='step',
            alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(H_cal_etotnorm,"h_cal"), bins=c.setbin(c.add_cut(H_cal_etotnorm,"h_cal"),200),label='no cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('H_cal_etotnorm')
    plt.ylabel('Count')

    ax = f.add_subplot(332)
    ax.hist(H_cer_npeSum,bins=c.setbin(H_cer_npeSum,200),label='no cut',histtype='step', alpha=0.5,
            stacked=True, fill=True)
    ax.hist(c.add_cut(H_cer_npeSum,"h_cer"), bins=c.setbin(c.add_cut(H_cer_npeSum,"h_cer"),200),label='no cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('H_cer_npeSum')
    plt.ylabel('Count')

    ax = f.add_subplot(333)
    ax.hist(P_cal_etotnorm,bins=c.setbin(P_cal_etotnorm,200),label='no cut',histtype='step',
            alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(P_cal_etotnorm,"p_cal"), bins=c.setbin(c.add_cut(P_cal_etotnorm,"p_cal"),200),label='no cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('P_cal_etotnorm')
    plt.ylabel('Count')

    ax = f.add_subplot(334)
    ax.hist(P_hgcer_npeSum,bins=c.setbin(P_hgcer_npeSum,200),label='no cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(P_hgcer_npeSum,"p_hgcer"), bins=c.setbin(c.add_cut(P_hgcer_npeSum,"p_hgcer"),200),label='no cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('P_hgcer_npeSum')
    plt.ylabel('Count')

    ax = f.add_subplot(335)
    ax.hist(P_aero_npeSum,bins=c.setbin(P_aero_npeSum,200),label='no cut',histtype='step',
            alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(P_aero_npeSum,"p_aero"), bins=c.setbin(c.add_cut(P_aero_npeSum,"p_aero"),200),label='no cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('P_aero_npeSum')
    plt.ylabel('Count')    
    
    ax = f.add_subplot(336)
    ax.hist(P_ngcer_npeSum,bins=c.setbin(P_ngcer_npeSum,200,0,250),label='no cut',histtype='step',alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(P_ngcer_npeSum,"p_ngcer_nt"), bins=c.setbin(P_ngcer_npeSum,200,0,250),label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('P_ngcer_npeSum')
    plt.ylabel('Count')

#    bins=c.setbin(c.add_cut(CTime_ePiCoinTime_ROC1,"coin_pid_notrack"),50, -2.5, 2.5)
#    print(bins)
#    ax = f.add_subplot(337)
#    ax.hist(CTime_ePiCoinTime_ROC1,bins=c.setbin(c.add_cut(CTime_ePiCoinTime_ROC1,"coin_pid_notrack"),50,-2.5,2.5),label='no cut',histtype='step',alpha=0.5, stacked=True, fill=True)
#    ax.hist(c.add_cut(CTime_ePiCoinTime_ROC1,"coin_pid_notrack"), bins=c.setbin(c.add_cut(CTime_ePiCoinTime_ROC1,"coin_pid_notrack"),200), label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
#    plt.yscale('log')
#    plt.xlabel('Coin Time Prompt Peak')
#    plt.ylabel('Count')

#    ax = f.add_subplot(338)
#    ax.hist(CTime_ePiCoinTime_ROC1, bins=c.setbin(c.add_cut(CTime_ePiCoinTime_ROC1,"coin_pid_notrack"),200),label='no cut',histtype='step',alpha=0.5, stacked=True, fill=True)
#    ax.hist(c.add_cut(CTime_ePiCoinTime_ROC1,"coin_pid_notrack_rand"), bins=c.setbin(c.add_cut(CTime_ePiCoinTime_ROC1,"coin_pid_notrack"),200), label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
#    plt.yscale('log')
#    plt.xlabel('Coin Time Randoms')
#    plt.ylabel('Count')


    ### 2D plots ###

    f = plt.figure(figsize=(19.20,8.00))
    f.suptitle("Run %s" % runNum)

    ax = f.add_subplot(241)
    ax.hist2d(H_cal_etotnorm,H_cer_npeSum,bins=[c.setbin(H_cal_etotnorm,400,0,2.0),c.setbin(H_cer_npeSum,400,0,30)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(H_cal_etotnorm,"h_ecut_lumi_nt"), c.add_cut(H_cer_npeSum,"h_ecut_lumi_nt"), bins=[c.setbin(H_cal_etotnorm,400,0,2.0),c.setbin(H_cer_npeSum,400,0,30)],cmin=1,label='cut', alpha=0.5)
    plt.xlabel('H_cal_etotnorm')
    plt.ylabel('H_cer_npeSum')

    ax = f.add_subplot(242)
    ax.hist2d(P_cal_etotnorm,P_hgcer_npeSum,bins=[c.setbin(P_cal_etotnorm,400,0,4),c.setbin(P_hgcer_npeSum,400,0,80)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(P_cal_etotnorm,"p_picut_lumi_nt"),c.add_cut(P_hgcer_npeSum,"p_picut_lumi_nt"),bins=[c.setbin(P_cal_etotnorm,400,0,4),c.setbin(P_hgcer_npeSum,400,0,80)],cmin=1,label='cut', alpha=0.5)
    plt.xlabel('P_cal_etotnorm')
    plt.ylabel('P_hgcer_npeSum')

    ax = f.add_subplot(243)
    ax.hist2d(P_cal_etotnorm,P_aero_npeSum,bins=[c.setbin(P_cal_etotnorm,400,0,4),c.setbin(P_aero_npeSum,400,0,100)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(P_cal_etotnorm,"p_picut_lumi_nt"),c.add_cut(P_aero_npeSum,"p_picut_lumi_nt"),bins=[c.setbin(P_cal_etotnorm,400,0,4),c.setbin(P_aero_npeSum,400,0,100)],cmin=1,label='cut', alpha=0.5)
    plt.xlabel('P_cal_etotnorm')
    plt.ylabel('P_aero_npeSum')

    ax = f.add_subplot(244)
    ax.hist2d(P_cal_etotnorm,P_ngcer_npeSum,bins=[c.setbin(P_cal_etotnorm,400,0,4),c.setbin(P_ngcer_npeSum,400,0,80)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(P_cal_etotnorm,"p_picut_lumi_nt"),c.add_cut(P_ngcer_npeSum,"p_picut_lumi_nt"),bins=[c.setbin(P_cal_etotnorm,400,0,4),c.setbin(P_ngcer_npeSum,400,0,80)],cmin=1,label='cut', alpha=0.5)
    plt.xlabel('P_cal_etotnorm')
    plt.ylabel('P_ngcer_npeSum')

    ax = f.add_subplot(245)
    ax.hist2d(P_aero_npeSum,P_hgcer_npeSum,bins=[c.setbin(P_aero_npeSum,400,0,100),c.setbin(P_hgcer_npeSum,400,0,80)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(P_aero_npeSum,"p_picut_lumi_nt"),c.add_cut(P_hgcer_npeSum,"p_picut_lumi_nt"),bins=[c.setbin(P_aero_npeSum,400,0,100),c.setbin(P_hgcer_npeSum,400,0,80)],cmin=1,label='cut', alpha=0.5)
    plt.xlabel('P_aero_npeSum')
    plt.ylabel('P_hgcer_npeSum')

    ax = f.add_subplot(246)
    ax.hist2d(P_ngcer_npeSum,P_hgcer_npeSum,bins=[c.setbin(P_ngcer_npeSum,400,0,80),c.setbin(P_hgcer_npeSum,400,0,80)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(P_ngcer_npeSum,"p_picut_lumi_nt"),c.add_cut(P_hgcer_npeSum,"p_picut_lumi_nt"),bins=[c.setbin(P_ngcer_npeSum,400,0,80),c.setbin(P_hgcer_npeSum,400,0,80)],cmin=1,label='cut', alpha=0.5)
    plt.xlabel('P_ngcer_npeSum')
    plt.ylabel('P_hgcer_npeSum')

    ax = f.add_subplot(247)
    ax.hist2d(P_aero_npeSum,P_ngcer_npeSum,bins=[c.setbin(P_aero_npeSum,400,0,100),c.setbin(P_ngcer_npeSum,400,0,80)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(P_aero_npeSum,"p_picut_lumi_nt"),c.add_cut(P_ngcer_npeSum,"p_picut_lumi_nt"),bins=[c.setbin(P_aero_npeSum,400,0,100),c.setbin(P_ngcer_npeSum,400,0,80)],cmin=1,label='cut', alpha=0.5)
    plt.xlabel('P_aero_npeSum')
    plt.ylabel('P_ngcer_npeSum')

    #ax = f.add_subplot(248)
    #plt.axis('off')
    #i=0
    #plt.text(-0.15,1.00,"HMS cuts...",fontsize=8)
    #for cut,val in zip(cuts,cutVals):
    #    if cut == "h_ecut_lumi_nt":
    #        for v in val:
    #            plt.text(-0.15,0.95-i/10," {}".format(v),fontsize=8)
    #            i+=1
    #plt.text(-0.15,0.95-((i)/10+0.05),"SHMS cuts...",fontsize=8)
    #for cut,val in zip(cuts,cutVals):
    #    if cut == "p_ecut_lumi_nt":
    #        for v in val:
    #            plt.text(-0.15,0.95-(i+1)/10," {}".format(v),fontsize=8)
    #            i+=1

    plt.tight_layout(rect=[0,0.03,1,0.95])
    plt.savefig(UTILPATH+'/scripts/luminosity/OUTPUTS/plots/pid/pid2D_%s.png' % (runNum))

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
    coin_etotnorm = c.add_cut(H_cal_etotnorm, "coin_pid_notrack")
    coin_etotnorm_rand = c.add_cut(H_cal_etotnorm, "coin_pid_notrack_rand")
    coin_etottracknorm = c.add_cut(H_cal_etottracknorm, "coin_pid_track")
    coin_etottracknorm_rand = c.add_cut(H_cal_etottracknorm, "coin_pid_track_rand")

    COIN_pid_noCT = c.add_cut(H_cal_etotnorm, "coin_pid_only")

    # Applies PID cuts, once integrated this will give the events (no track)
    h_etotnorm = c.add_cut(H_cal_etotnorm,"h_ecut_lumi_nt") 
    p_etotnorm = c.add_cut(P_cal_etotnorm,"p_picut_lumi_nt")
    #h_etotnorm = c.add_cut(H_cal_etotnorm,"h_ecut_lumi_eff")

    # Applies PID cuts, once integrated this will give the events (track)
    h_ecuts_etottracknorm = c.add_cut(H_cal_etottracknorm,"h_ecut_lumi_eff")
    p_pcuts_etottracknorm = c.add_cut(P_cal_etottracknorm,"p_picut_lumi_eff")
    
    #HMS_Track_Eff = scipy.integrate.simps(c.add_cut(H_cal_etotnorm,"h_etrack_lumi_after"))/scipy.integrate.simps(c.add_cut(H_cal_etotnorm,"h_etrack_lumi_before"))
    #SHMS_Track_Eff = scipy.integrate.simps(c.add_cut(P_cal_etotnorm,"p_pitrack_lumi_after"))/scipy.integrate.simps(c.add_cut(P_cal_etotnorm,"p_pitrack_lumi_before"))

    # Creates a dictionary for the calculated luminosity values 
    track_info = {
        
        "tot_events" : len(EventType),
        "h_int_etotnorm_evts" : scipy.integrate.simps(h_etotnorm),
        "p_int_etotnorm_evts" : scipy.integrate.simps(p_etotnorm),
        "h_int_etottracknorm_evts" : scipy.integrate.simps(h_ecuts_etottracknorm),
        "p_int_etottracknorm_evts" : scipy.integrate.simps(p_pcuts_etottracknorm),
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
            "coin_int_etotnorm_evts" : scipy.integrate.simps(coin_etotnorm),
            "coin_int_etottracknorm_evts" : scipy.integrate.simps(coin_etottracknorm),
            "coin_int_etotnorm_evts_rand" : scipy.integrate.simps(coin_etotnorm_rand),
            "coin_int_etottracknorm_evts_rand" : scipy.integrate.simps(coin_etottracknorm_rand),
            "coin_int_noct_notrack" : scipy.integrate.simps(COIN_pid_noCT),
            "coin_int_noct" : scipy.integrate.simps(COIN_pid_noCT),
            "h_int_etottracknorm_evts" : scipy.integrate.simps(h_ecuts_etottracknorm),
            "p_int_etottracknorm_evts" : scipy.integrate.simps(p_pcuts_etottracknorm),
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

    print("\nNumber of HMS good events: %.0f +/- %.0f " % ((HMS_PS*track_info["h_int_etottracknorm_evts"]), math.sqrt(HMS_PS*track_info["h_int_etottracknorm_evts"])))
    print("ReportFile HMS tracking efficiency: %f +/- %f\n" % ((track_info["HMS_track"]), (track_info["HMS_track_uncern"])))

    print("Number of SHMS good events: %.0f +/- %.0f " % ((SHMS_PS*track_info["p_int_etottracknorm_evts"]), math.sqrt(SHMS_PS*track_info["p_int_etottracknorm_evts"])))
    print("ReportFile SHMS tracking efficiency: %f +/- %f\n" % ((track_info["SHMS_track"]), (track_info["SHMS_track_uncern"])))
    
#    print("Calculated HMS Tracking efficiency:  %f +/- ?" % (HMS_Track_Eff))
#    print("Calculated SHMS Tracking efficiency:  %f +/- ?" % (SHMS_Track_Eff))

    print("Number of HMS good untrack events: %.0f +/- %.0f" % ((track_info["h_int_etotnorm_evts"]), math.sqrt(track_info["h_int_etotnorm_evts"])))
    print("Number of SHMS good untrack events: %.0f +/- %.0f" % ((track_info["p_int_etotnorm_evts"]), math.sqrt(track_info["p_int_etotnorm_evts"])))
    print("Number of COIN good no Coin time cuts events: %.0f +/- %.0f" % ((track_info["coin_int_noct"]), math.sqrt(track_info["coin_int_noct"])))
    print("Number of COIN good untrack events: %.0f +/- %.0f" % ((track_info["coin_int_etotnorm_evts"]), math.sqrt(track_info["coin_int_etotnorm_evts"])))
    print("Number of COIN good track events: %.0f +/- %.0f" % ((track_info["coin_int_etottracknorm_evts"]), math.sqrt(track_info["coin_int_etottracknorm_evts"])))
    print("Number of COIN rand untrack events: %.0f" % (track_info["coin_int_etotnorm_evts_rand"]))
    print("Number of COIN rand track events: %.0f" % (track_info["coin_int_etottracknorm_evts_rand"]))

    

    print("============================================================================\n\n")
          
    return track_info

def main():

    pid_cuts()
    plt.show()

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

#    # Updates csv file with luminosity calculated values for later analysis (see plot_yield.py)
#    if file_exists:
#        table.to_csv(out_f, index = False, header=False, mode='a',)
#    else:
#        table.to_csv(out_f, index = False, header=True, mode='a',)

if __name__ == '__main__':
    main()
