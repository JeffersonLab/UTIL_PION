#! /usr/bin/python
# Description:
# ================================================================
# Time-stamp: "2021-08-31 00:21:09 trottar"
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
import matplotlib.pyplot as plt
import sys, math, os, subprocess

ROOTPrefix = sys.argv[1]
runNum = sys.argv[2]
MaxEvent=sys.argv[3]

# Add this to all files for more dynamic pathing
USER = subprocess.getstatusoutput("whoami") # Grab user info for file finding
HOST = subprocess.getstatusoutput("hostname")

if ("farm" in HOST[1]):
    REPLAYPATH = "/group/c-kaonlt/USERS/%s/hallc_replay_lt" % USER[1]
elif ("lark" in HOST[1]):
    REPLAYPATH = "/home/%s/work/JLab/hallc_replay_lt" % USER[1]
elif ("cdaq" in HOST[1]):
    REPLAYPATH = "/home/cdaq/hallc-online/hallc_replay_lt"
elif ("trottar" in HOST[1]):
    REPLAYPATH = "/home/trottar/Analysis/hallc_replay_lt"

sys.path.insert(0, '%s/UTIL_PION/bin/python/' % REPLAYPATH)
import kaonlt as klt
import scaler

print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER[1], HOST[1], REPLAYPATH))

thres_curr = 2.5

out_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/lumi_data.csv" % REPLAYPATH

# Construct the name of the rootfile based upon the info we provided
OUTPATH = "%s/UTIL_PION/OUTPUT/Analysis/PionLT" % REPLAYPATH        # Output folder location
rootName = "%s/UTIL_PION/ROOTfiles/Analysis/Lumi/%s_%s_%s.root" % (REPLAYPATH,ROOTPrefix,runNum,MaxEvent)     # Input file location and variables taking
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

report = "%s/UTIL_PION/REPORT_OUTPUT/Analysis/Lumi/%s_%s_%s.report" % (REPLAYPATH,ROOTPrefix,runNum,MaxEvent)

f = open(report)
    
psList = ['SW_Ps1_factor','SW_Ps2_factor','SW_Ps3_factor','SW_Ps4_factor','SW_Ps5_factor','SW_Ps6_factor']
    
psActual = [-1,1,2,3,5,9,17,33,65,129,257,513,1025,2049,4097,8193,16385,32769]
psValue = [-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]

for line in f:
    data = line.split(':')
    curr_data = line.split(':')
    if ('SW_BCM4A_Beam_Cut_Current' in curr_data[0]) :
        report_current_tmp = curr_data[1].split("uA")[0].strip()
    if ('SW_SHMS_Electron_Singles_TRACK_EFF' in curr_data[0]):
        SHMS_track_info = curr_data[1].split("+-")
    if ('SW_HMS_Electron_Singles_TRACK_EFF' in curr_data[0]):
        HMS_track_info = curr_data[1].split("+-")
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
report_current = float(report_current_tmp)
SHMS_track_eff = float(SHMS_track_info[0])
SHMS_track_uncern = float(SHMS_track_info[1])
HMS_track_eff = float(HMS_track_info[0])
HMS_track_uncern = float(HMS_track_info[1])

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

PS_list = [["PS1",PS1],["PS2",PS2],["PS3",PS3],["PS4",PS4],["PS5",PS5],["PS6",PS6]]
PS_used = []

for val in PS_list:
    if val[1] != 0:
        PS_used.append(val)

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
s_branch = klt.pyBranch(s_tree)

P_BCM4A_scalerCharge = s_tree.array("P.BCM4A.scalerCharge")

'''
ANALYSIS TREE, T
'''

tree = up.open(rootName)["T"]
branch = klt.pyBranch(tree)

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
# T_coin_pEDTM_tdcTime = tree.array("T.coin.pEDTM_tdcTime")
T_coin_pEDTM_tdcTime = tree.array("T.coin.pEDTM_tdcTimeRaw")
EvtType = tree.array("fEvtHdr.fEvtType")

fout = REPLAYPATH+'/UTIL_PION/DB/CUTS/run_type/lumi.cuts'

# read in cuts file and make dictionary
c = klt.pyPlot(REPLAYPATH)
# apply RF cuts to timing cuts file
#c.cut_RF(runNum,MaxEvent)
readDict = c.read_dict(fout,runNum)

# This method calls several methods in kaonlt package. It is required to create properly formated
# dictionaries. The evaluation must be in the analysis script because the analysis variables (i.e. the
# leaves of interest) are not defined in the kaonlt package. This makes the system more flexible
# overall, but a bit more cumbersome in the analysis script. Perhaps one day a better solution will be
# implimented.
def make_cutDict(cut,inputDict=None):

    global c

    c = klt.pyPlot(REPLAYPATH,readDict)
    x = c.w_dict(cut)
    # print("%s" % cut)
    # print("x ", x)
    
    if inputDict == None:
        inputDict = {}
        
    for key,val in readDict.items():
        if key == cut:
            inputDict.update({key : {}})

    for i,val in enumerate(x):
        tmp = x[i]
        if tmp == "":
            continue
        else:
            inputDict[cut].update(eval(tmp))
        
    return inputDict

cutDict = make_cutDict("p_track_lumi_before")
cutDict = make_cutDict("p_hadtrack_lumi_before",cutDict)
cutDict = make_cutDict("p_pitrack_lumi_before",cutDict)
cutDict = make_cutDict("p_ktrack_lumi_before",cutDict)
cutDict = make_cutDict("p_ptrack_lumi_before",cutDict)
cutDict = make_cutDict("p_track_lumi_after",cutDict)
cutDict = make_cutDict("p_hadtrack_lumi_after",cutDict)
cutDict = make_cutDict("p_pitrack_lumi_after",cutDict)
cutDict = make_cutDict("p_ktrack_lumi_after",cutDict)
cutDict = make_cutDict("p_ptrack_lumi_after",cutDict)
cutDict = make_cutDict("p_etrack_lumi_before",cutDict)
cutDict = make_cutDict("p_etrack_lumi_after",cutDict)
cutDict = make_cutDict("p_pcut_lumi_eff",cutDict)
cutDict = make_cutDict("h_track_lumi_before",cutDict)
cutDict = make_cutDict("h_etrack_lumi_before",cutDict)
cutDict = make_cutDict("h_track_lumi_after",cutDict)
cutDict = make_cutDict("h_etrack_lumi_after",cutDict)
cutDict = make_cutDict("h_etrack_lumi_after",cutDict)
cutDict = make_cutDict("h_hadcut_lumi_eff",cutDict)
cutDict = make_cutDict("h_cal",cutDict)
cutDict = make_cutDict("h_cer",cutDict)
cutDict = make_cutDict("p_cal",cutDict)
cutDict = make_cutDict("p_hgcer",cutDict)
cutDict = make_cutDict("p_aero",cutDict)
cutDict = make_cutDict("c_noedtm",cutDict)
cutDict = make_cutDict("c_edtm",cutDict)
c = klt.pyPlot(REPLAYPATH,cutDict)

def pid_cuts():

    f = plt.figure(figsize=(11.69,8.27))
    ax = f.add_subplot(231)
    ax.hist(H_cal_etotnorm,bins=c.setbin(H_cal_etotnorm,200),label='no cut',histtype='step',
            alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(H_cal_etotnorm,"h_cal"),
             bins=c.setbin(c.add_cut(H_cal_etotnorm,"h_cal"),200),label='no cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('H_cal_etotnorm')
    plt.ylabel('Count')

    ax = f.add_subplot(232)
    ax.hist(H_cer_npeSum,bins=c.setbin(H_cer_npeSum,200),label='no cut',histtype='step', alpha=0.5,
            stacked=True, fill=True)
    ax.hist(c.add_cut(H_cer_npeSum,"h_cer"),
            bins=c.setbin(c.add_cut(H_cer_npeSum,"h_cer"),200),label='no cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('H_cer_npeSum')
    plt.ylabel('Count')

    ax = f.add_subplot(233)
    ax.hist(P_cal_etotnorm,bins=c.setbin(P_cal_etotnorm,200),label='no cut',histtype='step',
            alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(P_cal_etotnorm,"p_cal"),
             bins=c.setbin(c.add_cut(P_cal_etotnorm,"p_cal"),200),label='no cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('P_cal_etotnorm')
    plt.ylabel('Count')

    ax = f.add_subplot(234)
    ax.hist(P_hgcer_npeSum,bins=c.setbin(P_hgcer_npeSum,200),label='no cut',histtype='step',
            alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(P_hgcer_npeSum,"p_hgcer"),
             bins=c.setbin(c.add_cut(P_hgcer_npeSum,"p_hgcer"),200),label='no cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('P_hgcer_npeSum')
    plt.ylabel('Count')

    ax = f.add_subplot(235)
    ax.hist(P_aero_npeSum,bins=c.setbin(P_aero_npeSum,200),label='no cut',histtype='step',
            alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(P_aero_npeSum,"p_aero"),
             bins=c.setbin(c.add_cut(P_aero_npeSum,"p_aero"),200),label='no cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('P_aero_npeSum')
    plt.ylabel('Count')    
    
def analysis():
    
    bcm_before = H_bcm_bcm4a_AvgCurrent
    bcm_after = [x for x in H_bcm_bcm4a_AvgCurrent if x > thres_curr ]

    EDTM = c.add_cut(T_coin_pEDTM_tdcTime,"c_edtm")
    
    SHMSTRIG = [x
             for x,bcm in zip(T_coin_pTRIG_SHMS_ROC2_tdcTime,bcm_after)
             if bcm > thres_curr             
             if x != 0.0]
    HMSTRIG  = [x
             for x,bcm in zip(T_coin_pTRIG_HMS_ROC1_tdcTime,bcm_after)
             if bcm > thres_curr
             if x !=0.0]
    if len(PS_used) > 2:
        COINTRIG  = [x
                     for x,bcm in zip(T_coin_pTRIG_COIN_ROC1_tdcTime,bcm_after)
                     if bcm > thres_curr
                     if x !=0.0]
    
    EventType = [x
             for x,bcm in zip(EvtType,bcm_after)
             if bcm > thres_curr]

    SHMSTRIG_cut = [trig1
                    for (trig1,evt,bcm) in zip(T_coin_pTRIG_SHMS_ROC2_tdcTime,EvtType,bcm_after)
                    if bcm > thres_curr
                    if evt == 1]
    
    # p_track_lumi_before
    p_track_lumi_before = c.add_cut(P_dc_ntrack,"p_track_lumi_before")
    
    # p_hadtrack_lumi_before
    p_hadtrack_lumi_before = c.add_cut(P_dc_ntrack,"p_hadtrack_lumi_before")

    # p_pitrack_lumi_before
    p_pitrack_lumi_before = c.add_cut(P_dc_ntrack,"p_pitrack_lumi_before")

    # p_ktrack_lumi_before
    p_ktrack_lumi_before = c.add_cut(P_dc_ntrack,"p_ktrack_lumi_before")

    # p_ptrack_lumi_before
    p_ptrack_lumi_before = c.add_cut(P_dc_ntrack,"p_ptrack_lumi_before")

    # p_track_lumi_after
    p_track_lumi_after = c.add_cut(P_dc_ntrack,"p_track_lumi_after")

    # p_hadtrack_lumi_after
    p_hadtrack_lumi_after = c.add_cut(P_dc_ntrack,"p_hadtrack_lumi_after")

    # p_pitrack_lumi_after
    p_pitrack_lumi_after = c.add_cut(P_dc_ntrack,"p_pitrack_lumi_after")

    # p_ktrack_lumi_after
    p_ktrack_lumi_after = c.add_cut(P_dc_ntrack,"p_ktrack_lumi_after")

    # p_ptrack_lumi_after
    p_ptrack_lumi_after = c.add_cut(P_dc_ntrack,"p_ptrack_lumi_after")

    # p_etrack_lumi_before
    p_etrack_lumi_before = c.add_cut(P_hgcer_npeSum,"p_etrack_lumi_before")

    # p_show_before
    p_show_before = c.add_cut(P_cal_etotnorm,"p_etrack_lumi_before")

    # p_etrack_lumi_after
    p_etrack_lumi_after  = c.add_cut(P_hgcer_npeSum,"p_etrack_lumi_after")

    # p_pcut_lumi_eff
    p_pcut_lumi_eff  = c.add_cut(P_hgcer_npeSum,"p_pcut_lumi_eff")

    # p_show_after
    p_show_after  = c.add_cut(P_cal_etotnorm,"p_pcut_lumi_eff")

    HMSTRIG_cut = [ x
                    for (x, evt, bcm ) in zip(T_coin_pTRIG_HMS_ROC1_tdcTime, EvtType, bcm_after)
                    if bcm > thres_curr
                    if evt == 2]

    # h_track_lumi_before
    h_track_lumi_before = c.add_cut(H_dc_ntrack,"h_track_lumi_before")

    # h_etrack_lumi_before
    h_etrack_lumi_before = c.add_cut(H_dc_ntrack,"h_etrack_lumi_before")
    
    # h_track_lumi_after
    h_track_lumi_after = c.add_cut(H_dc_ntrack,"h_track_lumi_after")


    # h_etrack_lumi_after
    h_etrack_lumi_after = c.add_cut(H_dc_ntrack,"h_etrack_lumi_after")    

    # h_etrack_lumi_before
    h_etrack_lumi_before_iterate = [H_cer_npeSum, bcm_after]
    h_etrack_lumi_before = [cer
                      for (cer, bcm) in zip(*h_etrack_lumi_before_iterate)
                      if bcm > thres_curr]

    # h_dp_before
    h_dp_before_iterate = [H_gtr_dp, bcm_after]
    h_dp_before = [h_dp
                      for (h_dp, bcm) in zip(*h_dp_before_iterate)
                      if bcm > thres_curr]

    # h_th_before
    h_th_before_iterate = [H_tr_tg_th, bcm_after]
    h_th_before = [h_th
                      for (h_th, bcm) in zip(*h_th_before_iterate)
                      if bcm > thres_curr]

    # h_ph_before
    h_ph_before_iterate = [H_tr_tg_ph, bcm_after]
    h_ph_before = [h_ph
                      for (h_ph, bcm) in zip(*h_ph_before_iterate)
                      if bcm > thres_curr]

    # h_show_before
    h_show_before_iterate = [H_tr_tg_th, bcm_after]
    h_show_before = [h_caletot
                      for (h_caletot, bcm) in zip(*h_show_before_iterate)
                      if bcm > thres_curr]
    
    # h_etrack_lumi_after
    h_etrack_lumi_after = c.add_cut(H_cer_npeSum,"h_etrack_lumi_after")
    
    # h_dp_after
    h_dp_after = c.add_cut(H_gtr_dp,"h_etrack_lumi_after")
    
    # h_th_after
    h_th_after = c.add_cut(H_tr_tg_th,"h_etrack_lumi_after")
    
    # h_ph_after
    h_ph_after = c.add_cut(H_tr_tg_ph,"h_etrack_lumi_after")
    
    # h_show_after
    h_show_after = c.add_cut(H_cal_etotnorm,"h_etrack_lumi_after")
    
    # h_hadcut_lumi_eff
    h_hadcut_lumi_eff = c.add_cut(H_cal_etotnorm,"h_hadcut_lumi_eff")

    # goodscinhit cut
    h_hadcuts_goodscinhit = c.add_cut(H_hod_goodscinhit,"h_hadcut_lumi_eff")
    p_pcuts_goodscinhit = c.add_cut(P_hod_goodscinhit,"p_pcut_lumi_eff")
    
    track_info = {
        
        "HMS_evts_scaler" : len(h_hadcut_lumi_eff),
        "HMS_evts_scaler_uncern" : math.sqrt(len(h_hadcut_lumi_eff)),
        "SHMS_evts_scaler" : len(p_pcut_lumi_eff),
        "SHMS_evts_scaler_uncern" : math.sqrt(len(p_pcut_lumi_eff)),
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

    print("Terminate","Selection rules have been applied, plotting results")
    print("Using prescale factors: %s %.0f, %s %.0f\n" % (PS_names[0],SHMS_PS,PS_names[1],HMS_PS))
    print("Total number of events: %.0f" % (len(EventType)))
    print("Number of EDTM  Events: %.0f" % (len(EDTM)))
    print("Number of SHMSTRIG Events: %.0f" % (SHMS_PS*len(SHMSTRIG_cut)))
    print("Number of HMSTRIG Events: %.0f\n" % (HMS_PS*len(HMSTRIG_cut)))
    print("Number of HMS good events: %.0f +/- %.0f " % ((HMS_PS*len(h_hadcut_lumi_eff))
                                                         ,math.sqrt(HMS_PS*len(h_hadcut_lumi_eff))))
    print("Calculated HMS tracking efficiency: %f +/- %f\n" % ((HMS_track_eff), (HMS_track_uncern)))
    print("Number of SHMS good events: %.0f +/- %.0f" % ((SHMS_PS*len(p_pcut_lumi_eff)),
                                                         math.sqrt(SHMS_PS*len(p_pcut_lumi_eff))))
    print("Calculated SHMS tracking efficiency: %f +/- %f\n" % ((SHMS_track_eff),(SHMS_track_uncern)))
    print("============================================================================\n\n")
          
    return track_info

def main():

    #id_cuts()
    #plt.show()

    # combine dictionaries
    scalers = scaler.scaler(PS_names, SHMS_PS, SHMS_PS, thres_curr,report_current,REPLAYPATH,runNum,MaxEvent,s_tree,s_branch) 
    track_info = analysis()

    # lumi_data = {**scalers , **track_info} # only python 3.5+
    data = {}
    for d in (scalers, track_info): 
        data.update(d)
    lumi_data = {i : data[i] for i in sorted(data.keys())}

    table  = pd.DataFrame([lumi_data], columns=lumi_data.keys())
    table = table.reindex(sorted(table.columns), axis=1)
    
    file_exists = os.path.isfile(out_f)

    if file_exists:
        table.to_csv(out_f, index = False, header=False, mode='a',)
    else:
        table.to_csv(out_f, index = False, header=True, mode='a',)

if __name__ == '__main__':
    main()
