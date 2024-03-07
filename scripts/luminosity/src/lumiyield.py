#! /usr/bin/python
#
# Description: This is where the variables for the yield calculations are formulated.
# Variables calculated: tot_events, h_int_etottracknorm_evts, p_int_etottracknorm_evts, SHMSTRIG_cut, HMSTRIG_cut, HMS_track, HMS_track_uncern, SHMS_track, SHMS_track_uncern, accp_edtm
# ================================================================
# Time-stamp: "2023-06-01 11:24:41 trottar"
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
from ltsep import Root

lt=Root(os.path.realpath(__file__))

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
SCRIPTPATH=lt.SCRIPTPATH
ANATYPE=lt.ANATYPE

################################################################################################################################################
'''
Grab prescale values and tracking efficiencies from report file
'''

# Open report file to grab prescale values and tracking efficiency
report = UTILPATH+"/REPORT_OUTPUT/Analysis/Lumi/%s_%s_%s.report" % (ROOTPrefix,runNum,MaxEvent)
f = open(report)
psList = ['KLT_Ps1_factor','KLT_Ps2_factor','KLT_Ps3_factor','KLT_Ps4_factor','KLT_Ps5_factor','KLT_Ps6_factor']
    
# Prescale input value (psValue) to its actual DAQ understanding (psActual)
psActual = [-1,1,2,3,5,9,17,33,65,129,257,513,1025,2049,4097,8193,16385,32769]
psValue = [-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]

# Search root file for prescale values and tracking efficiency, then save as variables
for line in f:
    data = line.split(':')
    track_data = line.split(':')
    if (5149 <= int(runNum) <= 5303):
        if ('KLT_SHMS_Pion_SING_TRACK_EFF' in track_data[0]):
            SHMS_track_info = track_data[1].split("+-")
    else:
        #if ('KLT_SHMS_Elec_SING_TRACK_EFF' in track_data[0]):
        if ('KLT_SHMS_Elec_ALL_TRACK_EFF' in track_data[0]):
            SHMS_track_info = track_data[1].split("+-")
    if ('KLT_HMS_Elec_SING_TRACK_EFF' in track_data[0]):
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
    
try:
    ps1=int(ps1_tmp)
except NameError:
    ps1=-1
try:
    ps2=int(ps2_tmp)
except NameError:
    ps2=-1
try:
    ps3=int(ps3_tmp)
except NameError:
    ps3=-1
try:
    ps4=int(ps4_tmp)
except NameError:
    ps4=-1
try:
    ps5=int(ps5_tmp)
except NameError:
    ps5=-1
try:
    ps6=int(ps6_tmp)
except NameError:
    ps6=-1
SHMS_track_eff = float(SHMS_track_info[0]) # Also define below, I'll probably use the report for consistency's sake
SHMS_track_uncern = float(SHMS_track_info[1])
HMS_track_eff = float(HMS_track_info[0]) # Also define below, I'll probably use the report for consistency's sake
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

PS_list = [["PS1",PS1],["PS2",PS2],["PS3",PS3],["PS4",PS4],["PS5",PS5],["PS6",PS6]]
PS_names = []
PSDict = {}
for val in PS_list:
    PSDict.update({val[0] : val[1]})
    if val[0] == "PS1" or val[0] == "PS2":
        if val[1] != 0:
            SHMS_PS = val[1]
            PS_names.append(val[0])    
    if val[0] == "PS3" or val[0] == "PS4":        
        if val[1] != 0:
            HMS_PS = val[1] 
            PS_names.append(val[0])
    if val[0] == "PS5" or val[0] == "PS6":        
        if val[1] != 0:
            COIN_PS = val[1]
            PS_names.append(val[0])

try:
    SHMS_PS
except NameError:
    SHMS_PS = None

try:
    HMS_PS
except NameError:
    HMS_PS = None

try:
    COIN_PS
except NameError:
    COIN_PS = None
            
################################################################################################################################################
'''
Check PID of luminosity run
'''

from check_pid import check_pid

pid = check_pid(runNum,ANATYPE)

HMS_PID = pid[0]
SHMS_PID = pid[1]

################################################################################################################################################

cut_f = '/DB/CUTS/run_type/lumi.cuts'

if ANATYPE == "Pion":
    cuts = ["h_cal_nt","h_cer_nt","p_cal_nt","p_hgcer_nt","p_aero_nt","h_cal","h_cer","p_cal","p_hgcer","p_aero","p_ngcer_nt","p_%scut_lumi_nt" % SHMS_PID,"h_%scut_lumi_nt" % HMS_PID,"p_%scut_lumi" % SHMS_PID,"h_%scut_lumi" % HMS_PID,"c_noedtm","c_edtm","c_edtmHMS","c_edtmSHMS","c_curr","h_%strack_lumi_before" % HMS_PID,"h_%strack_lumi_after" % HMS_PID,"p_%strack_lumi_before" % SHMS_PID,"p_%strack_lumi_after" % SHMS_PID]
else:
    cuts = ["h_cal_nt","h_cer_nt","p_cal_nt","p_hgcer_nt","p_aero_nt","h_cal","h_cer","p_cal","p_hgcer","p_aero","p_%scut_lumi_nt" % SHMS_PID,"h_%scut_lumi_nt" % HMS_PID,"p_%scut_lumi" % SHMS_PID,"h_%scut_lumi" % HMS_PID,"c_noedtm","c_edtm","c_edtmHMS","c_edtmSHMS","c_curr","h_%strack_lumi_before" % HMS_PID,"h_%strack_lumi_after" % HMS_PID,"p_%strack_lumi_before" % SHMS_PID,"p_%strack_lumi_after" % SHMS_PID]
    
for ps in PS_names:
    if ps == "PS1" or ps == "PS2":
        cuts+=["c_ptrigSHMS%s" % ps.replace("PS","")]        
    if ps == "PS3" or ps == "PS4":
        cuts+=["c_ptrigHMS%s" % ps.replace("PS","")]
    if ps == "PS5" or ps == "PS6":
        cuts+=["c_ptrigCOIN%s" % ps.replace("PS","")]    

lt=Root(os.path.realpath(__file__),"Lumi",ROOTPrefix,runNum,MaxEvent,cut_f,cuts)

proc_root = lt.setup_ana()
c = proc_root[0] # Cut object
tree = proc_root[1] # Dictionary of branches
strDict = proc_root[2] # Dictionary of cuts as strings

################################################################################################################################################

print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))
# Output for luminosity table
out_f = UTILPATH+"/scripts/luminosity/OUTPUTS/lumi_data.csv"
rootName = UTILPATH+"/ROOTfiles/Analysis/Lumi/%s_%s_%s.root" % (ROOTPrefix,runNum,MaxEvent)     # Input file location and variables taking

################################################################################################################################################
'''
Import scaler script for use in luminosity analysis
'''

# Import scaler table
import scaler

################################################################################################################################################
'''
SCALER TREE, TSP
'''

if SHMS_PS == None:
    s_tree = up.open(rootName)["TSH"]
else:
    s_tree = up.open(rootName)["TSP"]

################################################################################################################################################

for ps in PS_names:
    if ps == "PS1" or ps == "PS2":
        T_coin_pTRIG_SHMS_ROC1_tdcTimeRaw = tree["T_coin_pTRIG%s_ROC1_tdcTimeRaw" % ps.replace("PS","")]
        T_coin_pTRIG_SHMS_ROC2_tdcTimeRaw = tree["T_coin_pTRIG%s_ROC2_tdcTimeRaw" % ps.replace("PS","")]
        T_coin_pTRIG_SHMS_ROC1_tdcTime = tree["T_coin_pTRIG%s_ROC1_tdcTime" % ps.replace("PS","")]
        T_coin_pTRIG_SHMS_ROC2_tdcTime = tree["T_coin_pTRIG%s_ROC2_tdcTime" % ps.replace("PS","")]
        
    if ps == "PS3" or ps == "PS4":
        T_coin_pTRIG_HMS_ROC1_tdcTimeRaw = tree["T_coin_pTRIG%s_ROC1_tdcTimeRaw" % ps.replace("PS","")]
        T_coin_pTRIG_HMS_ROC2_tdcTimeRaw = tree["T_coin_pTRIG%s_ROC2_tdcTimeRaw" % ps.replace("PS","")]
        T_coin_pTRIG_HMS_ROC1_tdcTime = tree["T_coin_pTRIG%s_ROC1_tdcTime" % ps.replace("PS","")]
        T_coin_pTRIG_HMS_ROC2_tdcTime = tree["T_coin_pTRIG%s_ROC2_tdcTime" % ps.replace("PS","")]
        
    if ps == "PS5" or ps == "PS6":
        T_coin_pTRIG_COIN_ROC1_tdcTimeRaw = tree["T_coin_pTRIG%s_ROC1_tdcTimeRaw" % ps.replace("PS","")]
        T_coin_pTRIG_COIN_ROC2_tdcTimeRaw = tree["T_coin_pTRIG%s_ROC2_tdcTimeRaw" % ps.replace("PS","")]
        T_coin_pTRIG_COIN_ROC1_tdcTime = tree["T_coin_pTRIG%s_ROC1_tdcTime" % ps.replace("PS","")]
        T_coin_pTRIG_COIN_ROC2_tdcTime = tree["T_coin_pTRIG%s_ROC2_tdcTime" % ps.replace("PS","")]

T_coin_pEDTM_tdcTimeRaw = tree["T_coin_pEDTM_tdcTimeRaw"]
T_coin_pEDTM_tdcTime = tree["T_coin_pEDTM_tdcTime"]

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
    ax.hist(tree['H_cal_etotnorm'],bins=c.setbin(tree['H_cal_etotnorm'],200,0,2.0),label='no cut',histtype='step',alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(tree['H_cal_etotnorm'],"h_cal_nt"),bins=c.setbin(tree['H_cal_etotnorm'],200,0,2.0),label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('H_cal_etotnorm')
    plt.ylabel('Count')

    ax = f.add_subplot(232)
    ax.hist(tree['H_cer_npeSum'],bins=c.setbin(tree['H_cer_npeSum'],200,0,60),label='no cut',histtype='step', alpha=0.5,stacked=True, fill=True)
    ax.hist(c.add_cut(tree['H_cer_npeSum'],"h_cer_nt"),bins=c.setbin(tree['H_cer_npeSum'],200,0,60),label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('H_cer_npeSum')
    plt.ylabel('Count')

    ax = f.add_subplot(233)
    ax.hist(tree['P_cal_etotnorm'],bins=c.setbin(tree['P_cal_etotnorm'],200,0,4),label='no cut',histtype='step',alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(tree['P_cal_etotnorm'],"p_cal_nt"),bins=c.setbin(tree['P_cal_etotnorm'],200,0,4),label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('P_cal_etotnorm')
    plt.ylabel('Count')

    ax = f.add_subplot(234)
    ax.hist(tree['P_hgcer_npeSum'],bins=c.setbin(tree['P_hgcer_npeSum'],200,0,200),label='no cut',histtype='step',alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(tree['P_hgcer_npeSum'],"p_hgcer_nt"),bins=c.setbin(tree['P_hgcer_npeSum'],200,0,200),label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('P_hgcer_npeSum')
    plt.ylabel('Count')

    ax = f.add_subplot(235)
    ax.hist(tree['P_aero_npeSum'],bins=c.setbin(tree['P_aero_npeSum'],200,0,400),label='no cut',histtype='step',alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(tree['P_aero_npeSum'],"p_aero_nt"),bins=c.setbin(tree['P_aero_npeSum'],200,0,400),label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('P_aero_npeSum')
    plt.ylabel('Count')   

    if ANATYPE == "Pion":
        ax = f.add_subplot(236)
        ax.hist(tree['P_ngcer_npeSum'],bins=c.setbin(tree['P_ngcer_npeSum'],200,0,250),label='no cut',histtype='step',alpha=0.5, stacked=True, fill=True)
        ax.hist(c.add_cut(tree['P_ngcer_npeSum'],"p_ngcer_nt"), bins=c.setbin(tree['P_ngcer_npeSum'],200,0,250),label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
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
    ax.hist2d(tree['H_cal_etotnorm'],tree['H_cer_npeSum'],bins=[c.setbin(tree['H_cal_etotnorm'],400,0,2.0),c.setbin(tree['H_cer_npeSum'],400,0,30)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(tree['H_cal_etotnorm'],"h_%scut_lumi_nt" % HMS_PID), c.add_cut(tree['H_cer_npeSum'],"h_%scut_lumi_nt" % HMS_PID), bins=[c.setbin(tree['H_cal_etotnorm'],400,0,2.0),c.setbin(tree['H_cer_npeSum'],400,0,30)],cmin=1,label='cut', alpha=1.0)
    plt.xlabel('H_cal_etotnorm')
    plt.ylabel('H_cer_npeSum')

    ax = f.add_subplot(242)
    ax.hist2d(tree['P_cal_etotnorm'],tree['P_hgcer_npeSum'],bins=[c.setbin(tree['P_cal_etotnorm'],400,0,4),c.setbin(tree['P_hgcer_npeSum'],400,0,80)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(tree['P_cal_etotnorm'],"p_%scut_lumi_nt" % SHMS_PID),c.add_cut(tree['P_hgcer_npeSum'],"p_%scut_lumi_nt" % SHMS_PID),bins=[c.setbin(tree['P_cal_etotnorm'],400,0,4),c.setbin(tree['P_hgcer_npeSum'],400,0,80)],cmin=1,label='cut', alpha=1.0)
    plt.xlabel('P_cal_etotnorm')
    plt.ylabel('P_hgcer_npeSum')

    ax = f.add_subplot(243)
    ax.hist2d(tree['P_cal_etotnorm'],tree['P_aero_npeSum'],bins=[c.setbin(tree['P_cal_etotnorm'],400,0,4),c.setbin(tree['P_aero_npeSum'],400,0,100)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(tree['P_cal_etotnorm'],"p_%scut_lumi_nt" % SHMS_PID),c.add_cut(tree['P_aero_npeSum'],"p_%scut_lumi_nt" % SHMS_PID),bins=[c.setbin(tree['P_cal_etotnorm'],400,0,4),c.setbin(tree['P_aero_npeSum'],400,0,100)],cmin=1,label='cut', alpha=1.0)
    plt.xlabel('P_cal_etotnorm')
    plt.ylabel('P_aero_npeSum')

    if ANATYPE == "Pion":
        ax = f.add_subplot(244)
        ax.hist2d(tree['P_cal_etotnorm'],tree['P_ngcer_npeSum'],bins=[c.setbin(tree['P_cal_etotnorm'],400,0,4),c.setbin(tree['P_ngcer_npeSum'],400,0,80)],cmin=1,label='no cut',alpha=0.5)
        ax.hist2d(c.add_cut(tree['P_cal_etotnorm'],"p_%scut_lumi_nt" % SHMS_PID),c.add_cut(tree['P_ngcer_npeSum'],"p_%scut_lumi_nt" % SHMS_PID),bins=[c.setbin(tree['P_cal_etotnorm'],400,0,4),c.setbin(tree['P_ngcer_npeSum'],400,0,80)],cmin=1,label='cut', alpha=1.0)
        plt.xlabel('P_cal_etotnorm')
        plt.ylabel('P_ngcer_npeSum')

    ax = f.add_subplot(245)
    ax.hist2d(tree['P_aero_npeSum'],tree['P_hgcer_npeSum'],bins=[c.setbin(tree['P_aero_npeSum'],400,0,100),c.setbin(tree['P_hgcer_npeSum'],400,0,80)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(tree['P_aero_npeSum'],"p_%scut_lumi_nt" % SHMS_PID),c.add_cut(tree['P_hgcer_npeSum'],"p_%scut_lumi_nt" % SHMS_PID),bins=[c.setbin(tree['P_aero_npeSum'],400,0,100),c.setbin(tree['P_hgcer_npeSum'],400,0,80)],cmin=1,label='cut', alpha=1.0)
    plt.xlabel('P_aero_npeSum')
    plt.ylabel('P_hgcer_npeSum')

    if ANATYPE == "Pion":
        ax = f.add_subplot(246)
        ax.hist2d(tree['P_ngcer_npeSum'],tree['P_hgcer_npeSum'],bins=[c.setbin(tree['P_ngcer_npeSum'],400,0,80),c.setbin(tree['P_hgcer_npeSum'],400,0,80)],cmin=1,label='no cut',alpha=0.5)
        ax.hist2d(c.add_cut(tree['P_ngcer_npeSum'],"p_%scut_lumi_nt" % SHMS_PID),c.add_cut(tree['P_hgcer_npeSum'],"p_%scut_lumi_nt" % SHMS_PID),bins=[c.setbin(tree['P_ngcer_npeSum'],400,0,80),c.setbin(tree['P_hgcer_npeSum'],400,0,80)],cmin=1,label='cut', alpha=1.0)
        plt.xlabel('P_ngcer_npeSum')
        plt.ylabel('P_hgcer_npeSum')

        ax = f.add_subplot(247)
        ax.hist2d(tree['P_aero_npeSum'],tree['P_ngcer_npeSum'],bins=[c.setbin(tree['P_aero_npeSum'],400,0,100),c.setbin(tree['P_ngcer_npeSum'],400,0,80)],cmin=1,label='no cut',alpha=0.5)
        ax.hist2d(c.add_cut(tree['P_aero_npeSum'],"p_%scut_lumi_nt" % SHMS_PID),c.add_cut(tree['P_ngcer_npeSum'],"p_%scut_lumi_nt" % SHMS_PID),bins=[c.setbin(tree['P_aero_npeSum'],400,0,100),c.setbin(tree['P_ngcer_npeSum'],400,0,80)],cmin=1,label='cut', alpha=1.0)
        plt.xlabel('P_aero_npeSum')
        plt.ylabel('P_ngcer_npeSum')

    ax = f.add_subplot(248)
    plt.axis('off')
    i=0
    plt.text(-0.15,1.00,"HMS cuts... (%s)" % HMS_PID,fontsize=8)
    for key,val in strDict.items():
        if key == "h_%scut_lumi_nt" % HMS_PID:
            for v in val:
                plt.text(-0.15,0.95-i/10," {}".format(v),fontsize=8)
                i+=1
    plt.text(-0.15,0.95-((i)/10+0.05),"SHMS cuts... (%s)" % SHMS_PID,fontsize=8)
    for key,val in strDict.items():
        if key == "p_%scut_lumi_nt" % SHMS_PID:
            for v in val:
                plt.text(-0.15,0.95-(i+1)/10," {}".format(v),fontsize=8)
                i+=1
        if key == "c_curr":
            global thres_curr, report_current
            # e.g. Grabbing threshold current (ie 2.5) from something like this [' {"H_bcm_bcm1_AvgCurrent" : (abs(H_bcm_bcm1_AvgCurrent-55) < 2.5)}']
            thres_curr = float(val[0].split(":")[1].split("<")[1].split(")")[0].strip())
            # e.g. Grabbing set current for run (ie 55) from something like this [' {"H_bcm_bcm1_AvgCurrent" : (abs(H_bcm_bcm1_AvgCurrent-55) < 2.5)}']
            report_current = float(val[0].split(":")[1].split("<")[0].split(")")[0].split("-")[1].strip())

    plt.tight_layout(rect=[0,0.03,1,0.95])   
    plt.savefig(UTILPATH+'/scripts/luminosity/OUTPUTS/plots/pid/pid2D_%s.png' % (runNum))

def track_pid_cuts():
    '''
    Plots of track pid cuts that will be applied to the event selection
    '''

    ###########################
    ######## 1D plots  ########
    ###########################

    f = plt.figure(figsize=(11.69,8.27))
    f.suptitle("Run %s" % runNum)

    ax = f.add_subplot(231)
    ax.hist(tree['H_cal_etottracknorm'],bins=c.setbin(tree['H_cal_etottracknorm'],200,0,2.0),label='no cut',histtype='step',alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(tree['H_cal_etottracknorm'],"h_cal"),bins=c.setbin(tree['H_cal_etottracknorm'],200,0,2.0),label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('H_cal_etottracknorm')
    plt.ylabel('Count')

    ax = f.add_subplot(232)
    ax.hist(tree['H_cer_npeSum'],bins=c.setbin(tree['H_cer_npeSum'],200,0,60),label='no cut',histtype='step', alpha=0.5,stacked=True, fill=True)
    ax.hist(c.add_cut(tree['H_cer_npeSum'],"h_cer"),bins=c.setbin(tree['H_cer_npeSum'],200,0,60),label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('H_cer_npeSum')
    plt.ylabel('Count')

    ax = f.add_subplot(233)
    ax.hist(tree['P_cal_etottracknorm'],bins=c.setbin(tree['P_cal_etottracknorm'],200,0,4),label='no cut',histtype='step',alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(tree['P_cal_etottracknorm'],"p_cal"),bins=c.setbin(tree['P_cal_etottracknorm'],200,0,4),label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('P_cal_etottracknorm')
    plt.ylabel('Count')

    ax = f.add_subplot(234)
    ax.hist(tree['P_hgcer_npeSum'],bins=c.setbin(tree['P_hgcer_npeSum'],200,0,200),label='no cut',histtype='step',alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(tree['P_hgcer_npeSum'],"p_hgcer"),bins=c.setbin(tree['P_hgcer_npeSum'],200,0,200),label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('P_hgcer_npeSum')
    plt.ylabel('Count')

    ax = f.add_subplot(235)
    ax.hist(tree['P_aero_npeSum'],bins=c.setbin(tree['P_aero_npeSum'],200,0,400),label='no cut',histtype='step',alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(tree['P_aero_npeSum'],"p_aero"),bins=c.setbin(tree['P_aero_npeSum'],200,0,400),label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('P_aero_npeSum')
    plt.ylabel('Count')   

    if ANATYPE == "Pion":
        ax = f.add_subplot(236)
        ax.hist(tree['P_ngcer_npeSum'],bins=c.setbin(tree['P_ngcer_npeSum'],200,0,250),label='no cut',histtype='step',alpha=0.5, stacked=True, fill=True)
        ax.hist(c.add_cut(tree['P_ngcer_npeSum'],"p_ngcer"), bins=c.setbin(tree['P_ngcer_npeSum'],200,0,250),label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
        plt.yscale('log')
        plt.xlabel('P_ngcer_npeSum')
        plt.ylabel('Count')

    plt.legend(loc="upper right")

    plt.tight_layout(rect=[0,0.03,1,0.95])   
    plt.savefig(UTILPATH+'/scripts/luminosity/OUTPUTS/plots/track_pid/track_pid_%s.png' % (runNum))

    ###########################
    ######## 2D plots  ########
    ###########################

    f = plt.figure(figsize=(19.20,8.00))
    f.suptitle("Run %s" % runNum)

    ax = f.add_subplot(241)
    ax.hist2d(tree['H_cal_etottracknorm'],tree['H_cer_npeSum'],bins=[c.setbin(tree['H_cal_etottracknorm'],400,0,2.0),c.setbin(tree['H_cer_npeSum'],400,0,30)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(tree['H_cal_etottracknorm'],"h_%scut_lumi" % HMS_PID), c.add_cut(tree['H_cer_npeSum'],"h_%scut_lumi" % HMS_PID), bins=[c.setbin(tree['H_cal_etottracknorm'],400,0,2.0),c.setbin(tree['H_cer_npeSum'],400,0,30)],cmin=1,label='cut', alpha=1.0)
    plt.xlabel('H_cal_etottracknorm')
    plt.ylabel('H_cer_npeSum')

    ax = f.add_subplot(242)
    ax.hist2d(tree['P_cal_etottracknorm'],tree['P_hgcer_npeSum'],bins=[c.setbin(tree['P_cal_etottracknorm'],400,0,4),c.setbin(tree['P_hgcer_npeSum'],400,0,80)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(tree['P_cal_etottracknorm'],"p_%scut_lumi" % SHMS_PID),c.add_cut(tree['P_hgcer_npeSum'],"p_%scut_lumi" % SHMS_PID),bins=[c.setbin(tree['P_cal_etottracknorm'],400,0,4),c.setbin(tree['P_hgcer_npeSum'],400,0,80)],cmin=1,label='cut', alpha=1.0)
    plt.xlabel('P_cal_etottracknorm')
    plt.ylabel('P_hgcer_npeSum')

    ax = f.add_subplot(243)
    ax.hist2d(tree['P_cal_etottracknorm'],tree['P_aero_npeSum'],bins=[c.setbin(tree['P_cal_etottracknorm'],400,0,4),c.setbin(tree['P_aero_npeSum'],400,0,100)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(tree['P_cal_etottracknorm'],"p_%scut_lumi" % SHMS_PID),c.add_cut(tree['P_aero_npeSum'],"p_%scut_lumi" % SHMS_PID),bins=[c.setbin(tree['P_cal_etottracknorm'],400,0,4),c.setbin(tree['P_aero_npeSum'],400,0,100)],cmin=1,label='cut', alpha=1.0)
    plt.xlabel('P_cal_etottracknorm')
    plt.ylabel('P_aero_npeSum')

    if ANATYPE == "Pion":
        ax = f.add_subplot(244)
        ax.hist2d(tree['P_cal_etottracknorm'],tree['P_ngcer_npeSum'],bins=[c.setbin(tree['P_cal_etottracknorm'],400,0,4),c.setbin(tree['P_ngcer_npeSum'],400,0,80)],cmin=1,label='no cut',alpha=0.5)
        ax.hist2d(c.add_cut(tree['P_cal_etottracknorm'],"p_%scut_lumi" % SHMS_PID),c.add_cut(tree['P_ngcer_npeSum'],"p_%scut_lumi" % SHMS_PID),bins=[c.setbin(tree['P_cal_etottracknorm'],400,0,4),c.setbin(tree['P_ngcer_npeSum'],400,0,80)],cmin=1,label='cut', alpha=1.0)
        plt.xlabel('P_cal_etottracknorm')
        plt.ylabel('P_ngcer_npeSum')

    ax = f.add_subplot(245)
    ax.hist2d(tree['P_aero_npeSum'],tree['P_hgcer_npeSum'],bins=[c.setbin(tree['P_aero_npeSum'],400,0,100),c.setbin(tree['P_hgcer_npeSum'],400,0,80)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(tree['P_aero_npeSum'],"p_%scut_lumi" % SHMS_PID),c.add_cut(tree['P_hgcer_npeSum'],"p_%scut_lumi" % SHMS_PID),bins=[c.setbin(tree['P_aero_npeSum'],400,0,100),c.setbin(tree['P_hgcer_npeSum'],400,0,80)],cmin=1,label='cut', alpha=1.0)
    plt.xlabel('P_aero_npeSum')
    plt.ylabel('P_hgcer_npeSum')

    if ANATYPE == "Pion":
        ax = f.add_subplot(246)
        ax.hist2d(tree['P_ngcer_npeSum'],tree['P_hgcer_npeSum'],bins=[c.setbin(tree['P_ngcer_npeSum'],400,0,80),c.setbin(tree['P_hgcer_npeSum'],400,0,80)],cmin=1,label='no cut',alpha=0.5)
        ax.hist2d(c.add_cut(tree['P_ngcer_npeSum'],"p_%scut_lumi" % SHMS_PID),c.add_cut(tree['P_hgcer_npeSum'],"p_%scut_lumi" % SHMS_PID),bins=[c.setbin(tree['P_ngcer_npeSum'],400,0,80),c.setbin(tree['P_hgcer_npeSum'],400,0,80)],cmin=1,label='cut', alpha=1.0)
        plt.xlabel('P_ngcer_npeSum')
        plt.ylabel('P_hgcer_npeSum')

        ax = f.add_subplot(247)
        ax.hist2d(tree['P_aero_npeSum'],tree['P_ngcer_npeSum'],bins=[c.setbin(tree['P_aero_npeSum'],400,0,100),c.setbin(tree['P_ngcer_npeSum'],400,0,80)],cmin=1,label='no cut',alpha=0.5)
        ax.hist2d(c.add_cut(tree['P_aero_npeSum'],"p_%scut_lumi" % SHMS_PID),c.add_cut(tree['P_ngcer_npeSum'],"p_%scut_lumi" % SHMS_PID),bins=[c.setbin(tree['P_aero_npeSum'],400,0,100),c.setbin(tree['P_ngcer_npeSum'],400,0,80)],cmin=1,label='cut', alpha=1.0)
        plt.xlabel('P_aero_npeSum')
        plt.ylabel('P_ngcer_npeSum')


    ax = f.add_subplot(248)
    plt.axis('off')
    i=0
    plt.text(-0.15,1.00,"HMS cuts... (%s)" % HMS_PID,fontsize=8)
    for key,val in strDict.items():
        if key == "h_%scut_lumi" % HMS_PID:
            for v in val:
                plt.text(-0.15,0.95-i/10," {}".format(v),fontsize=8)
                i+=1
    plt.text(-0.15,0.95-((i)/10+0.05),"SHMS cuts... (%s)" % SHMS_PID,fontsize=8)
    for key,val in strDict.items():
        if key == "p_%scut_lumi" % SHMS_PID:
            for v in val:
                plt.text(-0.15,0.95-(i+1)/10," {}".format(v),fontsize=8)
                i+=1
        if key == "c_curr":
            global thres_curr, report_current
            # e.g. Grabbing threshold current (ie 2.5) from something like this [' {"H_bcm_bcm1_AvgCurrent" : (abs(H_bcm_bcm1_AvgCurrent-55) < 2.5)}']
            thres_curr = float(val[0].split(":")[1].split("<")[1].split(")")[0].strip())
            # e.g. Grabbing set current for run (ie 55) from something like this [' {"H_bcm_bcm1_AvgCurrent" : (abs(H_bcm_bcm1_AvgCurrent-55) < 2.5)}']
            report_current = float(val[0].split(":")[1].split("<")[0].split(")")[0].split("-")[1].strip())

    plt.tight_layout(rect=[0,0.03,1,0.95])   
    plt.savefig(UTILPATH+'/scripts/luminosity/OUTPUTS/plots/track_pid/track_pid2D_%s.png' % (runNum))


    #######################
    ######## Beta  ########
    #######################

    f = plt.figure(figsize=(19.20,8.00))
    f.suptitle("Run %s" % runNum)

    ax = f.add_subplot(241)
    ax.hist2d(tree['H_cer_npeSum'],tree['P_gtr_beta'],bins=[c.setbin(tree['H_cer_npeSum'],400,0,80),c.setbin(tree['P_gtr_beta'],400,-2.0,2.0)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(tree['H_cer_npeSum'],"p_%scut_lumi" % SHMS_PID),c.add_cut(tree['P_gtr_beta'],"p_%scut_lumi" % SHMS_PID),bins=[c.setbin(tree['H_cer_npeSum'],400,0,80),c.setbin(tree['P_gtr_beta'],400,-2.0,2.0)],cmin=1,label='cut', alpha=1.0)
    plt.xlabel('H_cer_npeSum')
    plt.ylabel('P_gtr_beta')

    ax = f.add_subplot(242)
    ax.hist2d(tree['H_cal_etottracknorm'],tree['P_gtr_beta'],bins=[c.setbin(tree['H_cal_etottracknorm'],400,0,4),c.setbin(tree['P_gtr_beta'],400,-2.0,2.0)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(tree['H_cal_etottracknorm'],"p_%scut_lumi" % SHMS_PID),c.add_cut(tree['P_gtr_beta'],"p_%scut_lumi" % SHMS_PID),bins=[c.setbin(tree['H_cal_etottracknorm'],400,0,4),c.setbin(tree['P_gtr_beta'],400,-2.0,2.0)],cmin=1,label='cut', alpha=1.0)
    plt.xlabel('H_cal_etottracknorm')
    plt.ylabel('P_gtr_beta')

    ax = f.add_subplot(243)
    ax.hist(tree['P_gtr_beta'],bins=c.setbin(tree['P_gtr_beta'],200,-2.0,2.0),label='no cut',histtype='step',alpha=0.5, stacked=True, fill=True)
    ax.hist(c.add_cut(tree['P_gtr_beta'],"h_cal"),bins=c.setbin(tree['P_gtr_beta'],200,-2.0,2.0),label='cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.yscale('log')
    plt.xlabel('P_gtr_beta')
    plt.ylabel('Count')

    ax = f.add_subplot(244)
    ax.hist2d(tree['P_cal_etottracknorm'],tree['P_gtr_beta'],bins=[c.setbin(tree['P_cal_etottracknorm'],400,0,4),c.setbin(tree['P_gtr_beta'],400,-2.0,2.0)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(tree['P_cal_etottracknorm'],"p_%scut_lumi" % SHMS_PID),c.add_cut(tree['P_gtr_beta'],"p_%scut_lumi" % SHMS_PID),bins=[c.setbin(tree['P_cal_etottracknorm'],400,0,4),c.setbin(tree['P_gtr_beta'],400,-2.0,2.0)],cmin=1,label='cut', alpha=1.0)
    plt.xlabel('P_cal_etottracknorm')
    plt.ylabel('P_gtr_beta')

    ax = f.add_subplot(245)
    ax.hist2d(tree['P_hgcer_npeSum'],tree['P_gtr_beta'],bins=[c.setbin(tree['P_hgcer_npeSum'],400,0,80),c.setbin(tree['P_gtr_beta'],400,-2.0,2.0)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(tree['P_hgcer_npeSum'],"p_%scut_lumi" % SHMS_PID),c.add_cut(tree['P_gtr_beta'],"p_%scut_lumi" % SHMS_PID),bins=[c.setbin(tree['P_hgcer_npeSum'],400,0,80),c.setbin(tree['P_gtr_beta'],400,-2.0,2.0)],cmin=1,label='cut', alpha=1.0)
    plt.xlabel('P_hgcer_npeSum')
    plt.ylabel('P_gtr_beta')

    ax = f.add_subplot(246)
    ax.hist2d(tree['P_aero_npeSum'],tree['P_gtr_beta'],bins=[c.setbin(tree['P_aero_npeSum'],400,0,100),c.setbin(tree['P_gtr_beta'],400,-2.0,2.0)],cmin=1,label='no cut',alpha=0.5)
    ax.hist2d(c.add_cut(tree['P_aero_npeSum'],"p_%scut_lumi" % SHMS_PID),c.add_cut(tree['P_gtr_beta'],"p_%scut_lumi" % SHMS_PID),bins=[c.setbin(tree['P_aero_npeSum'],400,0,100),c.setbin(tree['P_gtr_beta'],400,-2.0,2.0)],cmin=1,label='cut', alpha=1.0)
    plt.xlabel('P_aero_npeSum')
    plt.ylabel('P_gtr_beta')

    ax = f.add_subplot(248)
    plt.axis('off')
    i=0
    plt.text(-0.15,1.00,"HMS cuts... (%s)" % HMS_PID,fontsize=8)
    for key,val in strDict.items():
        if key == "h_%scut_lumi_nt" % HMS_PID:
            for v in val:
                plt.text(-0.15,0.95-i/10," {}".format(v),fontsize=8)
                i+=1
    plt.text(-0.15,0.95-((i)/10+0.05),"SHMS cuts... (%s)" % SHMS_PID,fontsize=8)
    for key,val in strDict.items():
        if key == "p_%scut_lumi_nt" % SHMS_PID:
            for v in val:
                plt.text(-0.15,0.95-(i+1)/10," {}".format(v),fontsize=8)
                i+=1
        if key == "c_curr":
            #global thres_curr, report_current
            # e.g. Grabbing threshold current (ie 2.5) from something like this [' {"H_bcm_bcm1_AvgCurrent" : (abs(H_bcm_bcm1_AvgCurrent-55) < 2.5)}']
            thres_curr = float(val[0].split(":")[1].split("<")[1].split(")")[0].strip())
            # e.g. Grabbing set current for run (ie 55) from something like this [' {"H_bcm_bcm1_AvgCurrent" : (abs(H_bcm_bcm1_AvgCurrent-55) < 2.5)}']
            report_current = float(val[0].split(":")[1].split("<")[0].split(")")[0].split("-")[1].strip())


    plt.tight_layout(rect=[0,0.03,1,0.95])   
    plt.savefig(UTILPATH+'/scripts/luminosity/OUTPUTS/plots/track_pid/track_beta_%s.png' % (runNum))


################################################################################################################################################
    
def analysis():
    '''
    Calculate variables for table
    '''

    # Applies edtm window cuts to edtm to get accepted edtm events
    EDTM = c.add_cut(T_coin_pEDTM_tdcTimeRaw,"c_edtm")
    EDTM_SHMS = [x
                 for (x,evt) in zip(c.add_cut(T_coin_pEDTM_tdcTimeRaw,"c_edtm"),tree["EvtType"])
                 if evt == 1 or evt == 3]
    EDTM_HMS = [x
                for (x,evt) in zip(c.add_cut(T_coin_pEDTM_tdcTimeRaw,"c_edtm"),tree["EvtType"])
                if evt == 2 or evt == 3]

    for ps in PS_names:
        if ps == "PS1" or ps == "PS2":
            # Applies trigger window cuts to trigger to get accepted trigger events
            SHMSTRIG = [x
                        for x in c.add_cut(T_coin_pTRIG_SHMS_ROC2_tdcTime,"c_ptrigSHMS%s" % ps.replace("PS",""))
                        if x != 0.0]

        if ps == "PS3" or ps == "PS4":
            # Applies trigger window cuts to trigger to get accepted trigger events
            HMSTRIG  = [x
                        for x in c.add_cut(T_coin_pTRIG_HMS_ROC1_tdcTime,"c_ptrigHMS%s" % ps.replace("PS",""))
                        if x !=0.0]
            
        if ps == "PS5" or ps == "PS6":
            # Applies trigger window cuts to trigger to get accepted trigger events
            COINTRIG  = [x
                         for x in c.add_cut(T_coin_pTRIG_COIN_ROC1_tdcTime,"c_ptrigCOIN%s" % ps.replace("PS",""))
                         if x !=0.0]

    try:
        SHMSTRIG
    except NameError:
        SHMSTRIG = []

    try:
        HMSTRIG
    except NameError:
        HMSTRIG = []

    try:
        COINTRIG
    except NameError:
        COINTRIG = []

    # Cuts event type with current cut, probably pointless?
    EventType = c.add_cut(tree["EvtType"],"c_curr")

    for ps in PS_names:
        if ps == "PS1" or ps == "PS2":
            # Applies trigger window cuts to trigger to get accepted trigger events
            SHMSTRIG_cut = [trig1
                            for (trig1,evt) in zip(c.add_cut(T_coin_pTRIG_SHMS_ROC2_tdcTime,"c_ptrigSHMS%s" % ps.replace("PS","")),tree["EvtType"])
                            if evt == 1 or evt == 3]

        if ps == "PS3" or ps == "PS4":
            # Applies trigger window cuts to trigger to get accepted trigger events
            HMSTRIG_cut = [ x
                            for (x, evt) in zip(c.add_cut(T_coin_pTRIG_HMS_ROC1_tdcTime,"c_ptrigHMS%s" % ps.replace("PS","")), tree["EvtType"])
                            if evt == 2 or evt == 3] # 3/21/2023, added evttype 3 as suggested by Carlos
            
        if ps == "PS5" or ps == "PS6":
            # Applies trigger window cuts to trigger to get accepted trigger events
            COINTRIG_cut = [ x
                             for (x, evt) in zip(c.add_cut(T_coin_pTRIG_COIN_ROC1_tdcTime,"c_ptrigCOIN%s" % ps.replace("PS","")), tree["EvtType"])
                             if (evt == 1 or evt == 2)]

    try:
        SHMSTRIG_cut
    except NameError:
        SHMSTRIG_cut = []

    try:
        HMSTRIG_cut
    except NameError:
        HMSTRIG_cut = []

    try:
        COINTRIG_cut
    except NameError:
        COINTRIG_cut = []

    h_et_should = len(c.add_cut(tree["H_cal_etotnorm"],"h_%strack_lumi_before" % HMS_PID))
    h_et_did = len(c.add_cut(tree["H_cal_etotnorm"],"h_%strack_lumi_after" % HMS_PID))
    #HMS_track_eff = h_et_did/h_et_should

    p_et_should = len(c.add_cut(tree["P_cal_etotnorm"],"p_%strack_lumi_before" % SHMS_PID))
    p_et_did = len(c.add_cut(tree["P_cal_etotnorm"],"p_%strack_lumi_after" % SHMS_PID))
    #SHMS_track_eff = p_et_did/p_et_should

    # Applies PID cuts, once integrated this will give the events (no track)
    h_etotnorm = c.add_cut(tree["H_cal_etotnorm"],"h_%scut_lumi_nt" % HMS_PID) 
    p_etotnorm = c.add_cut(tree["P_cal_etotnorm"],"p_%scut_lumi_nt" % SHMS_PID)

    # Applies PID cuts, once integrated this will give the events (track)
    #h_goodscinhit = c.add_cut(tree["H_hod_goodscinhit"],"h_%scut_lumi" % HMS_PID)
    #p_goodscinhit = c.add_cut(tree["P_hod_goodscinhit"],"p_%scut_lumi" % SHMS_PID)
    h_etottracknorm = c.add_cut(tree["H_cal_etottracknorm"],"h_%scut_lumi" % HMS_PID)
    p_etottracknorm = c.add_cut(tree["P_cal_etottracknorm"],"p_%scut_lumi" % SHMS_PID)    
    
    # Creates a dictionary for the calculated luminosity values 
    track_info = {

        "tot_events" : len(EventType),
        "h_int_etotnorm_evts" : scipy.integrate.simps(h_etotnorm),
        "p_int_etotnorm_evts" : scipy.integrate.simps(p_etotnorm),
        "h_int_etottracknorm_evts" : scipy.integrate.simps(h_etottracknorm),
        "p_int_etottracknorm_evts" : scipy.integrate.simps(p_etottracknorm),                
        "accp_edtm" : (len(EDTM)),
        "paccp_edtm" : (len(EDTM_SHMS)),
        "haccp_edtm" : (len(EDTM_HMS)),

    }

    track_info.update(PSDict)
    track_info.update({"SHMSTRIG_cut" : len(SHMSTRIG_cut)})
    track_info.update({"SHMS_track" : SHMS_track_eff})
    track_info.update({"SHMS_track_uncern" : SHMS_track_uncern})
    track_info.update({"HMSTRIG_cut" : len(HMSTRIG_cut)})
    track_info.update({"HMS_track" : HMS_track_eff})
    track_info.update({"HMS_track_uncern" : HMS_track_uncern})
    track_info.update({"COINTRIG_cut" : len(COINTRIG_cut)})
    
    print("\nTerminate","Selection rules have been applied, plotting results")
    print("Total number of events: %.0f" % (track_info['tot_events']))
    print("Number of EDTM  Events: %.0f" % (track_info['accp_edtm']))
    for ps in PS_names:
        if ps == "PS1" or ps == "PS2":
            print("Number of SHMSTRIG Events: %.0f" % (SHMS_PS*track_info['SHMSTRIG_cut']))
            print("Number of SHMS EDTM  Events: %.0f" % (track_info['paccp_edtm']))
            print("Number of SHMS good events: %.0f +/- %.0f " % ((SHMS_PS*track_info['h_int_etottracknorm_evts']), math.sqrt(SHMS_PS*track_info['h_int_etottracknorm_evts'])))
            print("Calculated SHMS tracking efficiency: %f +/- %f\n" % ((track_info['SHMS_track']), (track_info['SHMS_track_uncern'])))            
        if ps == "PS3" or ps == "PS4":
            print("Number of HMSTRIG Events: %.0f" % (HMS_PS*track_info['HMSTRIG_cut']))
            print("Number of HMS EDTM  Events: %.0f" % (track_info['haccp_edtm']))
            print("\nNumber of HMS good events: %.0f +/- %.0f " % ((HMS_PS*track_info['h_int_etottracknorm_evts']), math.sqrt(HMS_PS*track_info['h_int_etottracknorm_evts'])))
            print("Calculated HMS tracking efficiency: %f +/- %f\n" % ((track_info['HMS_track']), (track_info['HMS_track_uncern'])))
        if ps == "PS5" or ps == "PS6":
            print("Number of COINTRIG Events: %.0f" % (COIN_PS*track_info['COINTRIG_cut']))

    print("============================================================================\n\n")
          
    return track_info

################################################################################################################################################

def main():

    pid_cuts()
    track_pid_cuts()
    #plt.show()

    # lumi_data = {**scalers , **track_info} # only python 3.5+

    # Import dictionaries
    scalers = scaler.scaler(PS_names, HMS_PS, SHMS_PS, COIN_PS, thres_curr, report_current, runNum, MaxEvent, s_tree) 
    track_info = analysis()

    # Merge and sort the two dictionaries of calculations
    data = {}
    for d in (scalers, track_info): 
        data.update(d)
    lumi_data = {i : data[i] for i in sorted(data.keys())}

    # Convert merged dictionary to a pandas dataframe then sort it
    table  = pd.DataFrame([lumi_data], columns=lumi_data.keys())
    table = table.reindex(sorted(table.columns), axis=1)
    
    # Replace zeros with NaN
    table.replace(0,np.nan,inplace=True)
    # Replace None with NaN
    table.replace([None],np.nan,inplace=True)

    file_exists = os.path.isfile(out_f)

    # Updates csv file with luminosity calculated values for later analysis (see plot_yield.py)
    if file_exists:
        try:
            out_data = pd.read_csv(out_f)
        except IOError:
            print("Error: %s does not appear to exist." % out_f)
        # Checks if run number is alread in csv and replaces it if it is there
        run_index = out_data.index[out_data['run number'] == int(runNum)].tolist()
        out_data.drop(run_index, inplace=True)
        out_data = out_data.append(table,ignore_index=True)
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print("Output luminosity values\n",out_data)
        out_data.to_csv(out_f, index = False, header=True, mode='w+',)
    else:
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print("Output luminosity values\n",table)
        table.to_csv(out_f, index = False, header=True, mode='a',)

if __name__ == '__main__':
    main()
