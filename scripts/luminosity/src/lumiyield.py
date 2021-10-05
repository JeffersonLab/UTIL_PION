#! /usr/bin/python
# Description:
# ================================================================
# Time-stamp: "2021-10-05 01:57:04 trottar"
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
    REPLAYPATH="/group/c-pionlt/online_analysis/hallc_replay_lt"
elif ("lark" in HOST[1]):
    REPLAYPATH = "/home/%s/work/JLab/hallc_replay_lt" % USER[1]
elif ("cdaq" in HOST[1]):
    REPLAYPATH = "/home/cdaq/hallc-online/hallc_replay_lt"
elif ("trottar" in HOST[1]):
    REPLAYPATH = "/home/trottar/Analysis/hallc_replay_lt"

sys.path.insert(0, '%s/UTIL_PION/bin/python/' % REPLAYPATH)
import kaonlt as klt
import scaler
#import scaler_nocut as scaler

print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER[1], HOST[1], REPLAYPATH))

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
T_coin_pEDTM_tdcTimeRaw = tree.array("T.coin.pEDTM_tdcTimeRaw")
EvtType = tree.array("fEvtHdr.fEvtType")

fout = REPLAYPATH+'/UTIL_PION/DB/CUTS/run_type/lumi.cuts'

# read in cuts file and make dictionary
c = klt.pyPlot(REPLAYPATH)
# apply RF cuts to timing cuts file
#c.cut_RF(runNum,MaxEvent)
readDict = c.read_dict(fout,runNum)

def make_cutDict(cut,inputDict=None):
    '''
    This method calls several methods in kaonlt package. It is required to create properly formated
    dictionaries. The evaluation must be in the analysis script because the analysis variables (i.e. the
    leaves of interest) are not defined in the kaonlt package. This makes the system more flexible
    overall, but a bit more cumbersome in the analysis script. Perhaps one day a better solution will be
    implimented.
    '''

    global c

    c = klt.pyPlot(REPLAYPATH,readDict)
    x = c.w_dict(cut)
    print("\n%s" % cut)
    print(x, "\n")

    # Threshold current
    if cut == "c_curr":
        global thres_curr, report_current
        # e.g. Grabbing threshold current (ie 2.5) from something like this [' {"H_bcm_bcm4a_AvgCurrent" : (abs(H_bcm_bcm4a_AvgCurrent-55) < 2.5)}']
        thres_curr = float(x[0].split(":")[1].split("<")[1].split(")")[0].strip())
        # e.g. Grabbing set current for run (ie 55) from something like this [' {"H_bcm_bcm4a_AvgCurrent" : (abs(H_bcm_bcm4a_AvgCurrent-55) < 2.5)}']
        report_current = float(x[0].split(":")[1].split("<")[0].split(")")[0].split("-")[1].strip())
    
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

cutDict = make_cutDict("h_cal")
cutDict = make_cutDict("h_cer",cutDict)
cutDict = make_cutDict("p_cal",cutDict)
cutDict = make_cutDict("p_hgcer",cutDict)
cutDict = make_cutDict("p_aero",cutDict)
cutDict = make_cutDict("p_ecut_lumi_eff",cutDict)
cutDict = make_cutDict("p_picut_lumi_eff",cutDict)
cutDict = make_cutDict("p_kcut_lumi_eff",cutDict)
cutDict = make_cutDict("p_pcut_lumi_eff",cutDict)
cutDict = make_cutDict("p_hadcut_lumi_eff",cutDict)
cutDict = make_cutDict("h_ecut_lumi_eff",cutDict)
cutDict = make_cutDict("h_picut_lumi_eff",cutDict)
cutDict = make_cutDict("h_hadcut_lumi_eff",cutDict)
cutDict = make_cutDict("c_noedtm",cutDict)
cutDict = make_cutDict("c_edtm",cutDict)
cutDict = make_cutDict("c_ptrigHMS",cutDict)
cutDict = make_cutDict("c_ptrigSHMS",cutDict)
if len(PS_used) > 2:
    cutDict = make_cutDict("c_ptrigCOIN",cutDict)
cutDict = make_cutDict("c_curr",cutDict)
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

    EDTM = c.add_cut(T_coin_pEDTM_tdcTimeRaw,"c_edtm")
    
    SHMSTRIG = [x
                for x in c.add_cut(T_coin_pTRIG_SHMS_ROC2_tdcTime,"c_ptrigSHMS")
                if x != 0.0]
    HMSTRIG  = [x
                for x in c.add_cut(T_coin_pTRIG_HMS_ROC1_tdcTime,"c_ptrigHMS")
                if x !=0.0]

    if len(PS_used) > 2:
        COINTRIG  = [x
                     for x in c.add_cut(T_coin_pTRIG_COIN_ROC1_tdcTime,"c_ptrigCOIN")
                     if x !=0.0]
    
    EventType = c.add_cut(EvtType,"c_curr")

    SHMSTRIG_cut = [trig1
                    for (trig1,evt) in zip(c.add_cut(T_coin_pTRIG_SHMS_ROC2_tdcTime,"c_ptrigSHMS"),EvtType)
                    if evt == 1]

    HMSTRIG_cut = [ x
                    for (x, evt) in zip(c.add_cut(T_coin_pTRIG_HMS_ROC1_tdcTime,"c_ptrigHMS"), EvtType)
                    if evt == 2]

    # goodscinhit cut
    h_hadcuts_goodscinhit = c.add_cut(H_hod_goodscinhit,"h_hadcut_lumi_eff")
    p_pcuts_goodscinhit = c.add_cut(P_hod_goodscinhit,"p_pcut_lumi_eff")
    
    track_info = {
        
        "tot_events" : len(EventType),
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

def main():

    #id_cuts()
    #plt.show()

    # combine dictionaries
    scalers = scaler.scaler(PS_names, HMS_PS, SHMS_PS, thres_curr, report_current, REPLAYPATH, runNum, MaxEvent, s_tree, s_branch) 
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
