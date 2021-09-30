#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2021-09-30 07:04:51 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from csv import DictReader
import sys, os, subprocess

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

inp_name = sys.argv[1]
if "1" in inp_name:
    if "LH2" in inp_name.upper():
        target = "LH2"
        inp_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_1/LH2/lumi_data_l1_lh2.csv" % str(REPLAYPATH)
        out_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_1/LH2/yield_data_l1_lh2.csv" % str(REPLAYPATH)
        print("\nGrabbing input...\n\n%s" % str(inp_f))
    if "LD2" in inp_name.upper():
        target = "LD2"
        inp_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_1/LD2/lumi_data_l1_ld2.csv" % str(REPLAYPATH)
        out_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_1/LD2/yield_data_l1_ld2.csv" % str(REPLAYPATH)
        print("\nGrabbing input...\n\n%s" % str(inp_f))
    if "C" in inp_name.upper():
        target = "carbon"
        inp_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_1/Carbon0p5/lumi_data_l1_c0p5.csv" % str(REPLAYPATH)
        out_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_1/Carbon0p5/yield_data_l1_c0p5.csv" % str(REPLAYPATH)
        print("\nGrabbing input...\n\n%s" % str(inp_f))
elif "2" in inp_name:
    if "LH2" in inp_name.upper():
        target = "LH2"
        inp_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_2/LH2/lumi_data_l2_lh2.csv" % str(REPLAYPATH)
        out_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_2/LH2/yield_data_l2_lh2.csv" % str(REPLAYPATH)
        print("\nGrabbing input...\n\n%s" % str(inp_f))
    if "LD2" in inp_name.upper():
        target = "LD2"
        inp_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_2/LD2/lumi_data_l2_ld2.csv" % str(REPLAYPATH)
        out_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_2/LD2/yield_data_l2_ld2.csv" % str(REPLAYPATH)
        print("\nGrabbing input...\n\n%s" % str(inp_f))
    if "C" in inp_name.upper():
        target = "carbon"
        inp_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_2/Carbon0p5/lumi_data_l2_c0p5.csv" % str(REPLAYPATH)
        out_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_2/Carbon0p5/yield_data_l2_c0p5.csv" % str(REPLAYPATH)
        print("\nGrabbing input...\n\n%s" % str(inp_f))
elif inp_name == None:
    inp_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/lumi_data.csv" % str(REPLAYPATH)
    out_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/yield_data.csv" % str(REPLAYPATH)
    print("\nError: Invalid input...\nGrabbing default input...\n\n%s" % str(inp_f))
else:
    inp_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/lumi_data.csv" % str(REPLAYPATH)
    out_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/yield_data.csv" % str(REPLAYPATH)
    print("\nGrabbing default input...\n\n%s" % str(inp_f))

print("\nRunning as %s on %s, hallc_replay_lt path assumed as %s" % (USER[1], HOST[1], REPLAYPATH))

try:
    lumi_data = dict(pd.read_csv(inp_f))
except IOError:
    print("Error: %s does not appear to exist." % inp_f)
print(lumi_data.keys())
    
# prints first instance of run number-> print(lumi_data["run number"][0])

if "PS1" in lumi_data.keys():
    SHMS_PS = lumi_data["PS1"]
if "PS2" in lumi_data.keys():
    SHMS_PS = lumi_data["PS2"]
if "PS3" in lumi_data.keys():
    HMS_PS = lumi_data["PS3"]
if "PS4" in lumi_data.keys():
    HMS_PS = lumi_data["PS4"]    

numRuns = len(lumi_data["run number"])

def makeList(lumi_input):
    new_lst = [lumi_data[lumi_input][i] for i,evts in enumerate(lumi_data["run number"])]
    new_lst = np.asarray(pd.Series(new_lst).fillna(0)) # changes NaN to zeros and convert to numpy
    return new_lst

def calc_yield():

    yield_dict = {
        "current" : makeList("charge")/makeList("time"),
    
        "rate_HMS" : makeList("HMSTRIG_scaler")/makeList("time"),
        "rate_SHMS" : makeList("SHMSTRIG_scaler")/makeList("time"),
        
        "TLT" : 1-(makeList("accp_edtm")/makeList("sent_edtm")),

        "uncern_HMS_evts_scaler" : np.sqrt(makeList("HMSTRIG_scaler"))/makeList("HMSTRIG_scaler"),

        "uncern_SHMS_evts_scaler" : np.sqrt(makeList("SHMSTRIG_scaler"))/makeList("SHMSTRIG_scaler"),

        "uncern_HMS_evts_notrack" : np.sqrt(makeList("h_int_goodscin_evts"))/makeList("h_int_goodscin_evts"),

        "uncern_SHMS_evts_notrack" : np.sqrt(makeList("p_int_goodscin_evts"))/makeList("p_int_goodscin_evts"),

        "uncern_HMS_evts_track" : np.sqrt(makeList("h_int_goodscin_evts"))/makeList("h_int_goodscin_evts"),

        "uncern_SHMS_evts_track" : np.sqrt(makeList("p_int_goodscin_evts"))/makeList("p_int_goodscin_evts"),

        "HMS_scaler_accp" : makeList("HMSTRIG_scaler")-makeList("sent_edtm"),

        "SHMS_scaler_accp" : makeList("SHMSTRIG_scaler")-makeList("sent_edtm"),
    }

    # Calculate yield values

    yield_HMS_scaler = (yield_dict["HMS_scaler_accp"])/(makeList("charge")*makeList("CPULT_scaler")*makeList("HMS_eLT"))
    yield_HMS_notrack = makeList("h_int_goodscin_evts")/(makeList("charge")*yield_dict["TLT"])
    yield_HMS_track = makeList("h_int_goodscin_evts")/(makeList("charge")*yield_dict["TLT"]*makeList("HMS_track"))
    yield_dict.update({"yield_HMS_scaler" : yield_HMS_scaler})
    yield_dict.update({"yield_HMS_notrack" : yield_HMS_notrack})
    yield_dict.update({"yield_HMS_track" : yield_HMS_track})

    yield_SHMS_scaler = (yield_dict["SHMS_scaler_accp"])/(makeList("charge")*makeList("CPULT_scaler")*makeList("SHMS_eLT"))
    yield_SHMS_notrack = makeList("p_int_goodscin_evts")/(makeList("charge")*yield_dict["TLT"])
    yield_SHMS_track = makeList("p_int_goodscin_evts")/(makeList("charge")*yield_dict["TLT"]*makeList("SHMS_track"))
    yield_dict.update({"yield_SHMS_scaler" : yield_SHMS_scaler})
    yield_dict.update({"yield_SHMS_notrack" : yield_SHMS_notrack})
    yield_dict.update({"yield_SHMS_track" : yield_SHMS_track})

    for i,curr in enumerate(yield_dict["current"]):
        if curr == max(yield_dict["current"]):
            max_yield_HMS_scaler = yield_dict["yield_HMS_scaler"][i]
            max_yield_SHMS_scaler = yield_dict["yield_SHMS_scaler"][i]
    yield_dict.update({"max_yield_HMS_scaler" : max_yield_HMS_scaler})
    yield_dict.update({"max_yield_SHMS_scaler" : max_yield_SHMS_scaler})
                
    for i,curr in enumerate(yield_dict["current"]):
        if curr == max(yield_dict["current"]):
            max_yield_HMS_notrack = yield_dict["yield_HMS_notrack"][i]
            max_yield_SHMS_notrack = yield_dict["yield_SHMS_notrack"][i]
    yield_dict.update({"max_yield_HMS_notrack" : max_yield_HMS_notrack})
    yield_dict.update({"max_yield_SHMS_notrack" : max_yield_SHMS_notrack})

    for i,curr in enumerate(yield_dict["current"]):
        if curr == max(yield_dict["current"]):
            max_yield_HMS_track = yield_dict["yield_HMS_track"][i]
            max_yield_SHMS_track = yield_dict["yield_SHMS_track"][i]
    yield_dict.update({"max_yield_HMS_track" : max_yield_HMS_track})
    yield_dict.update({"max_yield_SHMS_track" : max_yield_SHMS_track})

    yieldRel_HMS_scaler = yield_dict["yield_HMS_scaler"]/yield_dict["max_yield_HMS_scaler"]
    yieldRel_HMS_notrack = yield_dict["yield_HMS_notrack"]/yield_dict["max_yield_HMS_notrack"]
    yieldRel_HMS_track = yield_dict["yield_HMS_track"]/yield_dict["max_yield_HMS_track"]
    yield_dict.update({"yieldRel_HMS_scaler" : yieldRel_HMS_scaler})
    yield_dict.update({"yieldRel_HMS_notrack" : yieldRel_HMS_notrack})
    yield_dict.update({"yieldRel_HMS_track" : yieldRel_HMS_track})

    yieldRel_SHMS_scaler = yield_dict["yield_SHMS_scaler"]/yield_dict["max_yield_SHMS_scaler"]
    yieldRel_SHMS_notrack = yield_dict["yield_SHMS_notrack"]/yield_dict["max_yield_SHMS_notrack"]
    yieldRel_SHMS_track = yield_dict["yield_SHMS_track"]/yield_dict["max_yield_SHMS_track"]
    yield_dict.update({"yieldRel_SHMS_scaler" : yieldRel_SHMS_scaler})
    yield_dict.update({"yieldRel_SHMS_notrack" : yieldRel_SHMS_notrack})
    yield_dict.update({"yieldRel_SHMS_track" : yieldRel_SHMS_track})

    yield_table = pd.DataFrame(yield_dict, columns=yield_dict.keys())
    yield_table = yield_table.reindex(sorted(yield_table.columns), axis=1)

    return yield_table

def mergeDicts():
    yield_data = calc_yield()
    # data = {**lumi_data, **yield_data} # only python 3.5+
    
    for key, val in lumi_data.items():
        lumi_data[key] = val

    datadict = {}
    for d in (lumi_data, yield_data): 
        datadict.update(d)
    data = {i : datadict[i] for i in sorted(datadict.keys())}

    table  = pd.DataFrame(data, columns=data.keys())
    table = table.reindex(sorted(table.columns), axis=1)

    return table

def plot_yield():

    yield_data = mergeDicts()

    for i, val in enumerate(yield_data["run number"]):
        print("Run numbers:",yield_data["run number"][i],"Current Values:",yield_data["current"][i])
    
    relYieldPlot = plt.figure(figsize=(12,8))

    #HMS plot scaler
    plt.subplot(2,3,1)    
    plt.grid(zorder=1)
    plt.xlim(0,100)
    plt.ylim(0.9,1.1)
    plt.plot([0,100], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["current"],yield_data["yieldRel_HMS_scaler"],yerr=yield_data["uncern_HMS_evts_scaler"],color='black',linestyle='None',zorder=3)
    plt.scatter(yield_data["current"],yield_data["yieldRel_HMS_scaler"],color='blue',zorder=4)
    plt.ylabel('Rel. Yield Scaler', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    #HMS plot no track
    plt.subplot(2,3,2)    
    plt.grid(zorder=1)
    plt.xlim(0,100)
    #plt.ylim(0.9,1.1)
    plt.plot([0,100], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["current"],yield_data["yieldRel_HMS_notrack"],yerr=yield_data["uncern_HMS_evts_notrack"],color='black',linestyle='None',zorder=3)
    plt.scatter(yield_data["current"],yield_data["yieldRel_HMS_notrack"],color='blue',zorder=4)
    plt.ylabel('Rel. Yield no track', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    #HMS plot track
    plt.subplot(2,3,3)    
    plt.grid(zorder=1)
    plt.xlim(0,100)
    #plt.ylim(0.9,1.1)
    plt.plot([0,100], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["current"],yield_data["yieldRel_HMS_track"],yerr=yield_data["uncern_HMS_evts_track"],color='black',linestyle='None',zorder=3)
    plt.scatter(yield_data["current"],yield_data["yieldRel_HMS_track"],color='blue',zorder=4)
    plt.ylabel('Rel. Yield track', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

        
    #SHMS plot scaler
    plt.subplot(2,3,4)    
    plt.grid(zorder=1)
    plt.xlim(0,100)
    plt.ylim(0.9,1.1)
    plt.plot([0,100], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["current"],yield_data["yieldRel_SHMS_scaler"],yerr=yield_data["uncern_SHMS_evts_scaler"],color='black',linestyle='None',zorder=3)
    plt.scatter(yield_data["current"],yield_data["yieldRel_SHMS_scaler"],color='blue',zorder=4)
    plt.ylabel('Rel. Yield Scaler', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('SHMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('SHMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('SHMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    #SHMS plot no track
    plt.subplot(2,3,5)    
    plt.grid(zorder=1)
    plt.xlim(0,100)
    #plt.ylim(0.9,1.1)
    plt.plot([0,100], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["current"],yield_data["yieldRel_SHMS_notrack"],yerr=yield_data["uncern_SHMS_evts_notrack"],color='black',linestyle='None',zorder=3)
    plt.scatter(yield_data["current"],yield_data["yieldRel_SHMS_notrack"],color='blue',zorder=4)
    plt.ylabel('Rel. Yield no track', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('SHMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('SHMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('SHMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    #SHMS plot track
    plt.subplot(2,3,6)    
    plt.grid(zorder=1)
    plt.xlim(0,100)
    #plt.ylim(0.9,1.1)
    plt.plot([0,100], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["current"],yield_data["yieldRel_SHMS_track"],yerr=yield_data["uncern_SHMS_evts_track"],color='black',linestyle='None',zorder=3)
    plt.scatter(yield_data["current"],yield_data["yieldRel_SHMS_track"],color='blue',zorder=4)
    plt.ylabel('Rel. Yield track', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('SHMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('SHMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('SHMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    plt.tight_layout()      

    edtmPlot = plt.figure(figsize=(12,8))

    #EDTM vs Current
    plt.subplot(2,3,1)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.scatter(yield_data["current"],yield_data["accp_edtm"],color='blue',zorder=4)
    plt.ylabel('Accept EDTM', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    #HMS rate vs Current
    plt.subplot(2,3,2)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.scatter(yield_data["current"],yield_data["rate_HMS"]/1000,color='blue',zorder=4)
    plt.ylabel('HMS Rate [kHz]', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    #EDTM vs HMS Rate
    plt.subplot(2,3,3)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.scatter(yield_data["rate_HMS"]/1000,yield_data["accp_edtm"]/(yield_data["time"]*1000),color='blue',zorder=4)
    plt.ylabel('EDTM Rate [kHz]', fontsize=16)
    plt.xlabel('HMS Rate [kHz]', fontsize =16)
    if target == 'LD2' :
        plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    #TLT vs Current
    plt.subplot(2,3,4)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.scatter(yield_data["current"],yield_data["TLT"],color='blue',zorder=4)
    plt.ylabel('TLT', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    #SHMS rate vs Current
    plt.subplot(2,3,5)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.scatter(yield_data["current"],yield_data["rate_SHMS"]/1000,color='blue',zorder=4)
    plt.ylabel('SHMS Rate [kHz]', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('SHMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('SHMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('SHMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    #EDTM vs SHMS Rate
    plt.subplot(2,3,6)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.scatter(yield_data["rate_SHMS"]/1000,yield_data["accp_edtm"]/(yield_data["time"]*1000),color='blue',zorder=4)
    plt.ylabel('EDTM Rate [kHz]', fontsize=16)
    plt.xlabel('SHMS Rate [kHz]', fontsize =16)
    if target == 'LD2' :
        plt.title('SHMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('SHMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('SHMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    plt.tight_layout()      

    logPlot = plt.figure(figsize=(12,8))

    #HMS plot scaler
    plt.subplot(2,4,1)    
    plt.grid(zorder=1)
    plt.xlim(0,100)
    plt.ylim(0.9,1.1)
    plt.plot([0,100], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["current"],yield_data["yieldRel_HMS_scaler"],yerr=yield_data["uncern_HMS_evts_scaler"],color='black',linestyle='None',zorder=3)
    plt.scatter(yield_data["current"],yield_data["yieldRel_HMS_scaler"],color='blue',zorder=4)
    plt.ylabel('Rel. Yield Scaler', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =12)
    if target == 'LD2' :
        plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    elif target == 'LH2' :
        plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    else :
        plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)

    #HMS rate vs Current
    plt.subplot(2,4,2)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.scatter(yield_data["current"],yield_data["rate_HMS"]/1000,color='blue',zorder=4)
    plt.ylabel('HMS Rate [kHz]', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =12)
    if target == 'LD2' :
        plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    elif target == 'LH2' :
        plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    else :
        plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)

    #EDTM vs HMS Rate
    plt.subplot(2,4,3)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.scatter(yield_data["rate_HMS"]/1000,yield_data["accp_edtm"]/(yield_data["time"]*1000),color='blue',zorder=4)
    plt.ylabel('EDTM Rate [kHz]', fontsize=16)
    plt.xlabel('HMS Rate [kHz]', fontsize =12)
    if target == 'LD2' :
        plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    elif target == 'LH2' :
        plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    else :
        plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)

    #TLT vs Current
    plt.subplot(2,4,4)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.scatter(yield_data["current"],yield_data["TLT"],color='blue',zorder=4)
    plt.ylabel('TLT', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =12)
    if target == 'LD2' :
        plt.title('LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    elif target == 'LH2' :
        plt.title('LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    else :
        plt.title('Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)

    #Scaler CPULT vs Current
    plt.subplot(2,4,8)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.scatter(yield_data["current"],yield_data["CPULT_scaler"],color='blue',zorder=4)
    plt.ylabel('Scaler CPULT', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =12)
    if target == 'LD2' :
        plt.title('LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    elif target == 'LH2' :
        plt.title('LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    else :
        plt.title('Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)

    #SHMS plot scaler
    plt.subplot(2,4,5)    
    plt.grid(zorder=1)
    plt.xlim(0,100)
    plt.ylim(0.9,1.1)
    plt.plot([0,100], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["current"],yield_data["yieldRel_SHMS_scaler"],yerr=yield_data["uncern_SHMS_evts_scaler"],color='black',linestyle='None',zorder=3)
    plt.scatter(yield_data["current"],yield_data["yieldRel_SHMS_scaler"],color='blue',zorder=4)
    plt.ylabel('Rel. Yield Scaler', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =12)
    if target == 'LD2' :
        plt.title('SHMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    elif target == 'LH2' :
        plt.title('SHMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    else :
        plt.title('SHMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)

    #SHMS rate vs Current
    plt.subplot(2,4,6)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.scatter(yield_data["current"],yield_data["rate_SHMS"]/1000,color='blue',zorder=4)
    plt.ylabel('SHMS Rate [kHz]', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =12)
    if target == 'LD2' :
        plt.title('SHMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    elif target == 'LH2' :
        plt.title('SHMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    else :
        plt.title('SHMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)

    #EDTM vs SHMS Rate
    plt.subplot(2,4,7)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.scatter(yield_data["rate_SHMS"]/1000,yield_data["accp_edtm"]/(yield_data["time"]*1000),color='blue',zorder=4)
    plt.ylabel('EDTM Rate [kHz]', fontsize=16)
    plt.xlabel('SHMS Rate [kHz]', fontsize =12)
    if target == 'LD2' :
        plt.title('SHMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    elif target == 'LH2' :
        plt.title('SHMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)
    else :
        plt.title('SHMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =12)

    plt.tight_layout()            
        
    plt.show()

    print("\nYield info...\n",yield_data[["run number","yieldRel_HMS_scaler","yieldRel_SHMS_scaler","yieldRel_HMS_notrack","yieldRel_SHMS_notrack","yieldRel_HMS_track","yieldRel_SHMS_track"]])

    return yield_data

def debug():
    '''
    Prints various values to debug, customize to your heart's content
    '''
    data = mergeDicts()
    print("\n\n=======================")
    print("DEBUG data")
    print("=======================")
    ### Debug prints
    print(data[["run number","sent_edtm","TLT","CPULT_scaler","current","time"]])
    print("EDTM scaler rate: ", data["sent_edtm"]/data["time"])
    print("Accepted EDTM rate: ", data["accp_edtm"]/data["time"])
    print("Run numbers: ", data["run number"].sort_values())
    ###
    print("=======================\n\n")

def main():

    debug()

    # Plot yields
    yield_data = plot_yield()

    table = mergeDicts()
    
    file_exists = os.path.isfile(out_f)

    table.to_csv(out_f, index=False, header=True, mode='w+',)

if __name__ == '__main__':
    main()
