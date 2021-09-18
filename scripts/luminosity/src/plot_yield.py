#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2021-08-31 01:03:55 trottar"
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
    REPLAYPATH = "/group/c-kaonlt/USERS/%s/hallc_replay_lt" % USER[1]
elif ("lark" in HOST[1]):
    REPLAYPATH = "/home/%s/work/JLab/hallc_replay_lt" % USER[1]
elif ("cdaq" in HOST[1]):
    REPLAYPATH = "/home/cdaq/hallc-online/hallc_replay_lt"
elif ("trottar" in HOST[1]):
    REPLAYPATH = "/home/trottar/Analysis/hallc_replay_lt"


print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER[1], HOST[1], REPLAYPATH))

inp_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/lumi_data.csv" % str(REPLAYPATH)
out_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/yield_data.csv" % str(REPLAYPATH)

try:
    lumi_data = dict(pd.read_csv(inp_f))
except IOError:
    print("Error: %s does not appear to exist." % inp_f)
print(lumi_data.keys())
    
# prints first instance of run number-> print(lumi_data["run number"][0])
target = "carbon"

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
        "charge" : makeList("charge"),
        "current" : makeList("charge")/makeList("time"),
    
        "rate_HMS" : makeList("HMS_evts_scaler")/makeList("time"),
        "rate_SHMS" : makeList("HMS_evts_scaler")/makeList("time"),
        
        "cpuLT_scaler" : makeList("CPULT_scaler"),
        "cpuLT_scaler_uncern" : makeList("CPULT_scaler_uncern"),
        
        "TLT" : makeList("accp_edtm")/makeList("sent_edtm"),

        "HMS_evts_scaler" : makeList("HMSTRIG_scaler"),
        "uncern_HMS_evts_scaler" : np.sqrt(makeList("HMSTRIG_scaler"))/makeList("HMSTRIG_scaler"),

        "SHMS_evts_scaler" : makeList("SHMSTRIG_scaler"),
        "uncern_SHMS_evts_scaler" : np.sqrt(makeList("SHMSTRIG_scaler"))/makeList("SHMSTRIG_scaler"),

        "HMS_evts_notrack" : makeList("h_int_goodscin_evts"),
        "uncern_HMS_evts_notrack" : np.sqrt(makeList("h_int_goodscin_evts"))/makeList("h_int_goodscin_evts"),

        "SHMS_evts_notrack" : makeList("p_int_goodscin_evts"),
        "uncern_SHMS_evts_notrack" : np.sqrt(makeList("p_int_goodscin_evts"))/makeList("p_int_goodscin_evts"),

        "HMS_evts_track" : makeList("h_int_goodscin_evts"),
        "uncern_HMS_evts_track" : np.sqrt(makeList("h_int_goodscin_evts"))/makeList("h_int_goodscin_evts"),

        "SHMS_evts_track" : makeList("p_int_goodscin_evts"),
        "uncern_SHMS_evts_track" : np.sqrt(makeList("p_int_goodscin_evts"))/makeList("p_int_goodscin_evts"),

        "HMS_scaler_accp" : makeList("HMSTRIG_scaler")-makeList("sent_edtm"),
        "HMS_scaler" : makeList("HMSTRIG_scaler"),

        "SHMS_scaler_accp" : makeList("SHMSTRIG_scaler")-makeList("sent_edtm"),
        "SHMS_scaler" : makeList("SHMSTRIG_scaler"),

        "yield_HMS_scaler" : (makeList("HMSTRIG_scaler")-makeList("sent_edtm"))/makeList("charge"),
        "yield_HMS_notrack" : makeList("h_int_goodscin_evts")/(makeList("charge")*makeList("CPULT_scaler")),
        "yield_HMS_track" : makeList("h_int_goodscin_evts")/(makeList("charge")*makeList("CPULT_scaler")*makeList("HMS_track")),

        "yield_SHMS_scaler" : (makeList("SHMSTRIG_scaler")-makeList("sent_edtm"))/makeList("charge"),
        "yield_SHMS_notrack" : makeList("p_int_goodscin_evts")/(makeList("charge")*makeList("CPULT_scaler")),
        "yield_SHMS_track" : makeList("p_int_goodscin_evts")/(makeList("charge")*makeList("CPULT_scaler")*makeList("SHMS_track")),

        "HMS_track" : makeList("HMS_track"),

        "SHMS_track" : makeList("SHMS_track"),

    }

    for i,curr in enumerate(yield_dict["current"]):
        if curr == min(yield_dict["current"]):
            min_yield_HMS_scaler = yield_dict["yield_HMS_scaler"][i]
            min_yield_SHMS_scaler = yield_dict["yield_HMS_scaler"][i]
    yield_dict.update({"min_yield_HMS_scaler" : min_yield_HMS_scaler})
    yield_dict.update({"min_yield_SHMS_scaler" : min_yield_SHMS_scaler})
                
    for i,curr in enumerate(yield_dict["current"]):
        if curr == min(yield_dict["current"]):
            min_yield_HMS_notrack = yield_dict["yield_HMS_notrack"][i]
            min_yield_SHMS_notrack = yield_dict["yield_HMS_notrack"][i]
    yield_dict.update({"min_yield_HMS_notrack" : min_yield_HMS_notrack})
    yield_dict.update({"min_yield_SHMS_notrack" : min_yield_SHMS_notrack})

    for i,curr in enumerate(yield_dict["current"]):
        if curr == min(yield_dict["current"]):
            min_yield_HMS_track = yield_dict["yield_HMS_track"][i]
            min_yield_SHMS_track = yield_dict["yield_HMS_track"][i]
    yield_dict.update({"min_yield_HMS_track" : min_yield_HMS_track})
    yield_dict.update({"min_yield_SHMS_track" : min_yield_SHMS_track})

    yieldRel_HMS_scaler = yield_dict["yield_HMS_scaler"]/yield_dict["min_yield_HMS_scaler"]
    yieldRel_HMS_notrack = yield_dict["yield_HMS_notrack"]/yield_dict["min_yield_HMS_notrack"]
    yieldRel_HMS_track = yield_dict["yield_HMS_track"]/yield_dict["min_yield_HMS_track"]
    yield_dict.update({"yieldRel_HMS_scaler" : yieldRel_HMS_scaler})
    yield_dict.update({"yieldRel_HMS_notrack" : yieldRel_HMS_notrack})
    yield_dict.update({"yieldRel_HMS_track" : yieldRel_HMS_track})

    yieldRel_SHMS_scaler = yield_dict["yield_SHMS_scaler"]/yield_dict["min_yield_SHMS_scaler"]
    yieldRel_SHMS_notrack = yield_dict["yield_SHMS_notrack"]/yield_dict["min_yield_SHMS_notrack"]
    yieldRel_SHMS_track = yield_dict["yield_SHMS_track"]/yield_dict["min_yield_SHMS_track"]
    yield_dict.update({"yieldRel_SHMS_scaler" : yieldRel_SHMS_scaler})
    yield_dict.update({"yieldRel_SHMS_notrack" : yieldRel_SHMS_notrack})
    yield_dict.update({"yieldRel_SHMS_track" : yieldRel_SHMS_track})

    return yield_dict

def plot_yield():

    yield_data = calc_yield()

    for i, val in enumerate(lumi_data["run number"]):
        print("Run numbers:",lumi_data["run number"][i],"Current Values:",yield_data["current"][i])
    
    relYieldPlot = plt.figure(figsize=(12,8))

    #HMS plot scaler
    plt.subplot(2,3,1)    
    plt.grid(zorder=1)
    plt.xlim(0,100)
    plt.ylim(0.9,1.1)
    plt.plot([0,70], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["current"],yield_data["yieldRel_HMS_scaler"],yerr=yield_data["uncern_HMS_evts_scaler"],color='black',linestyle='None',zorder=3)
    plt.scatter(yield_data["current"],yield_data["yieldRel_HMS_scaler"],color='blue',zorder=4)
    plt.ylabel('Rel. Yield Scaler', fontsize=16)
    if target == 'LD2' :
        plt.title('HMS LD2 %s-%s' % (lumi_data["run number"][0],lumi_data["run number"][numRuns-1]), fontsize =16)
        plt.xlabel('Current [uA]', fontsize =16)
    elif target == 'LH2' :
        plt.title('HMS LH2 %s-%s' % (lumi_data["run number"][0],lumi_data["run number"][numRuns-1]), fontsize =16)
        plt.xlabel('Current [uA]', fontsize =16)
    else :
        plt.title('HMS Carbon %s-%s' % (lumi_data["run number"][0],lumi_data["run number"][numRuns-1]), fontsize =16)
        plt.xlabel('Current [uA]', fontsize =16)

    #HMS plot no track
    plt.subplot(2,3,2)    
    plt.grid(zorder=1)
    plt.xlim(0,100)
    # plt.ylim(0.9,1.1)
    plt.plot([0,70], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["current"],yield_data["yieldRel_HMS_notrack"],yerr=yield_data["uncern_HMS_evts_notrack"],color='black',linestyle='None',zorder=3)
    plt.scatter(yield_data["current"],yield_data["yieldRel_HMS_notrack"],color='blue',zorder=4)
    plt.ylabel('Rel. Yield no track', fontsize=16)
    if target == 'LD2' :
        plt.title('HMS LD2 %s-%s' % (lumi_data["run number"][0],lumi_data["run number"][numRuns-1]), fontsize =16)
        plt.xlabel('Current [uA]', fontsize =16)
    elif target == 'LH2' :
        plt.title('HMS LH2 %s-%s' % (lumi_data["run number"][0],lumi_data["run number"][numRuns-1]), fontsize =16)
        plt.xlabel('Current [uA]', fontsize =16)
    else :
        plt.title('HMS Carbon %s-%s' % (lumi_data["run number"][0],lumi_data["run number"][numRuns-1]), fontsize =16)
        plt.xlabel('Current [uA]', fontsize =16)

    #HMS plot track
    plt.subplot(2,3,3)    
    plt.grid(zorder=1)
    plt.xlim(0,100)
    # plt.ylim(0.9,1.1)
    plt.plot([0,70], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["current"],yield_data["yieldRel_HMS_track"],yerr=yield_data["uncern_HMS_evts_track"],color='black',linestyle='None',zorder=3)
    plt.scatter(yield_data["current"],yield_data["yieldRel_HMS_track"],color='blue',zorder=4)
    plt.ylabel('Rel. Yield track', fontsize=16)
    if target == 'LD2' :
        plt.title('HMS LD2 %s-%s' % (lumi_data["run number"][0],lumi_data["run number"][numRuns-1]), fontsize =16)
        plt.xlabel('Current [uA]', fontsize =16)
    elif target == 'LH2' :
        plt.title('HMS LH2 %s-%s' % (lumi_data["run number"][0],lumi_data["run number"][numRuns-1]), fontsize =16)
        plt.xlabel('Current [uA]', fontsize =16)
    else :
        plt.title('HMS Carbon %s-%s' % (lumi_data["run number"][0],lumi_data["run number"][numRuns-1]), fontsize =16)
        plt.xlabel('Current [uA]', fontsize =16)

        
    #SHMS plot scaler
    plt.subplot(2,3,4)    
    plt.grid(zorder=1)
    plt.xlim(0,100)
    plt.ylim(0.9,1.1)
    plt.plot([0,70], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["current"],yield_data["yieldRel_SHMS_scaler"],yerr=yield_data["uncern_SHMS_evts_scaler"],color='black',linestyle='None',zorder=3)
    plt.scatter(yield_data["current"],yield_data["yieldRel_SHMS_scaler"],color='blue',zorder=4)
    plt.ylabel('Rel. Yield Scaler', fontsize=16)
    if target == 'LD2' :
        plt.title('SHMS LD2 %s-%s' % (lumi_data["run number"][0],lumi_data["run number"][numRuns-1]), fontsize =16)
        plt.xlabel('Current [uA]', fontsize =16)
    elif target == 'LH2' :
        plt.title('SHMS LH2 %s-%s' % (lumi_data["run number"][0],lumi_data["run number"][numRuns-1]), fontsize =16)
        plt.xlabel('Current [uA]', fontsize =16)
    else :
        plt.title('SHMS Carbon %s-%s' % (lumi_data["run number"][0],lumi_data["run number"][numRuns-1]), fontsize =16)
        plt.xlabel('Current [uA]', fontsize =16)

    #SHMS plot no track
    plt.subplot(2,3,5)    
    plt.grid(zorder=1)
    plt.xlim(0,100)
    #plt.ylim(0.98,1.02)
    plt.plot([0,70], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["current"],yield_data["yieldRel_SHMS_notrack"],yerr=yield_data["uncern_SHMS_evts_notrack"],color='black',linestyle='None',zorder=3)
    plt.scatter(yield_data["current"],yield_data["yieldRel_SHMS_notrack"],color='blue',zorder=4)
    plt.ylabel('Rel. Yield no track', fontsize=16)
    if target == 'LD2' :
        plt.title('SHMS LD2 %s-%s' % (lumi_data["run number"][0],lumi_data["run number"][numRuns-1]), fontsize =16)
        plt.xlabel('Current [uA]', fontsize =16)
    elif target == 'LH2' :
        plt.title('SHMS LH2 %s-%s' % (lumi_data["run number"][0],lumi_data["run number"][numRuns-1]), fontsize =16)
        plt.xlabel('Current [uA]', fontsize =16)
    else :
        plt.title('SHMS Carbon %s-%s' % (lumi_data["run number"][0],lumi_data["run number"][numRuns-1]), fontsize =16)
        plt.xlabel('Current [uA]', fontsize =16)

    #SHMS plot track
    plt.subplot(2,3,6)    
    plt.grid(zorder=1)
    plt.xlim(0,100)
    #plt.ylim(0.98,1.02)
    plt.plot([0,70], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["current"],yield_data["yieldRel_SHMS_track"],yerr=yield_data["uncern_SHMS_evts_track"],color='black',linestyle='None',zorder=3)
    plt.scatter(yield_data["current"],yield_data["yieldRel_SHMS_track"],color='blue',zorder=4)
    plt.ylabel('Rel. Yield track', fontsize=16)
    if target == 'LD2' :
        plt.title('SHMS LD2 %s-%s' % (lumi_data["run number"][0],lumi_data["run number"][numRuns-1]), fontsize =16)
        plt.xlabel('Current [uA]', fontsize =16)
    elif target == 'LH2' :
        plt.title('SHMS LH2 %s-%s' % (lumi_data["run number"][0],lumi_data["run number"][numRuns-1]), fontsize =16)
        plt.xlabel('Current [uA]', fontsize =16)
    else :
        plt.title('SHMS Carbon %s-%s' % (lumi_data["run number"][0],lumi_data["run number"][numRuns-1]), fontsize =16)
        plt.xlabel('Current [uA]', fontsize =16)        

    plt.tight_layout()            
        
    plt.show()

    print("HMS Scaler Yield",yield_data["yieldRel_HMS_scaler"])
    print("SHMS Scaler Yield",yield_data["yieldRel_SHMS_scaler"])
    print("HMS No Track Yield",yield_data["yieldRel_HMS_notrack"])
    print("SHMS No Track Yield",yield_data["yieldRel_SHMS_notrack"])
    print("HMS Track Yield",yield_data["yieldRel_HMS_track"])
    print("SHMS Track Yield",yield_data["yieldRel_SHMS_track"])

    return yield_data

def main():
    
    yield_data = plot_yield()
    # data = {**lumi_data, **yield_data} # only python 3.5+
    
    for key, val in lumi_data.items():
        lumi_data[key] = val.tolist()

    datadict = {}
    for d in (lumi_data, yield_data): 
        datadict.update(d)
    data = {i : datadict[i] for i in sorted(datadict.keys())}
    
    # table  = pd.DataFrame([data], columns=data.keys(), index=False)
    table  = pd.DataFrame([data], columns=data.keys())
    table = table.reindex(sorted(table.columns), axis=1)
    
    file_exists = os.path.isfile(out_f)

    if file_exists:
        table.to_csv(out_f, index = True, header=False, mode='a',)
    else:
        table.to_csv(out_f, index = True, header=True, mode='a',)

if __name__ == '__main__':
    main()
