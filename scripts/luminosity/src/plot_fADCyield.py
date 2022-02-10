#! /usr/bin/python
#
# Description: Grabs lumi data from corresponding csv depending on run setting. Then plots the yields and creates a comprehensive table.
# Variables calculated: current, rate_HMS, rate_SHMS, sent_edtm_PS, uncern_HMS_evts_scaler, uncern_SHMS_evts_scaler, uncern_HMS_evts_notrack, uncern_SHMS_evts_notrack, uncern_HMS_evts_track, uncern_SHMS_evts_track
# ================================================================
# Time-stamp: "2021-11-04 02:01:18 trottar"
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
SCRIPTPATH = lt.SetPath(os.path.realpath(__file__)).getPath("SCRIPTPATH")

################################################################################################################################################

# Depending on input, the corresponding data setting csv data will be grabbed
inp_name = sys.argv[1]
if "fADC" in inp_name:
    target = "LH2"
    inp_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/fADC_data.csv" % str(REPLAYPATH)
    out_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/fADCyield_pt1_data.csv" % str(REPLAYPATH)
    print("\nError: Invalid input...\nGrabbing default input...\n\n%s" % str(inp_f))

print("\nRunning as %s on %s, hallc_replay_lt path assumed as %s" % (USER[1], HOST[1], REPLAYPATH))

# Converts csv data to dataframe
try:
    lumi_data = pd.read_csv(inp_f)
except IOError:
    print("Error: %s does not appear to exist." % inp_f)
print(lumi_data.keys())

def removeRun(runNum):
    '''
    Removes runs from DF and subsequently will not be plotted or included in yield csv output
    '''
    global lumi_data
    lumi_data = lumi_data[lumi_data["run number"] != runNum].reset_index(drop=True)
    return lumi_data

# Remove runs, removeRun(runNumber)

# Convert to dict for proper formatting when eventually merging dictionaries
lumi_data = dict(lumi_data)
print(lumi_data.keys())

# Define prescale variables
if "PS1" in lumi_data.keys():
    SHMS_PS = lumi_data["PS1"]
if "PS2" in lumi_data.keys():
    SHMS_PS = lumi_data["PS2"]
if "PS3" in lumi_data.keys():
    HMS_PS = lumi_data["PS3"]
if "PS4" in lumi_data.keys():
    HMS_PS = lumi_data["PS4"]
if "PS5" in lumi_data.keys():
    COIN_PS = lumi_data["PS5"]
if "PS6" in lumi_data.keys():
    COIN_PS = lumi_data["PS6"]

try:
    COIN_PS
except NameError:
    COIN_PS = None

# Define number of runs to analyze
numRuns = len(lumi_data["run number"])

def makeList(lumi_input):
    '''
    Takes input list and converts to numpy (with NaN converted to zeros) so it can be used mathematically
    '''
    new_lst = [lumi_data[lumi_input][i] for i,evts in enumerate(lumi_data["run number"])]
    new_lst = np.asarray(pd.Series(new_lst).fillna(0)) # changes NaN to zeros and convert to numpy
    return new_lst

def calc_yield():
    '''
    Creates a new dictionary with yield calculations. The relative yield is defined relative to the maximum current.
    '''
    # Create dictionary for calculations that were not calculated in previous scripts.
    yield_dict = {
        "current" : makeList("charge")/makeList("time"),
            
        "rate_HMS" : makeList("HMSTRIG_scaler")/makeList("time"),
        "rate_SHMS" : makeList("SHMSTRIG_scaler")/makeList("time"),
        "rate_COIN" : makeList("COINTRIG_scaler")/makeList("time"),
            
        "sent_edtm_PS" : makeList("sent_edtm"),
            
        "uncern_HMS_evts_scaler" : np.sqrt(makeList("HMSTRIG_scaler"))/makeList("HMSTRIG_scaler"),
            
        "uncern_SHMS_evts_scaler" : np.sqrt(makeList("SHMSTRIG_scaler"))/makeList("SHMSTRIG_scaler"),

        "uncern_COIN_evts_scaler" : np.sqrt(makeList("COINTRIG_scaler"))/makeList("COINTRIG_scaler"),
            
        "uncern_HMS_evts_notrack" : np.sqrt(makeList("h_int_etotnorm_evts"))/makeList("h_int_etotnorm_evts"),
            
        "uncern_SHMS_evts_notrack" : np.sqrt(makeList("p_int_etotnorm_evts"))/makeList("p_int_etotnorm_evts"),

        "uncern_COIN_evts_notrack" : (np.sqrt(makeList("coin_int_etotnorm_evts"))+ np.sqrt(makeList("coin_int_etotnorm_evts_rand")/6))/makeList("coin_int_etotnorm_evts"),
            
        "uncern_HMS_evts_track" : np.sqrt(makeList("h_int_etottracknorm_evts"))/makeList("h_int_etottracknorm_evts"),

        "uncern_SHMS_evts_track" : np.sqrt(makeList("p_int_etottracknorm_evts"))/makeList("p_int_etottracknorm_evts"),

        "uncern_COIN_evts_track" : (np.sqrt(makeList("coin_int_etottracknorm_evts"))+ np.sqrt(makeList("coin_int_etottracknorm_evts_rand")/6))/makeList("coin_int_etottracknorm_evts"),

        "HMS_track"  : makeList("HMS_track"),
        
        "SHMS_track" : makeList("SHMS_track"),

        "HMS_track_uncern"  : makeList("HMS_track_uncern"),
        
        "SHMS_track_uncern" : makeList("SHMS_track_uncern")
    }

    # Total livetime calculation
    TLT = makeList("accp_edtm")/yield_dict["sent_edtm_PS"]
    yield_dict.update({"TLT" : TLT})

    # Accepted scalers 
    HMS_scaler_accp = makeList("HMSTRIG_scaler")-yield_dict["sent_edtm_PS"]
    SHMS_scaler_accp = makeList("SHMSTRIG_scaler")-yield_dict["sent_edtm_PS"]
    yield_dict.update({"HMS_scaler_accp" : HMS_scaler_accp})
    yield_dict.update({"SHMS_scaler_accp" : SHMS_scaler_accp})
    COIN_scaler_accp = makeList("COINTRIG_scaler")-yield_dict["sent_edtm_PS"]
    yield_dict.update({"COIN_scaler_accp" : COIN_scaler_accp})

    # Calculate yield values

    yield_HMS_scaler = (yield_dict["HMS_scaler_accp"])/(makeList("charge")) #*makeList("CPULT_scaler")*makeList("HMS_eLT")) don't apply any efficency to scalers
    yield_HMS_notrack = ((makeList("h_int_etotnorm_evts")*makeList("PS4")))/(makeList("charge")*yield_dict["TLT"]) #+makeList("coin_int_hms_notrack")
    yield_HMS_track = ((makeList("h_int_etottracknorm_evts")*makeList("PS4")))/(makeList("charge")*yield_dict["TLT"]*makeList("HMS_track")) #+makeList("coin_int_hms")
    yield_dict.update({"yield_HMS_scaler" : yield_HMS_scaler})
    yield_dict.update({"yield_HMS_notrack" : yield_HMS_notrack})
    yield_dict.update({"yield_HMS_track" : yield_HMS_track})

    yield_SHMS_scaler = (yield_dict["SHMS_scaler_accp"])/(makeList("charge")) #*makeList("CPULT_scaler")*makeList("SHMS_eLT"))
    yield_SHMS_notrack = (makeList("p_int_etotnorm_evts")*makeList("PS1") )/(makeList("charge")*yield_dict["TLT"]) #+ makeList("coin_int_shms_notrack")
    yield_SHMS_track = (makeList("h_int_etottracknorm_evts")*makeList("PS1") )/(makeList("charge")*yield_dict["TLT"]*makeList("SHMS_track")) #+ makeList("coin_int_shms")
    yield_dict.update({"yield_SHMS_scaler" : yield_SHMS_scaler})
    yield_dict.update({"yield_SHMS_notrack" : yield_SHMS_notrack})
    yield_dict.update({"yield_SHMS_track" : yield_SHMS_track})

    yield_COIN_scaler = (yield_dict["COIN_scaler_accp"])/(makeList("charge"))
    COIN_notrack_Promt = (makeList("coin_int_etotnorm_evts"))
    COIN_notrack_Rand = (makeList("coin_int_etotnorm_evts_rand"))
    yield_COIN_notrack = (COIN_notrack_Promt - (COIN_notrack_Rand / 6))/(makeList("charge")*yield_dict["TLT"]*makeList("HMS_track")*makeList("SHMS_track"))
    COIN_track_Promt = (makeList("coin_int_etottracknorm_evts"))
    COIN_track_Rand = (makeList("coin_int_etottracknorm_evts_rand"))
    yield_COIN_track = (COIN_track_Promt - (COIN_track_Rand / 6))/(makeList("charge")*yield_dict["TLT"]*makeList("HMS_track")*makeList("SHMS_track"))
    yield_dict.update({"yield_COIN_scaler": yield_COIN_scaler})
    yield_dict.update({"yield_COIN_notrack": yield_COIN_notrack})
    yield_dict.update({"yield_COIN_track": yield_COIN_track}) 

    # Define relative yield relative to minimum current
    for i,curr in enumerate(yield_dict["current"]):
        if curr == min(yield_dict["current"]):
            max_yield_HMS_scaler = yield_dict["yield_HMS_scaler"][i]
            max_yield_SHMS_scaler = yield_dict["yield_SHMS_scaler"][i]
            max_yield_COIN_scaler = yield_dict["yield_COIN_scaler"][i]
    yield_dict.update({"max_yield_HMS_scaler" : max_yield_HMS_scaler})
    yield_dict.update({"max_yield_SHMS_scaler" : max_yield_SHMS_scaler})
    yield_dict.update({"max_yield_COIN_scaler" : max_yield_COIN_scaler})
    for i,curr in enumerate(yield_dict["current"]):
        if curr == max(yield_dict["current"]):
            max_yield_HMS_notrack = yield_dict["yield_HMS_notrack"][i]
            max_yield_SHMS_notrack = yield_dict["yield_SHMS_notrack"][i]
            max_yield_COIN_notrack = yield_dict["yield_COIN_notrack"][i]
    yield_dict.update({"max_yield_HMS_notrack" : max_yield_HMS_notrack})
    yield_dict.update({"max_yield_SHMS_notrack" : max_yield_SHMS_notrack})
    yield_dict.update({"max_yield_COIN_notrack" : max_yield_COIN_notrack})
    for i,curr in enumerate(yield_dict["current"]):
        if curr == max(yield_dict["current"]):
            max_yield_HMS_track = yield_dict["yield_HMS_track"][i]
            max_yield_SHMS_track = yield_dict["yield_SHMS_track"][i]
            max_yield_COIN_track = yield_dict["yield_COIN_track"][i]
    yield_dict.update({"max_yield_HMS_track" : max_yield_HMS_track})
    yield_dict.update({"max_yield_SHMS_track" : max_yield_SHMS_track})
    yield_dict.update({"max_yield_COIN_track" : max_yield_COIN_track})
    

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

    yieldRel_COIN_scaler = yield_dict["yield_COIN_scaler"]/yield_dict["max_yield_COIN_scaler"]
    yieldRel_COIN_notrack = yield_dict["yield_COIN_notrack"]/yield_dict["max_yield_COIN_notrack"]
    yieldRel_COIN_track = yield_dict["yield_COIN_track"]/yield_dict["max_yield_COIN_track"] 
    yield_dict.update({"yieldRel_COIN_scaler" : yieldRel_COIN_scaler})
    yield_dict.update({"yieldRel_COIN_notrack" : yieldRel_COIN_notrack})
    yield_dict.update({"yieldRel_COIN_track" : yieldRel_COIN_track})

    # Restructure dictionary to dataframe format so it matches lumi_data
    yield_table = pd.DataFrame(yield_dict, columns=yield_dict.keys())
    yield_table = yield_table.reindex(sorted(yield_table.columns), axis=1)

    return yield_table

def mergeDicts():
    '''
    Merge dictionaries/dataframes, convert to dataframe and sort
    '''
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
    '''
    Plot yields and various other analysis plots
    '''
    yield_data = mergeDicts()

    for i, val in enumerate(yield_data["run number"]):
        print("Run numbers:",yield_data["run number"][i],"Current Values:",yield_data["current"][i])
    
    relYieldPlot = plt.figure(figsize=(12,12))

    #HMS plot scaler
    plt.subplot(3,3,1)    
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
    plt.subplot(3,3,2)    
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
    plt.subplot(3,3,3)    
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
    plt.subplot(3,3,4)    
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
    plt.subplot(3,3,5)    
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
    plt.subplot(3,3,6)    
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

    #COIN plot scaler
    plt.subplot(3,3,7)    
    plt.grid(zorder=1)
    plt.xlim(0,100)
    #plt.ylim(0.9,1.1)
    plt.plot([0,100], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["current"],yield_data["yieldRel_COIN_scaler"],yerr=yield_data["uncern_COIN_evts_scaler"],color='black',linestyle='None',zorder=3)
    plt.scatter(yield_data["current"],yield_data["yieldRel_COIN_scaler"],color='blue',zorder=4)
    plt.ylabel('Rel. Yield Scaler', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('COIN LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('COIN LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('COIN Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    #COIN plot notrack
    plt.subplot(3,3,8)    
    plt.grid(zorder=1)
    plt.xlim(0,100)
    #plt.ylim(0.9,1.1)
    plt.plot([0,100], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["current"],yield_data["yieldRel_COIN_notrack"],yerr=yield_data["uncern_COIN_evts_notrack"],color='black',linestyle='None',zorder=3)
    plt.scatter(yield_data["current"],yield_data["yieldRel_COIN_notrack"],color='blue',zorder=4)
    plt.ylabel('Rel. Yield No-Track', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('COIN LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('COIN LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('COIN Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    #COIN plot track
    plt.subplot(3,3,9)    
    plt.grid(zorder=1)
    plt.xlim(0,100)
    #plt.ylim(0.9,1.1)
    plt.plot([0,100], [1,1], 'r-',zorder=2)
    plt.errorbar(yield_data["current"],yield_data["yieldRel_COIN_track"],yerr=yield_data["uncern_COIN_evts_track"],color='black',linestyle='None',zorder=3)
    plt.scatter(yield_data["current"],yield_data["yieldRel_COIN_track"],color='blue',zorder=4)
    plt.ylabel('Rel. Yield Track', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('COIN LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('COIN LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('COIN Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    plt.tight_layout()
    plt.savefig('%s/UTIL_PION/scripts/luminosity/OUTPUTS/plots/yield_%s.png' % (REPLAYPATH,"relYieldPlot_fADC"))     # Input file location and variables taking)


    #########################################################################################################################################################

    edtmPlot = plt.figure(figsize=(12,8))

    #Ratio accp/total scaler vs current
    plt.subplot(2,4,1)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.scatter(yield_data["current"],(yield_data["HMS_scaler_accp"]+yield_data["SHMS_scaler_accp"])/(yield_data["HMSTRIG_scaler"]+yield_data["SHMSTRIG_scaler"]),color='blue',zorder=4)
    plt.ylabel('Accp/Total scaler count', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    #Scaler EDTM rate vs current
    plt.subplot(2,4,2)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.scatter(yield_data["current"],yield_data["sent_edtm"]/yield_data["time"],color='blue',zorder=4)
    plt.scatter(yield_data["current"],yield_data["sent_edtm_PS"]/yield_data["time"],color='red',zorder=4)
    plt.ylabel('Scaler EDTM rate [Hz]', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    #EDTM vs HMS Rate
    plt.subplot(2,4,3)    
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
    plt.subplot(2,4,4)    
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

    #Time vs current
    plt.subplot(2,4,5)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.scatter(yield_data["current"],yield_data["time"]/60,color='blue',zorder=4)
    plt.ylabel('Time [min]', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)


    #EDTM vs SHMS Rate
    plt.subplot(2,4,6)    
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

    #HMS track vs current
    plt.subplot(2,4,7)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.errorbar(yield_data["current"],yield_data["HMS_track"],yerr=yield_data["HMS_track_uncern"],color='black',linestyle='None',zorder=3)    
    plt.scatter(yield_data["current"],yield_data["HMS_track"],color='blue',zorder=4)
    plt.ylabel('HMS Tracking Effiency', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('HMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('HMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('HMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    #SHMS track vs current
    plt.subplot(2,4,8)    
    plt.grid(zorder=1)
    #plt.xlim(0,100)
    plt.errorbar(yield_data["current"],yield_data["SHMS_track"],yerr=yield_data["SHMS_track_uncern"],color='black',linestyle='None',zorder=3)    
    plt.scatter(yield_data["current"],yield_data["SHMS_track"],color='blue',zorder=4)
    plt.ylabel('SHMS Tracking Effiency', fontsize=16)
    plt.xlabel('Current [uA]', fontsize =16)
    if target == 'LD2' :
        plt.title('SHMS LD2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    elif target == 'LH2' :
        plt.title('SHMS LH2 %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)
    else :
        plt.title('SHMS Carbon %s-%s' % (int(min(yield_data["run number"])),int(max(yield_data["run number"]))), fontsize =16)

    plt.tight_layout()      
    plt.savefig('%s/UTIL_PION/scripts/luminosity/OUTPUTS/plots/yield_%s.png' % (REPLAYPATH,"edtmPlot_fADC"))

    #########################################################################################################################################################

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
    plt.savefig('%s/UTIL_PION/scripts/luminosity/OUTPUTS/plots/yield_%s.png' % (REPLAYPATH,"logPlot_fADC"))
        
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
    print(data[["run number","PS1","PS4","sent_edtm","TLT","CPULT_scaler","current","time"]])
   # print("EDTM scaler rate: ", data["sent_edtm"]/data["time"])
   # print("Accepted EDTM rate: ", data["accp_edtm"]/data["time"])
   # print("Run numbers: ", data["run number"].sort_values())
   # print("HMS scaler ratio",data["HMS_scaler_accp"]/data["HMSTRIG_scaler"])
   # print("SHMS scaler ratio",data["SHMS_scaler_accp"]/data["SHMSTRIG_scaler"])
    print("HMS etotnorm events\n",data[["h_int_etotnorm_evts","current"]])
    print("SHMS etotnorm events\n",data[["p_int_etotnorm_evts","current"]])
    print("COIN etotnorm events\n", data[["coin_int_etotnorm_evts", "current"]])
    print("HMS yield\n",data["yield_HMS_notrack"])
    print("SHMS yield\n",data["yield_SHMS_notrack"])
    print("COIN yield\n", data["yield_COIN_notrack"])
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
