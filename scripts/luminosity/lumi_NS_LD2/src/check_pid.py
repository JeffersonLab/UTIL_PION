#! /usr/bin/python
#
# Description: Grabs lumi data from corresponding csv depending on run setting. Then plots the yields and creates a comprehensive table.
# Variables calculated: current, rate_HMS, rate_SHMS, sent_edtm_PS, uncern_HMS_evts_scaler, uncern_SHMS_evts_scaler, uncern_HMS_evts_notrack, uncern_SHMS_evts_notrack, uncern_HMS_evts_track, uncern_SHMS_evts_track
# ================================================================
# Time-stamp: "2022-06-29 03:25:33 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import pandas as pd

def check_pid(runNum,ANATYPE):

    inp_f = "%sLT/pid_type.csv" % ANATYPE

    # Converts csv data to dataframe
    try:
        pid_data = pd.read_csv(inp_f)
        #print(inp_f)
        #print(pid_data.keys())
    except IOError:
        print("Error: %s does not appear to exist." % inp_f)
        sys.exit(0)

    pid_data_runNum = pid_data[(pid_data['Run_Start'] <= int(runNum)) & (pid_data['Run_End'] >= int(runNum))]

    particleDict = {
        "pion" : "pi",
        "electron" : "e",
        "kaon" : "k",
        "proton" : "p",
        "hadron" : "had",
    } 

    for key,val in particleDict.items():
        if pid_data_runNum['HMS_PID'].values[0] == key:
            HMS_PID = val
        if pid_data_runNum['SHMS_PID'].values[0] == key:
            SHMS_PID = val

    print('''
    For run {0}...
    HMS_PID = {1}
    SHMS_PID = {2}
    '''.format(runNum,HMS_PID,SHMS_PID))

    return [HMS_PID,SHMS_PID]
