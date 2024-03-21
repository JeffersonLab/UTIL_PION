#! /usr/bin/python
#
# Description: Grabs lumi data from corresponding csv depending on run setting. Then plots the yields and creates a comprehensive table.
# Variables calculated: current, rate_HMS, rate_SHMS, sent_edtm_PS, uncern_HMS_evts_scaler, uncern_SHMS_evts_scaler, uncern_HMS_evts_notrack, uncern_SHMS_evts_notrack, uncern_HMS_evts_track, uncern_SHMS_evts_track
# ================================================================
# Time-stamp: "2022-06-27 02:10:13 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Update - Nathan Heinrich (heinricn)
# Reworked Paths so that Data for PionLT experiment would be read in instead
#

def get_file(inp_name,SCRIPTPATH):
    '''
    Grab proper lumi data file
    '''

    # Depending on input, the corresponding data setting csv data will be grabbed
    if "9-2" in inp_name:
        if "LH2" in inp_name.upper():
            target = "LH2"
            inp_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi9-2/lumi_data_LH2.csv"
            out_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi9-2/yield_data_LH2.csv"
            print("\nGrabbing input...\n\n%s" % str(inp_f))
        if "LD2" in inp_name.upper():
            target = "LD2"
            inp_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi9-2/lumi_data_LD2.csv"
            out_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi9-2/yield_data_LD2.csv"
            print("\nGrabbing input...\n\n%s" % str(inp_f))
        if "C" in inp_name.upper():
            target = "carbon"
            inp_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi9-2/lumi_data_Carbon.csv"
            out_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi9-2/yield_data_Carbon.csv"
            print("\nGrabbing input...\n\n%s" % str(inp_f))

    else:
        target = "carbon"
        inp_f = SCRIPTPATH+"/luminosity/OUTPUTS/lumi_data.csv"
        out_f = SCRIPTPATH+"/luminosity/OUTPUTS/yield_data.csv"
        print("\nError: Invalid input...\nGrabbing default input...\n\n%s" % str(inp_f))

    return [target,inp_f,out_f]
