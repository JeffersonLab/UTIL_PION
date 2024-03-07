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
# Copyright (c) trottar
#


def get_file(inp_name,SCRIPTPATH):
    '''
    Grab proper lumi data file
    '''

    # Depending on input, the corresponding data setting csv data will be grabbed
    if "1" in inp_name:
        if "LH2" in inp_name.upper():
            target = "LH2"
            inp_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_1/LH2/lumi_data_l1_lh2.csv"
            out_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_1/LH2/yield_data_l1_lh2.csv"
            print("\nGrabbing input...\n\n%s" % str(inp_f))
        if "LD2" in inp_name.upper():
            target = "LD2"
            inp_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_1/LD2/lumi_data_l1_ld2.csv"
            out_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_1/LD2/yield_data_l1_ld2.csv"
            print("\nGrabbing input...\n\n%s" % str(inp_f))
        if "C" in inp_name.upper():
            target = "carbon"
            inp_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_1/Carbon0p5/lumi_data_l1_c0p5.csv"
            out_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_1/Carbon0p5/yield_data_l1_c0p5.csv"
            print("\nGrabbing input...\n\n%s" % str(inp_f))

    elif "2" in inp_name:
        if "LH2" in inp_name.upper():
            target = "LH2"
            inp_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_2/LH2/lumi_data_l2_lh2.csv"
            out_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_2/LH2/yield_data_l2_lh2.csv"
            print("\nGrabbing input...\n\n%s" % str(inp_f))
        if "LD2" in inp_name.upper():
            target = "LD2"
            inp_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_2/LD2/lumi_data_l2_ld2.csv"
            out_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_2/LD2/yield_data_l2_ld2.csv"
            print("\nGrabbing input...\n\n%s" % str(inp_f))
        if "C" in inp_name.upper():
            target = "carbon"
            inp_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_2/Carbon0p5/lumi_data_l2_c0p5.csv"
            out_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_2/Carbon0p5/yield_data_l2_c0p5.csv"
            print("\nGrabbing input...\n\n%s" % str(inp_f))

    elif "3" in inp_name:
        if "LH2" in inp_name.upper():
            target = "LH2"
            inp_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_3/LH2/lumi_data_l3_lh2.csv"
            out_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_3/LH2/yield_data_l3_lh2.csv"
            print("\nGrabbing input...\n\n%s" % str(inp_f))
        if "LD2" in inp_name.upper():
            target = "LD2"
            inp_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_3/LD2/lumi_data_l3_ld2.csv"
            out_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_3/LD2/yield_data_l3_ld2.csv"
            print("\nGrabbing input...\n\n%s" % str(inp_f))
        if "C" in inp_name.upper():
            target = "carbon"
            inp_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_3/Carbon0p5/lumi_data_l3_c0p5.csv"
            out_f = SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_3/Carbon0p5/yield_data_l3_c0p5.csv"
            print("\nGrabbing input...\n\n%s" % str(inp_f))

    else:
        target = "carbon"
        inp_f = SCRIPTPATH+"/luminosity/OUTPUTS/lumi_data.csv"
        out_f = SCRIPTPATH+"/luminosity/OUTPUTS/yield_data.csv"
        print("\nError: Invalid input...\nGrabbing default input...\n\n%s" % str(inp_f))

    return [target,inp_f,out_f]
