#! /usr/bin/python
#
# Description: Script is used to reanalyze all lumi data or to organize lumi data values into subdirectories
# ================================================================
# Time-stamp: "2023-05-31 11:32:18 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import numpy as np
import pandas as pd
import sys, os, subprocess

################################################################################################################################################
'''
Define a flag to reanalyze lumi data
'''

# Defines a flag for python which will reanalyze all lumi data
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-reana","--reanalyze", help="Reanalyze lumi data",action="store_true")
args = parser.parse_args() 

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
Determine which runs to reanalyze
'''

if ANATYPE == "Pion":

    # Flag to chose which runs to plot (mainly for debugging, keep as "all")
    l_flag = "1"

    if l_flag == "all":
        # All lumi runs
        lumi_list = [12126,12128,12129,12130,12131,12132,12133,12135,12137,12140,12141,12142,12143,12144,12145,12147,12148,12150,
                     12151,12152,12153,12154,12156,12157,12158,12159,12160,12161,12162,12163,12164,12166,12167,12170,12171,12172,
                     12173,12174,12175,12178,12179,12181,12184,12185,12186,12187,12189,12190,12191,12192,12193,12194,12195,12196,
                     12197,12198,12199,13874,13876,13877,13878,13879,13880,13881,13882,13885,13886,13887,13891,13892,13894,13897,
                     13898,13899,13900,13901,13902,13903,13904,13905,13906,13907,13908,13910,13911,13912,13913,13914,13915,13916,13918,13919]
    elif l_flag == "1":
        # One run (change as needed)
        lumi_list = [13874]
    else:
        # Any number of runs
        #lumi_list = [12126,12128,12129,12130,12131,12132,12133,12135,12137,12140,12141,12142]
        lumi_list = [12126,12144]
else:
    # Flag to chose which runs to plot (mainly for debugging, keep as "all")
    l_flag = "1"

    if l_flag == "all":
        # All lumi runs
        lumi_list = [5149,5150,5151,5152,5153,5154,5155,5156,5157,5158,5159,5160,5161,5162,5163,5164,5165,5166,5167,5168,5169,5170,5171,5173,
                     5175,5176,5177,5178,5179,5180,5181,5295,5297,5298,5299,5300,5301,5302,5303,7841,7842,7843,7844,7845,7846,7847,7862,7859,
                     7863,7864,7862,7859,7863,7864,7865,7948,7949,7950,7951,7952,7954,7956,7957,7958,7959,7960,7961,8299,8300,8301,8302,8303]
    elif l_flag == "1":
        # One run (change as needed)
        lumi_list = [5154]
    else:
        # Any number of runs
        lumi_list = [5154,5155,5156,5157,5158,5298,5299]

################################################################################################################################################
'''
Reanalyze data
'''

# If reanalyze argument is called then all lumi data will be processed again
if args.reanalyze:
    for l in lumi_list:
        os.system("../replay_lumi.sh %s" % l) 

################################################################################################################################################
'''
Move lumi analyzed data from general csv to setting and target specific
'''

# Location of lumi data csv
inp_f = SCRIPTPATH+"/luminosity/OUTPUTS/lumi_data.csv"
try:
    lumi_data = pd.read_csv(inp_f)
    print(lumi_data.keys())
except IOError:
    print("Error: %s does not appear to exist." % inp_f)
    sys.exit(0)

def convertDFtoCSV(inp_data,out_f):
    '''
    Converts dataframe to csv file for arbitrary df and output location
    '''
    file_exists = os.path.isfile(out_f)

    table  = pd.DataFrame(inp_data, columns=inp_data.keys())
    table = table.reindex(sorted(table.columns), axis=1)

    if file_exists:
        table.to_csv(out_f, index = False, header=True, mode='w',)
    else:
        table.to_csv(out_f, index = False, header=True, mode='a',)


if ANATYPE == "Pion":
    # Redefine lumi data by run number
    l1_lh2 = dict(lumi_data.loc[(lumi_data["run number"] >= 12126) & (lumi_data["run number"] <= 12142)])
    # Convert to csv from dataframe
    convertDFtoCSV(l1_lh2,SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_1/LH2/lumi_data_l1_lh2.csv")
    print("\n\nLumi #1 LH2 runs {0} are now in {1}".format(list(l1_lh2["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_1/LH2/lumi_data_l1_lh2.csv"))
    l1_c = dict(lumi_data.loc[(lumi_data["run number"] >= 12144) & (lumi_data["run number"] <= 12154)])
    convertDFtoCSV(l1_c,SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_1/Carbon0p5/lumi_data_l1_c0p5.csv")
    print("Lumi #1 C runs {0} are now in {1}".format(list(l1_c["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_1/Carbon0p5/lumi_data_l1_c0p5.csv"))
    l1_ld2 = dict(lumi_data.loc[(lumi_data["run number"] >= 12157) & (lumi_data["run number"] <= 12166)])
    convertDFtoCSV(l1_ld2,SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_1/LD2/lumi_data_l1_ld2.csv")
    print("Lumi #1 LD2 runs {0} are now in {1}".format(list(l1_ld2["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_1/LD2/lumi_data_l1_ld2.csv"))
    l2_lh2 = dict(lumi_data.loc[(lumi_data["run number"] >= 12192) & (lumi_data["run number"] <= 12199)])
    convertDFtoCSV(l2_lh2,SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_2/LH2/lumi_data_l2_lh2.csv")
    print("Lumi #2 LH2 runs {0} are now in {1}".format(list(l2_lh2["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_2/LH2/lumi_data_l2_lh2.csv"))
    l2_c = dict(lumi_data.loc[(lumi_data["run number"] >= 12181) & (lumi_data["run number"] <= 12191)])
    convertDFtoCSV(l2_c,SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_2/Carbon0p5/lumi_data_l2_c0p5.csv")
    print("Lumi #2 C runs {0} are now in {1}".format(list(l2_c["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_2/Carbon0p5/lumi_data_l2_c0p5.csv"))
    l2_ld2 = dict(lumi_data.loc[(lumi_data["run number"] >= 12167) & (lumi_data["run number"] <= 12175)])
    convertDFtoCSV(l2_ld2,SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_2/LD2/lumi_data_l2_ld2.csv")
    print("Lumi #2 LD2 runs {0} are now in {1}".format(list(l2_ld2["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_2/LD2/lumi_data_l2_ld2.csv"))
    l3_lh2 = dict(lumi_data.loc[(lumi_data["run number"] >= 13874) & (lumi_data["run number"] <= 13881)])
    convertDFtoCSV(l3_lh2,SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_3/LH2/lumi_data_l3_lh2.csv")
    print("Lumi #3 LH2 runs {0} are now in {1}".format(list(l3_lh2["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_3/LH2/lumi_data_l3_lh2.csv"))
    l3_c = dict(lumi_data.loc[(lumi_data["run number"] >= 13897) & (lumi_data["run number"] <= 13910)])
    convertDFtoCSV(l3_c,SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_3/Carbon0p5/lumi_data_l3_c0p5.csv")
    print("Lumi #3 C runs {0} are now in {1}".format(list(l3_c["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_3/Carbon0p5/lumi_data_l3_c0p5.csv"))
    l3_ld2 = dict(lumi_data.loc[(lumi_data["run number"] >= 13885) & (lumi_data["run number"] <= 13894)])
    convertDFtoCSV(l3_ld2,SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_3/LD2/lumi_data_l3_ld2.csv")
    print("Lumi #3 LD2 runs {0} are now in {1}".format(list(l3_ld2["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_3/LD2/lumi_data_l3_ld2.csv"))
    l3_edtm = dict(lumi_data.loc[(lumi_data["run number"] >= 13911) & (lumi_data["run number"] <= 13919)])
    convertDFtoCSV(l3_edtm,SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_3/EDTM/lumi_data_l3_edtm.csv")
    print("Lumi #3 EDTM runs {0} are now in {1}".format(list(l3_edtm["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/Lumi_3/EDTM/lumi_data_l2_edtm.csv"))

else:
    # Redefine lumi data by run number
    lh2_l1_10p6 = dict(lumi_data.loc[((lumi_data["run number"] >= 5150) & (lumi_data["run number"] <= 5153)) 
                                     | ((lumi_data["run number"] >= 5160) & (lumi_data["run number"] <= 5166)) 
                                     | (lumi_data["run number"] == 5295)
                                     | (lumi_data["run number"] == 5297)])
    # Convert to csv from dataframe
    convertDFtoCSV(lh2_l1_10p6,SCRIPTPATH+"/luminosity/OUTPUTS/10p6/Lumi_1/LH2/lumi_data_lh2_l1_10p6.csv")
    print("\n\nLumi #1 LH2 runs {0} are now in {1}".format(list(lh2_l1_10p6["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/10p6/Lumi_1/LH2/lumi_data_lh2_l1_10p6.csv"))
    c_l1_10p6 = dict(lumi_data.loc[((lumi_data["run number"] >= 5154) & (lumi_data["run number"] <= 5158)) 
                                   | ((lumi_data["run number"] >= 5298) & (lumi_data["run number"] <= 5299))])
    convertDFtoCSV(c_l1_10p6,SCRIPTPATH+"/luminosity/OUTPUTS/10p6/Lumi_1/Carbon0p5/lumi_data_c_l1_10p6.csv")
    print("\n\nLumi #1 Carbon0p5 runs {0} are now in {1}".format(list(c_l1_10p6["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/10p6/Lumi_1/Carbon0p5/lumi_data_c_l1_10p6.csv"))
    lh2_l2_10p6 = dict(lumi_data.loc[((lumi_data["run number"] >= 5167) & (lumi_data["run number"] <= 5171)) 
                                     | ((lumi_data["run number"] >= 5302) & (lumi_data["run number"] <= 5303))])
    convertDFtoCSV(lh2_l2_10p6,SCRIPTPATH+"/luminosity/OUTPUTS/10p6/Lumi_2/LH2/lumi_data_lh2_l2_10p6.csv")
    print("\n\nLumi #2 LH2 runs {0} are now in {1}".format(list(lh2_l2_10p6["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/10p6/Lumi_2/LH2/lumi_data_lh2_l2_10p6.csv"))
    c_l2_10p6 = dict(lumi_data.loc[((lumi_data["run number"] >= 5175) & (lumi_data["run number"] <= 5176)) 
                                   | ((lumi_data["run number"] >= 5178) & (lumi_data["run number"] <= 5179)) 
                                   | (lumi_data["run number"] == 5181) 
                                   | ((lumi_data["run number"] >= 5300) & (lumi_data["run number"] <= 5301))])
    convertDFtoCSV(c_l2_10p6,SCRIPTPATH+"/luminosity/OUTPUTS/10p6/Lumi_2/Carbon0p5/lumi_data_c_l2_10p6.csv")
    print("\n\nLumi #2 Carbon0p5 runs {0} are now in {1}".format(list(c_l2_10p6["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/10p6/Lumi_2/Carbon0p5/lumi_data_c_l2_10p6.csv"))
    lh2_l3_10p6 = dict(lumi_data.loc[((lumi_data["run number"] >= 5342) & (lumi_data["run number"] <= 5350)) 
                                     | (lumi_data["run number"] == 5362)
                                     | (lumi_data["run number"] == 5366)])
    convertDFtoCSV(lh2_l3_10p6,SCRIPTPATH+"/luminosity/OUTPUTS/10p6/Lumi_3/LH2/lumi_data_lh2_l3_10p6.csv")
    print("\n\nLumi #3 LH2 runs {0} are now in {1}".format(list(lh2_l3_10p6["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/10p6/Lumi_3/LH2/lumi_data_lh2_l3_10p6.csv"))
    al_l3_10p6 = dict(lumi_data.loc[((lumi_data["run number"] >= 5359) & (lumi_data["run number"] <= 5361))])
    convertDFtoCSV(al_l3_10p6,SCRIPTPATH+"/luminosity/OUTPUTS/10p6/Lumi_3/Aluminum/lumi_data_al_l3_10p6.csv")
    print("\n\nLumi #3 Aluminum runs {0} are now in {1}".format(list(al_l3_10p6["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/10p6/Lumi_3/Aluminum/lumi_data_al_l3_10p6.csv"))
    c_l3_10p6 = dict(lumi_data.loc[((lumi_data["run number"] >= 5351) & (lumi_data["run number"] <= 5358))])
    convertDFtoCSV(c_l3_10p6,SCRIPTPATH+"/luminosity/OUTPUTS/10p6/Lumi_3/Carbon0p5/lumi_data_c_l3_10p6.csv")
    print("\n\nLumi #3 Carbon0p5 runs {0} are now in {1}".format(list(c_l3_10p6["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/10p6/Lumi_3/Carbon0p5/lumi_data_c_l3_10p6.csv"))
    lh2_l1_6p2 = dict(lumi_data.loc[((lumi_data["run number"] >= 7843) & (lumi_data["run number"] <= 7845)) 
                                    | ((lumi_data["run number"] >= 7862) & (lumi_data["run number"] <= 7863))])
    convertDFtoCSV(lh2_l1_6p2,SCRIPTPATH+"/luminosity/OUTPUTS/6p2/Lumi_1/LH2/lumi_data_lh2_l1_6p2.csv")
    print("\n\nLumi #1 LH2 runs {0} are now in {1}".format(list(lh2_l1_6p2["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/6p2/Lumi_1/LH2/lumi_data_lh2_l1_6p2.csv"))
    c_l1_6p2 = dict(lumi_data.loc[(lumi_data["run number"] == 7841) 
                                  | ((lumi_data["run number"] >= 7846) & (lumi_data["run number"] <= 7847)) 
                                  | ((lumi_data["run number"] >= 7864) & (lumi_data["run number"] <= 7865))])
    convertDFtoCSV(c_l1_6p2,SCRIPTPATH+"/luminosity/OUTPUTS/6p2/Lumi_1/Carbon0p5/lumi_data_c_l1_6p2.csv")
    print("\n\nLumi #1 Carbon0p5 runs {0} are now in {1}".format(list(c_l1_6p2["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/6p2/Lumi_1/Carbon0p5/lumi_data_c_l1_6p2.csv"))
    lh2_l1_8p2 = dict(lumi_data.loc[(lumi_data["run number"] == 7954)
                                    | ((lumi_data["run number"] >= 7956) & (lumi_data["run number"] <= 7960))])
    convertDFtoCSV(lh2_l1_8p2,SCRIPTPATH+"/luminosity/OUTPUTS/8p2/Lumi_1/LH2/lumi_data_lh2_l1_8p2.csv")
    print("\n\nLumi #1 LH2 runs {0} are now in {1}".format(list(lh2_l1_8p2["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/8p2/Lumi_1/LH2/lumi_data_lh2_l1_8p2.csv"))
    c_l1_8p2 = dict(lumi_data.loc[((lumi_data["run number"] >= 7948) & (lumi_data["run number"] <= 7952))])
    convertDFtoCSV(c_l1_8p2,SCRIPTPATH+"/luminosity/OUTPUTS/8p2/Lumi_1/Carbon0p5/lumi_data_c_l1_8p2.csv")
    print("\n\nLumi #1 Carbon0p5 runs {0} are now in {1}".format(list(c_l1_8p2["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/8p2/Lumi_1/Carbon0p5/lumi_data_c_l1_8p2.csv"))
    lh2_l2_8p2 = dict(lumi_data.loc[((lumi_data["run number"] >= 8038) & (lumi_data["run number"] <= 8085))])
    convertDFtoCSV(lh2_l2_8p2,SCRIPTPATH+"/luminosity/OUTPUTS/8p2/Lumi_1/LH2/lumi_data_lh2_l2_8p2.csv")
    print("\n\nLumi #1 LH2 runs {0} are now in {1}".format(list(lh2_l2_8p2["run number"]),SCRIPTPATH+"/luminosity/OUTPUTS/8p2/Lumi_1/LH2/lumi_data_lh2_l2_8p2.csv"))
