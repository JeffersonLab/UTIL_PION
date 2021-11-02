#! /usr/bin/python
#
# Description: Script is used to reanalyze all lumi data or to organize lumi data values into subdirectories
# ================================================================
# Time-stamp: "2021-11-01 01:29:30 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import numpy as np
import pandas as pd
import sys, os, subprocess

# Defines a flag for python which will reanalyze all lumi data
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-reana","--reanalyze", help="Reanalyze lumi data",action="store_true")
args = parser.parse_args() 

# Flag to chose which runs to plot (mainly for debugging, keep as "all")
l_flag = ""

if l_flag == "all":
    # All lumi runs
    lumi_list = [12126,12128,12129,12130,12131,12132,12133,12135,12137,12140,12141,12142,12143,12144,12145,12147,12148,12150,
                 12151,12152,12153,12154,12156,12157,12158,12159,12160,12161,12162,12163,12164,12166,12167,12170,12171,12172,
                 12173,12174,12175,12178,12179,12181,12184,12185,12186,12187,12189,12190,12191,12192,12193,12194,12195,12196,12197,12198,12199]
elif l_flag == "1":
    # One run (change as needed)
    lumi_list = [12126]
else:
    # Any number of runs
    #lumi_list = [12126,12128,12129,12130,12131,12132,12133,12135,12137,12140,12141,12142]
    lumi_list = [12126,12144]

# If reanalyze argument is called then all lumi data will be processed again
if args.reanalyze:
    for l in lumi_list:
        os.system("../replay_lumi.sh %s" % l) 

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

# Location of lumi data csv
inp_f = "%s/UTIL_PION/scripts/luminosity/OUTPUTS/lumi_data.csv" % str(REPLAYPATH)
try:
    lumi_data = pd.read_csv(inp_f)
except IOError:
    print("Error: %s does not appear to exist." % inp_f)
print(lumi_data.keys())

def convertDFtoCSV(inp_data,out_f):
    '''
    Converts dataframe to csv file for arbitrary df and output location
    '''
    table  = pd.DataFrame(inp_data, columns=inp_data.keys())
    table = table.reindex(sorted(table.columns), axis=1)
    
    file_exists = os.path.isfile(out_f)
    
    table.to_csv(out_f, index = False, header=True, mode='w+',)  

# Redefine lumi data by run number
l1_lh2 = dict(lumi_data.loc[(lumi_data["run number"] >= 12126) & (lumi_data["run number"] <= 12142)])
# Convert to csv from dataframe
convertDFtoCSV(l1_lh2,"%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_1/LH2/lumi_data_l1_lh2.csv" % str(REPLAYPATH))
print("\n\nLumi #1 LH2 runs {0} are now in {1}".format(list(l1_lh2["run number"]),"%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_1/LH2/lumi_data_l1_lh2.csv" % str(REPLAYPATH)))
l1_c = dict(lumi_data.loc[(lumi_data["run number"] >= 12144) & (lumi_data["run number"] <= 12154)])
convertDFtoCSV(l1_c,"%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_1/Carbon0p5/lumi_data_l1_c0p5.csv" % str(REPLAYPATH))
print("Lumi #1 C runs {0} are now in {1}".format(list(l1_c["run number"]),"%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_1/Carbon0p5/lumi_data_l1_c0p5.csv" % str(REPLAYPATH)))
l1_ld2 = dict(lumi_data.loc[(lumi_data["run number"] >= 12157) & (lumi_data["run number"] <= 12166)])
convertDFtoCSV(l1_ld2,"%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_1/LD2/lumi_data_l1_ld2.csv" % str(REPLAYPATH))
print("Lumi #1 LD2 runs {0} are now in {1}".format(list(l1_ld2["run number"]),"%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_1/LD2/lumi_data_l1_ld2.csv" % str(REPLAYPATH)))
l2_lh2 = dict(lumi_data.loc[(lumi_data["run number"] >= 12192) & (lumi_data["run number"] <= 12199)])
convertDFtoCSV(l2_lh2,"%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_2/LH2/lumi_data_l2_lh2.csv" % str(REPLAYPATH))
print("Lumi #2 LH2 runs {0} are now in {1}".format(list(l2_lh2["run number"]),"%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_2/LH2/lumi_data_l2_lh2.csv" % str(REPLAYPATH)))
l2_c = dict(lumi_data.loc[(lumi_data["run number"] >= 12181) & (lumi_data["run number"] <= 12191)])
convertDFtoCSV(l2_c,"%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_2/Carbon0p5/lumi_data_l2_c0p5.csv" % str(REPLAYPATH))
print("Lumi #2 C runs {0} are now in {1}".format(list(l2_c["run number"]),"%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_2/Carbon0p5/lumi_data_l2_c0p5.csv" % str(REPLAYPATH)))
l2_ld2 = dict(lumi_data.loc[(lumi_data["run number"] >= 12167) & (lumi_data["run number"] <= 12175)])
convertDFtoCSV(l2_ld2,"%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_2/LD2/lumi_data_l2_ld2.csv" % str(REPLAYPATH))
print("Lumi #2 LD2 runs {0} are now in {1}".format(list(l2_ld2["run number"]),"%s/UTIL_PION/scripts/luminosity/OUTPUTS/Lumi_2/LD2/lumi_data_l2_ld2.csv" % str(REPLAYPATH)))
