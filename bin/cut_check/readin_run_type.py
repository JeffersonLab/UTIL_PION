#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-03-24 16:52:50 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import numpy as np
import pandas as pd
# Required for up arrow returns
import readline
import os, sys

################################################################################################################################################
'''
User Inputs
'''

#ROOTPrefix = sys.argv[1]

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

run_type_dir = UTILPATH+"/DB/CUTS/run_type/"
general_dir = UTILPATH+"/DB/CUTS/general/"
param_dir = UTILPATH+"/DB/PARAM/"

runTypeDict = {

    "fADCdeadtime" : run_type_dir+"fADCdeadtime.cuts",
    "coin_heep" : run_type_dir+"coin_heep.cuts",
    "hSing_optics" : run_type_dir+"hSing_optics.cuts",
    "hSing_prod" : run_type_dir+"hSing_prod.cuts",
    "pSing_optics" : run_type_dir+"pSing_optics.cuts",
    "pSing_prod" : run_type_dir+"pSing_prod.cuts",
    "coinpeak" : run_type_dir+"coinpeak.cuts",
    "coin_prod" : run_type_dir+"coin_prod.cuts",
    "lumi" : run_type_dir+"lumi.cuts",
    "pid_eff" : run_type_dir+"pid_eff.cuts",
    "simc_coin_heep" : run_type_dir+"simc_coin_heep.cuts",
    "simc_sing_heep" : run_type_dir+"simc_sing_heep.cuts",

}

# Matches run type cuts with the general cuts (e.g pid, track, etc.)
generalDict = {
    "pid" : general_dir+"pid.cuts",
    "track" : general_dir+"track.cuts",
    "accept" : general_dir+"accept.cuts",
    "coin_time" : general_dir+"coin_time.cuts",
    "current" : general_dir+"current.cuts",
    "misc" : general_dir+"misc.cuts",
}

# Matches general cuts with run number specific cuts
paramDict = {
    "accept" : param_dir+"Acceptance_Parameters.csv",
    "track" : param_dir+"Tracking_Parameters.csv",
    "CT" : param_dir+"Timing_Parameters.csv",
    "pid" : param_dir+"PID_Parameters.csv",
    "misc" : param_dir+"Misc_Parameters.csv",
    "current" : param_dir+"Current_Parameters.csv"
}

for key, val in runTypeDict.items():
    print("{} -> {}\n".format(key,val),"-"*50)

def has_numbers(inputString):
    '''
    Check if string contains a number. Returns true if number is in string.
    '''
    return any(char.isdigit() for char in inputString)
    
def run_type_cut(cut):

    file_content = []
    with open(runTypeDict[cut], "r") as f:
        for line in f:
            if "#" not in line:
                file_content.append("\n"+line+"\n"+"-"*50)
        
    print(" ".join(file_content))

    return file_content

def general_cut(cuts, user_inp):

    cuts = cuts.split("=")
    cut_name = cuts[0].strip()
    cut_lst = cuts[1].split("+")

    out_cuts = ""
    
    file_content = []
    if user_inp in cut_name:
        for cut in cut_lst:
            for key, val in generalDict.items():
                if key in cut:
                    cut_key = cut.strip().split(".")[0]
                    cut_val = cut.strip().split(".")[1]
                    with open(generalDict[cut_key], "r") as f:
                        for line in f:
                            if "#" not in line:
                                if cut_val in line:
                                    if "==" in line:
                                        file_content.append(line.replace("==","&&").split("=")[1].replace("&&","=="))
                                    elif ">=" in line:
                                        file_content.append(line.replace(">=","&&").split("=")[1].replace("&&",">="))
                                    elif "<=" in line:
                                        file_content.append(line.replace("<=","&&").split("=")[1].replace("&&","<="))                                        
                                    else:
                                        file_content.append(line.split("=")[1])

        out_cuts = "\033[36m"+cut_name+"\033[0m = \033[32m"+",".join(file_content).replace("\n","")+"\033[0m"
    
    return out_cuts

def param_cut(cut, user_inp, runNum):

    cuts = general_cut(cut, user_cut_inp)
    print("?????",cuts)

    cuts = cuts.split(" = ")

    print("!!!!",cuts)

    out_cuts = ""
    
    if len(cuts) == 1:
        return out_cuts
    
    cut_name = cuts[0].strip()
    cut_lst = cuts[1].split(",")
    
    file_content = []
    for i, cut in enumerate(cut_lst):
        for key, val in paramDict.items():            
            if key in cut:
                param_tmp = cut.split(")")
                for param in param_tmp:
                    if "." in param:
                        print("~~~~",param)
                        paramVal = param.split(".")[1]
                        if has_numbers(paramVal):
                            cut = cut
                        else:
                            # Search param dictionary for values based off key
                            fout = paramDict[key]
                            try:
                                data = dict(pd.read_csv(fout))
                                data[paramVal]
                            except IOError:
                                print("ERROR 9: %s not found in %s" % (paramVal,fout))
                            except KeyError:
                                continue
                            for j,evt in enumerate(data['Run_Start']):
                                # Check if run number is defined in param file
                                if data['Run_Start'][j] <= np.int64(runNum) <= data['Run_End'][j]:
                                    cut = cut.replace(key+"."+paramVal,str(data[paramVal][j]))
                                    pass
                                else:
                                    # print("!!!!ERROR!!!!: Run %s not found in range %s-%s" % (np.int64(runNum),data['Run_Start'][i],data['Run_End'][i])) # Error 10
                                    continue
                                    
            else:
                cut = cut                
        file_content.append(cut)
        
    out_cuts = "\033[36m"+str(runNum)+"\n\n"+cut_name+"\033[0m = \033[32m"+",".join(file_content).replace("\n","")+"\033[0m"
    
    return out_cuts

user_run_type_inp =  input("\n\nPlease enter a run type cut (type 'exit' to end)...")

if user_run_type_inp[0:3] == "bye" or user_run_type_inp[0:4] == "exit":
    sys.exit(0)

cut_lst = run_type_cut(user_run_type_inp)

# Keep track of the last user input
last_user_cut_input = ""
last_user_check_input = ""
last_user_run_input = ""

while True:    

    user_cut_inp =  input("\n\nPlease enter a specific cut (type 'exit' to end and 'help' to relist run type cuts)...")

    # Returns previous prompt if up arrow is pressed
    if user_cut_inp == "\033[A":
        user_cut_inp = last_user_cut_input
        print(user_cut_inp)
    else:
        last_user_cut_input = user_cut_inp
    
    if user_cut_inp[0:3] == "bye" or user_cut_inp[0:4] == "exit":
        break

    if user_cut_inp[0:4] == "help":
        cut_lst = run_type_cut(user_run_type_inp)
        continue
    
    while True:

        user_check_inp =  input('\nWould you like to check cuts for specific run number? (yes or no)...')

        # Returns previous prompt if up arrow is pressed
        if user_check_inp == "\033[A":
            user_check_inp = last_user_check_input
            print(user_check_inp)
        else:
            last_user_check_input = user_check_inp        

        if "y" in user_check_inp:
            user_run_inp =  input('\nPlease enter run number (type exit to end)...')

            # Returns previous prompt if up arrow is pressed
            if user_run_inp == "\033[A":
                user_run_inp = last_user_run_input
                print(user_run_inp)
            else:
                last_user_run_input = user_run_inp
                
            if user_run_inp[0:3] == "bye" or user_run_inp[0:4] == "exit":
                break
            try:
                runNum = int(user_run_inp)
            except:
                print("Need a proper run number...")
                continue
            
            for cut in cut_lst:
                output = param_cut(cut.replace("-"*50,"").replace("\n",""), user_cut_inp, runNum)
                if output != "":
                    print("\n\n",output,"\n\n")
                    break
        elif "n" in user_check_inp:
            for cut in cut_lst:
                output = general_cut(cut.replace("-"*50,"").replace("\n",""), user_cut_inp)
                if output != "":
                    print("\n\n",output,"\n\n")
                    break
            break
        else:
            print("Please answer yes or no...")
            continue
