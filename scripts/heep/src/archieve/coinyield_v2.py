#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2022-06-30 01:28:15 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

# 15/01/21 - Stephen Kay, University of Regina
# 21/06/21 - Edited By - Muhammad Junaid, University of Regina, Canada

# Python version of the pion analysis script. Now utilises uproot to select event of each type and writes them to a root file
# Intention is to apply PID/selection cutting here and plot in a separate script
# Python should allow for easier reading of databases storing timing offsets e.t.c.
# 27/04/21 - Updated to use new hcana variables, old determinations removed

###################################################################################################################################################

# Import relevant packages
import uproot as up
import numpy as np
import root_numpy as rnp
import pandas as pd
import root_pandas as rpd
import ROOT
import scipy
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import sys, math, os, subprocess

##################################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=3:
    print("!!!!! ERROR !!!!!\n Expected 3 arguments\n Usage is with - ROOTfilePrefix RunNumber MaxEvents \n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################

# Input params - run number and max number of events
ROOTPrefix = sys.argv[1]
runNum = sys.argv[2]
MaxEvent = sys.argv[3]

################################################################################################################################################
'''
ltsep package import
'''

# Import package for cuts
from ltsep import Root

##############################################################################################################################################
'''
Define and set up cuts
'''

cut_f = '/DB/CUTS/run_type/coin_heep.cuts'

# defining Cuts
cuts = ["coin_ep_cut_all_RF", "coin_ep_cut_prompt_RF", "coin_ep_cut_rand_RF"]

lt=Root(os.path.realpath(__file__),"HeePCoin",ROOTPrefix,runNum,MaxEvent,cut_f,cuts)

OUTPATH=lt.OUTPATH

proc_root = lt.setup_ana()
c = proc_root[0] # Cut object
tree = proc_root[1] # Dictionary of branches
strDict = proc_root[2] # Dictionary of cuts as strings

#################################################################################################################################################################

def coin_protons():

    # Define the array of arrays containing the relevant HMS and SHMS info                              

    NoCut_COIN_Protons = [tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["CTime_epCoinTime_ROC1"], tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["MMp"], tree["H_RF_Dist"],tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"]]

    Uncut_COIN_Protons = [(tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["CTime_epCoinTime_ROC1"], tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["MMp"], tree["H_RF_Dist"], tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"]) for (tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["CTime_epCoinTime_ROC1"], tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["MMp"], tree["H_RF_Dist"], tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"]) in zip(*NoCut_COIN_Protons)
        ]

    # Create array of arrays of pions after cuts, all events, prompt and random          

    Cut_COIN_Protons_tmp = NoCut_COIN_Protons
    Cut_COIN_Protons_all_tmp = []
    Cut_COIN_Protons_prompt_tmp = []
    Cut_COIN_Protons_rand_tmp = []

    for arr in Cut_COIN_Protons_tmp:
        Cut_COIN_Protons_all_tmp.append(c.add_cut(arr, "coin_ep_cut_all_RF"))
        Cut_COIN_Protons_prompt_tmp.append(c.add_cut(arr, "coin_ep_cut_prompt_RF"))
        Cut_COIN_Protons_rand_tmp.append(c.add_cut(arr, "coin_ep_cut_rand_RF"))

    Cut_COIN_Protons_all = [(tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["CTime_epCoinTime_ROC1"], tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["MMp"], tree["H_RF_Dist"], tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"]) for (tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["CTime_epCoinTime_ROC1"], tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["MMp"], tree["H_RF_Dist"], tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"]) in zip(*Cut_COIN_Protons_all_tmp)
        ]

    Cut_COIN_Protons_prompt = [(tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["CTime_epCoinTime_ROC1"], tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["MMp"], tree["H_RF_Dist"], tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"]) for (tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["CTime_epCoinTime_ROC1"],tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["MMp"], tree["H_RF_Dist"], tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"]) in zip(*Cut_COIN_Protons_prompt_tmp)
        ]

    Cut_COIN_Protons_random = [(tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["CTime_epCoinTime_ROC1"], tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["MMp"], tree["H_RF_Dist"], tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"]) for (tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["CTime_epCoinTime_ROC1"],tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["MMp"], tree["H_RF_Dist"], tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"]) in zip(*Cut_COIN_Protons_rand_tmp)
        ]

    COIN_Protons = {
        "Uncut_Proton_Events" : Uncut_COIN_Protons,
        "Cut_Proton_Events_All" : Cut_COIN_Protons_all,
        "Cut_Proton_Events_Prompt" : Cut_COIN_Protons_prompt,
        "Cut_Proton_Events_Random" : Cut_COIN_Protons_random,
        }

    return COIN_Protons

##################################################################################################################################################################

def main():
        
    COIN_Proton_Data = coin_protons()

    # This is just the list of branches we use from the initial root file for each dict
    # I don't like re-defining this here as it's very prone to errors if you included (or removed something) earlier but didn't modify it here
    # Should base the branches to include based on some list and just repeat the list here (or call it again directly below)

    COIN_Proton_Data_Header = ["H_gtr_beta","H_gtr_xp","H_gtr_yp","H_gtr_dp", "H_gtr_p", "H_hod_goodscinhit","H_hod_goodstarttime","H_cal_etotnorm","H_cal_etottracknorm","H_cer_npeSum","CTime_epCoinTime_ROC1","P_gtr_beta","P_gtr_xp","P_gtr_yp","P_gtr_p","P_gtr_dp","P_hod_goodscinhit","P_hod_goodstarttime","P_cal_etotnorm","P_cal_etottracknorm","P_aero_npeSum","P_aero_xAtAero","P_aero_yAtAero","P_hgcer_npeSum","P_hgcer_xAtCer","P_hgcer_yAtCer","MMp","H_RF_Dist","P_RF_Dist", "Q2", "W", "epsilon", "ph_q", "MandelT", "pmiss", "pmiss_x", "pmiss_y", "pmiss_z"]

    # Need to create a dict for all the branches we grab                                                
    data = {}
    data.update(COIN_Proton_Data)
    data_keys = list(data.keys()) # Create a list of all the keys in all dicts added above, each is an array of data                                                                                       

    for i in range (0, len(data_keys)):
        if("Proton" in data_keys[i]):
            DFHeader=list(COIN_Proton_Data_Header)
        else:
            continue
            # Uncomment the line below if you want .csv file output, WARNING the files can be very large and take a long time to process!                                                                      
            #pd.DataFrame(data.get(data_keys[i])).to_csv("%s/%s_%s.csv" % (OUTPATH, data_keys[i], runNum), header=DFHeader, index=False) # Convert array to panda dataframe and write to csv with correct header                                                                                                      
        if (i == 0):
            pd.DataFrame(data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_Analysed_Data.root" % (OUTPATH, runNum, MaxEvent), key ="%s" % data_keys[i])
        elif (i != 0):
            pd.DataFrame(data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_Analysed_Data.root" % (OUTPATH, runNum, MaxEvent), key ="%s" % data_keys[i], mode ='a')

if __name__ == '__main__':
    main()
print ("Processing Complete")
