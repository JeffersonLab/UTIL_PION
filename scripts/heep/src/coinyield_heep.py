#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-03-15 01:28:15 junaid"
# ================================================================
#
# Author:  Muhammad Junaid <mjo147@uregina.ca>
#
# Copyright (c) junaid
#
###################################################################################################################################################

# Import relevant packages
import uproot
import uproot as up
import numpy as np

np.bool = bool
np.float = float

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
#cuts = ["coin_ep_cut_all_RF", "coin_ep_cut_prompt_RF", "coin_ep_cut_rand_RF"]
cuts = ["coin_ep_cut_all_noRF", "coin_ep_cut_prompt_noRF"]
lt=Root(os.path.realpath(__file__),"HeePCoin",ROOTPrefix,runNum,MaxEvent,cut_f,cuts)

OUTPATH=lt.OUTPATH

proc_root = lt.setup_ana()
c = proc_root[0] # Cut object
tree = proc_root[1] # Dictionary of branches
strDict = proc_root[2] # Dictionary of cuts as strings

#################################################################################################################################################################

def coin_protons():

    # Define the array of arrays containing the relevant HMS and SHMS info                              

    NoCut_COIN_Protons = [tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_dc_x_fp"], tree["H_dc_y_fp"], tree["H_dc_xp_fp"], tree["H_dc_yp_fp"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["CTime_epCoinTime_ROC1"], tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_dc_x_fp"], tree["P_dc_y_fp"], tree["P_dc_xp_fp"], tree["P_dc_yp_fp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["P_ngcer_npeSum"], tree["P_ngcer_xAtCer"], tree["P_ngcer_yAtCer"], tree["MMp"], tree["H_RF_Dist"],tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"], tree["emiss"], tree["Erecoil"], tree["Mrecoil"]]

    Uncut_COIN_Protons = [(tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_dc_x_fp"], tree["H_dc_y_fp"], tree["H_dc_xp_fp"], tree["H_dc_yp_fp"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["CTime_epCoinTime_ROC1"], tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_dc_x_fp"], tree["P_dc_y_fp"], tree["P_dc_xp_fp"], tree["P_dc_yp_fp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["P_ngcer_npeSum"], tree["P_ngcer_xAtCer"], tree["P_ngcer_yAtCer"], tree["MMp"], tree["H_RF_Dist"], tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"], tree["emiss"], tree["Erecoil"], tree["Mrecoil"]) for (tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_dc_x_fp"], tree["H_dc_y_fp"], tree["H_dc_xp_fp"], tree["H_dc_yp_fp"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["CTime_epCoinTime_ROC1"], tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_dc_x_fp"], tree["P_dc_y_fp"], tree["P_dc_xp_fp"], tree["P_dc_yp_fp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["P_ngcer_npeSum"], tree["P_ngcer_xAtCer"], tree["P_ngcer_yAtCer"], tree["MMp"], tree["H_RF_Dist"], tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"], tree["emiss"], tree["Erecoil"], tree["Mrecoil"]) in zip(*NoCut_COIN_Protons)
        ]

    # Create array of arrays of pions after cuts, all events, prompt and random          

    Cut_COIN_Protons_tmp = NoCut_COIN_Protons
    Cut_COIN_Protons_all_tmp = []
    Cut_COIN_Protons_prompt_tmp = []
#    Cut_COIN_Protons_rand_tmp = []

    for arr in Cut_COIN_Protons_tmp:
        Cut_COIN_Protons_all_tmp.append(c.add_cut(arr, "coin_ep_cut_all_noRF"))
        Cut_COIN_Protons_prompt_tmp.append(c.add_cut(arr, "coin_ep_cut_prompt_noRF"))
#        Cut_COIN_Protons_rand_tmp.append(c.add_cut(arr, "coin_ep_cut_rand_noRF"))

    Cut_COIN_Protons_all = [(tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_dc_x_fp"], tree["H_dc_y_fp"], tree["H_dc_xp_fp"], tree["H_dc_yp_fp"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["CTime_epCoinTime_ROC1"], tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_dc_x_fp"], tree["P_dc_y_fp"], tree["P_dc_xp_fp"], tree["P_dc_yp_fp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["P_ngcer_npeSum"], tree["P_ngcer_xAtCer"], tree["P_ngcer_yAtCer"], tree["MMp"], tree["H_RF_Dist"], tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"], tree["emiss"], tree["Erecoil"], tree["Mrecoil"]) for (tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_dc_x_fp"], tree["H_dc_y_fp"], tree["H_dc_xp_fp"], tree["H_dc_yp_fp"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["CTime_epCoinTime_ROC1"], tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_dc_x_fp"], tree["P_dc_y_fp"], tree["P_dc_xp_fp"], tree["P_dc_yp_fp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["P_ngcer_npeSum"], tree["P_ngcer_xAtCer"], tree["P_ngcer_yAtCer"], tree["MMp"], tree["H_RF_Dist"], tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"], tree["emiss"], tree["Erecoil"], tree["Mrecoil"]) in zip(*Cut_COIN_Protons_all_tmp)
        ]

    Cut_COIN_Protons_prompt = [(tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_dc_x_fp"], tree["H_dc_y_fp"], tree["H_dc_xp_fp"], tree["H_dc_yp_fp"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["CTime_epCoinTime_ROC1"], tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_dc_x_fp"], tree["P_dc_y_fp"], tree["P_dc_xp_fp"], tree["P_dc_yp_fp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["P_ngcer_npeSum"], tree["P_ngcer_xAtCer"], tree["P_ngcer_yAtCer"], tree["MMp"], tree["H_RF_Dist"], tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"], tree["emiss"], tree["Erecoil"], tree["Mrecoil"]) for (tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_dc_x_fp"], tree["H_dc_y_fp"], tree["H_dc_xp_fp"], tree["H_dc_yp_fp"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["CTime_epCoinTime_ROC1"],tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_dc_x_fp"], tree["P_dc_y_fp"], tree["P_dc_xp_fp"], tree["P_dc_yp_fp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["P_ngcer_npeSum"], tree["P_ngcer_xAtCer"], tree["P_ngcer_yAtCer"], tree["MMp"], tree["H_RF_Dist"], tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"], tree["emiss"], tree["Erecoil"], tree["Mrecoil"]) in zip(*Cut_COIN_Protons_prompt_tmp)
        ]

#    Cut_COIN_Protons_random = [(tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["CTime_epCoinTime_ROC1"], tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["P_ngcer_npeSum"], tree["P_ngcer_xAtCer"], tree["P_ngcer_yAtCer"], tree["MMp"], tree["H_RF_Dist"], tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"], tree["emiss"], tree["Erecoil"], tree["Mrecoil"]) for (tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["CTime_epCoinTime_ROC1"],tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["P_ngcer_npeSum"], tree["P_ngcer_xAtCer"], tree["P_ngcer_yAtCer"], tree["MMp"], tree["H_RF_Dist"], tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"], tree["emiss"], tree["Erecoil"], tree["Mrecoil"]) in zip(*Cut_COIN_Protons_rand_tmp)
#        ]

    COIN_Protons = {
        "Uncut_Proton_Events" : Uncut_COIN_Protons,
        "Cut_Proton_Events_All" : Cut_COIN_Protons_all,
        "Cut_Proton_Events_Prompt" : Cut_COIN_Protons_prompt,
#        "Cut_Proton_Events_Random" : Cut_COIN_Protons_random,
        }

    return COIN_Protons

##################################################################################################################################################################

#def dataframe_to_root(df, filename, key, mode='w'):
#    with uproot.recreate(filename) if mode == 'w' else uproot.update(filename) as file:
#        file[key] = {col: df[col].values for col in df.columns}

def main():
        
    COIN_Proton_Data = coin_protons()

    # This is just the list of branches we use from the initial root file for each dict
    # I don't like re-defining this here as it's very prone to errors if you included (or removed something) earlier but didn't modify it here
    # Should base the branches to include based on some list and just repeat the list here (or call it again directly below)

    COIN_Proton_Data_Header = ["H_gtr_beta","H_gtr_xp","H_gtr_yp","H_gtr_dp", "H_gtr_p", "H_dc_x_fp", "H_dc_y_fp", "H_dc_xp_fp", "H_dc_yp_fp", "H_hod_goodscinhit","H_hod_goodstarttime","H_cal_etotnorm","H_cal_etottracknorm","H_cer_npeSum","CTime_epCoinTime_ROC1","P_gtr_beta","P_gtr_xp","P_gtr_yp","P_gtr_p","P_gtr_dp", "P_dc_x_fp", "P_dc_y_fp", "P_dc_xp_fp", "P_dc_yp_fp", "P_hod_goodscinhit","P_hod_goodstarttime","P_cal_etotnorm","P_cal_etottracknorm","P_aero_npeSum","P_aero_xAtAero","P_aero_yAtAero","P_hgcer_npeSum","P_hgcer_xAtCer","P_hgcer_yAtCer","P_ngcer_npeSum","P_ngcer_xAtCer","P_ngcer_yAtCer","MMp","H_RF_Dist","P_RF_Dist", "Q2", "W", "epsilon", "ph_q", "MandelT", "pmiss", "pmiss_x", "pmiss_y", "pmiss_z","emiss", "Erecoil", "Mrecoil"]

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
#        df = pd.DataFrame(data.get(data_keys[i]), columns=DFHeader)
#        mode = 'w' if i == 0 else 'a'
#        dataframe_to_root(df, f"{OUTPATH}/{runNum}_{MaxEvent}_Analysed_Data.root", data_keys[i], mode)

        if (i == 0):
            pd.DataFrame(data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_Analysed_Data.root" % (OUTPATH, runNum, MaxEvent), key ="%s" % data_keys[i])
        elif (i != 0):
            pd.DataFrame(data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_Analysed_Data.root" % (OUTPATH, runNum, MaxEvent), key ="%s" % data_keys[i], mode ='a')

if __name__ == '__main__':
    main()
print ("Processing Complete")
