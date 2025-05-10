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
if len(sys.argv)-1!=4:
    print("!!!!! ERROR !!!!!\n Expected 3 arguments\n Usage is with - ROOTfilePrefix RunNumber MaxEvents \n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################

# Input params - run number and max number of events
ROOTPrefix = sys.argv[1]
runNum = sys.argv[2]
MaxEvent = sys.argv[3]
Spec = sys.argv[4]
Spec = Spec.upper()

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
if Spec == "HMS":
   cut_f = '/DB/CUTS/run_type/hSing_heep.cuts'
   # defining Cuts
   cuts = ["sing_ee_cut_all_noRF"]
   lt=Root(os.path.realpath(__file__),"HeePSing_HMS",ROOTPrefix,runNum,MaxEvent,cut_f,cuts)
if Spec == "SHMS":
   cut_f = '/DB/CUTS/run_type/pSing_heep.cuts'
   # defining Cuts
   cuts = ["sing_ee_cut_all_noRF"]
   lt=Root(os.path.realpath(__file__),"HeePSing_SHMS",ROOTPrefix,runNum,MaxEvent,cut_f,cuts)

OUTPATH=lt.OUTPATH

proc_root = lt.setup_ana()
c = proc_root[0] # Cut object
tree = proc_root[1] # Dictionary of branches
strDict = proc_root[2] # Dictionary of cuts as strings

#################################################################################################################################################################

def sing_electrons():
    if Spec == "HMS":
       # Define the array of arrays containing the relevant HMS and SHMS info                              

       NoCut_SING_Electrons = [tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_dc_x_fp"], tree["H_dc_y_fp"], tree["H_dc_xp_fp"], tree["H_dc_yp_fp"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["H_RF_Dist"], tree["Q2"], tree["H_W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"], tree["emiss"], tree["Erecoil"], tree["Mrecoil"]]

       Uncut_SING_Electrons = [(tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_dc_x_fp"], tree["H_dc_y_fp"], tree["H_dc_xp_fp"], tree["H_dc_yp_fp"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["H_RF_Dist"], tree["Q2"], tree["H_W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"], tree["emiss"], tree["Erecoil"], tree["Mrecoil"]) for (tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_dc_x_fp"], tree["H_dc_y_fp"], tree["H_dc_xp_fp"], tree["H_dc_yp_fp"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["H_RF_Dist"], tree["Q2"], tree["H_W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"], tree["emiss"], tree["Erecoil"], tree["Mrecoil"]) in zip(*NoCut_SING_Electrons)
            ]

       # Create array of arrays of pions after cuts, all events, prompt and random          
       Cut_SING_Electrons_tmp = NoCut_SING_Electrons
       Cut_SING_Electrons_all_tmp = []
       for arr in Cut_SING_Electrons_tmp:
           Cut_SING_Electrons_all_tmp.append(c.add_cut(arr, "sing_ee_cut_all_noRF"))

       Cut_SING_Electrons_all = [(tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_dc_x_fp"], tree["H_dc_y_fp"], tree["H_dc_xp_fp"], tree["H_dc_yp_fp"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["H_RF_Dist"], tree["Q2"], tree["H_W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"], tree["emiss"], tree["Erecoil"], tree["Mrecoil"]) for (tree["H_gtr_beta"], tree["H_gtr_xp"], tree["H_gtr_yp"], tree["H_gtr_dp"], tree["H_gtr_p"], tree["H_dc_x_fp"], tree["H_dc_y_fp"], tree["H_dc_xp_fp"], tree["H_dc_yp_fp"], tree["H_hod_goodscinhit"], tree["H_hod_goodstarttime"], tree["H_cal_etotnorm"], tree["H_cal_etottracknorm"], tree["H_cer_npeSum"], tree["H_RF_Dist"], tree["Q2"], tree["H_W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"], tree["emiss"], tree["Erecoil"], tree["Mrecoil"]) in zip(*Cut_SING_Electrons_all_tmp)
        ]

       SING_Electrons = {
           "Uncut_Electron_Events" : Uncut_SING_Electrons,
           "Cut_Electron_Events_All" : Cut_SING_Electrons_all,
           }

    if Spec == "SHMS":
       # Define the array of arrays containing the relevant HMS and SHMS info                              

       NoCut_SING_Electrons = [tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_dc_x_fp"], tree["P_dc_y_fp"], tree["P_dc_xp_fp"], tree["P_dc_yp_fp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["P_ngcer_npeSum"], tree["P_ngcer_xAtCer"], tree["P_ngcer_yAtCer"],tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"], tree["emiss"], tree["Erecoil"], tree["Mrecoil"]]

       Uncut_SING_Electrons = [(tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_dc_x_fp"], tree["P_dc_y_fp"], tree["P_dc_xp_fp"], tree["P_dc_yp_fp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["P_ngcer_npeSum"], tree["P_ngcer_xAtCer"], tree["P_ngcer_yAtCer"],tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"], tree["emiss"], tree["Erecoil"], tree["Mrecoil"]) for (tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_dc_x_fp"], tree["P_dc_y_fp"], tree["P_dc_xp_fp"], tree["P_dc_yp_fp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["P_ngcer_npeSum"], tree["P_ngcer_xAtCer"], tree["P_ngcer_yAtCer"],tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"], tree["emiss"], tree["Erecoil"], tree["Mrecoil"]) in zip(*NoCut_SING_Electrons)
        ]

       # Create array of arrays of pions after cuts, all events, prompt and random          
       Cut_SING_Electrons_tmp = NoCut_SING_Electrons
       Cut_SING_Electrons_all_tmp = []
       for arr in Cut_SING_Electrons_tmp:
           Cut_SING_Electrons_all_tmp.append(c.add_cut(arr, "sing_ee_cut_all_noRF"))

       Cut_SING_Electrons_all = [(tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_dc_x_fp"], tree["P_dc_y_fp"], tree["P_dc_xp_fp"], tree["P_dc_yp_fp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["P_ngcer_npeSum"], tree["P_ngcer_xAtCer"], tree["P_ngcer_yAtCer"],tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"], tree["emiss"], tree["Erecoil"], tree["Mrecoil"]) for (tree["P_gtr_beta"], tree["P_gtr_xp"], tree["P_gtr_yp"], tree["P_gtr_p"], tree["P_gtr_dp"], tree["P_dc_x_fp"], tree["P_dc_y_fp"], tree["P_dc_xp_fp"], tree["P_dc_yp_fp"], tree["P_hod_goodscinhit"], tree["P_hod_goodstarttime"], tree["P_cal_etotnorm"], tree["P_cal_etottracknorm"], tree["P_aero_npeSum"], tree["P_aero_xAtAero"], tree["P_aero_yAtAero"], tree["P_hgcer_npeSum"], tree["P_hgcer_xAtCer"], tree["P_hgcer_yAtCer"], tree["P_ngcer_npeSum"], tree["P_ngcer_xAtCer"], tree["P_ngcer_yAtCer"],tree["P_RF_Dist"], tree["Q2"], tree["W"], tree["epsilon"], tree["ph_q"], tree["MandelT"], tree["pmiss"], tree["pmiss_x"], tree["pmiss_y"], tree["pmiss_z"], tree["emiss"], tree["Erecoil"], tree["Mrecoil"]) in zip(*Cut_SING_Electrons_all_tmp)
        ]

       SING_Electrons = {
           "Uncut_Electron_Events" : Uncut_SING_Electrons,
           "Cut_Electron_Events_All" : Cut_SING_Electrons_all,
           }

    return SING_Electrons

##################################################################################################################################################################

def main():
        
    SING_Electron_Data = sing_electrons()

    if Spec == "HMS":
       SING_Electron_Data_Header = ["H_gtr_beta","H_gtr_xp","H_gtr_yp","H_gtr_dp", "H_gtr_p", "H_dc_x_fp", "H_dc_y_fp", "H_dc_xp_fp", "H_dc_yp_fp", "H_hod_goodscinhit","H_hod_goodstarttime","H_cal_etotnorm","H_cal_etottracknorm","H_cer_npeSum","H_RF_Dist","Q2", "H_W", "epsilon", "ph_q", "MandelT", "pmiss", "pmiss_x", "pmiss_y", "pmiss_z","emiss", "Erecoil", "Mrecoil"]
    if Spec == "SHMS":
       SING_Electron_Data_Header = ["P_gtr_beta","P_gtr_xp","P_gtr_yp","P_gtr_p","P_gtr_dp", "P_dc_x_fp", "P_dc_y_fp", "P_dc_xp_fp", "P_dc_yp_fp", "P_hod_goodscinhit","P_hod_goodstarttime","P_cal_etotnorm","P_cal_etottracknorm","P_aero_npeSum","P_aero_xAtAero","P_aero_yAtAero","P_hgcer_npeSum","P_hgcer_xAtCer","P_hgcer_yAtCer","P_ngcer_npeSum","P_ngcer_xAtCer","P_ngcer_yAtCer","P_RF_Dist", "Q2", "W", "epsilon", "ph_q", "MandelT", "pmiss", "pmiss_x", "pmiss_y", "pmiss_z","emiss", "Erecoil", "Mrecoil"]

    # Need to create a dict for all the branches we grab                                                
    data = {}
    data.update(SING_Electron_Data)
    data_keys = list(data.keys()) # Create a list of all the keys in all dicts added above, each is an array of data                                                                                       

    for i in range (0, len(data_keys)):
        if("Electron" in data_keys[i]):
            DFHeader=list(SING_Electron_Data_Header)
        else:
            continue
            # Uncomment the line below if you want .csv file output, WARNING the files can be very large and take a long time to process!                                                                      
            #pd.DataFrame(data.get(data_keys[i])).to_csv("%s/%s_%s.csv" % (OUTPATH, data_keys[i], runNum), header=DFHeader, index=False) # Convert array to panda dataframe and write to csv with correct header
        
        if (i == 0):
            pd.DataFrame(data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_%s_Analysed_Data.root" % (OUTPATH, runNum, MaxEvent, Spec), key ="%s" % data_keys[i])

        elif (i != 0):
            pd.DataFrame(data.get(data_keys[i]), columns = DFHeader, index = None).to_root("%s/%s_%s_%s_Analysed_Data.root" % (OUTPATH, runNum, MaxEvent, Spec), key ="%s" % data_keys[i], mode ='a')
if __name__ == '__main__':
    main()
print ("Processing Complete")
