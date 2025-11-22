#! /usr/bin/python
#
# Description:
# ================================================================
# Time-stamp: "2025-03-13 01:29:19 junaid"
# ================================================================
#
# Author:  Muhammad Junaid III <mjo147@uregina.ca>
#
# Copyright (c) junaid
#
###################################################################################################################################################

import re

def dictionary(UTILPATH,ROOTPrefix,runNum,MaxEvent, DEBUG=False):

    # Open report file to grab prescale values and tracking efficiency
#    report = UTILPATH+"/REPORT_OUTPUT/Analysis/HeeP/%s_%s_%s.report" % (ROOTPrefix,runNum,MaxEvent)
#    report = UTILPATH+"/REPORT_OUTPUT/Analysis/Lumi/%s_%s_%s.report" % (ROOTPrefix,runNum,MaxEvent)
#    report = UTILPATH+"/REPORT_OUTPUT/Analysis/pTRIG6/%s_%s_%s.report" % (ROOTPrefix,runNum,MaxEvent)
    report = UTILPATH+"/REPORT_OUTPUT/Analysis/PionLT/%s_%s_%s.report" % (ROOTPrefix,runNum,MaxEvent)

    with open(report) as f:
        effDict = {
            # General Info
            'Run_Number': None,
            'Target_Mass_(amu)': None ,
            'Beam_Energy' : None,
            'HMS_Particle_Mass' : None,
            'HMS_P_Central' : None,
            'HMS_Angle' : None,
            'SHMS_Particle_Mass' : None,
            'SHMS_P_Central' : None,
            'SHMS_Angle' : None,
            'BCM2_Charge': None,
#            'BCM2_Charge_ERROR': None,
            'BCM2_Beam_Cut_Charge': None,
            'BCM2_Current': None,
            'BCM2_Beam_Cut_Current': None,
            'SHMS_Run_Length': None,
            'HMS_Run_Length': None,
            'BCM_Cut_SHMS_Run_Length': None,
            'BCM_Cut_HMS_Run_Length': None,
            'Corrected_Current': None,
            'Corrected_Current_ERROR': None,
            'Corrected_Charge': None,
            'Corrected_Charge_ERROR': None,
            'Target_BoilingCorr': None,
            'Target_BoilingCorr_ERROR': None,
            'Target_BoilingCorr_SysUncer': None,
            'Ps1_factor': None,
            'Ps2_factor': None,
            'Ps3_factor': None,
            'Ps4_factor': None,
            'Ps5_factor': None,
            'Ps6_factor': None,
            'Total_SHMS_3/4_Triggers': None,
            '(current_cut)_Total_SHMS_3/4_Triggers': None ,
            'Pre-Scaled_SHMS_3/4_Triggers': None ,
            'SHMS_3/4_Trigger_Rate': None ,
            'Total_SHMS_EL-REAL_Triggers': None,
            '(current_cut)_Total_SHMS_EL-REAL_Triggers': None ,
            'Pre-Scaled_SHMS_EL-REAL_Triggers': None ,
            'SHMS_EL-REAL_Trigger_Rate': None ,
            'Accepted_SHMS_Triggers': None ,
            'Total_HMS_3/4_Triggers': None,
            '(current_cut)_Total_HMS_3/4_Triggers': None ,
            'Pre-Scaled_HMS_3/4_Triggers': None ,
            'HMS_3/4_Trigger_Rate': None ,
            'Total_HMS_EL-REAL_Triggers': None,
            '(current_cut)_Total_HMS_EL-REAL_Triggers': None ,
            'Pre-Scaled_HMS_EL-REAL_Triggers': None ,
            'HMS_EL-REAL_Trigger_Rate': None ,
            'Accepted_HMS_Triggers': None ,
            'Total_COIN_Triggers': None,
            '(current_cut)_Total_COIN_Triggers': None ,
            'Pre-Scaled_COIN_Triggers': None ,
            'Accepted_COIN_Triggers': None ,
            'COIN_Trigger_Rate': None ,
            # EDTM
            'EDTM_Accepted_Triggers': None ,
            'HMS_EDTM_Accepted_Triggers': None ,
            'SHMS_EDTM_Accepted_Triggers': None ,
            'Non_Scaler_EDTM_Live_Time': None,
            'Non_Scaler_EDTM_Live_Time_ERROR': None,
            'Non_Scaler_EDTM_Live_Time_Corr': None,
            'Non_Scaler_EDTM_Live_Time_Corr_ERROR': None,
            'Non_Scaler_EDTM_Live_Time_Syst_Uncer': None,
            # CPULT
            'HMS_CPULT': None,
            'HMS_CPULT_ERROR': None,
            'SHMS_CPULT': None,
            'SHMS_CPULT_ERROR': None,
            'COIN_CPULT': None,
            'COIN_CPULT_ERROR': None,
            # Coin Blocking Corrections
            'CoinBlocking': None,
            'CoinBlocking_ERROR': None,
            'CoinBlocking_SysUp': None,
            'CoinBlocking_SysDown': None,
            'CoinBlocking_Sys': None,
            # Tracking Efficiencies
            "SHMS_Elec_ALL_TRACK_EFF" : None,
            "SHMS_Elec_ALL_TRACK_EFF_ERROR" : None,
            "SHMS_Elec_COIN_TRACK_EFF" : None,
            "SHMS_Elec_COIN_TRACK_EFF_ERROR" : None,
            "SHMS_Elec_SING_TRACK_EFF" : None,
            "SHMS_Elec_SING_TRACK_EFF_ERROR" : None,
            "SHMS_Hadron_ALL_TRACK_EFF" : None,
            "SHMS_Hadron_ALL_TRACK_EFF_ERROR" : None,
            "SHMS_Prot_ALL_TRACK_EFF" : None,
            "SHMS_Prot_ALL_TRACK_EFF_ERROR" : None,
            "SHMS_Prot_COIN_TRACK_EFF" : None,
            "SHMS_Prot_COIN_TRACK_EFF_ERROR" : None,
            "SHMS_Prot_SING_TRACK_EFF" : None,
            "SHMS_Prot_SING_TRACK_EFF_ERROR" : None,
            "SHMS_Pion_ALL_TRACK_EFF" : None,
            "SHMS_Pion_ALL_TRACK_EFF_ERROR" : None,
            "SHMS_Pion_COIN_TRACK_EFF" : None,
            "SHMS_Pion_COIN_TRACK_EFF_ERROR" : None,
            "SHMS_Pion_SING_TRACK_EFF" : None,
            "SHMS_Pion_SING_TRACK_EFF_ERROR" : None,
            "HMS_Elec_ALL_TRACK_EFF" : None,
            "HMS_Elec_ALL_TRACK_EFF_ERROR" : None,
            "HMS_Elec_COIN_TRACK_EFF" : None,
            "HMS_Elec_COIN_TRACK_EFF_ERROR" : None,
            "HMS_Elec_SING_TRACK_EFF" : None,
            "HMS_Elec_SING_TRACK_EFF_ERROR" : None,
            # HGC (useless without cuts)
            'SHMS_HGC_ALL_Elec_Eff': None,
            'SHMS_HGC_ALL_Elec_Eff_ERROR': None,
            'SHMS_HGC_COIN_Elec_Eff': None,
            'SHMS_HGC_COIN_Elec_Eff_ERROR': None,
            'SHMS_HGC_SING_Elec_Eff': None,
            'SHMS_HGC_SING_Elec_Eff_ERROR': None,
            'SHMS_HGC_ALL_Pion_Eff': None,
            'SHMS_HGC_ALL_Pion_Eff_ERROR': None,
            'SHMS_HGC_COIN_Pion_Eff': None,
            'SHMS_HGC_COIN_Pion_Eff_ERROR': None,
            'SHMS_HGC_SING_Pion_Eff': None,
            'SHMS_HGC_SING_Pion_Eff_ERROR': None,
            'SHMS_HGC_ALL_Prot_Eff': None,
            'SHMS_HGC_ALL_Prot_Eff_ERROR': None,
            'SHMS_HGC_COIN_Prot_Eff': None,
            'SHMS_HGC_COIN_Prot_Eff_ERROR': None,
            'SHMS_HGC_SING_Prot_Eff': None,
            'SHMS_HGC_SING_Prot_Eff_ERROR': None,
            # Aerogel
            'SHMS_Aero_ALL_Elec_Eff': None,
            'SHMS_Aero_ALL_Elec_Eff_ERROR': None,
            'SHMS_Aero_COIN_Elec_Eff': None,
            'SHMS_Aero_COIN_Elec_Eff_ERROR': None,
            'SHMS_Aero_SING_Elec_Eff': None,
            'SHMS_Aero_SING_Elec_Eff_ERROR': None,
            'SHMS_Aero_ALL_Pion_Eff': None,
            'SHMS_Aero_ALL_Pion_Eff_ERROR': None,
            'SHMS_Aero_COIN_Pion_Eff': None,
            'SHMS_Aero_COIN_Pion_Eff_ERROR': None,
            'SHMS_Aero_SING_Pion_Eff': None,
            'SHMS_Aero_SING_Pion_Eff_ERROR': None,
            'SHMS_Aero_ALL_Prot_Eff': None,
            'SHMS_Aero_ALL_Prot_Eff_ERROR': None,
            'SHMS_Aero_COIN_Prot_Eff': None,
            'SHMS_Aero_COIN_Prot_Eff_ERROR': None,
            'SHMS_Aero_SING_Prot_Eff': None,
            'SHMS_Aero_SING_Prot_Eff_ERROR': None,
            # HMS Cer
            'HMS_Cer_ALL_Elec_Eff': None,
            'HMS_Cer_ALL_Elec_Eff_ERROR': None,
            'HMS_Cer_COIN_Elec_Eff': None,
            'HMS_Cer_COIN_Elec_Eff_ERROR': None,
            'HMS_Cer_SING_Elec_Eff': None,
            'HMS_Cer_SING_Elec_Eff_ERROR': None,
            # Hodoscope (no uncertainties...yet!!!!)
            'SHMS_Hodo_Plane_1' : None,
            'SHMS_Hodo_Plane_2' : None,
            'SHMS_Hodo_Plane_3' : None,
            'SHMS_Hodo_Plane_4' : None,
            'SHMS_Hodo_S1XY' : None,
            'SHMS_Hodo_S2XY' : None,
            'SHMS_Hodo_3_of_4_EFF' : None,
#            'SHMS_Hodo_3_of_4_EFF_ERROR' : None,
            'SHMS_Hodo_4_of_4_EFF' : None,
#            'SHMS_Hodo_4_of_4_EFF_ERROR' : None,
            'HMS_Hodo_Plane_1' : None,
            'HMS_Hodo_Plane_2' : None,
            'HMS_Hodo_Plane_3' : None,
            'HMS_Hodo_Plane_4' : None,
            'HMS_Hodo_S1XY' : None,
            'HMS_Hodo_S2XY' : None,
            'HMS_Hodo_3_of_4_EFF' : None,
#            'HMS_Hodo_3_of_4_EFF_ERROR' : None,
            'HMS_Hodo_4_of_4_EFF' : None,
#            'HMS_Hodo_4_of_4_EFF_ERROR' : None,
            'SHMS_Hodoscope_S1X_Triggers' : None,
            'HMS_Hodoscope_S1X_Triggers' : None,
            'SHMS_Hodoscope_S1X_Rate' : None,
            'HMS_Hodoscope_S1X_Rate' : None,
            # Calorimeter
	         'HMS_Cal_ALL_Elec_Eff' : None,            
            'HMS_Cal_ALL_Elec_Eff_ERROR' : None,
	        'HMS_Cal_COIN_Elec_Eff' : None,
            'HMS_Cal_COIN_Elec_Eff_ERROR' : None,
	        'HMS_Cal_SING_Elec_Eff' : None,
            'HMS_Cal_SING_Elec_Eff_ERROR' : None,
            'SHMS_Cal_ALL_Elec_Eff' : None,
            'SHMS_Cal_ALL_Elec_Eff_ERROR' : None,
            'SHMS_Cal_COIN_Elec_Eff' : None,
            'SHMS_Cal_COIN_Elec_Eff_ERROR' : None,
            'SHMS_Cal_SING_Elec_Eff' : None,
            'SHMS_Cal_SING_Elec_Eff_ERROR' : None,
        }

        # Search for keywords, then save as value in dictionary
        for line in f:
            data = line.split(':')
            if len(data) < 2:
                continue  # Skip invalid lines

            # Normalize key and value: strip PLT_, and unify ± to +-
            key_field = data[0].strip()
            if key_field.startswith("PLT_"):
                key_field = key_field[4:]
            value_field = data[1].replace("±", "+-")

            for key, val in effDict.items():
                # Extract error only from the line whose key exactly matches the base (nkey)
                if "ERROR" in key:
                    nkey = key.replace("_ERROR", "")
                    if key_field == nkey:
                        if DEBUG:
                            print("ERROR:", nkey, "\t", key)
                            print(data)
                        try:
                            effDict[key] = float(re.compile(r'[^\d.]+').sub("", value_field.split("+-")[1]))
                        except IndexError:
                            pass  # leave as None if no "+-"

                # Extract value only when the key matches exactly; prevents SysUncer overwrites
                if key_field == key:
                    if DEBUG:
                        print(key)
                        print(data)
                    if "+-" in value_field:
                        effDict[key] = float(re.compile(r'[^\d.]+').sub("", value_field.split("+-")[0]))
                    else:
                        effDict[key] = float(re.compile(r'[^\d.]+').sub("", value_field))

        if DEBUG:
            for key, val in effDict.items():
                if val is None:
                    print(key, "\t", val)
        # Remove any empty or None keys that are not used by this run type
        effDict = {key: val for key, val in effDict.items() if val is not None}

        if DEBUG:
            print("Final extracted efficiencies:\n", effDict) 

    return effDict
