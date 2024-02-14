#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2022-06-28 06:36:42 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from csv import DictReader
import sys, os

################################################################################################################################################
'''
User Inputs
'''
ROOTPrefix = sys.argv[1]
runType = sys.argv[2]
timestmp = sys.argv[3]

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
import ltsep as lt 

p=lt.SetPath(os.path.realpath(__file__))

# Add this to all files for more dynamic pathing
USER=p.getPath("USER") # Grab user info for file finding
HOST=p.getPath("HOST")
REPLAYPATH=p.getPath("REPLAYPATH")
UTILPATH=p.getPath("UTILPATH")
ANATYPE=p.getPath("ANATYPE")
SCRIPTPATH=p.getPath("SCRIPTPATH")

################################################################################################################################################

inp_f = UTILPATH+"/scripts/efficiency/OUTPUTS/%s_%s_efficiency_data_%s.csv"  % (ROOTPrefix.replace("replay_",""),runType,timestmp)

# Converts csv data to dataframe
try:
    efficiency_data = pd.read_csv(inp_f)
except IOError:
    print("Error: %s does not appear to exist." % inp_f)
print(efficiency_data.keys())

plt.figure(figsize=(12,8))

plt.subplot(141)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["Non_Scaler_EDTM_Live_Time"],yerr=efficiency_data["Non_Scaler_EDTM_Live_Time_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["Non_Scaler_EDTM_Live_Time"],color='blue',zorder=4)
plt.ylabel('EDTM', fontsize=12)
plt.xlabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(142)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Pion_SING_TRACK_EFF"],yerr=efficiency_data["SHMS_Pion_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Pion_SING_TRACK_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Pion_SING_TRACK_EFF', fontsize=12)
plt.xlabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(143)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Aero_SING_Pion_Eff"],yerr=efficiency_data["SHMS_Aero_SING_Pion_Eff_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Aero_SING_Pion_Eff"],color='blue',zorder=4)
plt.ylabel('SHMS_Aero_SING_Pion_Eff', fontsize=12)
plt.xlabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(144)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])   
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_3-4_%s.png' % (ROOTPrefix.replace("replay_","")))

plt.figure(figsize=(12,8))

plt.subplot(141)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["Non_Scaler_EDTM_Live_Time"],yerr=efficiency_data["Non_Scaler_EDTM_Live_Time_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["Non_Scaler_EDTM_Live_Time"],color='blue',zorder=4)
plt.ylabel('EDTM', fontsize=12)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(142)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Pion_SING_TRACK_EFF"],yerr=efficiency_data["SHMS_Pion_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Pion_SING_TRACK_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Pion_SING_TRACK_EFF', fontsize=12)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(143)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Aero_SING_Pion_Eff"],yerr=efficiency_data["SHMS_Aero_SING_Pion_Eff_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Aero_SING_Pion_Eff"],color='blue',zorder=4)
plt.ylabel('SHMS_Aero_SING_Pion_Eff', fontsize=12)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(144)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_S1X_%s.png' % (ROOTPrefix.replace("replay_","")))

plt.figure(figsize=(12,8))

plt.subplot(141)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["HMS_EL-REAL_Trigger_Rate"],efficiency_data["Non_Scaler_EDTM_Live_Time"],yerr=efficiency_data["Non_Scaler_EDTM_Live_Time_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["HMS_EL-REAL_Trigger_Rate"],efficiency_data["Non_Scaler_EDTM_Live_Time"],color='red',zorder=4)
plt.ylabel('EDTM', fontsize=12)
plt.xlabel('HMS EL-REAL Trigger Rate [kHz]', fontsize=12)
plt.title('HMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(142)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["HMS_EL-REAL_Trigger_Rate"],efficiency_data["HMS_Elec_SING_TRACK_EFF"],yerr=efficiency_data["HMS_Elec_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["HMS_EL-REAL_Trigger_Rate"],efficiency_data["HMS_Elec_SING_TRACK_EFF"],color='red',zorder=4)
plt.ylabel('HMS_Elec_SING_TRACK_EFF', fontsize=12)
plt.xlabel('HMS EL-REAL Trigger Rate [kHz]', fontsize=12)
plt.title('HMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(143)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["HMS_EL-REAL_Trigger_Rate"],efficiency_data["HMS_Cer_SING_Elec_Eff"],yerr=efficiency_data["HMS_Cer_SING_Elec_Eff_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["HMS_EL-REAL_Trigger_Rate"],efficiency_data["HMS_Cer_SING_Elec_Eff"],color='red',zorder=4)
plt.ylabel('HMS_Cer_SING_Elec_Eff', fontsize=12)
plt.xlabel('HMS EL-REAL Trigger Rate [kHz]', fontsize=12)
plt.title('HMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(144)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data["HMS_EL-REAL_Trigger_Rate"],efficiency_data["HMS_Hodo_3_of_4_EFF"],color='red',zorder=4)
plt.ylabel('HMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('HMS EL-REAL Trigger Rate [kHz]', fontsize=12)
plt.title('HMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])   
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/HMS_EL-REAL_%s.png' % (ROOTPrefix.replace("replay_","")))

plt.figure(figsize=(12,8))

plt.subplot(141)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["HMS_Hodoscope_S1X_Rate"],efficiency_data["Non_Scaler_EDTM_Live_Time"],yerr=efficiency_data["Non_Scaler_EDTM_Live_Time_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["HMS_Hodoscope_S1X_Rate"],efficiency_data["Non_Scaler_EDTM_Live_Time"],color='red',zorder=4)
plt.ylabel('EDTM', fontsize=12)
plt.xlabel('HMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('HMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(142)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["HMS_Hodoscope_S1X_Rate"],efficiency_data["HMS_Elec_SING_TRACK_EFF"],yerr=efficiency_data["HMS_Elec_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["HMS_Hodoscope_S1X_Rate"],efficiency_data["HMS_Elec_SING_TRACK_EFF"],color='red',zorder=4)
plt.ylabel('HMS_Elec_SING_TRACK_EFF', fontsize=12)
plt.xlabel('HMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('HMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(143)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["HMS_Hodoscope_S1X_Rate"],efficiency_data["HMS_Cer_SING_Elec_Eff"],yerr=efficiency_data["HMS_Cer_SING_Elec_Eff_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["HMS_Hodoscope_S1X_Rate"],efficiency_data["HMS_Cer_SING_Elec_Eff"],color='red',zorder=4)
plt.ylabel('HMS_Cer_SING_Elec_Eff', fontsize=12)
plt.xlabel('HMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('HMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(144)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data["HMS_Hodoscope_S1X_Rate"],efficiency_data["HMS_Hodo_3_of_4_EFF"],color='red',zorder=4)
plt.ylabel('HMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('HMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('HMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/HMS_S1X_%s.png' % (ROOTPrefix.replace("replay_","")))

plt.figure(figsize=(12,8))

plt.subplot(141)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["COIN_Trigger_Rate"],efficiency_data["Non_Scaler_EDTM_Live_Time"],yerr=efficiency_data["Non_Scaler_EDTM_Live_Time_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["COIN_Trigger_Rate"],efficiency_data["Non_Scaler_EDTM_Live_Time"],color='purple',zorder=4)
plt.ylabel('EDTM', fontsize=12)
plt.xlabel('COIN Trigger Rate [kHz]', fontsize=12)
plt.title('COIN %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(142)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["COIN_Trigger_Rate"],efficiency_data["COIN_CPULT"],yerr=efficiency_data["COIN_CPULT_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["COIN_Trigger_Rate"],efficiency_data["COIN_CPULT"],color='blue',zorder=4)
plt.ylabel('COIN CPULT', fontsize=12)
plt.xlabel('COIN Trigger Rate [kHz]', fontsize=12)
plt.title('COIN %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(143)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data["COIN_Trigger_Rate"],efficiency_data["SHMS_3/4_Trigger_Rate"],color='blue',zorder=4)
plt.ylabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
plt.xlabel('COIN Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(144)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data["COIN_Trigger_Rate"],efficiency_data["HMS_EL-REAL_Trigger_Rate"],color='red',zorder=4)
plt.ylabel('HMS EL-REAL Trigger Rate [kHz]', fontsize=12)
plt.xlabel('COIN Trigger Rate [kHz]', fontsize=12)
plt.title('HMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/COIN_%s.png' % (ROOTPrefix.replace("replay_","")))

#plt.show()
