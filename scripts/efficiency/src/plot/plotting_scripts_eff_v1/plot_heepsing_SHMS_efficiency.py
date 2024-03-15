#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-09 12:05:00 junaid"
# ================================================================
#
# Copied: Richard L. Trotta III <trotta@cua.edu>
# Created: Muhammad junaid  <mjo147@uregina.ca>
# Copyright (c) trottar & junaid
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
#print(efficiency_data.keys())

#############################################################################################################################################################################

plt.figure(figsize=(14,10))

plt.subplot(141)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["SHMS_EL-REAL_Trigger_Rate"],efficiency_data["Non_Scaler_EDTM_Live_Time"],yerr=efficiency_data["Non_Scaler_EDTM_Live_Time_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_EL-REAL_Trigger_Rate"],efficiency_data["Non_Scaler_EDTM_Live_Time"],color='blue',zorder=4)
plt.ylabel('EDTM', fontsize=12)
plt.xlabel('SHMS EL-REAL Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(142)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["SHMS_EL-REAL_Trigger_Rate"],efficiency_data["SHMS_Elec_SING_TRACK_EFF"],yerr=efficiency_data["SHMS_Elec_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_EL-REAL_Trigger_Rate"],efficiency_data["SHMS_Elec_SING_TRACK_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Elec_SING_TRACK_EFF', fontsize=12)
plt.xlabel('SHMS EL-REAL Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(143)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["SHMS_EL-REAL_Trigger_Rate"],efficiency_data["SHMS_Elec_ALL_TRACK_EFF"],yerr=efficiency_data["SHMS_Elec_ALL_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_EL-REAL_Trigger_Rate"],efficiency_data["SHMS_Elec_ALL_TRACK_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Elec_ALL_TRACK_EFF', fontsize=12)
plt.xlabel('SHMS EL-REAL Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(144)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data["SHMS_EL-REAL_Trigger_Rate"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('SHMS EL-REAL Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])   
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_EL-REAL_%s.png' % (ROOTPrefix.replace("replay_","")))

########################################################################################################################################################################################
'''
plt.figure(figsize=(10,14))

plt.subplot(411)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["Run_Number"],efficiency_data["Non_Scaler_EDTM_Live_Time"],yerr=efficiency_data["Non_Scaler_EDTM_Live_Time_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["Run_Number"],efficiency_data["Non_Scaler_EDTM_Live_Time"],color='blue',zorder=4)
plt.ylabel('EDTM', fontsize=12)
plt.xlabel('Run Number (HeePSingSHMS)', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(412)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["Run_Number"],efficiency_data["SHMS_Elec_SING_TRACK_EFF"],yerr=efficiency_data["SHMS_Elec_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["Run_Number"],efficiency_data["SHMS_Elec_SING_TRACK_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Elec_SING_TRACK_EFF', fontsize=12)
plt.xlabel('Run Number (HeePSingSHMS)', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(413)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["Run_Number"],efficiency_data["SHMS_Elec_ALL_TRACK_EFF"],yerr=efficiency_data["SHMS_Elec_ALL_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["Run_Number"],efficiency_data["SHMS_Elec_ALL_TRACK_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Elec_ALL_TRACK_EFF', fontsize=12)
plt.xlabel('Run Number (HeePSingSHMS)', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(414)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data["Run_Number"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('Run Number (HeePSingSHMS)', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/Run_Number_SHMS_%s.png' % (ROOTPrefix.replace("replay_","")))
'''

#########################################################################################################################################################################################

plt.figure(figsize=(14,10))

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
plt.errorbar(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Elec_SING_TRACK_EFF"],yerr=efficiency_data["SHMS_Elec_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Elec_SING_TRACK_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Elec_SING_TRACK_EFF', fontsize=12)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(143)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Elec_ALL_TRACK_EFF"],yerr=efficiency_data["SHMS_Elec_ALL_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Elec_ALL_TRACK_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Elec_ALL_TRACK_EFF', fontsize=12)
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

###################################################################################################################################################################################################
'''
plt.figure(figsize=(10,14))

plt.subplot(411)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.errorbar(efficiency_data["Run_Number"],efficiency_data["SHMS_CPULT"],yerr=efficiency_data["SHMS_CPULT_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["Run_Number"],efficiency_data["SHMS_CPULT"],color='green',zorder=4)
plt.ylabel('SHMS CPULT', fontsize=12)
plt.xlabel('Run Number (HeePSingSHMS)', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(412)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data["Run_Number"],efficiency_data["SHMS_Hodoscope_S1X_Rate"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodoscope_S1X_Rate [kHz]', fontsize=12)
plt.xlabel('Run Number (HeePSingSHMS)', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(413)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data["Run_Number"],efficiency_data["SHMS_EL-REAL_Trigger_Rate"],color='blue',zorder=4)
plt.ylabel('SHMS_EL-REAL_Trigger_Rate [kHz]', fontsize=12)
plt.xlabel('Run Number (HeePSingSHMS)', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/Run_Number_rates_%s.png' % (ROOTPrefix.replace("replay_","")))
'''
################################################################################################################################################################################################

#plt.show()

print("Plotting Complete")
