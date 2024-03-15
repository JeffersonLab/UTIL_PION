#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-09 12:05:00 junaid"
# ================================================================
#
# Created: Muhammad junaid  <mjo147@uregina.ca>
# Copyright (c) trottar & junaid
#

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from csv import DictReader
import sys, os
import matplotlib.gridspec as gridspec

################################################################################################################################################
'''
User Inputs
'''
'''
ROOTPrefix1 = sys.argv[1]
runType1 = sys.argv[2]
timestmp1 = sys.argv[3]

ROOTPrefix2 = sys.argv[4]
runType2 = sys.argv[5]
timestmp2 = sys.argv[6]

ROOTPrefix3 = sys.argv[7]
runType3 = sys.argv[8]
timestmp3 = sys.argv[9]
'''
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

#inp_f1 = UTILPATH+"/scripts/efficiency/OUTPUTS/%s_%s_efficiency_data_%s.csv"  % (ROOTPrefix1.replace("replay_",""),runType1,timestmp1)
#inp_f2 = UTILPATH+"/scripts/efficiency/OUTPUTS/%s_%s_efficiency_data_%s.csv"  % (ROOTPrefix2.replace("replay_",""),runType2,timestmp2)
#inp_f3 = UTILPATH+"/scripts/efficiency/OUTPUTS/%s_%s_efficiency_data_%s.csv"  % (ROOTPrefix3.replace("replay_",""),runType3,timestmp3)

inp_f1 = UTILPATH+"/scripts/efficiency/OUTPUTS/PionLT_SHMS_HeePSing_HeePSing_efficiency_data_2024_03_06.csv"
inp_f2 = UTILPATH+"/scripts/efficiency/OUTPUTS/PionLT_SHMS_Lumi_LumiSing_efficiency_data_2024_03_06.csv"
inp_f3 = UTILPATH+"/scripts/efficiency/OUTPUTS/PionLT_coin_production_Prod_efficiency_data_2024_03_02.csv"

efficiency_data1 = pd.read_csv(inp_f1)
efficiency_data2 = pd.read_csv(inp_f2)
efficiency_data3 = pd.read_csv(inp_f3)



#############################################################################################################################################################################
'''
plt.figure(figsize=(6,4))
plt.subplot(111)    
#plt.grid(zorder=1)
plt.xlim(0,1200)
plt.ylim(0.92,1.05)
plt.errorbar(efficiency_data1["SHMS_EL-REAL_Trigger_Rate"],efficiency_data1["SHMS_Elec_SING_TRACK_EFF"],yerr=efficiency_data1["SHMS_Elec_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data1["SHMS_EL-REAL_Trigger_Rate"],efficiency_data1["SHMS_Elec_SING_TRACK_EFF"],color='red',zorder=4, s=5.0)
plt.errorbar(efficiency_data2["SHMS_EL-REAL_Trigger_Rate"],efficiency_data2["SHMS_Elec_SING_TRACK_EFF"],yerr=efficiency_data2["SHMS_Elec_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data2["SHMS_EL-REAL_Trigger_Rate"],efficiency_data2["SHMS_Elec_SING_TRACK_EFF"],color='red',zorder=4, s=5.0)
plt.ylabel('Electron Tracking Efficiency', fontsize=12)
plt.xlabel('SHMS EL-REAL Trigger Rate [kHz]', fontsize=12)
plt.locator_params(axis='x', nbins=20) ### set number of bins for x axis only
plt.xticks(rotation=90)
plt.tight_layout(rect=[0,0.03,1,0.95])   
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_Electron_EFF_vs_3_4_Trigger_Rate.png')
'''
plt.figure(figsize=(6,4))
plt.subplot(111)
#plt.grid(zorder=1)
plt.xlim(0,3200)
plt.ylim(0.92,1.05)
plt.errorbar(efficiency_data1["SHMS_Hodoscope_S1X_Rate"],efficiency_data1["SHMS_Elec_SING_TRACK_EFF"],yerr=efficiency_data1["SHMS_Elec_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data1["SHMS_Hodoscope_S1X_Rate"],efficiency_data1["SHMS_Elec_SING_TRACK_EFF"],color='red',zorder=4, s=5.0)
plt.errorbar(efficiency_data2["SHMS_Hodoscope_S1X_Rate"],efficiency_data2["SHMS_Elec_SING_TRACK_EFF"],yerr=efficiency_data2["SHMS_Elec_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data2["SHMS_Hodoscope_S1X_Rate"],efficiency_data2["SHMS_Elec_SING_TRACK_EFF"],color='red',zorder=4, s=5.0)
plt.ylabel('Electron Tracking Efficiency', fontsize=14)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=14)
plt.locator_params(axis='x', nbins=20) ### set number of bins for x axis only
plt.xticks(rotation=90)
plt.tight_layout(rect=[0,0.03,1,0.95])   
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_Electron_EFF_vs_S1X_Rate.png', dpi=300)

########################################################################################################################################################################################
'''
plt.figure(figsize=(6,4))
plt.subplot(111)
#plt.grid(zorder=1)
plt.xlim(0,1200)
plt.ylim(0.92,1.02)
plt.errorbar(efficiency_data3["SHMS_3/4_Trigger_Rate"],efficiency_data3["SHMS_Pion_SING_TRACK_EFF"],yerr=efficiency_data3["SHMS_Pion_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data3["SHMS_3/4_Trigger_Rate"],efficiency_data3["SHMS_Pion_SING_TRACK_EFF"],color='blue',zorder=4, s=5.0)
plt.ylabel('Pion Tracking Efficiency', fontsize=12)
plt.xlabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
plt.locator_params(axis='x', nbins=20) ### set number of bins for x axis only
plt.xticks(rotation=90)
plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_Pion_EFF_vs_3_4_Trigger_Rate.png')
'''
plt.figure(figsize=(6, 4))
plt.subplot(111)
#plt.grid(zorder=1)
plt.xlim(0,3200)
plt.ylim(0.92,1.02)
plt.errorbar(efficiency_data3["SHMS_Hodoscope_S1X_Rate"],efficiency_data3["SHMS_Pion_SING_TRACK_EFF"],yerr=efficiency_data3["SHMS_Pion_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data3["SHMS_Hodoscope_S1X_Rate"],efficiency_data3["SHMS_Pion_SING_TRACK_EFF"],color='blue',zorder=4, s=5.0)
plt.ylabel('Pion Tracking Efficiency', fontsize=14)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=14)
plt.locator_params(axis='x', nbins=20) ### set number of bins for x axis only
plt.xticks(rotation=90)
plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_Pion_EFF_vs_S1X_Rate.png', dpi=300)

################################################################################################################################################################################################

fig, (ax1, ax2) = plt.subplots(2,figsize=(8, 7), sharex=True)
ax1.errorbar(efficiency_data1["SHMS_Hodoscope_S1X_Rate"],efficiency_data1["SHMS_Elec_SING_TRACK_EFF"],yerr=efficiency_data1["SHMS_Elec_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
ax1.scatter(efficiency_data1["SHMS_Hodoscope_S1X_Rate"],efficiency_data1["SHMS_Elec_SING_TRACK_EFF"],color='red',zorder=4, s=5.0)
ax1.errorbar(efficiency_data2["SHMS_Hodoscope_S1X_Rate"],efficiency_data2["SHMS_Elec_SING_TRACK_EFF"],yerr=efficiency_data2["SHMS_Elec_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
ax1.scatter(efficiency_data2["SHMS_Hodoscope_S1X_Rate"],efficiency_data2["SHMS_Elec_SING_TRACK_EFF"],color='red',zorder=4, s=5.0)
ax1.set_ylim(0.92,1.02)
ax1.tick_params(axis='y', labelsize=12)
ax1.set_ylabel('Electron Tracking Efficiency', fontsize=14)

ax2.errorbar(efficiency_data3["SHMS_Hodoscope_S1X_Rate"],efficiency_data3["SHMS_Pion_SING_TRACK_EFF"],yerr=efficiency_data3["SHMS_Pion_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
ax2.scatter(efficiency_data3["SHMS_Hodoscope_S1X_Rate"],efficiency_data3["SHMS_Pion_SING_TRACK_EFF"],color='blue',zorder=4, s=5.0)
ax2.set_xlim(0,3200)
ax2.set_ylabel('Pion Tracking Efficiency', fontsize=14)
ax2.set_xlabel('SHMS EL-REAL Trigger Rate [kHz]', fontsize=14)
ax2.tick_params(axis='x', labelsize=12)

ax2.set_ylim(0.92, 1.02)
ax2.tick_params(axis='y', labelsize=12)
ax1.locator_params(axis='x', nbins=20) ### set number of bins for x axis only
ax2.locator_params(axis='x', nbins=20) ### set number of bins for x axis only
ax2.tick_params(axis='x', rotation=90)
plt.subplots_adjust(hspace=0.01)
plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_Track_EFF_vs_S1X_Rate.png', dpi=300)

################################################################################################################################################################################################
#plt.show()

print("Plotting Complete")
