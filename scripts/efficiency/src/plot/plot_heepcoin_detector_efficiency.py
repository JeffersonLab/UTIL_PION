#! /usr/bin/python
#
# Description:
# ================================================================
# Time-stamp: "2024-03-13 01:29:19 junaid"
# ================================================================
#
# Author:  Muhammad Junaid III <mjo147@uregina.ca>
#
# Copyright (c) junaid
#
###################################################################################################################################################

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

#inp_f = UTILPATH+"/scripts/efficiency/OUTPUTS/%s_%s_efficiency_data_%s.csv"  % (ROOTPrefix.replace("replay_",""),runType,timestmp)
inp_f = UTILPATH+"/efficiencies/%s_%s_detectr_efficiency_data_%s.csv"  % (ROOTPrefix.replace("replay_",""),runType,timestmp)

# Converts csv data to dataframe
try:
    efficiency_data = pd.read_csv(inp_f)
    efficiency_data_wt = pd.read_csv(inp_f)
except IOError:
    print("Error: %s does not appear to exist." % inp_f)
#print(efficiency_data.keys())

# Replace zero errors with a very small number to avoid division by zero
#small_number = 1e-10
small_number = 1e-3
efficiency_data_wt["HMS_Cer_COIN_Elec_Eff_ERROR"] = efficiency_data["HMS_Cer_COIN_Elec_Eff_ERROR"].replace(0, small_number)
efficiency_data_wt["HMS_Cal_COIN_Elec_Eff_ERROR"] = efficiency_data["HMS_Cal_COIN_Elec_Eff_ERROR"].replace(0, small_number)

# Assuming efficiency_data is a pandas DataFrame
# Replace zero errors with NaN to filter them out later
#efficiency_data_wt["HMS_Cer_COIN_Elec_Eff_ERROR"] = efficiency_data["HMS_Cer_COIN_Elec_Eff_ERROR"].replace(0, pd.NA)
#efficiency_data_wt["HMS_Cal_COIN_Elec_Eff_ERROR"] = efficiency_data["HMS_Cal_COIN_Elec_Eff_ERROR"].replace(0, pd.NA)

# Drop rows with NaN values in the specified columns
#efficiency_data_wt = efficiency_data_wt.dropna(subset=["HMS_Cer_COIN_Elec_Eff_ERROR", "HMS_Cal_COIN_Elec_Eff_ERROR"])

#Average_HMS_Cer_SING_Elec_Eff = efficiency_data["HMS_Cer_SING_Elec_Eff"].mean()
#Average_HMS_Cer_SING_Elec_Eff_ERROR = efficiency_data["HMS_Cer_SING_Elec_Eff_ERROR"].mean()
#Average_HMS_Cal_SING_Elec_Eff = efficiency_data["HMS_Cal_SING_Elec_Eff"].mean()
#Average_HMS_Cal_SING_Elec_Eff_ERROR = efficiency_data["HMS_Cal_SING_Elec_Eff_ERROR"].mean()

#Average_HMS_Cer_COIN_Elec_Eff = efficiency_data["HMS_Cer_COIN_Elec_Eff"].mean()
#Average_HMS_Cer_COIN_Elec_Eff_ERROR = efficiency_data["HMS_Cer_COIN_Elec_Eff_ERROR"].mean()
#Average_HMS_Cal_COIN_Elec_Eff = efficiency_data["HMS_Cal_COIN_Elec_Eff"].mean()
#Average_HMS_Cal_COIN_Elec_Eff_ERROR = efficiency_data["HMS_Cal_COIN_Elec_Eff_ERROR"].mean()

#Average_HMS_Cer_ALL_Elec_Eff = efficiency_data["HMS_Cer_ALL_Elec_Eff"].mean()
#Average_HMS_Cer_ALL_Elec_Eff_ERROR = efficiency_data["HMS_Cer_ALL_Elec_Eff_ERROR"].mean()
#Average_HMS_Cal_ALL_Elec_Eff = efficiency_data["HMS_Cal_ALL_Elec_Eff"].mean()
#Average_HMS_Cal_ALL_Elec_Eff_ERROR = efficiency_data["HMS_Cal_ALL_Elec_Eff_ERROR"].mean()

# Calculate the weighted average efficiency and the associated error
#Weighted_Average_HMS_Cer_SING_Elec_Eff = sum(efficiency_data_wt["HMS_Cer_SING_Elec_Eff"] / efficiency_data_wt["HMS_Cer_SING_Elec_Eff_ERROR"] ** 2) / sum(1 / efficiency_data_wt["HMS_Cer_SING_Elec_Eff_ERROR"] ** 2)
#Weighted_Average_HMS_Cer_SING_Elec_Eff_ERROR = (1 / sum(1 / efficiency_data_wt["HMS_Cer_SING_Elec_Eff_ERROR"] ** 2)) ** 0.5
#Weighted_Average_HMS_Cal_SING_Elec_Eff = sum(efficiency_data_wt["HMS_Cal_SING_Elec_Eff"] / efficiency_data_wt["HMS_Cal_SING_Elec_Eff_ERROR"] ** 2) / sum(1 / efficiency_data_wt["HMS_Cal_SING_Elec_Eff_ERROR"] ** 2)
#Weighted_Average_HMS_Cal_SING_Elec_Eff_ERROR = (1 / sum(1 / efficiency_data_wt["HMS_Cal_SING_Elec_Eff_ERROR"] ** 2)) ** 0.5

Weighted_Average_HMS_Cer_COIN_Elec_Eff = sum(efficiency_data_wt["HMS_Cer_COIN_Elec_Eff"] / efficiency_data_wt["HMS_Cer_COIN_Elec_Eff_ERROR"] ** 2) / sum(1 / efficiency_data_wt["HMS_Cer_COIN_Elec_Eff_ERROR"] ** 2)
Weighted_Average_HMS_Cer_COIN_Elec_Eff_ERROR = (1 / sum(1 / efficiency_data_wt["HMS_Cer_COIN_Elec_Eff_ERROR"] ** 2)) ** 0.5
Weighted_Average_HMS_Cal_COIN_Elec_Eff = sum(efficiency_data_wt["HMS_Cal_COIN_Elec_Eff"] / efficiency_data_wt["HMS_Cal_COIN_Elec_Eff_ERROR"] ** 2) / sum(1 / efficiency_data_wt["HMS_Cal_COIN_Elec_Eff_ERROR"] ** 2)
Weighted_Average_HMS_Cal_COIN_Elec_Eff_ERROR = (1 / sum(1 / efficiency_data_wt["HMS_Cal_COIN_Elec_Eff_ERROR"] ** 2)) ** 0.5

#Weighted_Average_HMS_Cer_ALL_Elec_Eff = sum(efficiency_data_wt["HMS_Cer_ALL_Elec_Eff"] / efficiency_data_wt["HMS_Cer_ALL_Elec_Eff_ERROR"] ** 2) / sum(1 / efficiency_data_wt["HMS_Cer_ALL_Elec_Eff_ERROR"] ** 2)
#Weighted_Average_HMS_Cer_ALL_Elec_Eff_ERROR = (1 / sum(1 / efficiency_data_wt["HMS_Cer_ALL_Elec_Eff_ERROR"] ** 2)) ** 0.5
#Weighted_Average_HMS_Cal_ALL_Elec_Eff = sum(efficiency_data_wt["HMS_Cal_ALL_Elec_Eff"] / efficiency_data_wt["HMS_Cal_ALL_Elec_Eff_ERROR"] ** 2) / sum(1 / efficiency_data_wt["HMS_Cal_ALL_Elec_Eff_ERROR"] ** 2)
#Weighted_Average_HMS_Cal_ALL_Elec_Eff_ERROR = (1 / sum(1 / efficiency_data_wt["HMS_Cal_ALL_Elec_Eff_ERROR"] ** 2)) ** 0.5

################################################################################################################################################################################################
'''
plt.figure(figsize=(12,8))
plt.subplot(121)    
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.9,1.02)
plt.errorbar(efficiency_data["HMS_Hodoscope_S1X_Rate"],efficiency_data["HMS_Cer_SING_Elec_Eff"],yerr=efficiency_data["HMS_Cer_SING_Elec_Eff_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["HMS_Hodoscope_S1X_Rate"],efficiency_data["HMS_Cer_SING_Elec_Eff"],color='red',zorder=4)
plt.text(35, 0.95, f'Avg_HMS_Cer_Eff : {Weighted_Average_HMS_Cer_SING_Elec_Eff:.4f} +/- {Weighted_Average_HMS_Cer_SING_Elec_Eff_ERROR:.4f}',
         va='center', ha='left', color='blue', fontsize=14)
plt.ylabel('HMS_Cer_SING_Elec_Eff', fontsize=12)
plt.xlabel('HMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('HMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(122)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.9,1.02)
plt.errorbar(efficiency_data["HMS_Hodoscope_S1X_Rate"],efficiency_data["HMS_Cal_SING_Elec_Eff"],yerr=efficiency_data["HMS_Cal_SING_Elec_Eff_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["HMS_Hodoscope_S1X_Rate"],efficiency_data["HMS_Cal_SING_Elec_Eff"],color='red',zorder=4)
plt.text(35, 0.95, f'Avg_HMS_Cal_Eff : {Weighted_Average_HMS_Cal_SING_Elec_Eff:.4f} +/- {Weighted_Average_HMS_Cal_SING_Elec_Eff_ERROR:.4f}',
         va='center', ha='left', color='blue', fontsize=14)
plt.ylabel('HMS_Cal_SING_Elec_Eff', fontsize=12)
plt.xlabel('HMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('HMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])   
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/HMS_SING_ELec_%s.png' % (ROOTPrefix.replace("replay_","")))
'''
plt.figure(figsize=(12,8))
plt.subplot(121)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.96,1.01)
plt.errorbar(efficiency_data["HMS_Hodoscope_S1X_Rate"],efficiency_data["HMS_Cer_COIN_Elec_Eff"],yerr=efficiency_data["HMS_Cer_COIN_Elec_Eff_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["HMS_Hodoscope_S1X_Rate"],efficiency_data["HMS_Cer_COIN_Elec_Eff"],color='red',zorder=4)
plt.text(35, 0.975, f'Avg_HMS_Cer_Eff : {Weighted_Average_HMS_Cer_COIN_Elec_Eff:.4f} +/- {Weighted_Average_HMS_Cer_COIN_Elec_Eff_ERROR:.4f}',
         va='center', ha='left', color='blue', fontsize=14)
plt.ylabel('HMS_Cer_COIN_Elec_Eff', fontsize=12)
plt.xlabel('HMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('HMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(122)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.96,1.01)
plt.errorbar(efficiency_data["HMS_Hodoscope_S1X_Rate"],efficiency_data["HMS_Cal_COIN_Elec_Eff"],yerr=efficiency_data["HMS_Cal_COIN_Elec_Eff_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["HMS_Hodoscope_S1X_Rate"],efficiency_data["HMS_Cal_COIN_Elec_Eff"],color='red',zorder=4)
plt.text(35, 0.975, f'Avg_HMS_Cal_Eff : {Weighted_Average_HMS_Cal_COIN_Elec_Eff:.4f} +/- {Weighted_Average_HMS_Cal_COIN_Elec_Eff_ERROR:.4f}',
         va='center', ha='left', color='blue', fontsize=14)
plt.ylabel('HMS_Cal_COIN_Elec_Eff', fontsize=12)
plt.xlabel('HMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('HMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/HMS_COIN_ELec_%s.png' % (ROOTPrefix.replace("replay_","")))
'''
plt.figure(figsize=(12,8))
plt.subplot(121)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.9,1.02)
plt.errorbar(efficiency_data["HMS_Hodoscope_S1X_Rate"],efficiency_data["HMS_Cer_ALL_Elec_Eff"],yerr=efficiency_data["HMS_Cer_ALL_Elec_Eff_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["HMS_Hodoscope_S1X_Rate"],efficiency_data["HMS_Cer_ALL_Elec_Eff"],color='red',zorder=4)
plt.text(35, 0.95, f'Avg_HMS_Cer_Eff : {Weighted_Average_HMS_Cer_ALL_Elec_Eff:.4f} +/- {Weighted_Average_HMS_Cer_ALL_Elec_Eff_ERROR:.4f}',
         va='center', ha='left', color='blue', fontsize=14)
plt.ylabel('HMS_Cer_ALL_Elec_Eff', fontsize=12)
plt.xlabel('HMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('HMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(122)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.9,1.02)
plt.errorbar(efficiency_data["HMS_Hodoscope_S1X_Rate"],efficiency_data["HMS_Cal_ALL_Elec_Eff"],yerr=efficiency_data["HMS_Cal_ALL_Elec_Eff_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["HMS_Hodoscope_S1X_Rate"],efficiency_data["HMS_Cal_ALL_Elec_Eff"],color='red',zorder=4)
plt.text(35, 0.95, f'Avg_HMS_Cal_Eff : {Weighted_Average_HMS_Cal_ALL_Elec_Eff:.4f} +/- {Weighted_Average_HMS_Cal_ALL_Elec_Eff_ERROR:.4f}',
         va='center', ha='left', color='blue', fontsize=14)
plt.ylabel('HMS_Cal_ALL_Elec_Eff', fontsize=12)
plt.xlabel('HMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('HMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/HMS_ALL_ELec_%s.png' % (ROOTPrefix.replace("replay_","")))
'''
plt.figure(figsize=(12,8))
plt.subplot(121)
plt.grid(zorder=1)
plt.xlim(2.8,3.8)
plt.ylim(0.96,1.01)
plt.errorbar(efficiency_data["HMS_P_Central"],efficiency_data["HMS_Cer_COIN_Elec_Eff"],yerr=efficiency_data["HMS_Cer_COIN_Elec_Eff_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["HMS_P_Central"],efficiency_data["HMS_Cer_COIN_Elec_Eff"],color='red',zorder=4)
plt.text(2.9, 0.975, f'Avg_HMS_Cer_Eff : {Weighted_Average_HMS_Cer_COIN_Elec_Eff:.4f} +/- {Weighted_Average_HMS_Cer_COIN_Elec_Eff_ERROR:.4f}',
         va='center', ha='left', color='blue', fontsize=14)
plt.ylabel('HMS_Cer_COIN_Elec_Eff', fontsize=12)
plt.xlabel('HMS Central Momentum (MeV)', fontsize=12)
plt.title('HMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(122)
plt.grid(zorder=1)
plt.xlim(2.8,3.8)
plt.ylim(0.96,1.01)
plt.errorbar(efficiency_data["HMS_P_Central"],efficiency_data["HMS_Cal_COIN_Elec_Eff"],yerr=efficiency_data["HMS_Cal_COIN_Elec_Eff_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["HMS_P_Central"],efficiency_data["HMS_Cal_COIN_Elec_Eff"],color='red',zorder=4)
plt.text(2.9, 0.975, f'Avg_HMS_Cal_Eff : {Weighted_Average_HMS_Cal_COIN_Elec_Eff:.4f} +/- {Weighted_Average_HMS_Cal_COIN_Elec_Eff_ERROR:.4f}',
         va='center', ha='left', color='blue', fontsize=14)
plt.ylabel('HMS_Cal_COIN_Elec_Eff', fontsize=12)
plt.xlabel('HMS Central Momentum (MeV)', fontsize=12)
plt.title('HMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/HMS_momentum_vs_detctr_eff_%s.png' % (ROOTPrefix.replace("replay_","")))

########################################################################################################################################################################################

#plt.show()
#print(f"Avg_HMS_Cer_SING_Elec_Eff : {Average_HMS_Cer_SING_Elec_Eff:.4f} +/- {Average_HMS_Cer_SING_Elec_Eff_ERROR:.4f}")
#print(f"WtAvg_HMS_Cer_SING_Elec_Eff : {Weighted_Average_HMS_Cer_SING_Elec_Eff:.4f} +/- {Weighted_Average_HMS_Cer_SING_Elec_Eff_ERROR:.4f}")
#print(f"Avg_HMS_Cal_SING_Elec_Eff : {Average_HMS_Cal_SING_Elec_Eff:.4f} +/- {Average_HMS_Cal_SING_Elec_Eff_ERROR:.4f}")
#print(f"WtAvg_HMS_Cal_SING_Elec_Eff : {Weighted_Average_HMS_Cal_SING_Elec_Eff:.4f} +/- {Weighted_Average_HMS_Cal_SING_Elec_Eff_ERROR:.4f}")

#print(f"Avg_HMS_Cer_COIN_Elec_Eff : {Average_HMS_Cer_COIN_Elec_Eff:.4f} +/- {Average_HMS_Cer_COIN_Elec_Eff_ERROR:.4f}")
print(f"WtAvg_HMS_Cer_COIN_Elec_Eff : {Weighted_Average_HMS_Cer_COIN_Elec_Eff:.4f} +/- {Weighted_Average_HMS_Cer_COIN_Elec_Eff_ERROR:.4f}")
#print(f"Avg_HMS_Cal_COIN_Elec_Eff : {Average_HMS_Cal_COIN_Elec_Eff:.4f} +/- {Average_HMS_Cal_COIN_Elec_Eff_ERROR:.4f}")
print(f"WtAvg_HMS_Cal_COIN_Elec_Eff : {Weighted_Average_HMS_Cal_COIN_Elec_Eff:.4f} +/- {Weighted_Average_HMS_Cal_COIN_Elec_Eff_ERROR:.4f}")

#print(f"Avg_HMS_Cer_ALL_Elec_Eff : {Average_HMS_Cer_ALL_Elec_Eff:.4f} +/- {Average_HMS_Cer_ALL_Elec_Eff_ERROR:.4f}")
#print(f"WtAvg_HMS_Cer_ALL_Elec_Eff : {Weighted_Average_HMS_Cer_ALL_Elec_Eff:.4f} +/- {Weighted_Average_HMS_Cer_ALL_Elec_Eff_ERROR:.4f}")
#print(f"Avg_HMS_Cal_ALL_Elec_Eff : {Average_HMS_Cal_ALL_Elec_Eff:.4f} +/- {Average_HMS_Cal_ALL_Elec_Eff_ERROR:.4f}")
#print(f"WtAvg_HMS_Cal_ALL_Elec_Eff : {Weighted_Average_HMS_Cal_ALL_Elec_Eff:.4f} +/- {Weighted_Average_HMS_Cal_ALL_Elec_Eff_ERROR:.4f}")

print("Plotting Complete")
