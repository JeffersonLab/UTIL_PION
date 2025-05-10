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
RunList = sys.argv[4]

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

inp_f = UTILPATH+"/efficiencies/%s_%s_efficiency_data_%s.csv"  % (ROOTPrefix.replace("replay_",""),runType,timestmp)

# Read the run list file
try:
    with open(RunList, 'r') as file:
        run_numbers = [line.strip() for line in file if line.strip()]
except IOError:
    print("Error: %s does not appear to exist." % runlist_file)
    exit(1)

# Convert run numbers to float if necessary
try:
    run_numbers = [float(run) for run in run_numbers]
except ValueError:
    print("Error: One or more run numbers in the run list could not be converted to float.")
    exit(1)

# Converts csv data to dataframe
try:
    csv_data = pd.read_csv(inp_f)
except IOError:
    print("Error: %s does not appear to exist." % inp_f)
#print(efficiency_data.keys())

# Filter the data for the given RunNumber
efficiency_data = csv_data[csv_data['Run_Number'].astype(float).isin(run_numbers)]

print("Plotting for the Following Run Numbers")
print(run_numbers)

#for run in run_numbers:
#    print(run)

#############################################################################################################################################################################

plt.figure(figsize=(12,8))

plt.subplot(121)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.90,1.01)
plt.errorbar(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(122)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.90,1.01)
plt.errorbar(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_Hodo_3_of_4_EFF_%s.png' % (ROOTPrefix.replace("replay_","")))


plt.figure(figsize=(12,8))

plt.subplot(121)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.90,1.01)
plt.errorbar(efficiency_data["HMS_EL-REAL_Trigger_Rate"],efficiency_data["HMS_Hodo_3_of_4_EFF"],yerr=efficiency_data["HMS_Elec_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["HMS_EL-REAL_Trigger_Rate"],efficiency_data["HMS_Hodo_3_of_4_EFF"],color='red',zorder=4)
plt.ylabel('HMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('HMS EL-REAL Trigger Rate [kHz]', fontsize=12)
plt.title('HMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(122)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.90,1.01)
plt.errorbar(efficiency_data["HMS_Hodoscope_S1X_Rate"],efficiency_data["HMS_Hodo_3_of_4_EFF"],yerr=efficiency_data["HMS_Elec_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["HMS_Hodoscope_S1X_Rate"],efficiency_data["HMS_Hodo_3_of_4_EFF"],color='red',zorder=4)
plt.ylabel('HMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('HMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('HMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/HMS_Hodo_3_of_4_EFF_%s.png' % (ROOTPrefix.replace("replay_","")))

########################################################################################################################################################################################
#plt.show()

print("Plotting Complete")
