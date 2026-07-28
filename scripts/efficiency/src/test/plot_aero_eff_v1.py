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
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from csv import DictReader
import sys, os
import re

################################################################################################################################################
'''
User Inputs
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

inp_f = "/u/group/c-pionlt/USERS/junaid/hallc_replay_lt/UTIL_PION/efficiencies/test/PionLT_coin_production_Prod_efficiency_data_2025_11_21.csv"
runlist = "/u/group/c-pionlt/USERS/junaid/hallc_replay_lt/UTIL_BATCH/InputRunLists/eff_runlist/PionLT_Prod_runlist"

# Converts csv data to dataframe
try:
    efficiency_data = pd.read_csv(inp_f)
except IOError:
    print("Error: %s does not appear to exist." % inp_f)
#print(efficiency_data.keys())

#############################################################################################################################################################################

input_runs = []

try:
    with open(runlist, "r") as f:
        for line in f:
            line = line.strip()

            # Skip blank lines and comments
            if line == "" or line.startswith("#"):
                continue

            # Extract run numbers from the line
            # This works for simple lists like:
            # 12345
            # or lines with extra text containing run numbers
            nums = re.findall(r"\d+", line)

            for num in nums:
                input_runs.append(int(num))

except IOError:
    print("Error: %s does not appear to exist." % runlist)
    sys.exit(1)

input_runs = sorted(list(set(input_runs)))

if len(input_runs) == 0:
    print("Error: No run numbers found in runlist.")
    sys.exit(1)
    
efficiency_data = efficiency_data[
    efficiency_data["Run_Number"].astype(int).isin(input_runs)
].copy()

efficiency_data = efficiency_data.sort_values("Run_Number")

if efficiency_data.empty:
    print("Error: No matching run numbers found in efficiency CSV.")
    print("Runlist contains:", input_runs)
    sys.exit(1)

#print("Plotting the following runs:")
#print(efficiency_data["Run_Number"].astype(int).to_list())    
    
################################################################################################################################################################


plt.figure(figsize=(12,8))

plt.subplot(111)    
plt.grid(zorder=1)
plt.xlim(0,800)
plt.ylim(0.8,1.02)
plt.errorbar(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Aero_COIN_Pion_Eff"],yerr=efficiency_data["SHMS_Aero_COIN_Pion_Eff_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Aero_COIN_Pion_Eff"],color='blue',zorder=4)
plt.ylabel('Pion Efficiency (%)', fontsize=12)
plt.xlabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig("SHMS_aero1_eff.png")

plt.figure(figsize=(12,8))

plt.subplot(111)
plt.grid(zorder=1)
plt.xlim(0,3000)
plt.ylim(0.8,1.02)
plt.errorbar(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Aero_COIN_Pion_Eff"],yerr=efficiency_data["SHMS_Aero_COIN_Pion_Eff_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Aero_COIN_Pion_Eff"],color='blue',zorder=4)
plt.ylabel('Pion Efficiency (%)', fontsize=12)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])   
plt.savefig("SHMS_aero2_eff.png")

########################################################################################################################################################################################

plt.figure(figsize=(12,8))

plt.subplot(111)
plt.grid(zorder=1)
plt.xlim(0,2.5)
plt.ylim(0.8,1.02)
plt.errorbar(efficiency_data["COIN_Trigger_Rate"],efficiency_data["Non_Scaler_EDTM_Live_Time_Corr"],yerr=efficiency_data["Non_Scaler_EDTM_Live_Time_Corr_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["COIN_Trigger_Rate"],efficiency_data["Non_Scaler_EDTM_Live_Time_Corr"],color='blue',zorder=4)
plt.ylabel('EDTM LT (%)', fontsize=12)
plt.xlabel('COIN_Trigger_Rate [MHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig("edtm_eff.png")

###################################################################################################################################################################################################

#plt.show()

print("Plotting Complete")

