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
plt.ylim(0.9,1.01)
plt.errorbar(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Hodo_Plane_1"],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Hodo_Plane_1"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_Plane_1', fontsize=12)
plt.xlabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(122)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.9,1.01)
plt.errorbar(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Hodo_Plane_1"],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Hodo_Plane_1"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_Plane_1', fontsize=12)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])   
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_Hodo_Plane_1_%s.png' % (ROOTPrefix.replace("replay_","")))

################################################################################################################################################################################################

plt.figure(figsize=(12,8))

plt.subplot(121)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.8,1.01)
plt.errorbar(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Hodo_Plane_2"],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Hodo_Plane_2"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_Plane_2', fontsize=12)
plt.xlabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12/1)

plt.subplot(122)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.8,1.01)
plt.errorbar(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Hodo_Plane_2"],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Hodo_Plane_2"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_Plane_2', fontsize=12)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_Hodo_Plane_2_%s.png' % (ROOTPrefix.replace("replay_","")))

########################################################################################################################################################################################

plt.figure(figsize=(12,8))

plt.subplot(121)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.6,1.01)
plt.errorbar(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Hodo_Plane_3"],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Hodo_Plane_3"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_Plane_3', fontsize=12)
plt.xlabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(122)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.6,1.01)
plt.errorbar(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Hodo_Plane_3"],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Hodo_Plane_3"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_Plane_3', fontsize=12)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_Hodo_Plane_3_%s.png' % (ROOTPrefix.replace("replay_","")))

########################################################################################################################################################################################

plt.figure(figsize=(12,8))

plt.subplot(121)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.3,1.01)
plt.errorbar(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Hodo_Plane_4"],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Hodo_Plane_4"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_Plane_4', fontsize=12)
plt.xlabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(122)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.3,1.01)
plt.errorbar(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Hodo_Plane_4"],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Hodo_Plane_4"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_Plane_4', fontsize=12)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_Hodo_Plane_4_%s.png' % (ROOTPrefix.replace("replay_","")))

##########################################################################################################################################################################################

plt.figure(figsize=(12,8))

plt.subplot(121)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.75,1.01)
plt.errorbar(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(122)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.75,1.01)
plt.errorbar(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_Hodo_3_of_4_EFF_%s.png' % (ROOTPrefix.replace("replay_","")))

########################################################################################################################################################################################

plt.figure(figsize=(12,8))

plt.subplot(121)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.9,1.01)
plt.errorbar(efficiency_data["SHMS_P_Central"],efficiency_data["SHMS_Hodo_Plane_1"],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_P_Central"],efficiency_data["SHMS_Hodo_Plane_1"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_Plane_1', fontsize=12)
plt.xlabel('SHMS Central Momentum [GeV/c]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(122)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.8,1.01)
plt.errorbar(efficiency_data["SHMS_P_Central"],efficiency_data["SHMS_Hodo_Plane_2"],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_P_Central"],efficiency_data["SHMS_Hodo_Plane_2"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_Plane_2', fontsize=12)
plt.xlabel('SHMS Central Momentum [GeV/c]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12/1)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_Hodo_Plane_EFF_vs_mom_v1_%s.png' % (ROOTPrefix.replace("replay_","")))

########################################################################################################################################################################################

plt.figure(figsize=(12,8))

plt.subplot(121)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.6,1.01)
plt.errorbar(efficiency_data["SHMS_P_Central"],efficiency_data["SHMS_Hodo_Plane_3"],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_P_Central"],efficiency_data["SHMS_Hodo_Plane_3"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_Plane_3', fontsize=12)
plt.xlabel('SHMS Central Momentum [GeV/c]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(122)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.3,1.01)
plt.errorbar(efficiency_data["SHMS_P_Central"],efficiency_data["SHMS_Hodo_Plane_4"],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["SHMS_P_Central"],efficiency_data["SHMS_Hodo_Plane_4"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_Plane_4', fontsize=12)
plt.xlabel('SHMS Central Momentum [GeV/c]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_Hodo_Plane_EFF_vs_mom_v2_%s.png' % (ROOTPrefix.replace("replay_","")))

############################################################################################################################################################################################

plt.figure(figsize=(12,8))

mask_1 = (efficiency_data["SHMS_P_Central"] <= 2.5)
mask_2 = (efficiency_data["SHMS_P_Central"] > 2.5) & (efficiency_data["SHMS_P_Central"] <= 3.0)
mask_3 = (efficiency_data["SHMS_P_Central"] > 3.0) & (efficiency_data["SHMS_P_Central"] <= 4.5)
mask_4 = (efficiency_data["SHMS_P_Central"] > 4.5) & (efficiency_data["SHMS_P_Central"] <= 6.0)
mask_5 = (efficiency_data["SHMS_P_Central"] >= 6.0)

plt.subplot(111)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.85,1.01)
plt.scatter(efficiency_data["SHMS_P_Central"][mask_1], efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_1], color='brown', label='<2.5 GeV/c', zorder=4)
plt.scatter(efficiency_data["SHMS_P_Central"][mask_2], efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_2], color='black', label='2.5-3.0 GeV/c', zorder=4)
plt.scatter(efficiency_data["SHMS_P_Central"][mask_3], efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_3], color='blue', label='3.0-4.5  GeV/c', zorder=4)
plt.scatter(efficiency_data["SHMS_P_Central"][mask_4], efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_4], color='green', label='4.5-6.0 GeV/c', zorder=4)
plt.scatter(efficiency_data["SHMS_P_Central"][mask_5], efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_5], color='red', label='>6.0- GeV/c', zorder=4)
plt.errorbar(efficiency_data["SHMS_P_Central"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"],color='black',linestyle='None',zorder=3)
#plt.scatter(efficiency_data["SHMS_P_Central"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('SHMS Central Momentum [GeV/c]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)
plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_Hodo_3_of_4_EFF_vs_mom_%s.png' % (ROOTPrefix.replace("replay_","")))

############################################################################################################################################################################################

plt.figure(figsize=(12,8))

plt.subplot(121)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.75,1.01)
plt.scatter(efficiency_data["SHMS_3/4_Trigger_Rate"][mask_1], efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_1], color='brown', label='<2.5 GeV/c', zorder=4)
plt.errorbar(efficiency_data["SHMS_3/4_Trigger_Rate"][mask_1],efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_1],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"][mask_1],color='black',linestyle='None',zorder=3)
#plt.scatter(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(122)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.75,1.01)
plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"][mask_1], efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_1], color='brown', label='<2.5 GeV/c', zorder=4)
plt.errorbar(efficiency_data["SHMS_Hodoscope_S1X_Rate"][mask_1],efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_1],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"][mask_1],color='black',linestyle='None',zorder=3)
#plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_Hodo_3_of_4_EFF_%s_v1.png' % (ROOTPrefix.replace("replay_","")))


############################################################################################################################################################################################

plt.figure(figsize=(12,8))

plt.subplot(121)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.75,1.01)
plt.scatter(efficiency_data["SHMS_3/4_Trigger_Rate"][mask_2], efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_2], color='black', label='2.5-3.0 GeV/c', zorder=4)
plt.errorbar(efficiency_data["SHMS_3/4_Trigger_Rate"][mask_2],efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_2],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"][mask_2],color='black',linestyle='None',zorder=3)
#plt.scatter(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(122)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.75,1.01)
plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"][mask_2], efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_2], color='black', label='2.5-3.0 GeV/c', zorder=4)
plt.errorbar(efficiency_data["SHMS_Hodoscope_S1X_Rate"][mask_2],efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_2],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"][mask_2],color='black',linestyle='None',zorder=3)
#plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_Hodo_3_of_4_EFF_%s_v2.png' % (ROOTPrefix.replace("replay_","")))

##############################################################################################################################################################################################

plt.figure(figsize=(12,8))

plt.subplot(121)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.75,1.01)
plt.scatter(efficiency_data["SHMS_3/4_Trigger_Rate"][mask_3], efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_3], color='blue', label='3.0-4.5 GeV/c', zorder=4)
plt.errorbar(efficiency_data["SHMS_3/4_Trigger_Rate"][mask_3],efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_3],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"][mask_3],color='black',linestyle='None',zorder=3)
#plt.scatter(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(122)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.75,1.01)
plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"][mask_3], efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_3], color='blue', label='3.0-4.5 GeV/c', zorder=4)
plt.errorbar(efficiency_data["SHMS_Hodoscope_S1X_Rate"][mask_3],efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_3],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"][mask_3],color='black',linestyle='None',zorder=3)
#plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_Hodo_3_of_4_EFF_%s_v3.png' % (ROOTPrefix.replace("replay_","")))

##############################################################################################################################################################################################

plt.figure(figsize=(12,8))

plt.subplot(121)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.75,1.01)
plt.scatter(efficiency_data["SHMS_3/4_Trigger_Rate"][mask_4], efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_4], color='green', label='4.5-6.0 GeV/c', zorder=4)
plt.errorbar(efficiency_data["SHMS_3/4_Trigger_Rate"][mask_4],efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_4],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"][mask_4],color='black',linestyle='None',zorder=3)
#plt.scatter(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(122)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.75,1.01)
plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"][mask_4], efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_4], color='green', label='4.5-6.0 GeV/c', zorder=4)
plt.errorbar(efficiency_data["SHMS_Hodoscope_S1X_Rate"][mask_4],efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_4],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"][mask_4],color='black',linestyle='None',zorder=3)
#plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_Hodo_3_of_4_EFF_%s_v4.png' % (ROOTPrefix.replace("replay_","")))

#############################################################################################################################################################################################

plt.figure(figsize=(12,8))

plt.subplot(121)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.75,1.01)
plt.scatter(efficiency_data["SHMS_3/4_Trigger_Rate"][mask_5], efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_5], color='red', label='>6.0 GeV/c', zorder=4)
plt.errorbar(efficiency_data["SHMS_3/4_Trigger_Rate"][mask_5],efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_5],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"][mask_5],color='black',linestyle='None',zorder=3)
#plt.scatter(efficiency_data["SHMS_3/4_Trigger_Rate"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(122)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.75,1.01)
plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"][mask_5], efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_5], color='red', label='>6.0 GeV/c', zorder=4)
plt.errorbar(efficiency_data["SHMS_Hodoscope_S1X_Rate"][mask_5],efficiency_data["SHMS_Hodo_3_of_4_EFF"][mask_5],yerr=efficiency_data["SHMS_Prot_SING_TRACK_EFF_ERROR"][mask_5],color='black',linestyle='None',zorder=3)
#plt.scatter(efficiency_data["SHMS_Hodoscope_S1X_Rate"],efficiency_data["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_Hodo_3_of_4_EFF_%s_v5.png' % (ROOTPrefix.replace("replay_","")))

##############################################################################################################################################################################################

plt.figure(figsize=(12,8))

plt.subplot(121)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.9,1.01)
plt.errorbar(efficiency_data["COIN_Trigger_Rate"],efficiency_data["Non_Scaler_EDTM_Live_Time_Corr"],yerr=efficiency_data["Non_Scaler_EDTM_Live_Time_Corr_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["COIN_Trigger_Rate"],efficiency_data["Non_Scaler_EDTM_Live_Time_Corr"],color='blue',zorder=4)
plt.ylabel('EDTM Live Time', fontsize=12)
plt.xlabel('COIN Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(122)
plt.grid(zorder=1)
#plt.xlim(0,100)
plt.ylim(0.9,1.01)
plt.errorbar(efficiency_data["COIN_Trigger_Rate"],efficiency_data["COIN_CPULT"],yerr=efficiency_data["COIN_CPULT_ERROR"],color='black',linestyle='None',zorder=3)
plt.scatter(efficiency_data["COIN_Trigger_Rate"],efficiency_data["COIN_CPULT"],color='blue',zorder=4)
plt.ylabel('CPULT', fontsize=12)
plt.xlabel('COIN Trigger Rate [kHz]', fontsize=12)
plt.title('SHMS %s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_LT_vs_Rate_%s.png' % (ROOTPrefix.replace("replay_","")))

#############################################################################################################################################################################################
#plt.show()

print("Plotting Complete")
