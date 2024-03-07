#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-12-20 00:24:44 trottar"
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
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
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

# Including dummy
#efficiency_data_10p6 = efficiency_data[(efficiency_data['Run_Number'] >= 4865)  & (efficiency_data['Run_Number'] <= 5334)]
#efficiency_data_3p8 = efficiency_data[(efficiency_data['Run_Number'] >= 6638)  & (efficiency_data['Run_Number'] <= 6857)]
#efficiency_data_4p9 = efficiency_data[(efficiency_data['Run_Number'] >= 6885)  & (efficiency_data['Run_Number'] <= 7045)]
#efficiency_data_6p2 = efficiency_data[(efficiency_data['Run_Number'] >= 7871)  & (efficiency_data['Run_Number'] <= 7938)]
#efficiency_data_8p2 = efficiency_data[(efficiency_data['Run_Number'] >= 7978)  & (efficiency_data['Run_Number'] <= 8356)]

# Update 'your_file.txt' with the actual file path
with open(REPLAYPATH+'/UTIL_BATCH/InputRunLists/KaonLT_2018_2019/Prod_10p6_Autumn18', 'r') as file:
    # Assuming each line in the file contains a single Run_Number
    run_numbers = [int(line.strip()) for line in file]

# Assuming 'efficiency_data' is a DataFrame with a column named 'Run_Number'
efficiency_data_10p6 = efficiency_data[efficiency_data['Run_Number'].isin(run_numbers)]    

# Update 'your_file.txt' with the actual file path
with open(REPLAYPATH+'/UTIL_BATCH/InputRunLists/KaonLT_2018_2019/Prod_8p2_Spring19', 'r') as file:
    # Assuming each line in the file contains a single Run_Number
    run_numbers = [int(line.strip()) for line in file]

# Assuming 'efficiency_data' is a DataFrame with a column named 'Run_Number'
efficiency_data_8p2 = efficiency_data[efficiency_data['Run_Number'].isin(run_numbers)]    

# Update 'your_file.txt' with the actual file path
with open(REPLAYPATH+'/UTIL_BATCH/InputRunLists/KaonLT_2018_2019/Prod_6p2_Spring19', 'r') as file:
    # Assuming each line in the file contains a single Run_Number
    run_numbers = [int(line.strip()) for line in file]

# Assuming 'efficiency_data' is a DataFrame with a column named 'Run_Number'
efficiency_data_6p2 = efficiency_data[efficiency_data['Run_Number'].isin(run_numbers)]    

# Update 'your_file.txt' with the actual file path
with open(REPLAYPATH+'/UTIL_BATCH/InputRunLists/KaonLT_2018_2019/Prod_4p9_Autumn18', 'r') as file:
    # Assuming each line in the file contains a single Run_Number
    run_numbers = [int(line.strip()) for line in file]

# Assuming 'efficiency_data' is a DataFrame with a column named 'Run_Number'
efficiency_data_4p9 = efficiency_data[efficiency_data['Run_Number'].isin(run_numbers)]    

# Update 'your_file.txt' with the actual file path
with open(REPLAYPATH+'/UTIL_BATCH/InputRunLists/KaonLT_2018_2019/Prod_3p8_Autumn18', 'r') as file:
    # Assuming each line in the file contains a single Run_Number
    run_numbers = [int(line.strip()) for line in file]

# Assuming 'efficiency_data' is a DataFrame with a column named 'Run_Number'
efficiency_data_3p8 = efficiency_data[efficiency_data['Run_Number'].isin(run_numbers)]    

################################################################################################################################################

# Define the linear fit function
def linear_fit(x, m, b):
    return m * x + b

# Error weighted fit of data
def fit_data(plt, x_name, y_name):

    print("Plotting {} vs {}...".format(x_name, y_name))
    
    y_error_name = y_name+"_ERROR"
    
    # Make x data
    efficiency_xdata_10p6 = efficiency_data_10p6[x_name].copy()
    efficiency_xdata_3p8 = efficiency_data_3p8[x_name].copy()
    efficiency_xdata_4p9 = efficiency_data_4p9[x_name].copy()
    efficiency_xdata_6p2 = efficiency_data_6p2[x_name].copy()
    efficiency_xdata_8p2 = efficiency_data_8p2[x_name].copy()

    # Concatenate x data from different sources
    x_data = pd.concat([efficiency_xdata_10p6, efficiency_xdata_3p8, efficiency_xdata_4p9, efficiency_xdata_6p2, efficiency_xdata_8p2], ignore_index=True)

    # Make y data
    efficiency_ydata_10p6 = efficiency_data_10p6[y_name].copy()
    efficiency_ydata_3p8 = efficiency_data_3p8[y_name].copy()
    efficiency_ydata_4p9 = efficiency_data_4p9[y_name].copy()
    efficiency_ydata_6p2 = efficiency_data_6p2[y_name].copy()
    efficiency_ydata_8p2 = efficiency_data_8p2[y_name].copy()

    # Concatenate y data from different sources
    y_data = pd.concat([efficiency_ydata_10p6, efficiency_ydata_3p8, efficiency_ydata_4p9, efficiency_ydata_6p2, efficiency_ydata_8p2], ignore_index=True)
    
    # Make y error
    efficiency_error_10p6 = efficiency_data_10p6[y_error_name].copy()
    efficiency_error_3p8 = efficiency_data_3p8[y_error_name].copy()
    efficiency_error_4p9 = efficiency_data_4p9[y_error_name].copy()
    efficiency_error_6p2 = efficiency_data_6p2[y_error_name].copy()
    efficiency_error_8p2 = efficiency_data_8p2[y_error_name].copy()

    # Concatenate y error from different sources
    y_error = pd.concat([efficiency_error_10p6, efficiency_error_3p8, efficiency_error_4p9, efficiency_error_6p2, efficiency_error_8p2], ignore_index=True)
    y_error = y_error + 1e-10 # Prevent divide by zero
    
    # Perform the error-weighted linear fit
    params, covariance = curve_fit(linear_fit, x_data, y_data, sigma=y_error, absolute_sigma=True)

    # Extract the slope and intercept from the fit
    slope = params[0]
    intercept = params[1]

    # Calculate the standard deviations of the parameters
    slope_error = np.sqrt(covariance[0, 0])
    intercept_error = np.sqrt(covariance[1, 1])

    # Calculate the fitted values and residuals
    y_fit = linear_fit(x_data, slope, intercept)
    residuals = y_data - y_fit

    # Calculate the chi-square value
    chi_square = np.sum((residuals / y_error)**2)
    
    # Generate x values for the error band
    x_fit = np.linspace(min(x_data), max(x_data), 100)

    # Calculate y values for the error band
    y_fit = linear_fit(x_fit, slope, intercept)

    # Calculate upper and lower bounds for the error band
    y_upper = linear_fit(x_fit, slope + slope_error, intercept + intercept_error)
    y_lower = linear_fit(x_fit, slope - slope_error, intercept - intercept_error)

    # Plot the data and the fitted line
    plt.errorbar(x_data, y_data, yerr=y_error, label=None,color='black',linestyle='None',zorder=3)
    plt.plot(x_fit, y_fit, label='m={0:.2e}±{1:.2e}\nb={2:.2e}±{3:.2e}\nchisq={4:.2e}'.format(slope, slope_error, intercept, intercept_error, chi_square), color='limegreen', linewidth=2, zorder=6)
    plt.fill_between(x_fit, y_lower, y_upper, color='lightgreen', alpha=0.4,zorder=5)
    
    # Annotate the plot with the slope and intercept
    plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)

################################################################################################################################################

plt.figure(figsize=(12,8))

plt.subplot(221)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
fit_data(plt, "SHMS_3/4_Trigger_Rate", "Non_Scaler_EDTM_Live_Time")
plt.scatter(efficiency_data_10p6["SHMS_3/4_Trigger_Rate"],efficiency_data_10p6["Non_Scaler_EDTM_Live_Time"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["SHMS_3/4_Trigger_Rate"],efficiency_data_3p8["Non_Scaler_EDTM_Live_Time"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["SHMS_3/4_Trigger_Rate"],efficiency_data_4p9["Non_Scaler_EDTM_Live_Time"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["SHMS_3/4_Trigger_Rate"],efficiency_data_6p2["Non_Scaler_EDTM_Live_Time"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["SHMS_3/4_Trigger_Rate"],efficiency_data_8p2["Non_Scaler_EDTM_Live_Time"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('EDTM', fontsize=12)
plt.xlabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(222)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
fit_data(plt, "SHMS_3/4_Trigger_Rate", "SHMS_Pion_ALL_TRACK_EFF")
plt.scatter(efficiency_data_10p6["SHMS_3/4_Trigger_Rate"],efficiency_data_10p6["SHMS_Pion_ALL_TRACK_EFF"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["SHMS_3/4_Trigger_Rate"],efficiency_data_3p8["SHMS_Pion_ALL_TRACK_EFF"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["SHMS_3/4_Trigger_Rate"],efficiency_data_4p9["SHMS_Pion_ALL_TRACK_EFF"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["SHMS_3/4_Trigger_Rate"],efficiency_data_6p2["SHMS_Pion_ALL_TRACK_EFF"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["SHMS_3/4_Trigger_Rate"],efficiency_data_8p2["SHMS_Pion_ALL_TRACK_EFF"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('SHMS_Pion_ALL_TRACK_EFF', fontsize=12)
plt.xlabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(223)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
fit_data(plt, "SHMS_3/4_Trigger_Rate", "SHMS_Aero_ALL_Pion_Eff")
plt.scatter(efficiency_data_10p6["SHMS_3/4_Trigger_Rate"],efficiency_data_10p6["SHMS_Aero_ALL_Pion_Eff"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["SHMS_3/4_Trigger_Rate"],efficiency_data_3p8["SHMS_Aero_ALL_Pion_Eff"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["SHMS_3/4_Trigger_Rate"],efficiency_data_4p9["SHMS_Aero_ALL_Pion_Eff"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["SHMS_3/4_Trigger_Rate"],efficiency_data_6p2["SHMS_Aero_ALL_Pion_Eff"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["SHMS_3/4_Trigger_Rate"],efficiency_data_8p2["SHMS_Aero_ALL_Pion_Eff"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('SHMS_Aero_ALL_Pion_Eff', fontsize=12)
plt.xlabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(224)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
#fit_data(plt, "SHMS_3/4_Trigger_Rate", "SHMS_Hodo_3_of_4_EFF")
plt.scatter(efficiency_data_10p6["SHMS_3/4_Trigger_Rate"],efficiency_data_10p6["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["SHMS_3/4_Trigger_Rate"],efficiency_data_3p8["SHMS_Hodo_3_of_4_EFF"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["SHMS_3/4_Trigger_Rate"],efficiency_data_4p9["SHMS_Hodo_3_of_4_EFF"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["SHMS_3/4_Trigger_Rate"],efficiency_data_6p2["SHMS_Hodo_3_of_4_EFF"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["SHMS_3/4_Trigger_Rate"],efficiency_data_8p2["SHMS_Hodo_3_of_4_EFF"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])   
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_3-4_%s.png' % (ROOTPrefix.replace("replay_","")))

plt.figure(figsize=(12,8))

plt.subplot(221)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
fit_data(plt, "SHMS_Hodoscope_S1X_Rate", "Non_Scaler_EDTM_Live_Time")
plt.scatter(efficiency_data_10p6["SHMS_Hodoscope_S1X_Rate"],efficiency_data_10p6["Non_Scaler_EDTM_Live_Time"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["SHMS_Hodoscope_S1X_Rate"],efficiency_data_3p8["Non_Scaler_EDTM_Live_Time"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["SHMS_Hodoscope_S1X_Rate"],efficiency_data_4p9["Non_Scaler_EDTM_Live_Time"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["SHMS_Hodoscope_S1X_Rate"],efficiency_data_6p2["Non_Scaler_EDTM_Live_Time"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["SHMS_Hodoscope_S1X_Rate"],efficiency_data_8p2["Non_Scaler_EDTM_Live_Time"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('EDTM', fontsize=12)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(222)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
fit_data(plt, "SHMS_Hodoscope_S1X_Rate", "SHMS_Pion_ALL_TRACK_EFF")
plt.scatter(efficiency_data_10p6["SHMS_Hodoscope_S1X_Rate"],efficiency_data_10p6["SHMS_Pion_ALL_TRACK_EFF"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["SHMS_Hodoscope_S1X_Rate"],efficiency_data_3p8["SHMS_Pion_ALL_TRACK_EFF"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["SHMS_Hodoscope_S1X_Rate"],efficiency_data_4p9["SHMS_Pion_ALL_TRACK_EFF"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["SHMS_Hodoscope_S1X_Rate"],efficiency_data_6p2["SHMS_Pion_ALL_TRACK_EFF"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["SHMS_Hodoscope_S1X_Rate"],efficiency_data_8p2["SHMS_Pion_ALL_TRACK_EFF"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('SHMS_Pion_ALL_TRACK_EFF', fontsize=12)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(223)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
fit_data(plt, "SHMS_Hodoscope_S1X_Rate", "SHMS_Aero_ALL_Pion_Eff")
plt.scatter(efficiency_data_10p6["SHMS_Hodoscope_S1X_Rate"],efficiency_data_10p6["SHMS_Aero_ALL_Pion_Eff"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["SHMS_Hodoscope_S1X_Rate"],efficiency_data_3p8["SHMS_Aero_ALL_Pion_Eff"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["SHMS_Hodoscope_S1X_Rate"],efficiency_data_4p9["SHMS_Aero_ALL_Pion_Eff"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["SHMS_Hodoscope_S1X_Rate"],efficiency_data_6p2["SHMS_Aero_ALL_Pion_Eff"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["SHMS_Hodoscope_S1X_Rate"],efficiency_data_8p2["SHMS_Aero_ALL_Pion_Eff"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('SHMS_Aero_ALL_Pion_Eff', fontsize=12)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(224)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
#fit_data(plt, "SHMS_Hodoscope_S1X_Rate", "SHMS_Hodo_3_of_4_EFF")
plt.scatter(efficiency_data_10p6["SHMS_Hodoscope_S1X_Rate"],efficiency_data_10p6["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["SHMS_Hodoscope_S1X_Rate"],efficiency_data_3p8["SHMS_Hodo_3_of_4_EFF"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["SHMS_Hodoscope_S1X_Rate"],efficiency_data_4p9["SHMS_Hodo_3_of_4_EFF"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["SHMS_Hodoscope_S1X_Rate"],efficiency_data_6p2["SHMS_Hodo_3_of_4_EFF"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["SHMS_Hodoscope_S1X_Rate"],efficiency_data_8p2["SHMS_Hodo_3_of_4_EFF"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('SHMS S1X HODO Rate [kHz]', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_S1X_%s.png' % (ROOTPrefix.replace("replay_","")))

plt.figure(figsize=(12,8))

plt.subplot(221)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
fit_data(plt, "HMS_EL-REAL_Trigger_Rate", "Non_Scaler_EDTM_Live_Time")
plt.scatter(efficiency_data_10p6["HMS_EL-REAL_Trigger_Rate"],efficiency_data_10p6["Non_Scaler_EDTM_Live_Time"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["HMS_EL-REAL_Trigger_Rate"],efficiency_data_3p8["Non_Scaler_EDTM_Live_Time"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["HMS_EL-REAL_Trigger_Rate"],efficiency_data_4p9["Non_Scaler_EDTM_Live_Time"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["HMS_EL-REAL_Trigger_Rate"],efficiency_data_6p2["Non_Scaler_EDTM_Live_Time"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["HMS_EL-REAL_Trigger_Rate"],efficiency_data_8p2["Non_Scaler_EDTM_Live_Time"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('EDTM', fontsize=12)
plt.xlabel('HMS EL-REAL Trigger Rate [kHz]', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(222)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
fit_data(plt, "HMS_EL-REAL_Trigger_Rate", "HMS_Elec_ALL_TRACK_EFF")
plt.scatter(efficiency_data_10p6["HMS_EL-REAL_Trigger_Rate"],efficiency_data_10p6["HMS_Elec_ALL_TRACK_EFF"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["HMS_EL-REAL_Trigger_Rate"],efficiency_data_3p8["HMS_Elec_ALL_TRACK_EFF"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["HMS_EL-REAL_Trigger_Rate"],efficiency_data_4p9["HMS_Elec_ALL_TRACK_EFF"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["HMS_EL-REAL_Trigger_Rate"],efficiency_data_6p2["HMS_Elec_ALL_TRACK_EFF"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["HMS_EL-REAL_Trigger_Rate"],efficiency_data_8p2["HMS_Elec_ALL_TRACK_EFF"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('HMS_Elec_ALL_TRACK_EFF', fontsize=12)
plt.xlabel('HMS EL-REAL Trigger Rate [kHz]', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(223)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
fit_data(plt, "HMS_EL-REAL_Trigger_Rate", "HMS_Cer_ALL_Elec_Eff")
plt.scatter(efficiency_data_10p6["HMS_EL-REAL_Trigger_Rate"],efficiency_data_10p6["HMS_Cer_ALL_Elec_Eff"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["HMS_EL-REAL_Trigger_Rate"],efficiency_data_3p8["HMS_Cer_ALL_Elec_Eff"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["HMS_EL-REAL_Trigger_Rate"],efficiency_data_4p9["HMS_Cer_ALL_Elec_Eff"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["HMS_EL-REAL_Trigger_Rate"],efficiency_data_6p2["HMS_Cer_ALL_Elec_Eff"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["HMS_EL-REAL_Trigger_Rate"],efficiency_data_8p2["HMS_Cer_ALL_Elec_Eff"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('HMS_Cer_ALL_Elec_Eff', fontsize=12)
plt.xlabel('HMS EL-REAL Trigger Rate [kHz]', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(224)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
#fit_data(plt, "HMS_EL-REAL_Trigger_Rate", "HMS_Hodo_3_of_4_EFF")
plt.scatter(efficiency_data_10p6["HMS_EL-REAL_Trigger_Rate"],efficiency_data_10p6["HMS_Hodo_3_of_4_EFF"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["HMS_EL-REAL_Trigger_Rate"],efficiency_data_3p8["HMS_Hodo_3_of_4_EFF"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["HMS_EL-REAL_Trigger_Rate"],efficiency_data_4p9["HMS_Hodo_3_of_4_EFF"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["HMS_EL-REAL_Trigger_Rate"],efficiency_data_6p2["HMS_Hodo_3_of_4_EFF"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["HMS_EL-REAL_Trigger_Rate"],efficiency_data_8p2["HMS_Hodo_3_of_4_EFF"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('HMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('HMS EL-REAL Trigger Rate [kHz]', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])   
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/HMS_EL-REAL_%s.png' % (ROOTPrefix.replace("replay_","")))

plt.figure(figsize=(12,8))

plt.subplot(221)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
fit_data(plt, "HMS_Hodoscope_S1X_Rate", "Non_Scaler_EDTM_Live_Time")
plt.scatter(efficiency_data_10p6["HMS_Hodoscope_S1X_Rate"],efficiency_data_10p6["Non_Scaler_EDTM_Live_Time"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["HMS_Hodoscope_S1X_Rate"],efficiency_data_3p8["Non_Scaler_EDTM_Live_Time"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["HMS_Hodoscope_S1X_Rate"],efficiency_data_4p9["Non_Scaler_EDTM_Live_Time"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["HMS_Hodoscope_S1X_Rate"],efficiency_data_6p2["Non_Scaler_EDTM_Live_Time"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["HMS_Hodoscope_S1X_Rate"],efficiency_data_8p2["Non_Scaler_EDTM_Live_Time"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('EDTM', fontsize=12)
plt.xlabel('HMS S1X HODO Rate [kHz]', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(222)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
fit_data(plt, "HMS_Hodoscope_S1X_Rate", "HMS_Elec_ALL_TRACK_EFF")
plt.scatter(efficiency_data_10p6["HMS_Hodoscope_S1X_Rate"],efficiency_data_10p6["HMS_Elec_ALL_TRACK_EFF"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["HMS_Hodoscope_S1X_Rate"],efficiency_data_3p8["HMS_Elec_ALL_TRACK_EFF"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["HMS_Hodoscope_S1X_Rate"],efficiency_data_4p9["HMS_Elec_ALL_TRACK_EFF"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["HMS_Hodoscope_S1X_Rate"],efficiency_data_6p2["HMS_Elec_ALL_TRACK_EFF"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["HMS_Hodoscope_S1X_Rate"],efficiency_data_8p2["HMS_Elec_ALL_TRACK_EFF"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('HMS_Elec_ALL_TRACK_EFF', fontsize=12)
plt.xlabel('HMS S1X HODO Rate [kHz]', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(223)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
fit_data(plt, "HMS_Hodoscope_S1X_Rate", "HMS_Cer_ALL_Elec_Eff")
plt.scatter(efficiency_data_10p6["HMS_Hodoscope_S1X_Rate"],efficiency_data_10p6["HMS_Cer_ALL_Elec_Eff"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["HMS_Hodoscope_S1X_Rate"],efficiency_data_3p8["HMS_Cer_ALL_Elec_Eff"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["HMS_Hodoscope_S1X_Rate"],efficiency_data_4p9["HMS_Cer_ALL_Elec_Eff"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["HMS_Hodoscope_S1X_Rate"],efficiency_data_6p2["HMS_Cer_ALL_Elec_Eff"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["HMS_Hodoscope_S1X_Rate"],efficiency_data_8p2["HMS_Cer_ALL_Elec_Eff"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('HMS_Cer_ALL_Elec_Eff', fontsize=12)
plt.xlabel('HMS S1X HODO Rate [kHz]', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(224)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
#fit_data(plt, "HMS_Hodoscope_S1X_Rate", "HMS_Hodo_3_of_4_EFF")
plt.scatter(efficiency_data_10p6["HMS_Hodoscope_S1X_Rate"],efficiency_data_10p6["HMS_Hodo_3_of_4_EFF"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["HMS_Hodoscope_S1X_Rate"],efficiency_data_3p8["HMS_Hodo_3_of_4_EFF"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["HMS_Hodoscope_S1X_Rate"],efficiency_data_4p9["HMS_Hodo_3_of_4_EFF"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["HMS_Hodoscope_S1X_Rate"],efficiency_data_6p2["HMS_Hodo_3_of_4_EFF"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["HMS_Hodoscope_S1X_Rate"],efficiency_data_8p2["HMS_Hodo_3_of_4_EFF"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('HMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('HMS S1X HODO Rate [kHz]', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/HMS_S1X_%s.png' % (ROOTPrefix.replace("replay_","")))

plt.figure(figsize=(12,8))
'''
plt.subplot(221)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data_10p6["COIN_Trigger_Rate"],efficiency_data_10p6["Non_Scaler_EDTM_Live_Time"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["COIN_Trigger_Rate"],efficiency_data_3p8["Non_Scaler_EDTM_Live_Time"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["COIN_Trigger_Rate"],efficiency_data_4p9["Non_Scaler_EDTM_Live_Time"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["COIN_Trigger_Rate"],efficiency_data_6p2["Non_Scaler_EDTM_Live_Time"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["COIN_Trigger_Rate"],efficiency_data_8p2["Non_Scaler_EDTM_Live_Time"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('EDTM', fontsize=12)
plt.xlabel('COIN Trigger Rate [kHz]', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(222)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data_10p6["COIN_Trigger_Rate"],efficiency_data_10p6["COIN_CPULT"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["COIN_Trigger_Rate"],efficiency_data_3p8["COIN_CPULT"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["COIN_Trigger_Rate"],efficiency_data_4p9["COIN_CPULT"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["COIN_Trigger_Rate"],efficiency_data_6p2["COIN_CPULT"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["COIN_Trigger_Rate"],efficiency_data_8p2["COIN_CPULT"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)

plt.ylabel('COIN CPULT', fontsize=12)
plt.xlabel('COIN Trigger Rate [kHz]', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(223)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data_10p6["COIN_Trigger_Rate"],efficiency_data_10p6["SHMS_3/4_Trigger_Rate"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["COIN_Trigger_Rate"],efficiency_data_3p8["SHMS_3/4_Trigger_Rate"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["COIN_Trigger_Rate"],efficiency_data_4p9["SHMS_3/4_Trigger_Rate"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["COIN_Trigger_Rate"],efficiency_data_6p2["SHMS_3/4_Trigger_Rate"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["COIN_Trigger_Rate"],efficiency_data_8p2["SHMS_3/4_Trigger_Rate"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
plt.xlabel('COIN Trigger Rate [kHz]', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(224)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data_10p6["COIN_Trigger_Rate"],efficiency_data_10p6["HMS_EL-REAL_Trigger_Rate"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["COIN_Trigger_Rate"],efficiency_data_3p8["HMS_EL-REAL_Trigger_Rate"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["COIN_Trigger_Rate"],efficiency_data_4p9["HMS_EL-REAL_Trigger_Rate"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["COIN_Trigger_Rate"],efficiency_data_6p2["HMS_EL-REAL_Trigger_Rate"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["COIN_Trigger_Rate"],efficiency_data_8p2["HMS_EL-REAL_Trigger_Rate"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('HMS EL-REAL Trigger Rate [kHz]', fontsize=12)
plt.xlabel('COIN Trigger Rate [kHz]', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/COIN_%s.png' % (ROOTPrefix.replace("replay_","")))
'''

plt.figure(figsize=(12,8))

plt.subplot(221)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data_10p6["Run_Number"],efficiency_data_10p6["Non_Scaler_EDTM_Live_Time"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["Run_Number"],efficiency_data_3p8["Non_Scaler_EDTM_Live_Time"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["Run_Number"],efficiency_data_4p9["Non_Scaler_EDTM_Live_Time"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["Run_Number"],efficiency_data_6p2["Non_Scaler_EDTM_Live_Time"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["Run_Number"],efficiency_data_8p2["Non_Scaler_EDTM_Live_Time"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('EDTM', fontsize=12)
plt.xlabel('Run Number', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(222)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data_10p6["Run_Number"],efficiency_data_10p6["SHMS_Pion_ALL_TRACK_EFF"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["Run_Number"],efficiency_data_3p8["SHMS_Pion_ALL_TRACK_EFF"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["Run_Number"],efficiency_data_4p9["SHMS_Pion_ALL_TRACK_EFF"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["Run_Number"],efficiency_data_6p2["SHMS_Pion_ALL_TRACK_EFF"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["Run_Number"],efficiency_data_8p2["SHMS_Pion_ALL_TRACK_EFF"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('SHMS_Pion_ALL_TRACK_EFF', fontsize=12)
plt.xlabel('Run Number', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(223)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data_10p6["Run_Number"],efficiency_data_10p6["SHMS_Aero_ALL_Pion_Eff"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["Run_Number"],efficiency_data_3p8["SHMS_Aero_ALL_Pion_Eff"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["Run_Number"],efficiency_data_4p9["SHMS_Aero_ALL_Pion_Eff"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["Run_Number"],efficiency_data_6p2["SHMS_Aero_ALL_Pion_Eff"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["Run_Number"],efficiency_data_8p2["SHMS_Aero_ALL_Pion_Eff"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('SHMS_Aero_ALL_Pion_Eff', fontsize=12)
plt.xlabel('Run Number', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(224)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data_10p6["Run_Number"],efficiency_data_10p6["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["Run_Number"],efficiency_data_3p8["SHMS_Hodo_3_of_4_EFF"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["Run_Number"],efficiency_data_4p9["SHMS_Hodo_3_of_4_EFF"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["Run_Number"],efficiency_data_6p2["SHMS_Hodo_3_of_4_EFF"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["Run_Number"],efficiency_data_8p2["SHMS_Hodo_3_of_4_EFF"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('Run Number', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])   
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMS_3-4_%s_run.png' % (ROOTPrefix.replace("replay_","")))

plt.figure(figsize=(12,8))

plt.subplot(221)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data_10p6["Run_Number"],efficiency_data_10p6["Non_Scaler_EDTM_Live_Time"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["Run_Number"],efficiency_data_3p8["Non_Scaler_EDTM_Live_Time"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["Run_Number"],efficiency_data_4p9["Non_Scaler_EDTM_Live_Time"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["Run_Number"],efficiency_data_6p2["Non_Scaler_EDTM_Live_Time"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["Run_Number"],efficiency_data_8p2["Non_Scaler_EDTM_Live_Time"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('EDTM', fontsize=12)
plt.xlabel('Run Number', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(222)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data_10p6["Run_Number"],efficiency_data_10p6["HMS_Elec_ALL_TRACK_EFF"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["Run_Number"],efficiency_data_3p8["HMS_Elec_ALL_TRACK_EFF"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["Run_Number"],efficiency_data_4p9["HMS_Elec_ALL_TRACK_EFF"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["Run_Number"],efficiency_data_6p2["HMS_Elec_ALL_TRACK_EFF"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["Run_Number"],efficiency_data_8p2["HMS_Elec_ALL_TRACK_EFF"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('HMS_Elec_ALL_TRACK_EFF', fontsize=12)
plt.xlabel('Run Number', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(223)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data_10p6["Run_Number"],efficiency_data_10p6["HMS_Cer_ALL_Elec_Eff"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["Run_Number"],efficiency_data_3p8["HMS_Cer_ALL_Elec_Eff"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["Run_Number"],efficiency_data_4p9["HMS_Cer_ALL_Elec_Eff"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["Run_Number"],efficiency_data_6p2["HMS_Cer_ALL_Elec_Eff"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["Run_Number"],efficiency_data_8p2["HMS_Cer_ALL_Elec_Eff"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('HMS_Cer_ALL_Elec_Eff', fontsize=12)
plt.xlabel('Run Number', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(224)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data_10p6["Run_Number"],efficiency_data_10p6["HMS_Hodo_3_of_4_EFF"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["Run_Number"],efficiency_data_3p8["HMS_Hodo_3_of_4_EFF"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["Run_Number"],efficiency_data_4p9["HMS_Hodo_3_of_4_EFF"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["Run_Number"],efficiency_data_6p2["HMS_Hodo_3_of_4_EFF"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["Run_Number"],efficiency_data_8p2["HMS_Hodo_3_of_4_EFF"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('HMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('Run Number', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])   
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/HMS_EL-REAL_%s_run.png' % (ROOTPrefix.replace("replay_","")))

plt.figure(figsize=(12,8))

plt.subplot(221)    
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data_10p6["Run_Number"],efficiency_data_10p6["Non_Scaler_EDTM_Live_Time"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["Run_Number"],efficiency_data_3p8["Non_Scaler_EDTM_Live_Time"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["Run_Number"],efficiency_data_4p9["Non_Scaler_EDTM_Live_Time"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["Run_Number"],efficiency_data_6p2["Non_Scaler_EDTM_Live_Time"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["Run_Number"],efficiency_data_8p2["Non_Scaler_EDTM_Live_Time"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('EDTM', fontsize=12)
plt.xlabel('Run Number', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(222)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data_10p6["Run_Number"],efficiency_data_10p6["COIN_CPULT"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["Run_Number"],efficiency_data_3p8["COIN_CPULT"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["Run_Number"],efficiency_data_4p9["COIN_CPULT"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["Run_Number"],efficiency_data_6p2["COIN_CPULT"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["Run_Number"],efficiency_data_8p2["COIN_CPULT"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.subplot(223)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data_10p6["Run_Number"],efficiency_data_10p6["SHMS_3/4_Trigger_Rate"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["Run_Number"],efficiency_data_3p8["SHMS_3/4_Trigger_Rate"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["Run_Number"],efficiency_data_4p9["SHMS_3/4_Trigger_Rate"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["Run_Number"],efficiency_data_6p2["SHMS_3/4_Trigger_Rate"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["Run_Number"],efficiency_data_8p2["SHMS_3/4_Trigger_Rate"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('SHMS 3/4 Trigger Rate [kHz]', fontsize=12)
plt.xlabel('Run Number', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(224)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data_10p6["Run_Number"],efficiency_data_10p6["HMS_EL-REAL_Trigger_Rate"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["Run_Number"],efficiency_data_3p8["HMS_EL-REAL_Trigger_Rate"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["Run_Number"],efficiency_data_4p9["HMS_EL-REAL_Trigger_Rate"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["Run_Number"],efficiency_data_6p2["HMS_EL-REAL_Trigger_Rate"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["Run_Number"],efficiency_data_8p2["HMS_EL-REAL_Trigger_Rate"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('HMS EL-REAL Trigger Rate [kHz]', fontsize=12)
plt.xlabel('Run Number', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/COIN_%s_run.png' % (ROOTPrefix.replace("replay_","")))

plt.figure(figsize=(12,8))

plt.subplot(221)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
#fit_data(plt, "BCM1_Beam_Cut_Current", "HMS_Hodo_3_of_4_EFF")
plt.scatter(efficiency_data_10p6["BCM1_Beam_Cut_Current"],efficiency_data_10p6["HMS_Hodo_3_of_4_EFF"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["BCM1_Beam_Cut_Current"],efficiency_data_3p8["HMS_Hodo_3_of_4_EFF"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["BCM1_Beam_Cut_Current"],efficiency_data_4p9["HMS_Hodo_3_of_4_EFF"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["BCM1_Beam_Cut_Current"],efficiency_data_6p2["HMS_Hodo_3_of_4_EFF"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["BCM1_Beam_Cut_Current"],efficiency_data_8p2["HMS_Hodo_3_of_4_EFF"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('HMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('Current', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(222)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
fit_data(plt, "BCM1_Beam_Cut_Current", "HMS_Cal_ALL_Elec_Eff")
plt.scatter(efficiency_data_10p6["BCM1_Beam_Cut_Current"],efficiency_data_10p6["HMS_Cal_ALL_Elec_Eff"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["BCM1_Beam_Cut_Current"],efficiency_data_3p8["HMS_Cal_ALL_Elec_Eff"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["BCM1_Beam_Cut_Current"],efficiency_data_4p9["HMS_Cal_ALL_Elec_Eff"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["BCM1_Beam_Cut_Current"],efficiency_data_6p2["HMS_Cal_ALL_Elec_Eff"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["BCM1_Beam_Cut_Current"],efficiency_data_8p2["HMS_Cal_ALL_Elec_Eff"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('HMS_Cal_ALL_Elec_Eff', fontsize=12)
plt.xlabel('Current', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(223)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
fit_data(plt, "BCM1_Beam_Cut_Current", "HMS_Cer_ALL_Elec_Eff")
plt.scatter(efficiency_data_10p6["BCM1_Beam_Cut_Current"],efficiency_data_10p6["HMS_Cer_ALL_Elec_Eff"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["BCM1_Beam_Cut_Current"],efficiency_data_3p8["HMS_Cer_ALL_Elec_Eff"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["BCM1_Beam_Cut_Current"],efficiency_data_4p9["HMS_Cer_ALL_Elec_Eff"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["BCM1_Beam_Cut_Current"],efficiency_data_6p2["HMS_Cer_ALL_Elec_Eff"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["BCM1_Beam_Cut_Current"],efficiency_data_8p2["HMS_Cer_ALL_Elec_Eff"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('HMS_Cer_ALL_Elec_Eff', fontsize=12)
plt.xlabel('Current', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/HMSDet_%s_run.png' % (ROOTPrefix.replace("replay_","")))

plt.figure(figsize=(12,8))

plt.subplot(221)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
#fit_data(plt, "BCM1_Beam_Cut_Current", "SHMS_Hodo_3_of_4_EFF")
plt.scatter(efficiency_data_10p6["BCM1_Beam_Cut_Current"],efficiency_data_10p6["SHMS_Hodo_3_of_4_EFF"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["BCM1_Beam_Cut_Current"],efficiency_data_3p8["SHMS_Hodo_3_of_4_EFF"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["BCM1_Beam_Cut_Current"],efficiency_data_4p9["SHMS_Hodo_3_of_4_EFF"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["BCM1_Beam_Cut_Current"],efficiency_data_6p2["SHMS_Hodo_3_of_4_EFF"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["BCM1_Beam_Cut_Current"],efficiency_data_8p2["SHMS_Hodo_3_of_4_EFF"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('SHMS_Hodo_3_of_4_EFF', fontsize=12)
plt.xlabel('Current', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(222)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
fit_data(plt, "BCM1_Beam_Cut_Current", "SHMS_Aero_ALL_Pion_Eff")
plt.scatter(efficiency_data_10p6["BCM1_Beam_Cut_Current"],efficiency_data_10p6["SHMS_Aero_ALL_Pion_Eff"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["BCM1_Beam_Cut_Current"],efficiency_data_3p8["SHMS_Aero_ALL_Pion_Eff"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["BCM1_Beam_Cut_Current"],efficiency_data_4p9["SHMS_Aero_ALL_Pion_Eff"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["BCM1_Beam_Cut_Current"],efficiency_data_6p2["SHMS_Aero_ALL_Pion_Eff"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["BCM1_Beam_Cut_Current"],efficiency_data_8p2["SHMS_Aero_ALL_Pion_Eff"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('SHMS_Aero_ALL_Pion_Eff', fontsize=12)
plt.xlabel('Current', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/SHMSDet_%s_run.png' % (ROOTPrefix.replace("replay_","")))

plt.figure(figsize=(12,8))

plt.subplot(121)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
fit_data(plt, "BCM1_Beam_Cut_Current", "HMS_Elec_ALL_TRACK_EFF")
plt.scatter(efficiency_data_10p6["BCM1_Beam_Cut_Current"],efficiency_data_10p6["HMS_Elec_ALL_TRACK_EFF"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["BCM1_Beam_Cut_Current"],efficiency_data_3p8["HMS_Elec_ALL_TRACK_EFF"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["BCM1_Beam_Cut_Current"],efficiency_data_4p9["HMS_Elec_ALL_TRACK_EFF"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["BCM1_Beam_Cut_Current"],efficiency_data_6p2["HMS_Elec_ALL_TRACK_EFF"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["BCM1_Beam_Cut_Current"],efficiency_data_8p2["HMS_Elec_ALL_TRACK_EFF"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('HMS_Elec_ALL_TRACK_EFF', fontsize=12)
plt.xlabel('Current', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(122)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
fit_data(plt, "BCM1_Beam_Cut_Current", "SHMS_Pion_ALL_TRACK_EFF")
plt.scatter(efficiency_data_10p6["BCM1_Beam_Cut_Current"],efficiency_data_10p6["SHMS_Pion_ALL_TRACK_EFF"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["BCM1_Beam_Cut_Current"],efficiency_data_3p8["SHMS_Pion_ALL_TRACK_EFF"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["BCM1_Beam_Cut_Current"],efficiency_data_4p9["SHMS_Pion_ALL_TRACK_EFF"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["BCM1_Beam_Cut_Current"],efficiency_data_6p2["SHMS_Pion_ALL_TRACK_EFF"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["BCM1_Beam_Cut_Current"],efficiency_data_8p2["SHMS_Pion_ALL_TRACK_EFF"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('SHMS_Pion_ALL_TRACK_EFF', fontsize=12)
plt.xlabel('Current', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/Track_%s_run.png' % (ROOTPrefix.replace("replay_","")))

plt.figure(figsize=(12,8))

plt.subplot(221)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data_10p6["Run_Number"],efficiency_data_10p6["Non_Scaler_EDTM_Live_Time"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["Run_Number"],efficiency_data_3p8["Non_Scaler_EDTM_Live_Time"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["Run_Number"],efficiency_data_4p9["Non_Scaler_EDTM_Live_Time"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["Run_Number"],efficiency_data_6p2["Non_Scaler_EDTM_Live_Time"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["Run_Number"],efficiency_data_8p2["Non_Scaler_EDTM_Live_Time"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('EDTM', fontsize=12)
plt.xlabel('Run Number', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(222)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
plt.scatter(efficiency_data_10p6["Run_Number"],efficiency_data_10p6["BOIL_Eff"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["Run_Number"],efficiency_data_3p8["BOIL_Eff"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["Run_Number"],efficiency_data_4p9["BOIL_Eff"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["Run_Number"],efficiency_data_6p2["BOIL_Eff"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["Run_Number"],efficiency_data_8p2["BOIL_Eff"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('BOIL_Eff', fontsize=12)
plt.xlabel('Run Number', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(223)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
fit_data(plt, "BCM1_Beam_Cut_Current", "Non_Scaler_EDTM_Live_Time")
plt.scatter(efficiency_data_10p6["BCM1_Beam_Cut_Current"],efficiency_data_10p6["Non_Scaler_EDTM_Live_Time"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["BCM1_Beam_Cut_Current"],efficiency_data_3p8["Non_Scaler_EDTM_Live_Time"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["BCM1_Beam_Cut_Current"],efficiency_data_4p9["Non_Scaler_EDTM_Live_Time"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["BCM1_Beam_Cut_Current"],efficiency_data_6p2["Non_Scaler_EDTM_Live_Time"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["BCM1_Beam_Cut_Current"],efficiency_data_8p2["Non_Scaler_EDTM_Live_Time"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('EDTM', fontsize=12)
plt.xlabel('Current', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.subplot(224)
plt.grid(zorder=1)
#plt.xlim(0,100)
#plt.ylim(0.9,1.1)
fit_data(plt, "BCM1_Beam_Cut_Current", "BOIL_Eff")
plt.scatter(efficiency_data_10p6["BCM1_Beam_Cut_Current"],efficiency_data_10p6["BOIL_Eff"],color='blue',zorder=4,label='10p6')
plt.scatter(efficiency_data_3p8["BCM1_Beam_Cut_Current"],efficiency_data_3p8["BOIL_Eff"],color='red',zorder=4,label='3p8')
plt.scatter(efficiency_data_4p9["BCM1_Beam_Cut_Current"],efficiency_data_4p9["BOIL_Eff"],color='purple',zorder=4,label='4p9')
plt.scatter(efficiency_data_6p2["BCM1_Beam_Cut_Current"],efficiency_data_6p2["BOIL_Eff"],color='orange',zorder=4,label='6p2')
plt.scatter(efficiency_data_8p2["BCM1_Beam_Cut_Current"],efficiency_data_8p2["BOIL_Eff"],color='pink',zorder=4,label='8p2')
plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1, fontsize=10)
plt.ylabel('Boiling Correction', fontsize=12)
plt.xlabel('Current', fontsize=12)
#plt.title('%s-%s' % (int(min(efficiency_data["Run_Number"])),int(max(efficiency_data["Run_Number"]))), fontsize=12)

plt.tight_layout(rect=[0,0.03,1,0.95])
plt.savefig(UTILPATH+'/scripts/efficiency/OUTPUTS/plots/EDTM_%s_run.png' % (ROOTPrefix.replace("replay_","")))

#plt.show()
