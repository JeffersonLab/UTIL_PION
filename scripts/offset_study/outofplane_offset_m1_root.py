#! /usr/bin/pythoa
#
# Description:
# ================================================================
# Time-stamp: "2024-03-15 01:29:19 junaid"
# ================================================================
#
# Author:  Muhammad Junaid III <mjo147@uregina.ca>
#
# Copyright (c) junaid
#
###################################################################################################################################################

# Import relevant packages
import uproot as up
import numpy as np

np.bool = bool
np.float = float

import root_numpy as rnp
import pandas as pd
import root_pandas as rpd
import ROOT
import scipy
from scipy.optimize import curve_fit
import scipy.integrate as integrate
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from uncertainties import ufloat
import sys, math, os, subprocess
import array
import csv
import math
from ROOT import TCanvas, TPaveLabel, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TLegend, TGaxis, TLine, TMath, TLatex, TPaveText, TArc, TGraphPolar, TText
from ROOT import kBlack, kCyan, kRed, kGreen, kMagenta, kBlue
from functools import reduce

##################################################################################################################################################

'''
ltsep package import and pathing definitions
'''

# Import package for cuts
from ltsep import Root

lt=Root(os.path.realpath(__file__), "Plot_HeePCoin")

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH

#################################################################################################################################################

# Output PDF File Name
print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))
print("="*40)

# Defining Variables
E5p986_Pmy_simc = ufloat(0.00011174, 0.00002425)
E5p986_Pmy_data = ufloat(-0.00530417, 0.00007652)
E5p986_HMS_Pe = 3.271
E5p986_SHMS_Pp = 3.493
E5p986_y = (E5p986_Pmy_simc.nominal_value - E5p986_Pmy_data.nominal_value)/E5p986_HMS_Pe
E5p986_x = (E5p986_SHMS_Pp/E5p986_HMS_Pe)
E5p986_yerr = math.sqrt((E5p986_Pmy_simc.std_dev / E5p986_HMS_Pe)**2 + (E5p986_Pmy_data.std_dev / E5p986_HMS_Pe)**2 + (0.00376 / E5p986_HMS_Pe)**2)
E5p986_y = "{:.5f}".format(E5p986_y)
E5p986_yerr = "{:.5f}".format(E5p986_yerr)
E5p986_x = "{:.5f}".format(E5p986_x)
print('\nE5p986_y =', E5p986_y, '\n')
print('E5p986_yerr =', E5p986_yerr, '\n')
print('E5p986_x =', E5p986_x, '\n')
print("="*40)

E6p395s1_Pmy_simc = ufloat(0.00031749, 0.00003540)
E6p395s1_Pmy_data = ufloat(-0.01096507, 0.00005202)
E6p395s1_HMS_Pe = 4.752
E6p395s1_SHMS_Pp = 2.412
E6p395s1_y = (E6p395s1_Pmy_simc.nominal_value - E6p395s1_Pmy_data.nominal_value)/E6p395s1_HMS_Pe
E6p395s1_x = (E6p395s1_SHMS_Pp/E6p395s1_HMS_Pe)
E6p395s1_yerr = math.sqrt((E6p395s1_Pmy_simc.std_dev / E6p395s1_HMS_Pe)**2 + (E6p395s1_Pmy_data.std_dev / E6p395s1_HMS_Pe)**2 + (0.00376 / E6p395s1_HMS_Pe)**2)
E6p395s1_y = "{:.5f}".format(E6p395s1_y)
E6p395s1_yerr = "{:.5f}".format(E6p395s1_yerr)
E6p395s1_x = "{:.5f}".format(E6p395s1_x)
print('\nE6p395s1_y =', E6p395s1_y, '\n')
print('E6p395s1_yerr =', E6p395s1_yerr, '\n')
print('E6p395s1_x =', E6p395s1_x, '\n')
print("="*40)

E6p395s2_Pmy_simc = ufloat(0.00017365, 0.00003205)
E6p395s2_Pmy_data = ufloat(-0.00942472, 0.00007130)
E6p395s2_HMS_Pe = 4.391
E6p395s2_SHMS_Pp = 2.792
E6p395s2_y = (E6p395s2_Pmy_simc.nominal_value - E6p395s2_Pmy_data.nominal_value)/E6p395s2_HMS_Pe
E6p395s2_x = (E6p395s2_SHMS_Pp/E6p395s2_HMS_Pe)
E6p395s2_yerr = math.sqrt((E6p395s2_Pmy_simc.std_dev / E6p395s2_HMS_Pe)**2 + (E6p395s2_Pmy_data.std_dev / E6p395s2_HMS_Pe)**2 + (0.00376 / E6p395s2_HMS_Pe)**2)
E6p395s2_y = "{:.5f}".format(E6p395s2_y)
E6p395s2_yerr = "{:.5f}".format(E6p395s2_yerr)
E6p395s2_x = "{:.5f}".format(E6p395s2_x)
print('\nE6p395s2_y =', E6p395s2_y, '\n')
print('E6p395s2_yerr =', E6p395s2_yerr, '\n')
print('E6p395s2_x =', E6p395s2_x, '\n')
print("="*40)

E6p395s3_Pmy_simc = ufloat(0.00044029, 0.00002564)
E6p395s3_Pmy_data = ufloat(-0.00718003, 0.00011148)
E6p395s3_HMS_Pe = 3.014
E6p395s3_SHMS_Pp = 4.220
E6p395s3_y = (E6p395s3_Pmy_simc.nominal_value - E6p395s3_Pmy_data.nominal_value)/E6p395s3_HMS_Pe
E6p395s3_x = (E6p395s3_SHMS_Pp/E6p395s3_HMS_Pe)
E6p395s3_yerr = math.sqrt((E6p395s3_Pmy_simc.std_dev / E6p395s3_HMS_Pe)**2 + (E6p395s3_Pmy_data.std_dev / E6p395s3_HMS_Pe)**2 + (0.00376 / E6p395s3_HMS_Pe)**2)
E6p395s3_y = "{:.5f}".format(E6p395s3_y)
E6p395s3_yerr = "{:.5f}".format(E6p395s3_yerr)
E6p395s3_x = "{:.5f}".format(E6p395s3_x)
print('\nE6p395s3_y =', E6p395s3_y, '\n')
print('E6p395s3_yerr =', E6p395s3_yerr, '\n')
print('E6p395s3_x =', E6p395s3_x, '\n')
print("="*40)

E7p937_Pmy_simc = ufloat(0.00051255, 0.00002849)
E7p937_Pmy_data = ufloat(-0.00738244, 0.00017781)
E7p937_HMS_Pe = 3.283
E7p937_SHMS_Pp = 5.512
E7p937_y = (E7p937_Pmy_simc.nominal_value - E7p937_Pmy_data.nominal_value)/E7p937_HMS_Pe
E7p937_x = (E7p937_SHMS_Pp/E7p937_HMS_Pe)
E7p937_yerr = math.sqrt((E7p937_Pmy_simc.std_dev / E7p937_HMS_Pe)**2 + (E7p937_Pmy_data.std_dev / E7p937_HMS_Pe)**2 + (0.00376 / E7p937_HMS_Pe)**2)
E7p937_y = "{:.5f}".format(E7p937_y)
E7p937_yerr = "{:.5f}".format(E7p937_yerr)
E7p937_x = "{:.5f}".format(E7p937_x)
print('\nE7p937_y =', E7p937_y, '\n')
print('E7p937_yerr =', E7p937_yerr, '\n')
print('E7p937_x =', E7p937_x, '\n')
print("="*40)

E8p479_Pmy_simc = ufloat(0.00026799, 0.00003908)
E8p479_Pmy_data = ufloat(-0.00944413, 0.00013390)
E8p479_HMS_Pe = 5.587
E8p479_SHMS_Pp = 3.731
E8p479_y = (E8p479_Pmy_simc.nominal_value - E8p479_Pmy_data.nominal_value)/E8p479_HMS_Pe
E8p479_x = (E8p479_SHMS_Pp/E8p479_HMS_Pe)
E8p479_yerr = math.sqrt((E8p479_Pmy_simc.std_dev / E8p479_HMS_Pe)**2 + (E8p479_Pmy_data.std_dev / E8p479_HMS_Pe)**2 + (0.00376 / E8p479_HMS_Pe)**2)
E8p479_y = "{:.5f}".format(E8p479_y)
E8p479_yerr = "{:.5f}".format(E8p479_yerr)
E8p479_x = "{:.5f}".format(E8p479_x)
print('\nE8p479_y =', E8p479_y, '\n')
print('E8p479_yerr =', E8p479_yerr, '\n')
print('E8p479_x =', E8p479_x, '\n')
print("="*40)

E9p177_Pmy_simc = ufloat(0.00062312, 0.00003286)
E9p177_Pmy_data = ufloat(-0.00726819, 0.00027508)
E9p177_HMS_Pe = 3.738
E9p177_SHMS_Pp = 6.265
E9p177_y = (E9p177_Pmy_simc.nominal_value - E9p177_Pmy_data.nominal_value)/E9p177_HMS_Pe
E9p177_x = (E9p177_SHMS_Pp/E9p177_HMS_Pe)
E9p177_yerr = math.sqrt((E9p177_Pmy_simc.std_dev / E9p177_HMS_Pe)**2 + (E9p177_Pmy_data.std_dev / E9p177_HMS_Pe)**2 + (0.00376 / E9p177_HMS_Pe)**2)
E9p177_y = "{:.5f}".format(E9p177_y)
E9p177_yerr = "{:.5f}".format(E9p177_yerr)
E9p177_x = "{:.5f}".format(E9p177_x)
print('\nE9p177_y =', E9p177_y, '\n')
print('E9p177_yerr =', E9p177_yerr, '\n')
print('E9p177_x =', E9p177_x, '\n')
print("="*40)

E9p876_Pmy_simc = ufloat(0.00041417, 0.00004005)
E9p876_Pmy_data = ufloat(-0.01088800, 0.00017774)
E9p876_HMS_Pe = 5.366
E9p876_SHMS_Pp = 5.422
E9p876_y = (E9p876_Pmy_simc.nominal_value - E9p876_Pmy_data.nominal_value)/E9p876_HMS_Pe
E9p876_x = (E9p876_SHMS_Pp/E9p876_HMS_Pe)
E9p876_yerr = math.sqrt((E9p876_Pmy_simc.std_dev / E9p876_HMS_Pe)**2 + (E9p876_Pmy_data.std_dev / E9p876_HMS_Pe)**2 + (0.00376 / E9p876_HMS_Pe)**2)
E9p876_y = "{:.5f}".format(E9p876_y)
E9p876_yerr = "{:.5f}".format(E9p876_yerr)
E9p876_x = "{:.5f}".format(E9p876_x)
print('\nE9p876_y =', E9p876_y, '\n')
print('E9p876_yerr =', E9p876_yerr, '\n')
print('E9p876_x =', E9p876_x, '\n')
print("="*40)

E10p549_Pmy_simc = ufloat(0.00035501, 0.00004392)
E10p549_Pmy_data = ufloat(-0.00971681, 0.00017990)
E10p549_HMS_Pe = 5.878
E10p549_SHMS_Pp = 5.530
E10p549_y = (E10p549_Pmy_simc.nominal_value - E10p549_Pmy_data.nominal_value)/E10p549_HMS_Pe
E10p549_x = (E10p549_SHMS_Pp/E10p549_HMS_Pe)
E10p549_yerr = math.sqrt((E10p549_Pmy_simc.std_dev / E10p549_HMS_Pe)**2 + (E10p549_Pmy_data.std_dev / E10p549_HMS_Pe)**2 + (0.00376 / E10p549_HMS_Pe)**2)
E10p549_y = "{:.5f}".format(E10p549_y)
E10p549_yerr = "{:.5f}".format(E10p549_yerr)
E10p549_x = "{:.5f}".format(E10p549_x)
print('\nE10p549_y =', E10p549_y, '\n')
print('E10p549_yerr =', E10p549_yerr, '\n')
print('E10p549_x =', E10p549_x, '\n')
print("="*40)

#################################################################################################################################################

# Start batch mode
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# Create lists for x and y values along with their uncertainties
x_values = [E10p549_x, E5p986_x, E6p395s1_x, E6p395s2_x, E6p395s3_x, E7p937_x, E8p479_x, E9p177_x, E9p876_x]
y_values = [E10p549_y, E5p986_y, E6p395s1_y, E6p395s2_y, E6p395s3_y, E7p937_y, E8p479_y, E9p177_y, E9p876_y]
y_errors = [E10p549_yerr, E5p986_yerr, E6p395s1_yerr, E6p395s2_yerr, E6p395s3_yerr, E7p937_yerr, E8p479_yerr, E9p177_yerr, E9p876_yerr]

# Convert values to floats
x_min = float(min(x_values))
x_max = float(max(x_values))

# Create a TGraphErrors object
graph = ROOT.TGraphErrors(len(x_values))

# Fill the TGraphErrors with x, y values, and errors
for i in range(len(x_values)):
    graph.SetPoint(i, float(x_values[i]), float(y_values[i]))  # Ensure float conversion here
    graph.SetPointError(i, 0, float(y_errors[i]))  # Ensure float conversion here as well

# Define a fitting function
fit_func = ROOT.TF1("fit_func", "[0]*x + [1]", x_min, x_max)

# Set initial parameter values if needed
fit_func.SetParameters(1.0, 0.0)

# Perform the fit
result = graph.Fit(fit_func, "S")  # "S" option includes error on Y values

# Print fit result
result.Print()

# Access fit parameters
fit_parameters = result.GetParams()
fit_errors = result.GetErrors()
print("Fit parameters:", fit_parameters)

# Create a canvas
canvas = ROOT.TCanvas("canvas", "Graph with Error Bars", 1000, 600)

# Set titles and axis labels
graph.SetTitle("Out of Plane Offset")
graph.GetXaxis().SetTitle("$(P_p/P_{e^\prime})$")
graph.GetYaxis().SetTitle("$(PMY_\mathrm{SIMC} - PMY_\mathrm{DATA})/P_{e^\prime}$")

# Set different marker styles and colors for data points
graph.SetMarkerStyle(20)  # Circle
graph.SetMarkerColor(ROOT.kBlue)  # Blue
graph.GetXaxis().SetRangeUser(0.0,2.25)
graph.GetYaxis().SetRangeUser(0.001,0.003)

# Draw the graph
graph.Draw("AP")
fit_func.Draw("same")

# Add fit results to the plot
latex = ROOT.TLatex()
latex.SetTextSize(0.03)

# Write slope and intercept values
slope = fit_parameters[0]
intercept = fit_parameters[1]
slope_error = fit_errors[0]
intercept_error = fit_errors[1]

text_slope = "Slope: {:.2f}".format(slope*1000)
text_intercept = "Intercept: {:.2f}".format(intercept*1000)

latex.DrawLatexNDC(0.2, 0.85, text_slope)
latex.DrawLatexNDC(0.2, 0.80, text_intercept)

# Save the graph as PNG
canvas.SaveAs("Out_of_plane_offset_m1_root.png")

############################################################################################################################################

print ("Processing Complete")
