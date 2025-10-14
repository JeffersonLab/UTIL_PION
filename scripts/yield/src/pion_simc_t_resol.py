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

# Import relevant packages
import uproot
import uproot as up
import numpy as np

np.bool = bool
np.float = float

import root_numpy as rnp
import pandas as pd
import root_pandas as rpd
import ROOT
import scipy
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import sys, math, os, subprocess
import array
import csv
from ROOT import TCanvas, TList, TPaveLabel, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TLegend, TGaxis, TLine, TMath, TLatex, TPaveText, TArc, TGraphPolar, TText, TString, TCutG
from ROOT import kBlack, kCyan, kRed, kGreen, kMagenta, kBlue
from functools import reduce
import math as ma

##################################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=2:
    print("!!!!! ERROR !!!!!\n Expected 2 arguments\n Usage is with - Physics and SIMC Setting\n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################

# Input params - run number and max number of events
PHY_SETTING = sys.argv[1]
SIMC_SETTING = sys.argv[2]
#SIMC_Suffix_lowepsright = "Prod_Coin_{SIMC_SETTING}_lowepsright"
SIMC_Suffix_lowepscenter = "Prod_Coin_{}_lowepscenter".format(SIMC_SETTING)
SIMC_Suffix_lowepsleft = "Prod_Coin_{}_lowepsleft".format(SIMC_SETTING)
SIMC_Suffix_highepsright = "Prod_Coin_{}_highepsright".format(SIMC_SETTING)
SIMC_Suffix_highepscenter = "Prod_Coin_{}_highepscenter".format(SIMC_SETTING)
SIMC_Suffix_highepsleft = "Prod_Coin_{}_highepsleft".format(SIMC_SETTING)

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
from ltsep import Root

lt=Root(os.path.realpath(__file__), "Plot_ProdCoin")

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH
MMCUT_CSV   = "/u/group/c-pionlt/USERS/%s/hallc_replay_lt/UTIL_PION/LTSep_CSVs/mm_offset_cut_csv" % (USER)
DCUT_CSV    = "/u/group/c-pionlt/USERS/%s/hallc_replay_lt/UTIL_PION/LTSep_CSVs/diamond_cut_csv" % (USER)

# Extract the first three words from PHY_SETTING for the CSV file name
setting_name = "_".join(PHY_SETTING.split("_")[:3])
physet_dir_name = "%s_std" % (setting_name)
SIMCPATH = "/volatile/hallc/c-pionlt/%s/worksim/%s/" % (USER, physet_dir_name)

#################################################################################################################################################

# Output PDF File Name
print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))
Pion_Analysis_Distributions = "%s/%s_SIMC_Pion_Analysis_tresolution_Distributions.pdf" % (OUTPATH, PHY_SETTING)

# Extract the first three words from PHY_SETTING for the CSV file name
setting_name = "_".join(PHY_SETTING.split("_")[:3])
mmcut_csv_file = "%s/%s/%s_mm_offsets_cuts_parameters.csv" % (MMCUT_CSV, physet_dir_name, setting_name)
dcut_csv_file = "%s/%s/%s_diamond_cut_parameters.csv" % (DCUT_CSV, physet_dir_name, setting_name)

# Input file location and variables taking
#rootFile_SIMC_lowepsright = "%s/%s.root" % (SIMCPATH, SIMC_Suffix_lowepsright)
rootFile_SIMC_lowepscenter = "%s/%s.root" % (SIMCPATH, SIMC_Suffix_lowepscenter)
rootFile_SIMC_lowepsleft = "%s/%s.root" % (SIMCPATH, SIMC_Suffix_lowepsleft)
rootFile_SIMC_highepsright = "%s/%s.root" % (SIMCPATH, SIMC_Suffix_highepsright)
rootFile_SIMC_highepscenter = "%s/%s.root" % (SIMCPATH, SIMC_Suffix_highepscenter)
rootFile_SIMC_highepsleft = "%s/%s.root" % (SIMCPATH, SIMC_Suffix_highepsleft)

###################################################################################################################################################

# Cuts for SIMC
# Read the vertices from the CSV file
vertices = {}
try:
    with open(dcut_csv_file, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            vertex_name = row["Vertex"]
            x_value = round(float(row["X"]), 3)  # Round to 3 decimal places
            y_value = round(float(row["Y"]), 3)  # Round to 3 decimal places
            vertices[vertex_name] = [x_value, y_value]
except (FileNotFoundError, ValueError, KeyError) as e:
    print(f"Error reading diamond cut vertices from {dcut_csv_file}: {e}")
    sys.exit(1)

# Assign the vertices
vertex1 = vertices["vertex1"]  # bottom-left
vertex2 = vertices["vertex2"]  # top-left
vertex3 = vertices["vertex3"]  # top-right
vertex4 = vertices["vertex4"]  # bottom-right

# Print the vertex values rounded to 3 decimals
#print("Diamond Cut Vertices (rounded to 3 decimals):")
#print(f"Vertex 1 (bottom-left): [{vertex1[0]:.3f}, {vertex1[1]:.3f}]")
#print(f"Vertex 2 (top-left): [{vertex2[0]:.3f}, {vertex2[1]:.3f}]")
#print(f"Vertex 3 (top-right): [{vertex3[0]:.3f}, {vertex3[1]:.3f}]")
#print(f"Vertex 4 (bottom-right): [{vertex4[0]:.3f}, {vertex4[1]:.3f}]")

# Define the diamond cut
cutg_diamond_simc = ROOT.TCutG("cutg_diamond_simc", 5)
cutg_diamond_simc.SetVarX("Q2")
cutg_diamond_simc.SetVarY("W")
cutg_diamond_simc.SetPoint(0, vertex1[0], vertex1[1])  # bottom-left
cutg_diamond_simc.SetPoint(1, vertex2[0], vertex2[1])  # top-left
cutg_diamond_simc.SetPoint(2, vertex3[0], vertex3[1])  # top-right
cutg_diamond_simc.SetPoint(3, vertex4[0], vertex4[1])  # bottom-right
cutg_diamond_simc.SetPoint(4, vertex1[0], vertex1[1])  # bottom-left again to close the loop

# SIMC Cuts for Pions Selection
HMS_Acceptance = lambda event: (event.hsdelta >= -8.0) & (event.hsdelta <= 8.0) & (event.hsxpfp >= -0.08) & (event.hsxpfp <= 0.08) & (event.hsypfp >= -0.045) & (event.hsypfp <= 0.045)
SHMS_Acceptance = lambda event: (event.ssdelta >= -10.0) & (event.ssdelta <= 20.0) & (event.ssxpfp >= -0.06) & (event.ssxpfp <= 0.06) & (event.ssypfp >= -0.04) & (event.ssypfp <= 0.04)

# Read the MMpi cut values from the CSV file
try:
    with open(mmcut_csv_file, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        first_row = next(csv_reader)  # Get the first row
        if first_row:
            MM_Cut_lowvalue = float(first_row["MM_Cut_low"])  # Assign MM_Cut_low
            MM_Cut_highvalue = float(first_row["MM_Cut_high"])  # Assign MM_Cut_high
        else:
            raise ValueError(f"The CSV file {mmcut_csv_file} is empty or missing required columns.")
except (FileNotFoundError, ValueError, KeyError) as e:
    print(f"Error: {e}")
    sys.exit(1)

# Print the assigned values
#print(f"MMpi_Cut_lowvalue = {MMpi_Cut_lowvalue}")
#print(f"MMpi_Cut_highvalue = {MMpi_Cut_highvalue}")

# Missing Mass Cut
SIMC_MMpi_Cut = lambda event: (MM_Cut_lowvalue <= event.missmass <= MM_Cut_highvalue)
Diamond_Cut = lambda event: (cutg_diamond_simc.IsInside(event.Q2, event.W))

###############################################################################################################################################

# Grabs simc number of events and normalizaton factor
# Function to extract simc number of events and normalization factor
def extract_simc_factors(simc_suffix):
    simc_hist = "%s/%s.hist" % (SIMCPATH, simc_suffix)
    try:
        with open(simc_hist) as f_simc:
            for line in f_simc:
                if "Ncontribute" in line:
                    val = line.split("=")
                    simc_nevents = int(val[1])
                if "normfac" in line:
                    val = line.split("=")
                    simc_normfactor = float(val[1])
        if 'simc_nevents' in locals() and 'simc_normfactor' in locals():
            print('\nsimc_nevents = ', simc_nevents, '\nsimc_normfactor = ', simc_normfactor, '\n')
            return simc_nevents, simc_normfactor
        else:
            print("ERROR: Invalid simc hist file %s" % simc_hist)
            sys.exit(1)
    except FileNotFoundError:
        print("ERROR: File not found %s" % simc_hist)
        sys.exit(1)

# Extract factors for all SIMC files
simc_suffixes = {
    "lowepscenter": SIMC_Suffix_lowepscenter,
    "lowepsleft": SIMC_Suffix_lowepsleft,
    "highepsright": SIMC_Suffix_highepsright,
    "highepscenter": SIMC_Suffix_highepscenter,
    "highepsleft": SIMC_Suffix_highepsleft
}

normfac_simc_dict = {}

for key, suffix in simc_suffixes.items():
    simc_nevents, simc_normfactor = extract_simc_factors(suffix)
    normfac_simc = simc_normfactor / simc_nevents
    normfac_simc_dict[key] = normfac_simc
    print("normfac_simc for %s: " % key, normfac_simc)
    print("-" * 40)

###################################################################################################################################################
nbins = 400
min = -0.2
max = 0.2

# Defining Histograms for Pions
# Histograms for SIMC
t_ti_resol_pions_simc_lowepscenter_cut_all = ROOT.TH1D("t_ti_resol_pions_simc_lowepscenter_cut_all", "t - ti Distribution loweps_center; t - ti; Counts", nbins, min, max)
t_ti_resol_pions_simc_lowepsleft_cut_all = ROOT.TH1D("t_ti_resol_pions_simc_lowepsleft_cut_all", "t - ti Distribution loweps_left; t - ti; Counts", nbins, min, max)
t_ti_resol_pions_simc_highepsright_cut_all = ROOT.TH1D("t_ti_resol_pions_simc_highepsright_cut_all", "t - ti Distribution higheps_right; t - ti; Counts", nbins, min, max)
t_ti_resol_pions_simc_highepscenter_cut_all = ROOT.TH1D("t_ti_resol_pions_simc_highepscenter_cut_all", "t - ti Distribution higheps_center; t - ti; Counts", nbins, min, max)
t_ti_resol_pions_simc_highepsleft_cut_all = ROOT.TH1D("t_ti_resol_pions_simc_highepsleft_cut_all", "t - ti Distribution higheps_left; t - ti; Counts", nbins, min, max)

t_ti_resol_pions_simc_loweps_cut_all = ROOT.TH1D("t_ti_resol_pions_simc_loweps_cut_all", "t - ti Distribution loweps; t - ti; Counts", nbins, min, max)
t_ti_resol_pions_simc_higheps_cut_all = ROOT.TH1D("t_ti_resol_pions_simc_higheps_cut_all", "t - ti Distribution higheps; t - ti; Counts", nbins, min, max)

##########################################################################################################################################################################################################

# Read stuff from the main event tree
#infile_SIMC_lowepsright = ROOT.TFile.Open(rootFile_SIMC_lowepsright, "READ")
infile_SIMC_lowepscenter = ROOT.TFile.Open(rootFile_SIMC_lowepscenter, "READ")
infile_SIMC_lowepsleft = ROOT.TFile.Open(rootFile_SIMC_lowepsleft, "READ")
infile_SIMC_highepsright = ROOT.TFile.Open(rootFile_SIMC_highepsright, "READ")
infile_SIMC_highepscenter = ROOT.TFile.Open(rootFile_SIMC_highepscenter, "READ")
infile_SIMC_highepsleft = ROOT.TFile.Open(rootFile_SIMC_highepsleft, "READ")

Uncut_Pion_Events_SIMC_lowepscenter_tree = infile_SIMC_lowepscenter.Get("h10")
Uncut_Pion_Events_SIMC_lowepsleft_tree = infile_SIMC_lowepsleft.Get("h10")
Uncut_Pion_Events_SIMC_highepsright_tree = infile_SIMC_highepsright.Get("h10")
Uncut_Pion_Events_SIMC_highepscenter_tree = infile_SIMC_highepscenter.Get("h10")
Uncut_Pion_Events_SIMC_highepsleft_tree = infile_SIMC_highepsleft.Get("h10")

###################################################################################################################################################


# Fill histograms from SIMC ROOT File
for event in Uncut_Pion_Events_SIMC_lowepscenter_tree:
    # Define the acceptance cuts
    if HMS_Acceptance(event) & SHMS_Acceptance(event) & SIMC_MMpi_Cut(event) & Diamond_Cut(event):        
        t_ti_resol_pions_simc_lowepscenter_cut_all.Fill(event.t - event.ti, event.Weight)

for event in Uncut_Pion_Events_SIMC_lowepsleft_tree:
    # Define the acceptance cuts
    if HMS_Acceptance(event) & SHMS_Acceptance(event) & SIMC_MMpi_Cut(event) & Diamond_Cut(event):        
        t_ti_resol_pions_simc_lowepsleft_cut_all.Fill(event.t - event.ti, event.Weight)

for event in Uncut_Pion_Events_SIMC_highepsright_tree:
    # Define the acceptance cuts
    if HMS_Acceptance(event) & SHMS_Acceptance(event) & SIMC_MMpi_Cut(event) & Diamond_Cut(event):        
        t_ti_resol_pions_simc_highepsright_cut_all.Fill(event.t - event.ti, event.Weight)

for event in Uncut_Pion_Events_SIMC_highepscenter_tree:
    # Define the acceptance cuts
    if HMS_Acceptance(event) & SHMS_Acceptance(event) & SIMC_MMpi_Cut(event) & Diamond_Cut(event):        
        t_ti_resol_pions_simc_highepscenter_cut_all.Fill(event.t - event.ti, event.Weight)

for event in Uncut_Pion_Events_SIMC_highepsleft_tree:
    # Define the acceptance cuts
    if HMS_Acceptance(event) & SHMS_Acceptance(event) & SIMC_MMpi_Cut(event) & Diamond_Cut(event):        
        t_ti_resol_pions_simc_highepsleft_cut_all.Fill(event.t - event.ti, event.Weight)

print("####################################")
print("###### Histogram filling done ######")
print("####################################\n")

############################################################################################################################################

# SIMC Normalization
normfac_simc_lowepscenter = normfac_simc_dict["lowepscenter"]
normfac_simc_lowepsleft = normfac_simc_dict["lowepsleft"]
normfac_simc_highepsright = normfac_simc_dict["highepsright"]
normfac_simc_highepscenter = normfac_simc_dict["highepscenter"]
normfac_simc_highepsleft = normfac_simc_dict["highepsleft"]

# Normalizing histograms
t_ti_resol_pions_simc_lowepscenter_cut_all.Scale(normfac_simc_lowepscenter)
t_ti_resol_pions_simc_lowepsleft_cut_all.Scale(normfac_simc_lowepsleft)
t_ti_resol_pions_simc_highepsright_cut_all.Scale(normfac_simc_highepsright)
t_ti_resol_pions_simc_highepscenter_cut_all.Scale(normfac_simc_highepscenter)
t_ti_resol_pions_simc_highepsleft_cut_all.Scale(normfac_simc_highepsleft)

#############################################################################################################################################

# Resolution Histograms
# Adding left right and center settings for both low and high epsilon
t_ti_resol_pions_simc_loweps_cut_all.Add(t_ti_resol_pions_simc_lowepscenter_cut_all, t_ti_resol_pions_simc_lowepsleft_cut_all, 1, 1)
t_ti_resol_pions_simc_higheps_cut_all.Add(t_ti_resol_pions_simc_highepsright_cut_all, t_ti_resol_pions_simc_highepscenter_cut_all, 1, 1)
t_ti_resol_pions_simc_higheps_cut_all.Add(t_ti_resol_pions_simc_highepsleft_cut_all, 1)

#############################################################################################################################################

# Removes stat box
ROOT.gStyle.SetOptStat(0)

# Fit ranges for the t-t_i resolution histograms
fit_range_low = -0.008
fit_range_high = 0.008

# Saving histograms in PDF
c1_delta = TCanvas("c1_delta", "Variables Distributions", 100, 0, 1400, 1400)
c1_delta.Divide(1,2)
# Add PHY_SETTING to the left of the plots
phy_setting_text = ROOT.TLatex()
phy_setting_text.SetTextSize(0.05)
phy_setting_text.SetTextAlign(13)
phy_setting_text.SetNDC()
c1_delta.cd(1)
t_ti_resol_pions_simc_loweps_cut_all.Draw("hist")
t_ti_resol_pions_simc_loweps_cut_all.Fit("gaus", "Q", "", fit_range_low, fit_range_high)
fit_t_ti_loweps = t_ti_resol_pions_simc_loweps_cut_all.GetFunction('gaus')
fit_t_ti_loweps.SetLineColor(ROOT.kRed)
fit_t_ti_loweps.Draw("same")
loweps_sigma = ROOT.TLatex()
loweps_sigma.SetTextSize(0.05)
loweps_sigma.SetTextAlign(13)
loweps_sigma.SetNDC()
loweps_sigma.DrawLatex(0.6, 0.5, f'#color[2]{{#sigma}} = {fit_t_ti_loweps.GetParameter(2):.5f}#pm{fit_t_ti_loweps.GetParError(2):.5f}')
phy_setting_text.DrawLatex(0.13, 0.7, f"#color[1]{{Setting: {PHY_SETTING}}}")
loweps_t_ti_resol = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
loweps_t_ti_resol.AddEntry(t_ti_resol_pions_simc_loweps_cut_all, "low #epsilon", "l")
loweps_t_ti_resol.Draw()
c1_delta.cd(2)
t_ti_resol_pions_simc_higheps_cut_all.Draw("hist")
t_ti_resol_pions_simc_higheps_cut_all.Fit("gaus", "Q", "", fit_range_low, fit_range_high)
fit_t_ti_higheps = t_ti_resol_pions_simc_higheps_cut_all.GetFunction('gaus')
fit_t_ti_higheps.Draw("same")
higheps_sigma = ROOT.TLatex()
higheps_sigma.SetTextSize(0.05)
higheps_sigma.SetTextAlign(13)
higheps_sigma.SetNDC()
higheps_sigma.DrawLatex(0.6, 0.5, f'#color[2]{{#sigma}} = {fit_t_ti_higheps.GetParameter(2):.5f}#pm{fit_t_ti_higheps.GetParError(2):.5f}')
phy_setting_text.DrawLatex(0.13, 0.7, f"#color[1]{{Setting: {PHY_SETTING}}}")
higheps_t_ti_resol = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
higheps_t_ti_resol.AddEntry(t_ti_resol_pions_simc_higheps_cut_all, "high #epsilon", "l")
higheps_t_ti_resol.Draw()
c1_delta.Print(Pion_Analysis_Distributions)

#############################################################################################################################################

# Writing the vertices to the CSV file
csv_output_path = "%s/LTSep_CSVs/t_resolution_csv/%s/%s_simc_t_resolution_parameters.csv" % (UTILPATH, physet_dir_name, setting_name)

# Extract the fit parameters for low epsilon
loweps_value = fit_t_ti_loweps.GetParameter(2)
loweps_error = fit_t_ti_loweps.GetParError(2)

# Extract the fit parameters for high epsilon
higheps_value = fit_t_ti_higheps.GetParameter(2)
higheps_error = fit_t_ti_higheps.GetParError(2)

# Write the values to the CSV file
with open(csv_output_path, mode='w', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)
    # Write the header
    csv_writer.writerow(["Setting", "Value", "Error"])
    # Write the data for low epsilon
    csv_writer.writerow(["loweps", f"{loweps_value:.5f}", f"{loweps_error:.5f}"])
    # Write the data for high epsilon
    csv_writer.writerow(["higheps", f"{higheps_value:.5f}", f"{higheps_error:.5f}"])

print(f"t-resolution data saved to {csv_output_path}")

#############################################################################################################################################

# Making directories in output file
outHistFile = ROOT.TFile.Open("%s/%s_tresolution_Data.root" % (OUTPATH, PHY_SETTING) , "RECREATE")
d_Cut_Pion_Events_SIMC_Cut_All = outHistFile.mkdir("Cut_Pion_Events_SIMC_Cut_All")

# Writing Histograms for pions
d_Cut_Pion_Events_SIMC_Cut_All.cd()
t_ti_resol_pions_simc_lowepscenter_cut_all.Write()
t_ti_resol_pions_simc_lowepsleft_cut_all.Write()
t_ti_resol_pions_simc_highepsright_cut_all.Write()
t_ti_resol_pions_simc_highepscenter_cut_all.Write()
t_ti_resol_pions_simc_highepsleft_cut_all.Write()
t_ti_resol_pions_simc_loweps_cut_all.Write()
t_ti_resol_pions_simc_higheps_cut_all.Write()

print("####################################")
print("###### Histogram writing done ######")
print("####################################\n")

infile_SIMC_lowepscenter.Close()
infile_SIMC_lowepsleft.Close()
infile_SIMC_highepsright.Close()
infile_SIMC_highepscenter.Close()
infile_SIMC_highepsleft.Close()
outHistFile.Close()

print ("Processing Complete")