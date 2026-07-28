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
from ROOT import TCanvas, TList, TPaveLabel, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TLegend, TGaxis, TLine, TMath, TLatex, TPaveText, TArc, TGraphPolar, TText, TString
from ROOT import kBlack, kCyan, kRed, kGreen, kMagenta, kBlue
from functools import reduce
from array import array
import math as ma
import uncertainties as u
from ctypes import c_double

##################################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=1:
    print("!!!!! ERROR !!!!!\n Expected 2 arguments\n Usage is with - PHY_SETTING DATA_Suffix\n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################

# Input params - run number and max number of events
PHY_SETTING = sys.argv[1]
#DATA_Suffix = sys.argv[2]
DATA_Suffix = "-1_ProdCoin_Analysed_Data"

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

#################################################################################################################################################

# Output PDF File Name
print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))

# Input file location and variables taking
rootFile_DATA = "%s/%s_%s.root" % (OUTPATH, PHY_SETTING, DATA_Suffix)
mmcut_csv_file = "%s/%s/%s_mm_offsets_cuts_parameters.csv" % (MMCUT_CSV, physet_dir_name, setting_name)
dcut_csv_file = "%s/%s/%s_diamond_cut_parameters.csv" % (DCUT_CSV, physet_dir_name, setting_name)

###################################################################################################################################################

print ('\nPhysics Setting = ',PHY_SETTING, '\n')
print("-"*40)

# Cuts for Pions Selection - Change according to your data
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
cutg_diamond = ROOT.TCutG("cutg_diamond", 5)
cutg_diamond.SetVarX("Q2")
cutg_diamond.SetVarY("W")
cutg_diamond.SetPoint(0, vertex1[0], vertex1[1])  # bottom-left
cutg_diamond.SetPoint(1, vertex2[0], vertex2[1])  # top-left
cutg_diamond.SetPoint(2, vertex3[0], vertex3[1])  # top-right
cutg_diamond.SetPoint(3, vertex4[0], vertex4[1])  # bottom-right
cutg_diamond.SetPoint(4, vertex1[0], vertex1[1])  # bottom-left again to close the loop
Diamond_Cut = lambda event: (cutg_diamond.IsInside(event.Q2, event.W))

#------------------------------------------------------------------------------------------------------------------------------------

# Read the MMpi cut values from the CSV file
try:
    with open(mmcut_csv_file, mode='r', newline='') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            if row[csv_reader.fieldnames[0]].strip() == PHY_SETTING:  # Match first column (Physics_Setting)
                MM_Offset = float(row["MM_Offset"].strip())
                tshift_Offset = float(row["tshift_Offset"].strip())
                MM_Cut_lowvalue = float(row["MM_Cut_low"].strip())
                MM_Cut_highvalue = float(row["MM_Cut_high"].strip())
                break
        else:
            raise ValueError(f"No matching Physics_Setting '{PHY_SETTING}' found in {mmcut_csv_file}.")

except (FileNotFoundError, ValueError, KeyError) as e:
    print(f"Error: {e}")
    sys.exit(1)

# Print the assigned values
print(f"MMpi_Offset = {MM_Offset:.6f}")
print(f"tshift_Offset = {tshift_Offset:.6f}")
print(f"MMpi_Cut_lowvalue = {MM_Cut_lowvalue}")
print(f"MMpi_Cut_highvalue = {MM_Cut_highvalue}")
DATA_MMpi_Cut = lambda event: (MM_Cut_lowvalue <= (event.MMpi + MM_Offset) <= MM_Cut_highvalue)
print("-"*40)
# CoinTime nWindows
nWindows = 6

##########################################################################################################################################################################################################

# tmin Calculation
m2 = 0.9382720813
m3 = 0.13957039
m4 = 0.9395654133

def calc_tmin(q2, w, m2, m3, m4):
    m22 = m2 * m2
    m32 = m3 * m3
    m42 = m4 * m4

    s = w * w
    m12 = -q2

    omega = (s + q2 - m22) / (2.0 * m2)
    q = math.sqrt(max(q2 + omega**2, 0.0))

    e1cm = (s + m12 - m22) / (2.0 * w)
    e3cm = (s + m32 - m42) / (2.0 * w)

    p1lab = q
    p1cm = p1lab * m2 / w
    p3cm = math.sqrt(max(e3cm * e3cm - m32, 0.0))
    tmin = -((e1cm - e3cm)**2 - (p1cm - p3cm)**2)
    return tmin

###################################################################################################################################################

# Read stuff from the main event tree
infile_DATA = ROOT.TFile.Open(rootFile_DATA, "READ")

#Uncut_Pion_Events_Data_tree = infile_DATA.Get("Uncut_Pion_Events")
#Cut_Pion_Events_Accpt_Data_tree = infile_DATA.Get("Cut_Pion_Events_Accpt")
#Cut_Pion_Events_All_Data_tree = infile_DATA.Get("Cut_Pion_Events_All")
Cut_Pion_Events_Prompt_Data_tree = infile_DATA.Get("Cut_Pion_Events_Prompt")
Cut_Pion_Events_Random_Data_tree = infile_DATA.Get("Cut_Pion_Events_Random")
#nEntries_TBRANCH_DATA  = Cut_Pion_Events_Prompt_Data_tree.GetEntries()

# Making directories in output file
outHistFile = ROOT.TFile.Open("%s/tcheck/%s_tcheck_Data.root" % (OUTPATH, PHY_SETTING) , "RECREATE")

# Create output trees
Cut_Pion_Events_Prompt = ROOT.TTree("Cut_Pion_Events_Prompt", "Tree with Q2, W, MandelT, and calculated tmin")
Cut_Pion_Events_Random = ROOT.TTree("Cut_Pion_Events_Random", "Tree with Q2, W, MandelT, and calculated tmin")

# Prompt branches
Q2_prompt      = array('d', [0.0])
W_prompt       = array('d', [0.0])
MandelT_prompt = array('d', [0.0])
tmin_prompt    = array('d', [0.0])

# Random branches
Q2_random      = array('d', [0.0])
W_random       = array('d', [0.0])
MandelT_random = array('d', [0.0])
tmin_random    = array('d', [0.0])

# Prompt tree branches
Cut_Pion_Events_Prompt.Branch("Q2", Q2_prompt, "Q2/D")
Cut_Pion_Events_Prompt.Branch("W", W_prompt, "W/D")
Cut_Pion_Events_Prompt.Branch("MandelT", MandelT_prompt, "MandelT/D")
Cut_Pion_Events_Prompt.Branch("tmin", tmin_prompt, "tmin/D")

# Random tree branches
Cut_Pion_Events_Random.Branch("Q2", Q2_random, "Q2/D")
Cut_Pion_Events_Random.Branch("W", W_random, "W/D")
Cut_Pion_Events_Random.Branch("MandelT", MandelT_random, "MandelT/D")
Cut_Pion_Events_Random.Branch("tmin", tmin_random, "tmin/D")

# Fill prompt tree
for event in Cut_Pion_Events_Prompt_Data_tree:
    if DATA_MMpi_Cut(event) and Diamond_Cut(event):
        q2 = float(event.Q2)
        w = float(event.W)
        t = float(event.MandelT)
        Q2_prompt[0] = q2
        W_prompt[0] = w
        MandelT_prompt[0] = t
        tmin_prompt[0] = calc_tmin(q2, w, m2, m3, m4)
        Cut_Pion_Events_Prompt.Fill()

# Fill random tree
for event in Cut_Pion_Events_Random_Data_tree:
    if DATA_MMpi_Cut(event) and Diamond_Cut(event):
        q2 = float(event.Q2)
        w = float(event.W)
        t = float(event.MandelT)
        Q2_random[0] = q2
        W_random[0] = w
        MandelT_random[0] = t
        tmin_random[0] = calc_tmin(q2, w, m2, m3, m4)
        Cut_Pion_Events_Random.Fill()

# Write trees to file
outHistFile.cd()
Cut_Pion_Events_Prompt.Write()
Cut_Pion_Events_Random.Write()

###################################################################################################################################################

nbins = 200
t_min = 0.1
t_max = 0.7
tm_min = 0.1
tm_max = 0.4
# Define histograms
tmin_vs_t_pions_data_prompt_notshift_cut_all = ROOT.TH2D("tmin_vs_t_pions_data_prompt_notshift_cut_all", "tmin vs t Distribution (no tshift); tmin; t", nbins, tm_min, tm_max, nbins, t_min, t_max)
tmin_vs_t_pions_data_random_notshift_cut_all = ROOT.TH2D("tmin_vs_t_pions_data_random_notshift_cut_all", "tmin vs t Distribution (no tshift); tmin; t", nbins, tm_min, tm_max, nbins, t_min, t_max)
tmin_vs_t_pions_data_randsub_notshift_cut_all = ROOT.TH2D("tmin_vs_t_pions_data_randsub_notshift_cut_all", "tmin vs t Distribution (no tshift); tmin; t", nbins, tm_min, tm_max, nbins, t_min, t_max)

tmin_vs_t_pions_data_prompt_tshift_cut_all = ROOT.TH2D("tmin_vs_t_pions_data_prompt_tshift_cut_all", "tmin vs t Distribution (tshift); tmin; t", nbins, tm_min, tm_max, nbins, t_min, t_max)
tmin_vs_t_pions_data_random_tshift_cut_all = ROOT.TH2D("tmin_vs_t_pions_data_random_tshift_cut_all", "tmin vs t Distribution (tshift); tmin; t", nbins, tm_min, tm_max, nbins, t_min, t_max)
tmin_vs_t_pions_data_randsub_tshift_cut_all = ROOT.TH2D("tmin_vs_t_pions_data_randsub_tshift_cut_all", "tmin vs t Distribution (tshift); tmin; t", nbins, tm_min, tm_max, nbins, t_min, t_max)

# Fill prompt and random directly from trees
Cut_Pion_Events_Prompt.Draw("(-MandelT):tmin>>tmin_vs_t_pions_data_prompt_notshift_cut_all", "", "goff")
Cut_Pion_Events_Random.Draw("(-MandelT):tmin>>tmin_vs_t_pions_data_random_notshift_cut_all", "", "goff")
Cut_Pion_Events_Prompt.Draw(f"((-MandelT)+{tshift_Offset}):(tmin)>>tmin_vs_t_pions_data_prompt_tshift_cut_all","","goff")
Cut_Pion_Events_Random.Draw(f"((-MandelT)+{tshift_Offset}):(tmin)>>tmin_vs_t_pions_data_random_tshift_cut_all","","goff")

# Normalized random subtraction
tmin_vs_t_pions_data_random_notshift_cut_all.Scale(1.0/nWindows)
tmin_vs_t_pions_data_random_tshift_cut_all.Scale(1.0/nWindows)

# Random subtraction
tmin_vs_t_pions_data_randsub_notshift_cut_all.Add(tmin_vs_t_pions_data_prompt_notshift_cut_all, tmin_vs_t_pions_data_random_notshift_cut_all, 1.0, -1.0)
tmin_vs_t_pions_data_randsub_tshift_cut_all.Add(tmin_vs_t_pions_data_prompt_tshift_cut_all, tmin_vs_t_pions_data_random_tshift_cut_all, 1.0, -1.0)

###################################################################################################################################################

# saving in PDF format
Pion_Analysis_Distributions = "%s/tcheck/%s_ProdCoin_Pion_Analysis_Ratio_Comparison_Distributions.pdf" % (OUTPATH, PHY_SETTING)

c_tmin_compare = ROOT.TCanvas("c_tmin_compare", "Random Subtracted tmin vs t", 1600, 700)
c_tmin_compare.Divide(2, 1)

c_tmin_compare.cd(1)
ROOT.gPad.SetLogz()
#ROOT.gPad.SetRightMargin(0.15)
#tmin_vs_t_pions_data_randsub_notshift_cut_all.SetTitle("Random Subtracted (no tshift); -tmin; -t")
tmin_vs_t_pions_data_randsub_notshift_cut_all.Draw("COLZ")
line_t_eq_tmin = ROOT.TLine(t_min, t_min, t_max, t_max)
line_t_eq_tmin.SetLineColor(ROOT.kRed)
line_t_eq_tmin.SetLineWidth(2)
line_t_eq_tmin.SetLineStyle(2)
line_t_eq_tmin.Draw("SAME")

c_tmin_compare.cd(2)
ROOT.gPad.SetLogz()
#ROOT.gPad.SetRightMargin(0.15)
#tmin_vs_t_pions_data_randsub_tshift_cut_all.SetTitle("Random Subtracted (tshift); -tmin; -t")
tmin_vs_t_pions_data_randsub_tshift_cut_all.Draw("COLZ")
line_t_eq_tmin.Draw("SAME")
c_tmin_compare.Print(Pion_Analysis_Distributions)

###################################################################################################################################################

print("####################################")
print("###### Histogram writing done ######")
print("####################################\n")

infile_DATA.Close() 
outHistFile.Close()

print ("Processing Complete")