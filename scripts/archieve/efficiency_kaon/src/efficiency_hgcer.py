#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-12-18 13:10:47 trottar"
# ================================================================
# 
# Author:  Richard L. Trotta III <trotta@cua.edu>
# 
# Copyright (c) trottar
#

# 19/10/20 - Stephen Kay, University of Regina

# Import relevant packages
import uproot as up
import numpy as np
import ROOT
from ROOT import TCanvas, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TGaxis, TLine, TMath, TPaveText, TArc, TGraphPolar 
from ROOT import kBlack, kBlue, kRed
import os

# Input should be the input root file name (including suffix) and an output file name string (without any suffix)
def dictionary(UTILPATH,runNum,MaxEvent):

    ROOTPrefix = "efficiency"
    OutFilename = "%s_%s_%s" % (ROOTPrefix,runNum,MaxEvent)

    ################################################################################################################################################
    '''
    ltsep package import and pathing definitions
    '''

    # Import package for cuts
    from ltsep import Root

    lt=Root(os.path.realpath(__file__),"Plot_Prod_HGCer")

    # Add this to all files for more dynamic pathing
    USER=lt.USER # Grab user info for file finding
    HOST=lt.HOST
    REPLAYPATH=lt.REPLAYPATH
    UTILPATH=lt.UTILPATH
    ANATYPE=lt.ANATYPE
    OUTPATH=lt.OUTPATH

    ################################################################################################################################################

    rootName = "%s/%s_%s_%s.root" % (OUTPATH, runNum, MaxEvent, ROOTPrefix)

    #################################################################################################################################################

    InFile = ROOT.TFile.Open(rootName, "READ")
    TOutFilename = OutFilename
    # Establish the names of our output files quickly
    foutname = OUTPATH+"/" + TOutFilename + ".root"
    foutpdf = OUTPATH+"/" + TOutFilename + ".pdf"

    #Events_no_cal_hgc_aero_cuts  = InFile.Get("SHMS_cut_no_Cal_HGC_Aero")  
    Events_no_cal_hgc_aero_cuts  = InFile.Get("SHMS_Kaons_Without_HGC_Cuts")
    nEntries_Events_no_cal_hgc_aero_cuts  = Events_no_cal_hgc_aero_cuts.GetEntries()

    # Particles information no cuts
    SHMS_Events = InFile.Get("SHMS_Events")    
    nEntries_SHMS_Events = SHMS_Events.GetEntries()

    #################################################################################################################################################
    ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
    #ROOT.gStyle.SetOptStat(0) # Set style to hide stat box
    #ROOT.gROOT.ForceStyle() # Force style change
    #################################################################################################################################################

    # Defined Geomatrical cuts
    cutg = ROOT.TCutG("cutg",21)
    cutg.SetVarX("P_hgcer_yAtCer")
    cutg.SetVarY("P_hgcer_xAtCer")

    if runNum > 0 and runNum < 10000:
        cutg.SetPoint(0,-25,2)
        cutg.SetPoint(1,-2,2)
        cutg.SetPoint(2,-1,2.5)
        cutg.SetPoint(3,0,3)
        cutg.SetPoint(4,1,3)
        cutg.SetPoint(5,2,3.3)
        cutg.SetPoint(6,3,3.0)
        cutg.SetPoint(7,4,2.5)
        cutg.SetPoint(8,5,2)
        cutg.SetPoint(9,25,2)
        cutg.SetPoint(10,25,0.5)
        cutg.SetPoint(11,5,0.5)
        cutg.SetPoint(12,4,1)
        cutg.SetPoint(13,3,-1)
        cutg.SetPoint(14,2,-2)
        cutg.SetPoint(15,1,-2.3)
        cutg.SetPoint(16,0,-1.5)
        cutg.SetPoint(17,-1,-1)
        cutg.SetPoint(18,-2,0.5)
        cutg.SetPoint(19,-25,0.5)
        cutg.SetPoint(20,-25,2)

    elif runNum > 10000 and runNum < 20000:
        cutg.SetPoint(0,-10,2)
        cutg.SetPoint(1,-2,2)
        cutg.SetPoint(2,-1,2.5)
        cutg.SetPoint(3,0,3)
        cutg.SetPoint(4,1,3)
        cutg.SetPoint(5,2,3.3)
        cutg.SetPoint(6,3,3.0)
        cutg.SetPoint(7,4,2.5)
        cutg.SetPoint(8,5,2)
        cutg.SetPoint(9,10,2)
        cutg.SetPoint(10,10,1)
        cutg.SetPoint(11,5,1)
        cutg.SetPoint(12,4,1)
        cutg.SetPoint(13,3,-1)
        cutg.SetPoint(14,2,-2)
        cutg.SetPoint(15,1,-2.3)
        cutg.SetPoint(16,0,-1.5)
        cutg.SetPoint(17,-1,-1)
        cutg.SetPoint(18,-2,1)
        cutg.SetPoint(19,-10,1)
        cutg.SetPoint(20,-10,2)
        
    else:
        cutg.SetPoint(0,-25,2)
        cutg.SetPoint(1,-2,2)
        cutg.SetPoint(2,-1,2.5)
        cutg.SetPoint(3,0,3)
        cutg.SetPoint(4,1,3)
        cutg.SetPoint(5,2,3.3)
        cutg.SetPoint(6,3,3.0)
        cutg.SetPoint(7,4,2.5)
        cutg.SetPoint(8,5,2)
        cutg.SetPoint(9,25,2)
        cutg.SetPoint(10,25,0.5)
        cutg.SetPoint(11,5,0.5)
        cutg.SetPoint(12,4,1)
        cutg.SetPoint(13,3,-1)
        cutg.SetPoint(14,2,-2)
        cutg.SetPoint(15,1,-2.3)
        cutg.SetPoint(16,0,-1.5)
        cutg.SetPoint(17,-1,-1)
        cutg.SetPoint(18,-2,0.5)
        cutg.SetPoint(19,-25,0.5)
        cutg.SetPoint(20,-25,2)

    cutg.SetLineColor(kRed)
    cutg.SetLineWidth(2)

    #################################################################################################################################################

    h_hgcer_xAtCer_v_yAtCer  = ROOT.TH2D("hgcer_xAtCer_v_yAtCer","HGC; X; Y;" ,300, -50 ,50, 300, -50, 50)
    h3_Kaons_hgcer_XyAtCer_NPE = ROOT.TH3D("h3_Kaons_hgcer_XyAtCer_NPE","HGC; NPE Sum; X; Y", 300, -40, 40, 300, -40, 40, 300, 0.0, 40)

    h_aero_xAtCer_v_yAtCer  = ROOT.TH2D("aero_xAtCer_v_yAtCer","AERO; X; Y;" ,300, -60 ,60, 300, -60, 60)
    h3_Kaons_aero_XyAtCer_NPE = ROOT.TH3D("h3_Kaons_aero_XyAtCer_NPE","AERO; NPE Sum; X; Y", 300, -60, 60, 300, -60, 60, 300, 0.0, 40)

    h_hgcer_npeSum_v_aero_npeSum  = ROOT.TH2D("hgcer_npeSum_v_aero_npeSum","HGC vs AERO; hgcer; aero;" ,300,0,30, 300, 0, 40)

    for evt in Events_no_cal_hgc_aero_cuts:
        h_hgcer_xAtCer_v_yAtCer.Fill(evt.P_hgcer_xAtCer,evt.P_hgcer_yAtCer)
        h3_Kaons_hgcer_XyAtCer_NPE.Fill(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer, evt.P_hgcer_npeSum)
        h_aero_xAtCer_v_yAtCer.Fill(evt.P_aero_xAtCer,evt.P_aero_yAtCer)
        h3_Kaons_aero_XyAtCer_NPE.Fill(evt.P_aero_xAtCer, evt.P_aero_yAtCer, evt.P_aero_npeSum)
        h_hgcer_npeSum_v_aero_npeSum.Fill(evt.P_hgcer_npeSum,evt.P_aero_npeSum)


    h3_Kaons_hgcer_XyAtCer_NPE_pxy = ROOT.TProfile2D("h3_Kaons_hgcer_XyAtCer_NPE_pxy","HGC (vs NPE) NPE Sum; X; Y",300,-40,40, 300,-40,40,0.0,40)
    h3_Kaons_hgcer_XyAtCer_NPE.Project3DProfile("xy")

    h3_Kaons_aero_XyAtCer_NPE_pxy = ROOT.TProfile2D("h3_Kaons_aero_XyAtCer_NPE_pxy","AERO (vs NPE) NPE Sum; X; Y",300,-60,60, 300,-60,60,0.0,40)
    h3_Kaons_aero_XyAtCer_NPE.Project3DProfile("xy")

    c_hgcer_XY = TCanvas("c_hgcer_XY", "HGC XY")  
    c_hgcer_XY.Divide(2,2)   
    c_hgcer_XY.cd(1)
    h_hgcer_xAtCer_v_yAtCer.Draw("colz")
    cutg.Draw("same")
    c_hgcer_XY.cd(2)
    h3_Kaons_hgcer_XyAtCer_NPE_pxy.Draw("colz")
    cutg.Draw("same")
    c_hgcer_XY.cd(3)
    h_hgcer_xAtCer_v_yAtCer.Draw("[cutg],colz")
    c_hgcer_XY.cd(4)
    h3_Kaons_hgcer_XyAtCer_NPE_pxy.Draw("[cutg],colz")
    c_hgcer_XY.Print(foutpdf+'(')

    c_aero_XY = TCanvas("c_aero_XY", "AERO XY")  
    c_aero_XY.Divide(2,2)   
    c_aero_XY.cd(1)
    h3_Kaons_aero_XyAtCer_NPE_pxy.Draw("colz")
    c_aero_XY.cd(2)
    Events_no_cal_hgc_aero_cuts.Draw("P_aero_npeSum>>h(300,0.3,40)", "!cutg",  "colz")
    c_aero_XY.cd(3)
    h3_Kaons_aero_XyAtCer_NPE_pxy.Draw("[cutg],colz")
    c_aero_XY.cd(4)
    Events_no_cal_hgc_aero_cuts.Draw("P_aero_npeSum>>h2(300,0.3,40)", "cutg",  "colz")
    c_aero_XY.Print(foutpdf)

    c_hgcer_NPE = TCanvas("c_hgcer_NPE", "HGC (with TCutG)")  
    c_hgcer_NPE.Divide(2,2)   
    c_hgcer_NPE.cd(1)
    h_hgcer_npeSum_v_aero_npeSum.Draw("colz")
    c_hgcer_NPE.cd(2)
    Events_no_cal_hgc_aero_cuts.Draw("P_hgcer_npeSum>>h(300,0.3,40)", "!cutg",  "colz")
    c_hgcer_NPE.cd(3)
    h_hgcer_npeSum_v_aero_npeSum.Draw("[cutg],colz")
    c_hgcer_NPE.cd(4)
    Events_no_cal_hgc_aero_cuts.Draw("P_hgcer_npeSum>>h2(300,0.3,40)", "cutg",  "colz")
    c_hgcer_NPE.Print(foutpdf+')')

    #hgcer_did = cutg.IntegralHist(h_hgcer_npeSum_v_aero_npeSum)
    #hgcer_should = h_hgcer_npeSum_v_aero_npeSum.Integral(0,30,0,30)
    hgcer_did = cutg.IntegralHist(h3_Kaons_hgcer_XyAtCer_NPE_pxy)
    hgcer_should = h3_Kaons_hgcer_XyAtCer_NPE_pxy.Integral()

    hgcer_eff = hgcer_did/hgcer_should
    hgcer_error = np.sqrt(((hgcer_did*hgcer_should)+(hgcer_did*hgcer_did))/(hgcer_should*hgcer_should*hgcer_should))

    print("hgcer_did : ",hgcer_did,"\nhgcer_should : ",hgcer_should)
    print("HGCer efficiency : ",hgcer_eff,"+-",hgcer_error)

    effDict = {
        "SHMS_HGC_Kaon_Eff" : hgcer_eff,
        "SHMS_HGC_Kaon_Eff_ERROR" : hgcer_error,
    }

    return effDict
