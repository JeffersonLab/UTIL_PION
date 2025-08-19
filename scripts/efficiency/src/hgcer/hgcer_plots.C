/*
 * Description:
 * ================================================================
 * Time-stamp: "2022-06-15 17:07:30 trottar"
 * ================================================================
 *
 * Author:  Richard L. Trotta III <trotta@cua.edu>
 *
 * Copyright (c) trottar
 */

// 19/10/20 - Stephen Kay, University of Regina

// root .c macro plotting script, reads in desired trees from analysed root file and plots some stuff
// Saves  pdf file with plots and a .root file
#define Cherenkov_plots_cxx

// Include relevant stuff
#include <TStyle.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TGaxis.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <TSystem.h>
#include <TTree.h>
#include <TArc.h>
#include <TCutG.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>

// Input should be the input root file name (including suffix) and an output file name string (without any suffix)
void hgcer_plots(string InFilename = "", string OutFilename = "")
{

  TString Hostname = gSystem->HostName();
  TString User = (gSystem->GetUserInfo())->fUser;
  TString Replaypath;
  TString Outpath;
  TString Outpath1;
  TString rootFile;
  gStyle->SetPalette(55);

  // Set paths depending on system you're running on
  if(Hostname.Contains("farm")){
    Replaypath = "/group/c-kaonlt/USERS/"+User+"/hallc_replay_lt";
    // Output path for root file
    Outpath = Replaypath+"/UTIL_KAONLT/OUTPUT/Analysis/KaonLT";
    // Output path for output file
    Outpath1 = Replaypath+"/UTIL_KAONLT/scripts/efficiency/src/OUTPUTS";
  }
  else if(Hostname.Contains("qcd")){
    Replaypath = "/group/c-kaonlt/USERS/"+User+"/hallc_replay_lt";
    Outpath = Replaypath+"/UTIL_KAONLT/scripts/CoinTimePeak/OUTPUT";
  }
  else if (Hostname.Contains("phys.uregina.ca")){
    Replaypath = "/home/"+User+"/work/JLab/hallc_replay_lt";
    Outpath = Replaypath+"/UTIL_KAONLT/scripts/CoinTimePeak/OUTPUT";
  }

  // Add more as needed for your own envrionment
  if(InFilename == "") {
    cout << "Enter a Filename to analyse: ";
    cin >> InFilename;
  }  
  if(OutFilename == "") {
    cout << "Enter a Filename to output to: ";
    cin >> OutFilename;
  }

  TString TInFilename = InFilename;
  rootFile = Outpath+"/"+TInFilename;
  // Complain and exit if your file doesn't exist
  if (gSystem->AccessPathName(rootFile) == kTRUE){
    cerr << "!!!!! ERROR !!!!! " << endl << rootFile <<  " not found" << endl <<  "!!!!! ERRROR !!!!!" << endl;
    exit;
  }

  //################################################################################################################################################

  TFile *InFile = new TFile(rootFile, "READ");
  TString TOutFilename = OutFilename;
  // Establish the names of our output files quickly
  TString foutname = Outpath1+"/" + TOutFilename + ".root";
  TString foutpdf = Outpath1+"/" + TOutFilename + ".pdf";

  // Particles information with HGC cuts
  
  TTree* Pions         = (TTree*)InFile->Get("SHMS_Pions")  ; Long64_t nEntries_Pions = (Long64_t)Pions->GetEntries();
  TTree* Kaons         = (TTree*)InFile->Get("SHMS_Kaons")  ; Long64_t nEntries_Kaons = (Long64_t)Kaons->GetEntries();
  TTree* Protons       = (TTree*)InFile->Get("SHMS_Protons"); Long64_t nEntries_Protons = (Long64_t)Protons->GetEntries();
  TTree* Positrons     = (TTree*)InFile->Get("SHMS_Positrons")  ; Long64_t nEntries_Positrons = (Long64_t)Positrons->GetEntries();

 // Particles information no HGC cuts
  
  TTree* Pions_No_HGC_Cuts         = (TTree*)InFile->Get("SHMS_Pions_Without_HGC_Cuts")      ; Long64_t nEntries_Pions_No_HGC_Cuts       = (Long64_t)Pions_No_HGC_Cuts->GetEntries();
  TTree* Positrons_No_HGC_Cuts     = (TTree*)InFile->Get("SHMS_Positrons_Without_HGC_Cuts")  ; Long64_t nEntries_Positrons_No_HGC_Cuts   = (Long64_t)Positrons_No_HGC_Cuts->GetEntries();
  TTree* Kaons_No_HGC_Cuts         = (TTree*)InFile->Get("SHMS_Kaons_Without_HGC_Cuts")      ; Long64_t nEntries_Kaons_No_HGC_Cuts       = (Long64_t)Kaons_No_HGC_Cuts->GetEntries();
  TTree* Protons_No_HGC_Cuts       = (TTree*)InFile->Get("SHMS_Protons_Without_HGC_Cuts")    ; Long64_t nEntries_Protons_No_HGC_Cuts     = (Long64_t)Protons_No_HGC_Cuts->GetEntries();
  TTree* Pions_No_Aero_Cuts        = (TTree*)InFile->Get("SHMS_Pions_Aero_Without_Aero_Cuts"); Long64_t nEntries_Pions_No_Aero_Cuts  = (Long64_t)Pions_No_Aero_Cuts->GetEntries();
  TTree* Pions_No_Cal_Cuts         = (TTree*)InFile->Get("SHMS_Pions_Cal_Without_Cal_Cuts")  ; Long64_t nEntries_Pions_No_Cal_Cuts     = (Long64_t)Pions_No_Cal_Cuts->GetEntries();

 // Particles information no Cal HGC and Aero cuts

  TTree* Events_no_cal_hgc_cuts       = (TTree*)InFile->Get("SHMS_cut_no_Cal_HGC")       ; Long64_t nEntries_Events_no_cal_hgc_cuts       = (Long64_t)Events_no_cal_hgc_cuts->GetEntries();
  TTree* Events_no_cal_hgc_aero_cuts  = (TTree*)InFile->Get("SHMS_cut_no_Cal_HGC_Aero")  ; Long64_t nEntries_Events_no_cal_hgc_aero_cuts  = (Long64_t)Events_no_cal_hgc_aero_cuts->GetEntries();

 // Particles information no cuts
  
  TTree* SHMS_Events = (TTree*)InFile->Get("SHMS_Events")  ; Long64_t nEntries_SHMS_Events = (Long64_t)SHMS_Events->GetEntries();


  Double_t aero_npeSum; Pions_No_Aero_Cuts->SetBranchAddress("P_aero_npeSum", &aero_npeSum);
  Double_t cal_etotnorm; Pions_No_Cal_Cuts->SetBranchAddress("P_cal_etotnorm", &cal_etotnorm);

  //################################################################################################################################################

  // Defined Geomatrical cuts
  TCutG *cutg = new TCutG("cutg",21);
  cutg->SetVarX("P_hgcer_yAtCer");
  cutg->SetVarY("P_hgcer_xAtCer");
  cutg->SetPoint(0,-25,2);
  cutg->SetPoint(1,-2,2);
  cutg->SetPoint(2,-1,2.5);
  cutg->SetPoint(3,0,3);
  cutg->SetPoint(4,1,3);
  cutg->SetPoint(5,2,3.3);
  cutg->SetPoint(6,3,3.0);
  cutg->SetPoint(7,4,2.5);
  cutg->SetPoint(8,5,2);
  cutg->SetPoint(9,25,2);
  cutg->SetPoint(10,25,0.5);
  cutg->SetPoint(11,5,0.5);
  cutg->SetPoint(12,4,1);
  cutg->SetPoint(13,3,-1);
  cutg->SetPoint(14,2,-2);
  cutg->SetPoint(15,1,-2.3);
  cutg->SetPoint(16,0,-1.5);
  cutg->SetPoint(17,-1,-1);
  cutg->SetPoint(18,-2,0.5);
  cutg->SetPoint(19,-25,0.5);
  cutg->SetPoint(20,-25,2);

   cutg->SetLineColor(kRed);
  cutg->SetLineWidth(5);
  // cut for npe

  TCutG *cutg1 = new TCutG("cutg1",21);
  cutg1->SetVarX("P_hgcer_npeSum");
  cutg1->SetVarY("P_aero_npeSum");
  cutg1->SetPoint(0,-10,2);
  cutg1->SetPoint(1,-2,2);
  cutg1->SetPoint(2,-1,2.5);
  cutg1->SetPoint(3,0,3);
  cutg1->SetPoint(4,1,3);
  cutg1->SetPoint(5,2,3.3);
  cutg1->SetPoint(6,3,3.0);
  cutg1->SetPoint(7,4,2.5);
  cutg1->SetPoint(8,5,2);
  cutg1->SetPoint(9,10,2);
  cutg1->SetPoint(10,10,1);
  cutg1->SetPoint(11,5,1);
  cutg1->SetPoint(12,4,1);
  cutg1->SetPoint(13,3,-1);
  cutg1->SetPoint(14,2,-2);
  cutg1->SetPoint(15,1,-2.3);
  cutg1->SetPoint(16,0,-1.5);
  cutg1->SetPoint(17,-1,-1);
  cutg1->SetPoint(18,-2,1);
  cutg1->SetPoint(19,-10,1);
  cutg1->SetPoint(20,-10,2);

  cutg1->SetLineColor(kRed);
  cutg1->SetLineWidth(5);

  //################################################################################################################################################

  // Set branch address for no cuts
 
  Double_t XAtCer; SHMS_Events->SetBranchAddress("P_hgcer_xAtCer", &XAtCer);
  Double_t YAtCer; SHMS_Events->SetBranchAddress("P_hgcer_yAtCer", &YAtCer);
  Double_t Prshower; SHMS_Events->SetBranchAddress("P_cal_pr_eplane", &Prshower);
  Double_t Cal; SHMS_Events->SetBranchAddress("P_cal_fly_earray", &Cal);
  Double_t gtr_x; SHMS_Events->SetBranchAddress("P_gtr_x", &gtr_x);
  Double_t gtr_y; SHMS_Events->SetBranchAddress("P_gtr_y", &gtr_y);
  Double_t gtr_beta; SHMS_Events->SetBranchAddress("P_gtr_beta", &gtr_beta);
  Double_t gtr_p; SHMS_Events->SetBranchAddress("P_gtr_p", &gtr_p);
  Double_t gtr_dp; SHMS_Events->SetBranchAddress("P_gtr_dp", &gtr_dp);
  Double_t gtr_xp; SHMS_Events->SetBranchAddress("P_gtr_xp", &gtr_xp);
  Double_t gtr_yp; SHMS_Events->SetBranchAddress("P_gtr_yp", &gtr_yp);
  Double_t hgcer_npeSum; SHMS_Events->SetBranchAddress("P_hgcer_npeSum", &hgcer_npeSum);
  Double_t P_cal_etot; SHMS_Events->SetBranchAddress("P_cal_etotnorm", &P_cal_etot);
  Double_t CTime_ePi_ROC1; SHMS_Events->SetBranchAddress("CTime_ePiCoinTime_ROC1", &CTime_ePi_ROC1);
  Double_t CTime_eK_ROC1; SHMS_Events->SetBranchAddress("CTime_eKCoinTime_ROC1", &CTime_eK_ROC1);
  Double_t CTime_eP_ROC1; SHMS_Events->SetBranchAddress("CTime_epCoinTime_ROC1", &CTime_eP_ROC1);

  // Set branch address for no cal and  hgc cuts
  Double_t P_no_cal_hgc_cuts_hgcer_npeSum; Events_no_cal_hgc_cuts->SetBranchAddress("P_hgcer_npeSum", &P_no_cal_hgc_cuts_hgcer_npeSum);
  Double_t P_no_cal_hgc_etot; Events_no_cal_hgc_cuts->SetBranchAddress("P_cal_etotnorm", &P_no_cal_hgc_etot);

 // Set branch address for no cal, hgc and aero cuts
  Double_t P_no_cal_hgc_aero_cuts_hgcer_npeSum; Events_no_cal_hgc_aero_cuts->SetBranchAddress("P_hgcer_npeSum", &P_no_cal_hgc_aero_cuts_hgcer_npeSum);
  Double_t P_no_cal_hgc_aero_cuts_etot; Events_no_cal_hgc_aero_cuts->SetBranchAddress("P_cal_etotnorm", &P_no_cal_hgc_aero_cuts_etot);
  Double_t P_no_cal_hgc_aero_cuts_aero_npeSum; Events_no_cal_hgc_aero_cuts->SetBranchAddress("P_aero_npeSum", &P_no_cal_hgc_aero_cuts_aero_npeSum);
  Double_t P_no_cal_hgc_aero_cuts_aero_xAtCer; Events_no_cal_hgc_aero_cuts->SetBranchAddress("P_aero_xAtCer", &P_no_cal_hgc_aero_cuts_aero_xAtCer);
  Double_t P_no_cal_hgc_aero_cuts_aero_yAtCer; Events_no_cal_hgc_aero_cuts->SetBranchAddress("P_aero_yAtCer", &P_no_cal_hgc_aero_cuts_aero_yAtCer);
  Double_t P_no_cal_hgc_aero_cuts_hgcer_xAtCer; Events_no_cal_hgc_aero_cuts->SetBranchAddress("P_hgcer_xAtCer", &P_no_cal_hgc_aero_cuts_hgcer_xAtCer);
  Double_t P_no_cal_hgc_aero_cuts_hgcer_yAtCer; Events_no_cal_hgc_aero_cuts->SetBranchAddress("P_hgcer_yAtCer", &P_no_cal_hgc_aero_cuts_hgcer_yAtCer);

  //################################################################################################################################################
 
  TCutG *cutg2 = new TCutG("cutg2",21);
  cutg2->SetVarX("P_no_cal_hgc_aero_cuts_hgcer_yAtCer");
  cutg2->SetVarY("P_no_cal_hgc_aero_cuts_hgcer_XAtCer");
  cutg2->SetPoint(0,-25,2);
  cutg2->SetPoint(1,-2,2);
  cutg2->SetPoint(2,-1,2.5);
  cutg2->SetPoint(3,0,3);
  cutg2->SetPoint(4,1,3);
  cutg2->SetPoint(5,2,3.3);
  cutg2->SetPoint(6,3,3.0);
  cutg2->SetPoint(7,4,2.5);
  cutg2->SetPoint(8,5,2);
  cutg2->SetPoint(9,25,2);
  cutg2->SetPoint(10,25,0.5);
  cutg2->SetPoint(11,5,0.5);
  cutg2->SetPoint(12,4,1);
  cutg2->SetPoint(13,3,-1);
  cutg2->SetPoint(14,2,-2);
  cutg2->SetPoint(15,1,-2.3);
  cutg2->SetPoint(16,0,-1.5);
  cutg2->SetPoint(17,-1,-1);
  cutg2->SetPoint(18,-2,0.5);
  cutg2->SetPoint(19,-25,0.5);
  cutg2->SetPoint(20,-25,2);

  cutg2->SetLineColor(kRed);
  cutg2->SetLineWidth(5);

  //################################################################################################################################################
 
 // Set branch address for Pions with HGC cuts

 //Pions with HGC
  Double_t Pions_XAtCer; Pions->SetBranchAddress("P_hgcer_xAtCer", &Pions_XAtCer);
  Double_t Pions_YAtCer; Pions->SetBranchAddress("P_hgcer_yAtCer", &Pions_YAtCer);
  Double_t Pions_Cal_Pr_Shower; Pions->SetBranchAddress("P_cal_pr_eplane", &Pions_Cal_Pr_Shower);
  Double_t Pions_Cal_Shower; Pions->SetBranchAddress("P_cal_fly_earray", &Pions_Cal_Shower);
  Double_t Pions_SHMS_gtr_x; Pions->SetBranchAddress("P_gtr_x", &Pions_SHMS_gtr_x);
  Double_t Pions_SHMS_gtr_y; Pions->SetBranchAddress("P_gtr_y", &Pions_SHMS_gtr_y);
  Double_t Pions_SHMS_hgcer_npeSum; Pions->SetBranchAddress("P_hgcer_npeSum", &Pions_SHMS_hgcer_npeSum);
  Double_t Pions_SHMS_aero_npeSum; Pions->SetBranchAddress("P_aero_npeSum", &Pions_SHMS_aero_npeSum);
  Double_t Pions_gtr_beta; Pions->SetBranchAddress("P_gtr_beta", &Pions_gtr_beta);
  Double_t Pions_gtr_p; Pions->SetBranchAddress("P_gtr_p", &Pions_gtr_p);
  Double_t Pions_gtr_dp; Pions->SetBranchAddress("P_gtr_dp", &Pions_gtr_dp);
  Double_t Pions_gtr_xp; Pions->SetBranchAddress("P_gtr_xp", &Pions_gtr_xp);
  Double_t Pions_gtr_yp; Pions->SetBranchAddress("P_gtr_yp", &Pions_gtr_yp); 
  Double_t Pions_P_cal_etot; Pions->SetBranchAddress("P_cal_etotnorm", &Pions_P_cal_etot);
  Double_t Pions_CTime_ePi_ROC1; Pions->SetBranchAddress("CTime_ePiCoinTime_ROC1", &Pions_CTime_ePi_ROC1);

//Positrons with HGC
  Double_t Positrons_XAtCer; Positrons->SetBranchAddress("P_hgcer_xAtCer", &Positrons_XAtCer);
  Double_t Positrons_YAtCer; Positrons->SetBranchAddress("P_hgcer_yAtCer", &Positrons_YAtCer);
  Double_t Positrons_Cal_Pr_Shower; Positrons->SetBranchAddress("P_cal_pr_eplane", &Positrons_Cal_Pr_Shower);
  Double_t Positrons_Cal_Shower; Positrons->SetBranchAddress("P_cal_fly_earray", &Positrons_Cal_Shower);
  Double_t Positrons_SHMS_gtr_x; Positrons->SetBranchAddress("P_gtr_x", &Positrons_SHMS_gtr_x);
  Double_t Positrons_SHMS_gtr_y; Positrons->SetBranchAddress("P_gtr_y", &Positrons_SHMS_gtr_y);
  Double_t Positrons_SHMS_hgcer_npeSum; Positrons->SetBranchAddress("P_hgcer_npeSum", &Positrons_SHMS_hgcer_npeSum);
  Double_t Positrons_SHMS_aero_npeSum; Positrons->SetBranchAddress("P_aero_npeSum", &Positrons_SHMS_aero_npeSum);
  Double_t Positrons_gtr_beta; Positrons->SetBranchAddress("P_gtr_beta", &Positrons_gtr_beta);
  Double_t Positrons_gtr_p; Positrons->SetBranchAddress("P_gtr_p", &Positrons_gtr_p);
  Double_t Positrons_gtr_dp; Positrons->SetBranchAddress("P_gtr_dp", &Positrons_gtr_dp);
  Double_t Positrons_gtr_xp; Positrons->SetBranchAddress("P_gtr_xp", &Positrons_gtr_xp);
  Double_t Positrons_gtr_yp; Positrons->SetBranchAddress("P_gtr_yp", &Positrons_gtr_yp); 
  Double_t Positrons_P_cal_etot; Positrons->SetBranchAddress("P_cal_etotnorm", &Positrons_P_cal_etot);
  // Double_t Positrons_CTime_ePi_ROC1; Pions->SetBranchAddress("CTime_ePiCoinTime_ROC1", &_CTime_ePi_ROC1);

 // Set branch address for Kaons with HGC cuts

  //Kaons
  Double_t Kaons_XAtCer; Kaons->SetBranchAddress("P_hgcer_xAtCer", &Kaons_XAtCer);
  Double_t Kaons_YAtCer; Kaons->SetBranchAddress("P_hgcer_yAtCer", &Kaons_YAtCer);
  Double_t Kaons_Cal_Pr_Shower; Kaons->SetBranchAddress("P_cal_pr_eplane", &Kaons_Cal_Pr_Shower);
  Double_t Kaons_Cal_Shower; Kaons->SetBranchAddress("P_cal_fly_earray", &Kaons_Cal_Shower);
  Double_t Kaons_SHMS_gtr_x; Kaons->SetBranchAddress("P_gtr_x", &Kaons_SHMS_gtr_x);
  Double_t Kaons_SHMS_gtr_y; Kaons->SetBranchAddress("P_gtr_y", &Kaons_SHMS_gtr_y);
  Double_t Kaons_SHMS_hgcer_npeSum; Kaons->SetBranchAddress("P_hgcer_npeSum", &Kaons_SHMS_hgcer_npeSum);

 // Set branch address for Protons with HGC cuts
 
 //Protons
  Double_t Protons_XAtCer; Protons->SetBranchAddress("P_hgcer_xAtCer", &Protons_XAtCer);
  Double_t Protons_YAtCer; Protons->SetBranchAddress("P_hgcer_yAtCer", &Protons_YAtCer);
  Double_t Protons_Cal_Pr_Shower; Protons->SetBranchAddress("P_cal_pr_eplane", &Protons_Cal_Pr_Shower);
  Double_t Protons_Cal_Shower; Protons->SetBranchAddress("P_cal_fly_earray", &Protons_Cal_Shower);
  Double_t Protons_SHMS_gtr_x; Protons->SetBranchAddress("P_gtr_x", &Protons_SHMS_gtr_x);
  Double_t Protons_SHMS_gtr_y; Protons->SetBranchAddress("P_gtr_y", &Protons_SHMS_gtr_y);
  Double_t Protons_SHMS_hgcer_npeSum; Protons->SetBranchAddress("P_hgcer_npeSum", &Protons_SHMS_hgcer_npeSum);
  
  // Set branch address for Pions with no HGC cuts

  //Pions
  Double_t Pions_No_HGC_Cuts_XAtCer; Pions_No_HGC_Cuts->SetBranchAddress("P_hgcer_xAtCer", &Pions_No_HGC_Cuts_XAtCer);
  Double_t Pions_No_HGC_Cuts_YAtCer; Pions_No_HGC_Cuts->SetBranchAddress("P_hgcer_yAtCer", &Pions_No_HGC_Cuts_YAtCer);
  Double_t Pions_No_HGC_Cuts_Cal_Pr_Shower; Pions_No_HGC_Cuts->SetBranchAddress("P_cal_pr_eplane", &Pions_No_HGC_Cuts_Cal_Pr_Shower);
  Double_t Pions_No_HGC_Cuts_Cal_Shower; Pions_No_HGC_Cuts->SetBranchAddress("P_cal_fly_earray", &Pions_No_HGC_Cuts_Cal_Shower);
  Double_t Pions_No_HGC_Cuts_SHMS_gtr_x; Pions_No_HGC_Cuts->SetBranchAddress("P_gtr_x", &Pions_No_HGC_Cuts_SHMS_gtr_x);
  Double_t Pions_No_HGC_Cuts_SHMS_gtr_y; Pions_No_HGC_Cuts->SetBranchAddress("P_gtr_y", &Pions_No_HGC_Cuts_SHMS_gtr_y);
  Double_t Pions_No_HGC_Cuts_SHMS_hgcer_npeSum; Pions_No_HGC_Cuts->SetBranchAddress("P_hgcer_npeSum", &Pions_No_HGC_Cuts_SHMS_hgcer_npeSum);
  Double_t Pions_No_HGC_Cuts_SHMS_aero_npeSum; Pions_No_HGC_Cuts->SetBranchAddress("P_aero_npeSum", &Pions_No_HGC_Cuts_SHMS_aero_npeSum);
  Double_t Pions_No_HGC_Cuts_gtr_beta; Pions_No_HGC_Cuts->SetBranchAddress("P_gtr_beta", &Pions_No_HGC_Cuts_gtr_beta);
  Double_t Pions_No_HGC_Cuts_gtr_p; Pions_No_HGC_Cuts->SetBranchAddress("P_gtr_p", &Pions_No_HGC_Cuts_gtr_p);
  Double_t Pions_No_HGC_Cuts_gtr_dp; Pions_No_HGC_Cuts->SetBranchAddress("P_gtr_dp", &Pions_No_HGC_Cuts_gtr_dp);
  Double_t Pions_No_HGC_Cuts_gtr_xp; Pions_No_HGC_Cuts->SetBranchAddress("P_gtr_xp", &Pions_No_HGC_Cuts_gtr_xp);
  Double_t Pions_No_HGC_Cuts_gtr_yp; Pions_No_HGC_Cuts->SetBranchAddress("P_gtr_yp", &Pions_No_HGC_Cuts_gtr_yp); 
  Double_t Pions_No_HGC_Cuts_P_cal_etot; Pions_No_HGC_Cuts->SetBranchAddress("P_cal_etotnorm", &Pions_No_HGC_Cuts_P_cal_etot);
  Double_t Pions_No_HGC_Cuts_CTime_ePi_ROC1; Pions_No_HGC_Cuts->SetBranchAddress("CTime_ePiCoinTime_ROC1", &Pions_No_HGC_Cuts_CTime_ePi_ROC1);
  //Double_t CTime_eK_ROC1_no_HGC_cut; Pions_No_HGC_Cuts->SetBranchAddress("CTime_eKCoinTime_ROC1", &CTime_eK_ROC1_no_HGC_cut);
  // Double_t CTime_eP_ROC1_no_HGC_cut; Pions_No_HGC_Cuts->SetBranchAddress("CTime_epCoinTime_ROC1", &CTime_eP_ROC1_no_HGC_cut);

  //Positrons
  Double_t Positrons_No_HGC_Cuts_XAtCer; Positrons_No_HGC_Cuts->SetBranchAddress("P_hgcer_xAtCer", &Positrons_No_HGC_Cuts_XAtCer);
  Double_t Positrons_No_HGC_Cuts_YAtCer; Positrons_No_HGC_Cuts->SetBranchAddress("P_hgcer_yAtCer", &Positrons_No_HGC_Cuts_YAtCer);
  Double_t Positrons_No_HGC_Cuts_Cal_Pr_Shower; Positrons_No_HGC_Cuts->SetBranchAddress("P_cal_pr_eplane", &Positrons_No_HGC_Cuts_Cal_Pr_Shower);
  Double_t Positrons_No_HGC_Cuts_Cal_Shower; Positrons_No_HGC_Cuts->SetBranchAddress("P_cal_fly_earray", &Positrons_No_HGC_Cuts_Cal_Shower);
  Double_t Positrons_No_HGC_Cuts_SHMS_gtr_x; Positrons_No_HGC_Cuts->SetBranchAddress("P_gtr_x", &Positrons_No_HGC_Cuts_SHMS_gtr_x);
  Double_t Positrons_No_HGC_Cuts_SHMS_gtr_y; Positrons_No_HGC_Cuts->SetBranchAddress("P_gtr_y", &Positrons_No_HGC_Cuts_SHMS_gtr_y);
  Double_t Positrons_No_HGC_Cuts_SHMS_hgcer_npeSum; Positrons_No_HGC_Cuts->SetBranchAddress("P_hgcer_npeSum", &Positrons_No_HGC_Cuts_SHMS_hgcer_npeSum);
  Double_t Positrons_No_HGC_Cuts_SHMS_aero_npeSum; Positrons_No_HGC_Cuts->SetBranchAddress("P_aero_npeSum", &Positrons_No_HGC_Cuts_SHMS_aero_npeSum);
  Double_t Positrons_No_HGC_Cuts_gtr_beta; Positrons_No_HGC_Cuts->SetBranchAddress("P_gtr_beta", &Positrons_No_HGC_Cuts_gtr_beta);
  Double_t Positrons_No_HGC_Cuts_gtr_p; Positrons_No_HGC_Cuts->SetBranchAddress("P_gtr_p", &Positrons_No_HGC_Cuts_gtr_p);
  Double_t Positrons_No_HGC_Cuts_gtr_dp; Positrons_No_HGC_Cuts->SetBranchAddress("P_gtr_dp", &Positrons_No_HGC_Cuts_gtr_dp);
  Double_t Positrons_No_HGC_Cuts_gtr_xp; Positrons_No_HGC_Cuts->SetBranchAddress("P_gtr_xp", &Positrons_No_HGC_Cuts_gtr_xp);
  Double_t Positrons_No_HGC_Cuts_gtr_yp; Positrons_No_HGC_Cuts->SetBranchAddress("P_gtr_yp", &Positrons_No_HGC_Cuts_gtr_yp); 
  Double_t Positrons_No_HGC_Cuts_P_cal_etot; Positrons_No_HGC_Cuts->SetBranchAddress("P_cal_etotnorm", &Positrons_No_HGC_Cuts_P_cal_etot);
  // Double_t Pions_No_HGC_Cuts_CTime_ePi_ROC1; Pions_No_HGC_Cuts->SetBranchAddress("CTime_ePiCoinTime_ROC1", &Pions_No_HGC_Cuts_CTime_ePi_ROC1);
 
 // Set branch address for Kaons with no HGC cuts

  //Kaons
  Double_t Kaons_No_HGC_Cuts_XAtCer; Kaons_No_HGC_Cuts->SetBranchAddress("P_hgcer_xAtCer", &Kaons_No_HGC_Cuts_XAtCer);
  Double_t Kaons_No_HGC_Cuts_YAtCer; Kaons_No_HGC_Cuts->SetBranchAddress("P_hgcer_yAtCer", &Kaons_No_HGC_Cuts_YAtCer);
  Double_t Kaons_No_HGC_Cuts_Cal_Pr_Shower; Kaons_No_HGC_Cuts->SetBranchAddress("P_cal_pr_eplane", &Kaons_No_HGC_Cuts_Cal_Pr_Shower);
  Double_t Kaons_No_HGC_Cuts_Cal_Shower; Kaons_No_HGC_Cuts->SetBranchAddress("P_cal_fly_earray", &Kaons_No_HGC_Cuts_Cal_Shower);
  Double_t Kaons_No_HGC_Cuts_SHMS_gtr_x; Kaons_No_HGC_Cuts->SetBranchAddress("P_gtr_x", &Kaons_No_HGC_Cuts_SHMS_gtr_x);
  Double_t Kaons_No_HGC_Cuts_SHMS_gtr_y; Kaons_No_HGC_Cuts->SetBranchAddress("P_gtr_y", &Kaons_No_HGC_Cuts_SHMS_gtr_y);
  Double_t Kaons_No_HGC_Cuts_SHMS_hgcer_npeSum; Kaons_No_HGC_Cuts->SetBranchAddress("P_hgcer_npeSum", &Kaons_No_HGC_Cuts_SHMS_hgcer_npeSum);
  
 // Set branch address for Protons with no HGC cuts

 //Protons
  Double_t Protons_No_HGC_Cuts_XAtCer; Protons_No_HGC_Cuts->SetBranchAddress("P_hgcer_xAtCer", &Protons_No_HGC_Cuts_XAtCer);
  Double_t Protons_No_HGC_Cuts_YAtCer; Protons_No_HGC_Cuts->SetBranchAddress("P_hgcer_yAtCer", &Protons_No_HGC_Cuts_YAtCer);
  Double_t Protons_No_HGC_Cuts_Cal_Pr_Shower; Protons_No_HGC_Cuts->SetBranchAddress("P_cal_pr_eplane", &Protons_No_HGC_Cuts_Cal_Pr_Shower);
  Double_t Protons_No_HGC_Cuts_Cal_Shower; Protons_No_HGC_Cuts->SetBranchAddress("P_cal_fly_earray", &Protons_No_HGC_Cuts_Cal_Shower);
  Double_t Protons_No_HGC_Cuts_SHMS_gtr_x; Protons_No_HGC_Cuts->SetBranchAddress("P_gtr_x", &Protons_No_HGC_Cuts_SHMS_gtr_x);
  Double_t Protons_No_HGC_Cuts_SHMS_gtr_y; Protons_No_HGC_Cuts->SetBranchAddress("P_gtr_y", &Protons_No_HGC_Cuts_SHMS_gtr_y);
  Double_t Protons_No_HGC_Cuts_SHMS_hgcer_npeSum; Protons_No_HGC_Cuts->SetBranchAddress("P_hgcer_npeSum", &Protons_No_HGC_Cuts_SHMS_hgcer_npeSum);
 
  //################################################################################################################################################

  // Define Histograms for with HGC cuts
  
  // For 1D histos, can easily create directly from the corresponding branch
  //Pions->Draw("P_hgcer_npeSum >> h1_CT_Pions", "", "goff"); 
  //Kaons->Draw("CTime_eKCoinTime_ROC1 >> h1_CT_Kaons", "", "goff"); 
  //Protons->Draw("CTime_eKCoinTime_ROC1 >> h1_CT_Protons", "", "goff"); 

 //Histograms for aerogel and Cal
  TH1D *h1_aero_npeSum = new TH1D("h1_aero_npeSum","Aerogel; P_aero_npeSum; Events;", 300, 0.0, 30);

  for(Long64_t i = 0; i <nEntries_Pions_No_Aero_Cuts; i++){
    Pions_No_Aero_Cuts->GetEntry(i);
    h1_aero_npeSum->Fill(aero_npeSum);
  }

    TH1D *h1_cal_etot = new TH1D("h1_cal_etot","Calorimeter; P_cal_etotnorm; Events;", 300, 0.0, 5);
    //x TH2D *h2_cal_vs_hgc = new TH2D("h2_cal_vs_hgc","Calorimeter VS HGC; P_cal_etotnorm; P_hgcer_npeSum;", 300, 0.0, 5, 300, 0.0, 30);
      
    for(Long64_t i = 0; i <nEntries_Pions_No_Cal_Cuts; i++){
      Pions_No_Cal_Cuts->GetEntry(i);
      h1_cal_etot->Fill(cal_etotnorm); 
    }

  //Histograms for not  cuts

  TH1D *h1_hgcer_npeSum = new TH1D("h1_hgcer_npeSum","NPE vs Events; NPE; Events;", 300, 0.0, 30);
  TH1D *h1_XAtCer = new TH1D("h1_XAtCer","XAtCer; XAtCer; Events;", 300, -40, 40);
  TH1D *h1_YAtCer = new TH1D("h1_YAtCer","YAtCer; YAtCer; Events;", 300, -40, 40);
  TH1D *h1_Cal = new TH1D("h1_Cal","Calorimeter; Cal; Events;", 300, 0.0, 10.0);
  TH1D *h1_Prshower = new TH1D("h1_Prshower","Preshower; Prshower; Events;", 300, 0.0, 10.0);
  TH1D *h1_Xgtr = new TH1D("h1_Xgtr","X gtr; X gtr ; Events;", 300, -3.0, 3.0);
  TH1D *h1_Ygtr = new TH1D("h1_Ygtr","Y gtr; Y gtr ; Events;", 300, -3.0, 3.0);
  TH1D *h1_P_etot = new TH1D("h1_P_etot","P_cal_etotnorm; P_cal_etotnorm; Events;", 300, 0.0, 10.0);
  TH1D *h1_gtr_beta = new TH1D("h1_gtr_beta","P_gtr_beta; P_gtr_beta; Events;", 300, 0.0, 10.0);
  TH1D *h1_gtr_p = new TH1D("h1_gtr_p","P_gtr_p; P_gtr_p; Events;", 300, -10.0, 10.0);
  TH1D *h1_gtr_dp = new TH1D("h1_gtr_dp","P_gtr_dp; P_gtr_dp; Events;", 300, -30.0, 30.0);
  TH1D *h1_gtr_xp = new TH1D("h1_gtr_xp","P_gtr_xp; P_gtr_xp; Events;", 300, -40.0, 40.0);
  TH1D *h1_gtr_yp = new TH1D("h1_gtr_yp","P_gtr_yp; P_gtr_yp; Events;", 300, -40.0, 40.0);
  TH1D *h1_CTime_ePi_ROC1 = new TH1D("h1_CTime_ePi_ROC1","CTime_ePiCoinTime_ROC1; CTime_ePiCoinTime_ROC1; Events;", 300, -10.0, 100.0);
  TH1D *h1_CTime_eK_ROC1 = new TH1D("h1_CTime_eK_ROC1","CTime_eKCoinTime_ROC1; CTime_eKCoinTime_ROC1; Events;", 300, -10.0, 100.0);
  TH1D *h1_CTime_eP_ROC1 = new TH1D("h1_CTime_eP_ROC1","CTime_epCoinTime_ROC1; CTime_epCoinTime_ROC1; Events;", 300, -10.0, 100.0);
  
    for(Long64_t i = 0; i <nEntries_SHMS_Events; i++){
      SHMS_Events->GetEntry(i);
      h1_hgcer_npeSum->Fill(hgcer_npeSum);
      h1_XAtCer->Fill(XAtCer);
      h1_YAtCer->Fill(YAtCer);
      h1_Cal->Fill(Cal/gtr_p);
      h1_Prshower->Fill(Prshower/gtr_p);
      h1_Xgtr->Fill(gtr_x);
      h1_Ygtr->Fill(gtr_y);
      h1_P_etot->Fill(P_cal_etot);  
      h1_gtr_beta->Fill(gtr_beta);  
      h1_gtr_p->Fill(gtr_p);
      h1_gtr_dp->Fill(gtr_dp);
      h1_gtr_xp->Fill(gtr_xp);
      h1_gtr_yp->Fill(gtr_yp);
      h1_CTime_ePi_ROC1->Fill(CTime_ePi_ROC1);
      h1_CTime_eK_ROC1->Fill(CTime_eK_ROC1);
      h1_CTime_eP_ROC1->Fill(CTime_eP_ROC1);
  }
    //Histograms for cuts + no Cal, HGC and Aero cuts

    TH2D *h2_events_no_cal_hgc_cuts = new TH2D("h2_events_no_cal_hgc_cuts","Calorimeter VS HGC; P_cal_etotnorm; P_hgcer_npeSum;", 300, 0.0, 3.0, 300, 0.0, 30.0);
    for(Long64_t i = 0; i < nEntries_Events_no_cal_hgc_cuts; i++){
      Events_no_cal_hgc_cuts->GetEntry(i);
      h2_events_no_cal_hgc_cuts->Fill(P_no_cal_hgc_etot, P_no_cal_hgc_cuts_hgcer_npeSum);
    }

    TH2D *h2_events_no_cal_hgc_aero_cuts = new TH2D("h2_events_no_cal_hgc_aero_cuts","HGC vs Aero; P_hgcer_npeSum; P_aero_npeSum;", 300, 0.0, 30.0, 300, 0.0, 30.0);
    TH2D *h2_events_no_cal_aero_cuts = new TH2D("h2_events_no_cal_aero_cuts","Calorimeter VS Aero; P_cal_etotnorm; P_aero_npeSum;", 300, 0.0, 3.0, 300, 0.0, 30.0);
    TH3D *h3_events_no_cal_aero_cuts = new TH3D("h3_events_no_cal_aero_cuts","Aero; P_aero_xAtCer; P_aero_yAtCer;  P_aero_npeSum;", 300, -50, 50, 300, -50, 50, 300, 0, 30);
    for(Long64_t i = 0; i < nEntries_Events_no_cal_hgc_aero_cuts; i++){
      Events_no_cal_hgc_aero_cuts->GetEntry(i);
      h2_events_no_cal_hgc_aero_cuts->Fill(P_no_cal_hgc_aero_cuts_hgcer_npeSum, P_no_cal_hgc_aero_cuts_aero_npeSum);
      h2_events_no_cal_aero_cuts->Fill(P_no_cal_hgc_aero_cuts_etot, P_no_cal_hgc_aero_cuts_aero_npeSum);
      h3_events_no_cal_aero_cuts->Fill(P_no_cal_hgc_aero_cuts_aero_xAtCer, P_no_cal_hgc_aero_cuts_aero_yAtCer, P_no_cal_hgc_aero_cuts_aero_npeSum);   

 }

    //Histograms for cuts + HGC cuts
    //Pions
    //1-D Histograms
    TH1D *h1_Pions_XAtCer = new TH1D("h1_Pions_XAtCer","HGC; P_hgcer_xAtCer; Events;", 300, -40, 40);
    TH1D *h1_Pions_YAtCer = new TH1D("h1_Pions_YAtCer","HGC; P_hgcer_yAtCer; Events;", 300, -40, 40);
    TH1D *h1_Pions_Prshower = new TH1D("h1_Pions_Prshower","Calorimeter; P_cal_pr_eplane; Events;", 300, 0.0, 10.0);
    TH1D *h1_Pions_Shower = new TH1D("h1_Pions_Shower","Calorimeter; P_cal_fly_earray; Events;", 300, 0.0, 10.0);
    TH1D *h1_Pions_Xgtr = new TH1D("h1_Pions_Xgtr","HGC; P_gtr_x ; Events;", 300, -3.0, 3.0);
    TH1D *h1_Pions_Ygtr = new TH1D("h1_Pions_Ygtr","HGC; P_gtr_y ; Events;", 300, -3.0, 3.0);
    TH1D *h1_Pions_P_etot = new TH1D("h1_Pions_P_etot","Calorimeter; P_cal_etotnorm; Events;", 300, 0.0, 10.0);
    TH1D *h1_Pions_gtr_beta = new TH1D("h1_Pions_gtr_beta","HGC; P_gtr_beta; Events;", 300, 0.0, 10.0);
    TH1D *h1_Pions_gtr_p = new TH1D("h1_Pions_gtr_p","HGC; P_gtr_p; Events;", 300, -10.0, 10.0);
    TH1D *h1_Pions_gtr_dp = new TH1D("h1_Pions_gtr_dp","HGC; P_gtr_dp; Events;", 300, -30.0, 30.0);
    TH1D *h1_Pions_gtr_xp = new TH1D("h1_Pions_gtr_xp","HGC; P_gtr_xp; Events;", 300, -40.0, 40.0);
    TH1D *h1_Pions_gtr_yp = new TH1D("h1_Pions_gtr_yp","HGC; P_gtr_yp; Events;", 300, -40.0, 40.0);
    TH1D *h1_Pions_CTime_ePi_ROC1 = new TH1D("h1_Pions_CTime_ePi_ROC1","CTime_ePiCoinTime_ROC1; CTime_ePiCoinTime_ROC1; Events;", 300, -10.0, 100.0);
    TH1D *h1_Pions_hgcer_npeSum = new TH1D("h1_Pions_hgcer_npeSum","HGC; P_hgcer_npeSum; Events;", 300, 0.0, 40);
    TH1D *h1_Pions_aero_npeSum = new TH1D("h1_Pions_aero_npeSum","Aero; P_aero_npeSum; Events;", 300, 0.0, 40);

     //2-D Histograms
    TH3D *h3_Pions_XYAtCer_NPE = new TH3D("h3_Pions_XYAtCer_NPE","HGC; P_hgcer_xAtCer; P_hgcer_yAtCer; P_hgcer_npeSum", 300, -40, 40, 300, -40, 40, 300, 0.0, 30);
    TH2D *h2_Pions_XYAtCer = new TH2D("h2_Pions_XYAtCer","HGC; P_hgcer_yAtCer; P_hgcer_xAtCer;", 300, -40, 40, 300, -40, 40);
    TH2D *h2_Pions_Cal_Showers = new TH2D("h2_Pions_Cal_Showers","Calorimeter; P_cal_fly_earray; P_cal_pr_eplane;", 250, 0.0, 0.8, 250, 0.0, 0.7);
    TH2D *h2_Pions_XYgtr = new TH2D("h2_Pions_XYgtr","HGC; P_gtr_x ; P_gtr_y;", 250, -0.5, 0.8, 250, -2.5, 2.5);
    TH2D *h2_Pions_npeSum = new TH2D("h2_Pions_npeSum","HGC vs Aero; P_hgcer_npeSum ; P_aero_Sum;", 300, 0, 40, 300, 0, 40);
  
    for(Long64_t i = 0; i < nEntries_Pions; i++){
      Pions->GetEntry(i);
      //2-D Histograms
      h3_Pions_XYAtCer_NPE->Fill(Pions_XAtCer, Pions_YAtCer, Pions_SHMS_hgcer_npeSum);
      h2_Pions_Cal_Showers->Fill(Pions_Cal_Shower/Pions_gtr_p, Pions_Cal_Pr_Shower/Pions_gtr_p);
      h2_Pions_XYgtr->Fill(Pions_SHMS_gtr_x, Pions_SHMS_gtr_y);
      h2_Pions_XYAtCer->Fill(Pions_YAtCer, Pions_XAtCer);
      h2_Pions_npeSum->Fill(Pions_SHMS_hgcer_npeSum, Pions_SHMS_aero_npeSum);
      // 1-D Histograms
      h1_Pions_hgcer_npeSum->Fill(Pions_SHMS_hgcer_npeSum);
      h1_Pions_aero_npeSum->Fill(Pions_SHMS_aero_npeSum);     
      h1_Pions_XAtCer->Fill(Pions_XAtCer);
      h1_Pions_YAtCer->Fill(Pions_YAtCer);
      h1_Pions_Prshower->Fill(Pions_Cal_Pr_Shower/Pions_gtr_p);
      h1_Pions_Shower->Fill(Pions_Cal_Shower/Pions_gtr_p);
      h1_Pions_Xgtr->Fill(Pions_SHMS_gtr_x);
      h1_Pions_Ygtr->Fill(Pions_SHMS_gtr_y);
      h1_Pions_gtr_beta->Fill(Pions_gtr_beta);    
      h1_Pions_gtr_p->Fill(Pions_gtr_p);
      h1_Pions_gtr_dp->Fill(Pions_gtr_dp);
      h1_Pions_gtr_xp->Fill(Pions_gtr_xp);
      h1_Pions_gtr_yp->Fill(Pions_gtr_yp);
      h1_Pions_P_etot->Fill(Pions_P_cal_etot);  
      h1_Pions_CTime_ePi_ROC1->Fill(Pions_CTime_ePi_ROC1);

    }

    //Positrons
    //1-D Histograms
    TH1D *h1_Positrons_XAtCer = new TH1D("h1_Positrons_XAtCer","HGC; P_hgcer_xAtCer; Events;", 300, -40, 40);
    TH1D *h1_Positrons_YAtCer = new TH1D("h1_Positrons_YAtCer","HGC; P_hgcer_yAtCer; Events;", 300, -40, 40);
    TH1D *h1_Positrons_Prshower = new TH1D("h1_Positrons_Prshower","Calorimeter; P_cal_pr_eplane; Events;", 300, 0.0, 10.0);
    TH1D *h1_Positrons_Shower = new TH1D("h1_Positrons_Shower","Calorimeter; P_cal_fly_earray; Events;", 300, 0.0, 10.0);
    TH1D *h1_Positrons_Xgtr = new TH1D("h1_Positrons_Xgtr","HGC; P_gtr_x ; Events;", 300, -3.0, 3.0);
    TH1D *h1_Positrons_Ygtr = new TH1D("h1_Positrons_Ygtr","HGC; P_gtr_y ; Events;", 300, -3.0, 3.0);
    TH1D *h1_Positrons_P_etot = new TH1D("h1_Positrons_P_etot","Calorimeter; P_cal_etotnorm; Events;", 300, 0.0, 10.0);
    TH1D *h1_Positrons_gtr_beta = new TH1D("h1_Positrons_gtr_beta","HGC; P_gtr_beta; Events;", 300, 0.0, 10.0);
    TH1D *h1_Positrons_gtr_p = new TH1D("h1_Positrons_gtr_p","HGC; P_gtr_p; Events;", 300, -10.0, 10.0);
    TH1D *h1_Positrons_gtr_dp = new TH1D("h1_Positrons_gtr_dp","HGC; P_gtr_dp; Events;", 300, -30.0, 30.0);
    TH1D *h1_Positrons_gtr_xp = new TH1D("h1_Positrons_gtr_xp","HGC; P_gtr_xp; Events;", 300, -40.0, 40.0);
    TH1D *h1_Positrons_gtr_yp = new TH1D("h1_Positrons_gtr_yp","HGC; P_gtr_yp; Events;", 300, -40.0, 40.0);
    // TH1D *h1_Pions_CTime_ePi_ROC1 = new TH1D("h1_Pions_CTime_ePi_ROC1","CTime_ePiCoinTime_ROC1; CTime_ePiCoinTime_ROC1; Events;", 300, -10.0, 100.0);
    TH1D *h1_Positrons_hgcer_npeSum = new TH1D("h1_Positrons_hgcer_npeSum","HGC; P_hgcer_npeSum; Events;", 300, 0.0, 40);
    TH1D *h1_Positrons_aero_npeSum = new TH1D("h1_Positrons_aero_npeSum","Aero; P_aero_npeSum; Events;", 300, 0.0, 40);

     //2-D Histograms
    TH2D *h2_Positrons_XYAtCer = new TH2D("h2_Positrons_XYAtCer","HGC; P_hgcer_yAtCer; P_hgcer_xAtCer;", 300, -40, 40, 300, -40, 40);
    TH2D *h2_Positrons_Cal_Showers = new TH2D("h2_Positrons_Cal_Showers","Calorimeter; P_cal_fly_earray; P_cal_pr_eplane;", 250, 0.0, 1.5, 250, 0.0, 0.7);
    TH2D *h2_Positrons_XYgtr = new TH2D("h2_Positrons_XYgtr","HGC; P_gtr_x ; P_gtr_y;", 250, -0.5, 0.8, 250, -2.5, 2.5);
  
    for(Long64_t i = 0; i < nEntries_Positrons; i++){
      Positrons->GetEntry(i);
      //2-D Histograms
      h2_Positrons_XYAtCer->Fill(Positrons_YAtCer, Positrons_XAtCer);
      h2_Positrons_Cal_Showers->Fill(Positrons_Cal_Shower/Positrons_gtr_p, Positrons_Cal_Pr_Shower/Positrons_gtr_p);
      h2_Positrons_XYgtr->Fill(Positrons_SHMS_gtr_x, Positrons_SHMS_gtr_y);

      // 1-D Histograms
      h1_Positrons_hgcer_npeSum->Fill(Positrons_SHMS_hgcer_npeSum);
      h1_Positrons_aero_npeSum->Fill(Positrons_SHMS_aero_npeSum);     
      h1_Positrons_XAtCer->Fill(Positrons_XAtCer);
      h1_Positrons_YAtCer->Fill(Positrons_YAtCer);
      h1_Positrons_Prshower->Fill(Positrons_Cal_Pr_Shower/Positrons_gtr_p);
      h1_Positrons_Shower->Fill(Positrons_Cal_Shower/Positrons_gtr_p);
      h1_Positrons_Xgtr->Fill(Positrons_SHMS_gtr_x);
      h1_Positrons_Ygtr->Fill(Positrons_SHMS_gtr_y);
      h1_Positrons_gtr_beta->Fill(Positrons_gtr_beta);    
      h1_Positrons_gtr_p->Fill(Positrons_gtr_p);
      h1_Positrons_gtr_dp->Fill(Positrons_gtr_dp);
      h1_Positrons_gtr_xp->Fill(Positrons_gtr_xp);
      h1_Positrons_gtr_yp->Fill(Positrons_gtr_yp);
      h1_Positrons_P_etot->Fill(Positrons_P_cal_etot);  
      // h1_Pions_CTime_ePi_ROC1->Fill(Pions_CTime_ePi_ROC1);

    }

    //Histograms for cuts + HGC cuts
    //Kaons
    TH1D *h1_Kaons_hgcer_npeSum = new TH1D("h1_Kaons_hgcer_npeSum","NPE vs Events; NPE; Events;", 300, 0.0, 40);
    TH2D *h2_Kaons_XYAtCer = new TH2D("h2_Kaons_XYAtCer","YAtCer vs XAtCer; YAtCer; XAtCer;", 200, -40, 40, 200, -40, 40);
    TH2D *h2_Kaons_Cal_Showers = new TH2D("h2_Kaons_Cal_Showers","Cal vs Pr-Shower; Cal; Pr-Shower;", 300, 0.0, 2.0, 300, 0.0, 1.0);
    TH2D *h2_Kaons_XYgtr = new TH2D("h2_Kaons_XYgtr","XY gtr; X gtr ; Y gtr;", 300, -3.0, 3.0, 300, -10, 10);

    for(Long64_t i = 0; i < nEntries_Kaons; i++){
      Kaons->GetEntry(i);
      h2_Kaons_XYAtCer->Fill(Kaons_YAtCer, Kaons_XAtCer);
      h2_Kaons_Cal_Showers->Fill(Kaons_Cal_Shower/6.053, Kaons_Cal_Pr_Shower/6.053);
      h2_Kaons_XYgtr->Fill(Kaons_SHMS_gtr_x, Kaons_SHMS_gtr_y);
      h1_Kaons_hgcer_npeSum->Fill(Kaons_SHMS_hgcer_npeSum);
    }

    //Histograms for cuts + HGC cuts
    //Protons
    TH1D *h1_Protons_hgcer_npeSum = new TH1D("h1_Protons_hgcer_npeSum","NPE vs Events; NPE; Events;", 300, 0.0, 40);
    TH2D *h2_Protons_XYAtCer = new TH2D("h2_Protons_XYAtCer","YAtCer vs XAtCer; YAtCer; XAtCer;", 200, -40, 40, 200, -40, 40);
    TH2D *h2_Protons_Cal_Showers = new TH2D("h2_Protons_Cal_Showers","Cal vs Pr-Shower; Cal; Pr-Shower;", 300, 0.0, 2.0, 300, 0.0, 1.0);
    TH2D *h2_Protons_XYgtr = new TH2D("h2_Protons_XYgtr","XY gtr; X gtr ; Y gtr;", 300, -3.0, 3.0, 300, -10, 10);

    for(Long64_t i = 0; i < nEntries_Protons; i++){
      Protons->GetEntry(i);
      h2_Protons_XYAtCer->Fill(Protons_YAtCer, Protons_XAtCer);
      h2_Protons_Cal_Showers->Fill(Protons_Cal_Shower, Protons_Cal_Pr_Shower);
      h2_Protons_XYgtr->Fill(Protons_SHMS_gtr_x, Protons_SHMS_gtr_y);
      h1_Protons_hgcer_npeSum->Fill(Protons_SHMS_hgcer_npeSum);

    }

    //Histograms for cuts + no HGC cuts
    //Pions
    //1-D Histograms
    TH1D *h1_Pi_XAtCer = new TH1D("h1_Pi_XAtCer","HGC; P_hgcer_xAtCer; Events;", 300, -40, 40);
    TH1D *h1_Pi_YAtCer = new TH1D("h1_Pi_YAtCer","HGC; P_hgcer_yAtCer; Events;", 300, -40, 40);
    TH1D *h1_Pi_Prshower = new TH1D("h1_Pi_Prshower","Calorimeter; P_cal_pr_eplane; Events;", 300, 0.0, 10.0);
    TH1D *h1_Pi_Shower = new TH1D("h1_Pi_Shower","Calorimeter; P_cal_fly_earray; Events;", 300, 0.0, 10.0);
    TH1D *h1_Pi_Xgtr = new TH1D("h1_Pi_Xgtr","HGC; P_gtr_x ; Events;", 300, -3.0, 3.0);
    TH1D *h1_Pi_Ygtr = new TH1D("h1_Pi_Ygtr","HGC; P_gtr_y ; Events;", 300, -3.0, 3.0);
    TH1D *h1_Pi_P_etot = new TH1D("h1_Pi_P_etot","Calorimeter; P_cal_etotnorm; Events;", 300, 0.0, 10.0);
    TH1D *h1_Pi_gtr_beta = new TH1D("h1_Pi_gtr_beta","HGC; P_gtr_beta; Events;", 300, 0.0, 10.0);
    TH1D *h1_Pi_gtr_p = new TH1D("h1_Pi_gtr_p","HGC; P_gtr_p; Events;", 300, -10.0, 10.0);
    TH1D *h1_Pi_gtr_dp = new TH1D("h1_Pi_gtr_dp","HGC; P_gtr_dp; Events;", 300, -30.0, 30.0);
    TH1D *h1_Pi_gtr_xp = new TH1D("h1_Pi_gtr_xp","HGC; P_gtr_xp; Events;", 300, -40.0, 40.0);
    TH1D *h1_Pi_gtr_yp = new TH1D("h1_Pi_gtr_yp","HGC; P_gtr_yp; Events;", 300, -40.0, 40.0);
    TH1D *h1_Pi_CTime_ePi_ROC1 = new TH1D("h1_Pi_CTime_ePi_ROC1","CTime_ePiCoinTime_ROC1; CTime_ePiCoinTime_ROC1; Events;", 300, -10.0, 100.0);
    TH1D *h1_Pi_hgcer_npeSum = new TH1D("h1_Pi_hgcer_npeSum","HGC; P_hgcer_npeSum; Events;", 300, 0.0, 40);
    TH1D *h1_Pi_aero_npeSum = new TH1D("h1_Pi_aero_npeSum","Aero; P_aero_npeSum; Events;", 300, 0.0, 40);

     //2-D Histograms
    TH3D *h3_Pi_XYAtCer_NPE = new TH3D("h3_Pi_XYAtCer_NPE","HGC; P_hgcer_xAtCer; P_hgcer_yAtCer; P_hgcer_npeSum", 300, -40, 40, 300, -40, 40, 300, 0.0, 30);    
    TH2D *h2_Pi_XYAtCer = new TH2D("h2_Pi_XYAtCer","HGC; P_hgcer_yAtCer; P_hgcer_xAtCer;", 300, -40, 40, 300, -40, 40);
    TH2D *h2_Pi_Cal_Showers = new TH2D("h2_Pi_Cal_Showers","Calorimeter; P_cal_fly_earray; P_cal_pr_eplane;", 250, 0.0, 0.8, 250, 0.0, 0.7);
    TH2D *h2_Pi_XYgtr = new TH2D("h2_Pi_XYgtr","HGC; P_gtr_x ; P_gtr_y;", 250, -0.5, 0.8, 250, -2.5, 2.5);
  
    for(Long64_t i = 0; i < nEntries_Pions_No_HGC_Cuts; i++){
      Pions_No_HGC_Cuts->GetEntry(i);
      //2-D Histograms
      h3_Pi_XYAtCer_NPE->Fill(Pions_No_HGC_Cuts_XAtCer, Pions_No_HGC_Cuts_YAtCer,  Pions_No_HGC_Cuts_SHMS_hgcer_npeSum);
      h2_Pi_Cal_Showers->Fill(Pions_No_HGC_Cuts_Cal_Shower/Pions_No_HGC_Cuts_gtr_p, Pions_No_HGC_Cuts_Cal_Pr_Shower/Pions_No_HGC_Cuts_gtr_p);
      h2_Pi_XYgtr->Fill(Pions_No_HGC_Cuts_SHMS_gtr_x, Pions_No_HGC_Cuts_SHMS_gtr_y);
      h2_Pi_XYAtCer->Fill(Pions_No_HGC_Cuts_YAtCer, Pions_No_HGC_Cuts_XAtCer);
      // 1-D Histograms
      h1_Pi_hgcer_npeSum->Fill(Pions_No_HGC_Cuts_SHMS_hgcer_npeSum);
      h1_Pi_aero_npeSum->Fill(Pions_No_HGC_Cuts_SHMS_aero_npeSum);     
      h1_Pi_XAtCer->Fill(Pions_No_HGC_Cuts_XAtCer);
      h1_Pi_YAtCer->Fill(Pions_No_HGC_Cuts_YAtCer);
      h1_Pi_Prshower->Fill(Pions_No_HGC_Cuts_Cal_Pr_Shower/Pions_No_HGC_Cuts_gtr_p);
      h1_Pi_Shower->Fill(Pions_No_HGC_Cuts_Cal_Shower/Pions_No_HGC_Cuts_gtr_p);
      h1_Pi_Xgtr->Fill(Pions_No_HGC_Cuts_SHMS_gtr_x);
      h1_Pi_Ygtr->Fill(Pions_No_HGC_Cuts_SHMS_gtr_y);
      h1_Pi_gtr_beta->Fill(Pions_No_HGC_Cuts_gtr_beta);    
      h1_Pi_gtr_p->Fill(Pions_No_HGC_Cuts_gtr_p);
      h1_Pi_gtr_dp->Fill(Pions_No_HGC_Cuts_gtr_dp);
      h1_Pi_gtr_xp->Fill(Pions_No_HGC_Cuts_gtr_xp);
      h1_Pi_gtr_yp->Fill(Pions_No_HGC_Cuts_gtr_yp);
      h1_Pi_P_etot->Fill(Pions_No_HGC_Cuts_P_cal_etot);  
      h1_Pi_CTime_ePi_ROC1->Fill(Pions_No_HGC_Cuts_CTime_ePi_ROC1);

    }

    //Positrons
    //1-D Histograms
    TH1D *h1_Pos_XAtCer = new TH1D("h1_Pos_XAtCer","HGC; P_hgcer_xAtCer; Events;", 300, -40, 40);
    TH1D *h1_Pos_YAtCer = new TH1D("h1_Pos_YAtCer","HGC; P_hgcer_yAtCer; Events;", 300, -40, 40);
    TH1D *h1_Pos_Prshower = new TH1D("h1_Pos_Prshower","Calorimeter; P_cal_pr_eplane; Events;", 300, 0.0, 10.0);
    TH1D *h1_Pos_Shower = new TH1D("h1_Pos_Shower","Calorimeter; P_cal_fly_earray; Events;", 300, 0.0, 10.0);
    TH1D *h1_Pos_Xgtr = new TH1D("h1_Pos_Xgtr","HGC; P_gtr_x ; Events;", 300, -3.0, 3.0);
    TH1D *h1_Pos_Ygtr = new TH1D("h1_Pos_Ygtr","HGC; P_gtr_y ; Events;", 300, -3.0, 3.0);
    TH1D *h1_Pos_P_etot = new TH1D("h1_Pos_P_etot","Calorimeter; P_cal_etotnorm; Events;", 300, 0.0, 10.0);
    TH1D *h1_Pos_gtr_beta = new TH1D("h1_Pos_gtr_beta","HGC; P_gtr_beta; Events;", 300, 0.0, 10.0);
    TH1D *h1_Pos_gtr_p = new TH1D("h1_Pos_gtr_p","HGC; P_gtr_p; Events;", 300, -10.0, 10.0);
    TH1D *h1_Pos_gtr_dp = new TH1D("h1_Pos_gtr_dp","HGC; P_gtr_dp; Events;", 300, -30.0, 30.0);
    TH1D *h1_Pos_gtr_xp = new TH1D("h1_Pos_gtr_xp","HGC; P_gtr_xp; Events;", 300, -40.0, 40.0);
    TH1D *h1_Pos_gtr_yp = new TH1D("h1_Pos_gtr_yp","HGC; P_gtr_yp; Events;", 300, -40.0, 40.0);
    // TH1D *h1_Pos_CTime_ePi_ROC1 = new TH1D("h1_Pi_CTime_ePi_ROC1","CTime_ePiCoinTime_ROC1; CTime_ePiCoinTime_ROC1; Events;", 300, -10.0, 100.0);
    TH1D *h1_Pos_hgcer_npeSum = new TH1D("h1_Pos_hgcer_npeSum","HGC; P_hgcer_npeSum; Events;", 300, 0.0, 40);
    TH1D *h1_Pos_aero_npeSum = new TH1D("h1_Pos_aero_npeSum","Aero; P_aero_npeSum; Events;", 300, 0.0, 40);

     //2-D Histograms
    TH2D *h2_Pos_XYAtCer = new TH2D("h2_Pos_XYAtCer","HGC; P_hgcer_yAtCer; P_hgcer_xAtCer;", 300, -40, 40, 300, -40, 40);
    TH2D *h2_Pos_Cal_Showers = new TH2D("h2_Pos_Cal_Showers","Calorimeter; P_cal_fly_earray; P_cal_pr_eplane;", 250, 0.0, 1.5, 250, 0.0, 0.7);
    TH2D *h2_Pos_XYgtr = new TH2D("h2_Pos_XYgtr","HGC; P_gtr_x ; P_gtr_y;", 250, -0.5, 0.8, 250, -2.5, 2.5);
  
    for(Long64_t i = 0; i < nEntries_Positrons_No_HGC_Cuts; i++){
      Positrons_No_HGC_Cuts->GetEntry(i);
      //2-D Histograms
      h2_Pos_XYAtCer->Fill(Positrons_No_HGC_Cuts_YAtCer, Positrons_No_HGC_Cuts_XAtCer);
      h2_Pos_Cal_Showers->Fill(Positrons_No_HGC_Cuts_Cal_Shower/Positrons_No_HGC_Cuts_gtr_p, Positrons_No_HGC_Cuts_Cal_Pr_Shower/Positrons_No_HGC_Cuts_gtr_p);
      h2_Pos_XYgtr->Fill(Positrons_No_HGC_Cuts_SHMS_gtr_x, Positrons_No_HGC_Cuts_SHMS_gtr_y);

      // 1-D Histograms
      h1_Pos_hgcer_npeSum->Fill(Positrons_No_HGC_Cuts_SHMS_hgcer_npeSum);
      h1_Pos_aero_npeSum->Fill(Positrons_No_HGC_Cuts_SHMS_aero_npeSum);     
      h1_Pos_XAtCer->Fill(Positrons_No_HGC_Cuts_XAtCer);
      h1_Pos_YAtCer->Fill(Positrons_No_HGC_Cuts_YAtCer);
      h1_Pos_Prshower->Fill(Positrons_No_HGC_Cuts_Cal_Pr_Shower/Positrons_No_HGC_Cuts_gtr_p);
      h1_Pos_Shower->Fill(Positrons_No_HGC_Cuts_Cal_Shower/Positrons_No_HGC_Cuts_gtr_p);
      h1_Pos_Xgtr->Fill(Positrons_No_HGC_Cuts_SHMS_gtr_x);
      h1_Pos_Ygtr->Fill(Positrons_No_HGC_Cuts_SHMS_gtr_y);
      h1_Pos_gtr_beta->Fill(Positrons_No_HGC_Cuts_gtr_beta);    
      h1_Pos_gtr_p->Fill(Positrons_No_HGC_Cuts_gtr_p);
      h1_Pos_gtr_dp->Fill(Positrons_No_HGC_Cuts_gtr_dp);
      h1_Pos_gtr_xp->Fill(Positrons_No_HGC_Cuts_gtr_xp);
      h1_Pos_gtr_yp->Fill(Positrons_No_HGC_Cuts_gtr_yp);
      h1_Pos_P_etot->Fill(Positrons_No_HGC_Cuts_P_cal_etot);  
      //h1_Pi_CTime_ePi_ROC1->Fill(Pions_No_HGC_Cuts_CTime_ePi_ROC1);

    }

    TH1D *h1_Kaons_No_HGC_Cuts_hgcer_npeSum = new TH1D("h1_Kaons_No_HGC_Cuts_hgcer_npeSum","NPE vs Events; NPE; Events;", 300, 0.0, 40);
    TH2D *h2_Kaons_No_HGC_Cuts_XYAtCer = new TH2D("h2_Kaons_No_HGC_Cuts_XYAtCer","YAtCer vs XAtCer; YAtCer; XAtCer;", 200, -40, 40, 200, -40, 40);
    TH2D *h2_Kaons_No_HGC_Cuts_Cal_Showers = new TH2D("h2_Kaons_No_HGC_Cuts_Cal_Showers","Cal vs Pr-Shower; Cal; Pr-Shower;", 300, 0.0, 2.0, 300, 0.0, 1.0);
    TH2D *h2_Kaons_No_HGC_Cuts_XYgtr = new TH2D("h2_Kaons_No_HGC_Cuts_XYgtr","XY gtr; X gtr ; Y gtr;", 300, -3.0, 3.0, 300, -10, 10);
    
    for(Long64_t i = 0; i < nEntries_Kaons_No_HGC_Cuts; i++){
      Kaons_No_HGC_Cuts->GetEntry(i);
      h2_Kaons_No_HGC_Cuts_XYAtCer->Fill(Kaons_No_HGC_Cuts_YAtCer, Kaons_No_HGC_Cuts_XAtCer);
      h2_Kaons_No_HGC_Cuts_Cal_Showers->Fill(Kaons_No_HGC_Cuts_Cal_Shower/6.053, Kaons_No_HGC_Cuts_Cal_Pr_Shower/6.053);
      h2_Kaons_No_HGC_Cuts_XYgtr->Fill(Kaons_No_HGC_Cuts_SHMS_gtr_x, Kaons_No_HGC_Cuts_SHMS_gtr_y);
      h1_Kaons_No_HGC_Cuts_hgcer_npeSum->Fill(Kaons_No_HGC_Cuts_SHMS_hgcer_npeSum);

    }
    TH1D *h1_Protons_No_HGC_Cuts_hgcer_npeSum = new TH1D("h1_Protons_No_HGC_Cuts_hgcer_npeSum","NPE vs Events; NPE; Events;", 300, 0.0, 40);
    TH2D *h2_Protons_No_HGC_Cuts_XYAtCer = new TH2D("h2_Protons_No_HGC_Cuts_XYAtCer","YAtCer vs XAtCer; YAtCer; XAtCer;", 200, -40, 40, 200, -40, 40);
    TH2D *h2_Protons_No_HGC_Cuts_Cal_Showers = new TH2D("h2_Protons_No_HGC_Cuts_Cal_Showers","Cal vs Pr-Shower; Cal; Pr-Shower;", 300, 0.0, 2.0, 300, 0.0, 1.0);
    TH2D *h2_Protons_No_HGC_Cuts_XYgtr = new TH2D("h2_Protons_No_HGC_Cuts_XYgtr","XY gtr; X gtr ; Y gtr;", 300, -3.0, 3.0, 300, -10, 10);

    for(Long64_t i = 0; i < nEntries_Protons_No_HGC_Cuts; i++){
      Protons_No_HGC_Cuts->GetEntry(i);
      h2_Protons_No_HGC_Cuts_XYAtCer->Fill(Protons_No_HGC_Cuts_YAtCer, Protons_No_HGC_Cuts_XAtCer);
      h2_Protons_No_HGC_Cuts_Cal_Showers->Fill(Protons_No_HGC_Cuts_Cal_Shower/6.053, Protons_No_HGC_Cuts_Cal_Pr_Shower/6.053);
      h2_Protons_No_HGC_Cuts_XYgtr->Fill(Protons_No_HGC_Cuts_SHMS_gtr_x, Protons_No_HGC_Cuts_SHMS_gtr_y);
      h1_Protons_No_HGC_Cuts_hgcer_npeSum->Fill(Kaons_No_HGC_Cuts_SHMS_hgcer_npeSum);

    }
    TH2D *h = new TH2D("h","XY gtr; X gtr ; Y gtr;", 300, 0.0, 30.0, 300, 0, 30);
    for(Long64_t i = 0; i < nEntries_Events_no_cal_hgc_aero_cuts; i++){
      Events_no_cal_hgc_aero_cuts->GetEntry(i);
      if (!cutg2->IsInside(P_no_cal_hgc_aero_cuts_hgcer_yAtCer, P_no_cal_hgc_aero_cuts_hgcer_xAtCer)) continue;
      h->Fill(P_no_cal_hgc_aero_cuts_hgcer_npeSum, P_no_cal_hgc_aero_cuts_aero_npeSum);
    }
      
    //################################################################################################################################################
    gROOT->SetBatch(kTRUE); // Set ROOT to batch mode explicitly, does not splash anything to screen
    //################################################################################################################################################

    TCanvas *c_CT = new TCanvas("c_CT", "HGC (with TCutG)");  
    c_CT->Divide(2,2);   
    c_CT->cd(1);
    Events_no_cal_hgc_aero_cuts->Draw("P_hgcer_npeSum:P_aero_npeSum>>h1(300,0.0,30,300,0,30)", "cutg",  "colz");
    c_CT->cd(2);
    Events_no_cal_hgc_aero_cuts->Draw("P_hgcer_npeSum>>h2(300,0.3,30)", "cutg",  "colz");
    c_CT->cd(3);
    Events_no_cal_hgc_aero_cuts->Draw("P_hgcer_npeSum:P_aero_npeSum>>h3(300,0,30, 300, 0, 30)", "!cutg",  "colz"); 
    c_CT->cd(4);
    Events_no_cal_hgc_aero_cuts->Draw("P_hgcer_npeSum>>h4(300,0.3,30)", "!cutg",  "colz"); 
    c_CT->Print(foutpdf);

    //################################################################################################################################################

    //  TProfile2D h3_Pions_XYAtCer_NPE = Project3DProfile("xy") 
    TFile *OutHisto_file = new TFile(foutname,"RECREATE");
    TDirectory *Pions_info = OutHisto_file->mkdir("Pions_info");
    Pions_info->cd();

    // TH3D *Pions_XYAtCer_NPE;
    // Pions_XYAtCer_NPE = dynamic_cast<TH3D*> (GetOutputList()->FindObject("h3_Pions_XYAtCer_NPE"));
    TProfile2D *h3_Pions_XYAtCer_NPE_pxy = new TProfile2D("h3_Pions_XYAtCer_NPE_pxy","NPE vs X vs Y; X ; Y ",300,-40,40, 300,-40,40,0.0,40);
    h3_Pions_XYAtCer_NPE->Project3DProfile("xy");
    
    TProfile2D *h3_Pi_XYAtCer_NPE_pxy = new TProfile2D("h3_Pi_XYAtCer_NPE_pxy","NPE vs X vs Y; X ; Y ",300,-40,40, 300,-40,40,0.0,40);
    h3_Pi_XYAtCer_NPE->Project3DProfile("xy");
   
    TProfile2D *h3_events_no_cal_aero_cuts_pxy = new TProfile2D("h3_events_no_cal_aero_cuts_pxy","NPE vs X vs Y; X ; Y ",300,-50,50, 300,-50,50,0.0,30);
    h3_events_no_cal_aero_cuts->Project3DProfile("xy");

    //2-D Histograms
    h3_Pions_XYAtCer_NPE_pxy->GetListOfFunctions()->Add(cutg,"L"); 
    h3_Pions_XYAtCer_NPE_pxy->Write();
    //Pions_XYAtCer_NPE->Write();
    // h3_Pions_XYAtCer_NPE->Write();
    h2_Pions_Cal_Showers->Write();
    h2_Pions_XYgtr->Write();
    h2_Pions_XYAtCer->GetListOfFunctions()->Add(cutg,"L"); 
    h2_Pions_XYAtCer->Write();
    //    h2_Pions_npeSum->GetListOfFunctions()->Add(cutg1,"L");
    h2_Pions_npeSum->Write();
    //1-D Histograms
    h1_Pions_hgcer_npeSum->Write();
    h1_Pions_aero_npeSum->Write();   
    h1_Pions_XAtCer->Write();
    h1_Pions_YAtCer->Write();
    h1_Pions_Shower->Write();
    h1_Pions_Prshower->Write();
    h1_Pions_Xgtr->Write();
    h1_Pions_Ygtr->Write();
    h1_Pions_P_etot->Write();
    h1_Pions_gtr_beta->Write();
    h1_Pions_gtr_p->Write();
    h1_Pions_gtr_dp->Write();
    h1_Pions_gtr_xp->Write();  
    h1_Pions_gtr_yp->Write();  
    h1_Pions_CTime_ePi_ROC1->Write();  

    TDirectory *Positrons_info = OutHisto_file->mkdir("Positrons_info");

    Positrons_info->cd();
    //2-D Histograms
    h2_Positrons_XYAtCer->Write();
    h2_Positrons_Cal_Showers->Write();
    h2_Positrons_XYgtr->Write();

    //1-D Histograms
    h1_Positrons_hgcer_npeSum->Write();
    h1_Positrons_aero_npeSum->Write();   
    h1_Positrons_XAtCer->Write();
    h1_Positrons_YAtCer->Write();
    h1_Positrons_Shower->Write();
    h1_Positrons_Prshower->Write();
    h1_Positrons_Xgtr->Write();
    h1_Positrons_Ygtr->Write();
    h1_Positrons_P_etot->Write();
    h1_Positrons_gtr_beta->Write();
    h1_Positrons_gtr_p->Write();
    h1_Positrons_gtr_dp->Write();
    h1_Positrons_gtr_xp->Write();  
    h1_Positrons_gtr_yp->Write();  
    //h1_Pions_CTime_ePi_ROC1->Write();  

    
    TDirectory *Kaons_info = OutHisto_file->mkdir("Kaons_info");
    
    Kaons_info->cd();
    h2_Kaons_XYAtCer->Write();
    h2_Kaons_Cal_Showers->Write();
    h2_Kaons_XYgtr->Write();
    h1_Kaons_hgcer_npeSum->Write();

    TDirectory *Protons_info = OutHisto_file->mkdir("Protons_info");
    
    Protons_info->cd();
    h2_Protons_XYAtCer->Write();
    h2_Protons_Cal_Showers->Write();
    h2_Protons_XYgtr->Write();
    h1_Protons_hgcer_npeSum->Write();

    TDirectory *Pions_No_HGC_Cuts_info = OutHisto_file->mkdir("Pions_No_HGC_Cuts_info");
    
    Pions_No_HGC_Cuts_info->cd();

    //2-D Histograms
    h3_Pi_XYAtCer_NPE_pxy->GetListOfFunctions()->Add(cutg,"L");       
    h3_Pi_XYAtCer_NPE_pxy->Write();
    h2_Pi_Cal_Showers->Write();
    h2_Pi_XYgtr->Write();
    h2_Pi_XYAtCer->GetListOfFunctions()->Add(cutg,"L"); 
    h2_Pi_XYAtCer->Write();
    //1-D Histograms
    h1_Pi_hgcer_npeSum->Write();
    h1_Pi_aero_npeSum->Write();   
    h1_Pi_XAtCer->Write();
    h1_Pi_YAtCer->Write();
    h1_Pi_Shower->Write();
    h1_Pi_Prshower->Write();
    h1_Pi_Xgtr->Write();
    h1_Pi_Ygtr->Write();
    h1_Pi_P_etot->Write();
    h1_Pi_gtr_beta->Write();
    h1_Pi_gtr_p->Write();
    h1_Pi_gtr_dp->Write();
    h1_Pi_gtr_xp->Write();  
    h1_Pi_gtr_yp->Write();  
    h1_Pi_CTime_ePi_ROC1->Write();  

    TDirectory *Positrons_No_HGC_Cuts_info = OutHisto_file->mkdir("Positrons_No_HGC_Cuts_info");
    
    Positrons_No_HGC_Cuts_info->cd();

    //2-D Histograms
    h2_Pos_XYAtCer->Write();
    h2_Pos_Cal_Showers->Write();
    h2_Pos_XYgtr->Write();

    //1-D Histograms
    h1_Pos_hgcer_npeSum->Write();
    h1_Pos_aero_npeSum->Write();   
    h1_Pos_XAtCer->Write();
    h1_Pos_YAtCer->Write();
    h1_Pos_Shower->Write();
    h1_Pos_Prshower->Write();
    h1_Pos_Xgtr->Write();
    h1_Pos_Ygtr->Write();
    h1_Pos_P_etot->Write();
    h1_Pos_gtr_beta->Write();
    h1_Pos_gtr_p->Write();
    h1_Pos_gtr_dp->Write();
    h1_Pos_gtr_xp->Write();  
    h1_Pos_gtr_yp->Write();  
    //h1_Pos_CTime_ePi_ROC1->Write();  

   
    TDirectory *Kaons_No_HGC_Cuts_info = OutHisto_file->mkdir("Kaons_No_HGC_Cuts_info");
    
    Kaons_No_HGC_Cuts_info->cd();
    h2_Kaons_No_HGC_Cuts_XYAtCer->Write();
    h2_Kaons_No_HGC_Cuts_Cal_Showers->Write();
    h2_Kaons_No_HGC_Cuts_XYgtr->Write();
    h1_Kaons_No_HGC_Cuts_hgcer_npeSum->Write();
    TDirectory *Protons_No_HGC_Cuts_info = OutHisto_file->mkdir("Protons_No_HGC_Cuts_info");
    
    Protons_No_HGC_Cuts_info->cd();
    h2_Protons_No_HGC_Cuts_XYAtCer->Write();
    h2_Protons_No_HGC_Cuts_Cal_Showers->Write();
    h2_Protons_No_HGC_Cuts_XYgtr->Write();
    h1_Protons_No_HGC_Cuts_hgcer_npeSum->Write();

 TDirectory *SHMS_Events_No_Cuts = OutHisto_file->mkdir("SHMS_Events_No_Cuts");
    
    SHMS_Events_No_Cuts->cd();
    h1_XAtCer->Write();
    h1_YAtCer->Write();
    h1_Cal->Write();
    h1_Prshower->Write();
    h1_Xgtr->Write();
    h1_Ygtr->Write();
    h1_hgcer_npeSum->Write();
    h1_P_etot->Write();
    h1_gtr_beta->Write();
    h1_gtr_p->Write();
    h1_gtr_dp->Write();
    h1_gtr_xp->Write();  
    h1_gtr_yp->Write();  
    h1_CTime_ePi_ROC1->Write();  
    h1_CTime_eK_ROC1->Write();  
    h1_CTime_eP_ROC1->Write();  

    TDirectory *SHMS_Events_No_Aero_Cuts = OutHisto_file->mkdir("SHMS_Events_No_Aero_Cuts");
    
    SHMS_Events_No_Aero_Cuts->cd();
     
    h1_aero_npeSum->Write();

    TDirectory *SHMS_Events_No_Cal_HGC_Aero_Cuts = OutHisto_file->mkdir("SHMS_Events_No_Cal_HGC_Aero_Cuts");
    SHMS_Events_No_Cal_HGC_Aero_Cuts->cd();
     
    h1_cal_etot->Write();
    h2_events_no_cal_hgc_cuts->Write();
    // h2_events_no_cal_hgc_aero_cuts->GetListOfFunctions()->Add(cutg,"L"); 
    //h2_events_no_cal_hgc_aero_cuts->Write();
    h2_events_no_cal_aero_cuts->Write();
    h3_events_no_cal_aero_cuts_pxy->Write();
    h->Write(); 
    OutHisto_file->Close();

    /* TCanvas *c_CT = new TCanvas("c_CT", "HGC vs Aero (with TCutG)");  
    c_CT->Divide(2);   
    c_CT->cd(1);
    Events_no_cal_hgc_aero_cuts->Draw("P_hgcer_npeSum:P_aero_npeSum>>h1(300,0,30,300,0,30)", "cutg",  "colz");
    c_CT->Print(foutpdf);
    // c_CT->SaveAs("hgc_aer_TCutG.png");
    c_CT->cd(2);
    Events_no_cal_hgc_aero_cuts->Draw("P_hgcer_npeSum:P_aero_npeSum>>h2(300,0,30,300,0,30)", "!cutg",  "colz"); 
    c_CT->Print(foutpdf);*/
    

    //TString RunNumStr = TInFilename(0,4); Int_t RunNum=(RunNumStr.Atoi());
    //TString OutputStr = Form("%i,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f", RunNum, PionFit->GetParameter(1), PionFit->GetParError(1), PionFWHM, PionFWHMErr, KaonFit->GetParameter(1), KaonFit->GetParError(1), KaonFWHM, KaonFWHMErr, ProtonFit->GetParameter(1), ProtonFit->GetParError(1), ProtonFWHM, ProtonFWHMErr);
    //cout << OutputStr << endl;
}
