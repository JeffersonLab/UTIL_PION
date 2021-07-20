// 15/01/21 - Stephen Kay, University of Regina
// 03/06/21 - Edited by Muhammad Junaid, University of Regina. Canada

// root .c macro plotting script, reads in desired trees from analysed root file and plots some stuff
// Saves  pdf file with plots and a .root file

#define PlotPionPhysics_cxx

//////////////////////////////////////////////////////////////////////////////////////////////////

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

/////////////////////////////////////////////////////////////////////////////////////////////////

// Input should be the input root file name (including suffix) and an output file name string (without any suffix)
void PlotPionPhysics(string InFilename = "", string OutFilename = "")
{
  TString Hostname = gSystem->HostName();
  TString User = (gSystem->GetUserInfo())->fUser;
  TString Replaypath;
  TString Outpath;
  TString rootFile;
  Double_t nWindows = 6;
  gStyle->SetPalette(55);

/////////////////////////////////////////////////////////////////////////////////////////////////

  // Set paths depending on system you're running on
  if(Hostname.Contains("farm")){
    Replaypath = "/group/c-kaonlt/USERS/"+User+"/hallc_replay_lt";
    Outpath = Replaypath+"/UTIL_PION/scripts/demo/OUTPUT";
  }
  else if(Hostname.Contains("qcd")){
    Replaypath = "/group/c-kaonlt/USERS/"+User+"/hallc_replay_lt";
    Outpath = Replaypath+"/UTIL_PION/scripts/pionyield/OUTPUT";
  }
  else if (Hostname.Contains("phys.uregina.ca")){
    Replaypath = "/home/"+User+"/work/JLab/hallc_replay_lt";
    Outpath = Replaypath+"/UTIL_PION/scripts/pionyield/OUTPUT";
  }

////////////////////////////////////////////////////////////////////////////////////////////////

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
  TFile *InFile = new TFile(rootFile, "READ");

////////////////////////////////////////////////////////////////////////////////////////////////

  // This assumes a 4 digit run number! Should work for now, saves it as an additional input
  // This ALSO assumes the Filename is XXXX_..... too, which may not be true, edit as needed
  TString TOutFilename = OutFilename;
  // Establish the names of our output files quickly
  TString foutname = Outpath + "/" + TOutFilename + ".root";
//  TString foutpdf = Outpath + "/" + TOutFilename + ".pdf";
  
//////////////////////////////////////////////////////////////////////////////////////////////

  TTree* Uncut_Pion_Events = (TTree*)InFile->Get("Uncut_Pion_Events"); Long64_t nEntries_Uncut_Pion_Events = (Long64_t)Uncut_Pion_Events->GetEntries();
  TTree* Cut_Pion_Events_All = (TTree*)InFile->Get("Cut_Pion_Events_All"); Long64_t nEntries_Cut_Pion_Events_All = (Long64_t)Cut_Pion_Events_All->GetEntries();
  TTree* Cut_Pion_Events_Prompt = (TTree*)InFile->Get("Cut_Pion_Events_Prompt"); Long64_t nEntries_Cut_Pion_Events_Prompt = (Long64_t)Cut_Pion_Events_Prompt->GetEntries();
  TTree* Cut_Pion_Events_Random = (TTree*)InFile->Get("Cut_Pion_Events_Random"); Long64_t nEntries_Cut_Pion_Events_Random = (Long64_t)Cut_Pion_Events_Random->GetEntries();

  TTree* Uncut_Kaon_Events = (TTree*)InFile->Get("Uncut_Kaon_Events"); Long64_t nEntries_Uncut_Kaon_Events = (Long64_t)Uncut_Kaon_Events->GetEntries();
  TTree* Cut_Kaon_Events_All = (TTree*)InFile->Get("Cut_Kaon_Events_All"); Long64_t nEntries_Cut_Kaon_Events_All = (Long64_t)Cut_Kaon_Events_All->GetEntries();
  TTree* Cut_Kaon_Events_Prompt = (TTree*)InFile->Get("Cut_Kaon_Events_Prompt"); Long64_t nEntries_Cut_Kaon_Events_Prompt = (Long64_t)Cut_Kaon_Events_Prompt->GetEntries();
  TTree* Cut_Kaon_Events_Random = (TTree*)InFile->Get("Cut_Kaon_Events_Random"); Long64_t nEntries_Cut_Kaon_Events_Random = (Long64_t)Cut_Kaon_Events_Random->GetEntries();

  TTree* Uncut_Proton_Events = (TTree*)InFile->Get("Uncut_Proton_Events"); Long64_t nEntries_Uncut_Proton_Events = (Long64_t)Uncut_Proton_Events->GetEntries();
  TTree* Cut_Proton_Events_All = (TTree*)InFile->Get("Cut_Proton_Events_All"); Long64_t nEntries_Cut_Proton_Events_All = (Long64_t)Cut_Proton_Events_All->GetEntries();
  TTree* Cut_Proton_Events_Prompt = (TTree*)InFile->Get("Cut_Proton_Events_Prompt"); Long64_t nEntries_Cut_Proton_Events_Prompt = (Long64_t)Cut_Proton_Events_Prompt->GetEntries();
  TTree* Cut_Proton_Events_Random = (TTree*)InFile->Get("Cut_Proton_Events_Random"); Long64_t nEntries_Cut_Proton_Events_Random = (Long64_t)Cut_Proton_Events_Random->GetEntries();

///////////////////////////////////////////////////////////////////////////////////////////////

// Set branch address -> Need this to ensure event info is entangled correctly for 2D plots

// For Pions 
  Double_t H_gtr_beta_pions_uncut; Uncut_Pion_Events->SetBranchAddress("H_gtr_beta", &H_gtr_beta_pions_uncut);
  Double_t H_gtr_xp_pions_uncut; Uncut_Pion_Events->SetBranchAddress("H_gtr_xp", &H_gtr_xp_pions_uncut);
  Double_t H_gtr_yp_pions_uncut; Uncut_Pion_Events->SetBranchAddress("H_gtr_yp", &H_gtr_yp_pions_uncut);
  Double_t H_gtr_dp_pions_uncut; Uncut_Pion_Events->SetBranchAddress("H_gtr_dp", &H_gtr_dp_pions_uncut);
  Double_t H_hod_goodscinhit_pions_uncut; Uncut_Pion_Events->SetBranchAddress("H_hod_goodscinhit", &H_hod_goodscinhit_pions_uncut);
  Double_t H_hod_goodstarttime_pions_uncut; Uncut_Pion_Events->SetBranchAddress("H_hod_goodstarttime", &H_hod_goodstarttime_pions_uncut);
  Double_t H_cal_etotnorm_pions_uncut; Uncut_Pion_Events->SetBranchAddress("H_cal_etotnorm", &H_cal_etotnorm_pions_uncut);
  Double_t H_cal_etottracknorm_pions_uncut; Uncut_Pion_Events->SetBranchAddress("H_cal_etottracknorm", &H_cal_etottracknorm_pions_uncut);
  Double_t H_cer_npeSum_pions_uncut; Uncut_Pion_Events->SetBranchAddress("H_cer_npeSum", &H_cer_npeSum_pions_uncut);
  Double_t H_RFTime_Dist_pions_uncut; Uncut_Pion_Events->SetBranchAddress("H_RFTime_Dist", &H_RFTime_Dist_pions_uncut);
  Double_t P_gtr_beta_pions_uncut; Uncut_Pion_Events->SetBranchAddress("P_gtr_beta", &P_gtr_beta_pions_uncut);
  Double_t P_gtr_xp_pions_uncut; Uncut_Pion_Events->SetBranchAddress("P_gtr_xp", &P_gtr_xp_pions_uncut);
  Double_t P_gtr_yp_pions_uncut; Uncut_Pion_Events->SetBranchAddress("P_gtr_yp", &P_gtr_yp_pions_uncut);
  Double_t P_gtr_dp_pions_uncut; Uncut_Pion_Events->SetBranchAddress("P_gtr_dp", &P_gtr_dp_pions_uncut);
  Double_t P_gtr_p_pions_uncut; Uncut_Pion_Events->SetBranchAddress("P_gtr_p", &P_gtr_p_pions_uncut);
  Double_t P_hod_goodscinhit_pions_uncut; Uncut_Pion_Events->SetBranchAddress("P_hod_goodscinhit", &P_hod_goodscinhit_pions_uncut);
  Double_t P_hod_goodstarttime_pions_uncut; Uncut_Pion_Events->SetBranchAddress("P_hod_goodstarttime", &P_hod_goodstarttime_pions_uncut);
  Double_t P_cal_etotnorm_pions_uncut; Uncut_Pion_Events->SetBranchAddress("P_cal_etotnorm", &P_cal_etotnorm_pions_uncut);
  Double_t P_cal_etottracknorm_pions_uncut; Uncut_Pion_Events->SetBranchAddress("P_cal_etottracknorm", &P_cal_etottracknorm_pions_uncut);
  Double_t P_hgcer_npeSum_pions_uncut; Uncut_Pion_Events->SetBranchAddress("P_hgcer_npeSum", &P_hgcer_npeSum_pions_uncut);
  Double_t P_hgcer_xAtCer_pions_uncut; Uncut_Pion_Events->SetBranchAddress("P_hgcer_xAtCer", &P_hgcer_xAtCer_pions_uncut);
  Double_t P_hgcer_yAtCer_pions_uncut; Uncut_Pion_Events->SetBranchAddress("P_hgcer_yAtCer", &P_hgcer_yAtCer_pions_uncut);
  Double_t P_aero_npeSum_pions_uncut; Uncut_Pion_Events->SetBranchAddress("P_aero_npeSum", &P_aero_npeSum_pions_uncut);
  Double_t P_aero_xAtAero_pions_uncut; Uncut_Pion_Events->SetBranchAddress("P_aero_xAtAero", &P_aero_xAtAero_pions_uncut);
  Double_t P_aero_yAtAero_pions_uncut; Uncut_Pion_Events->SetBranchAddress("P_aero_yAtAero", &P_aero_yAtAero_pions_uncut);
  Double_t P_kin_MMpi_pions_uncut; Uncut_Pion_Events->SetBranchAddress("MMpi", &P_kin_MMpi_pions_uncut);
  Double_t P_RFTime_Dist_pions_uncut; Uncut_Pion_Events->SetBranchAddress("P_RF_Dist", &P_RFTime_Dist_pions_uncut);
  Double_t CTime_ePiCoinTime_ROC1_pions_uncut; Uncut_Pion_Events->SetBranchAddress("CTime_ePiCoinTime_ROC1", &CTime_ePiCoinTime_ROC1_pions_uncut);
  Double_t Q2_pions_uncut; Uncut_Pion_Events->SetBranchAddress("Q2", &Q2_pions_uncut);
  Double_t W_pions_uncut; Uncut_Pion_Events->SetBranchAddress("W", &W_pions_uncut);
  Double_t epsilon_pions_uncut; Uncut_Pion_Events->SetBranchAddress("Q2", &epsilon_pions_uncut);

  Double_t H_gtr_beta_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("H_gtr_beta", &H_gtr_beta_pions_cut_all);
  Double_t H_gtr_xp_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("H_gtr_xp", &H_gtr_xp_pions_cut_all);
  Double_t H_gtr_yp_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("H_gtr_yp", &H_gtr_yp_pions_cut_all);
  Double_t H_gtr_dp_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("H_gtr_dp", &H_gtr_dp_pions_cut_all);
  Double_t H_hod_goodscinhit_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("H_hod_goodscinhit", &H_hod_goodscinhit_pions_cut_all);
  Double_t H_hod_goodstarttime_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("H_hod_goodstarttime", &H_hod_goodstarttime_pions_cut_all);
  Double_t H_cal_etotnorm_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("H_cal_etotnorm", &H_cal_etotnorm_pions_cut_all);
  Double_t H_cal_etottracknorm_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("H_cal_etottracknorm", &H_cal_etottracknorm_pions_cut_all);
  Double_t H_cer_npeSum_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("H_cer_npeSum", &H_cer_npeSum_pions_cut_all);
  Double_t H_RFTime_Dist_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("H_RFTime_Dist", &H_RFTime_Dist_pions_cut_all);
  Double_t P_gtr_beta_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("P_gtr_beta", &P_gtr_beta_pions_cut_all);
  Double_t P_gtr_xp_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("P_gtr_xp", &P_gtr_xp_pions_cut_all);
  Double_t P_gtr_yp_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("P_gtr_yp", &P_gtr_yp_pions_cut_all);
  Double_t P_gtr_dp_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("P_gtr_dp", &P_gtr_dp_pions_cut_all);
  Double_t P_gtr_p_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("P_gtr_p", &P_gtr_p_pions_cut_all);
  Double_t P_hod_goodscinhit_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("P_hod_goodscinhit", &P_hod_goodscinhit_pions_cut_all);
  Double_t P_hod_goodstarttime_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("P_hod_goodstarttime", &P_hod_goodstarttime_pions_cut_all);
  Double_t P_cal_etotnorm_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("P_cal_etotnorm", &P_cal_etotnorm_pions_cut_all);
  Double_t P_cal_etottracknorm_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("P_cal_etottracknorm", &P_cal_etottracknorm_pions_cut_all);
  Double_t P_hgcer_npeSum_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("P_hgcer_npeSum", &P_hgcer_npeSum_pions_cut_all);
  Double_t P_hgcer_xAtCer_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("P_hgcer_xAtCer", &P_hgcer_xAtCer_pions_cut_all);
  Double_t P_hgcer_yAtCer_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("P_hgcer_yAtCer", &P_hgcer_yAtCer_pions_cut_all);
  Double_t P_aero_npeSum_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("P_aero_npeSum", &P_aero_npeSum_pions_cut_all);
  Double_t P_aero_xAtAero_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("P_aero_xAtAero", &P_aero_xAtAero_pions_cut_all);
  Double_t P_aero_yAtAero_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("P_aero_yAtAero", &P_aero_yAtAero_pions_cut_all);
  Double_t P_kin_MMpi_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("MMpi", &P_kin_MMpi_pions_cut_all);
  Double_t P_RFTime_Dist_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("P_RF_Dist", &P_RFTime_Dist_pions_cut_all);
  Double_t CTime_ePiCoinTime_ROC1_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("CTime_ePiCoinTime_ROC1", &CTime_ePiCoinTime_ROC1_pions_cut_all);
  Double_t Q2_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("Q2", &Q2_pions_cut_all);
  Double_t W_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("W", &W_pions_cut_all);
  Double_t epsilon_pions_cut_all; Cut_Pion_Events_All->SetBranchAddress("Q2", &epsilon_pions_cut_all);


  Double_t P_gtr_beta_pions_cut_pr; Cut_Pion_Events_Prompt->SetBranchAddress("P_gtr_beta", &P_gtr_beta_pions_cut_pr);
  Double_t P_RFTime_Dist_pions_cut_pr; Cut_Pion_Events_Prompt->SetBranchAddress("P_RF_Dist", &P_RFTime_Dist_pions_cut_pr);
  Double_t CTime_ePiCoinTime_ROC1_pions_cut_pr; Cut_Pion_Events_Prompt->SetBranchAddress("CTime_ePiCoinTime_ROC1", &CTime_ePiCoinTime_ROC1_pions_cut_pr);
  Double_t P_kin_MMpi_pions_cut_pr; Cut_Pion_Events_Prompt->SetBranchAddress("MMpi", &P_kin_MMpi_pions_cut_pr);

  Double_t P_gtr_beta_pions_cut_rn; Cut_Pion_Events_Random->SetBranchAddress("P_gtr_beta", &P_gtr_beta_pions_cut_rn);
  Double_t P_RFTime_Dist_pions_cut_rn; Cut_Pion_Events_Random->SetBranchAddress("P_RF_Dist", &P_RFTime_Dist_pions_cut_rn);
  Double_t CTime_ePiCoinTime_ROC1_pions_cut_rn; Cut_Pion_Events_Random->SetBranchAddress("CTime_ePiCoinTime_ROC1", &CTime_ePiCoinTime_ROC1_pions_cut_rn);
  Double_t P_kin_MMpi_pions_cut_rn; Cut_Pion_Events_Random->SetBranchAddress("MMpi", &P_kin_MMpi_pions_cut_rn);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// For Kaons
  Double_t H_gtr_beta_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("H_gtr_beta", &H_gtr_beta_kaons_uncut);
  Double_t H_gtr_xp_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("H_gtr_xp", &H_gtr_xp_kaons_uncut);
  Double_t H_gtr_yp_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("H_gtr_yp", &H_gtr_yp_kaons_uncut);
  Double_t H_gtr_dp_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("H_gtr_dp", &H_gtr_dp_kaons_uncut);
  Double_t H_hod_goodscinhit_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("H_hod_goodscinhit", &H_hod_goodscinhit_kaons_uncut);
  Double_t H_hod_goodstarttime_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("H_hod_goodstarttime", &H_hod_goodstarttime_kaons_uncut);
  Double_t H_cal_etotnorm_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("H_cal_etotnorm", &H_cal_etotnorm_kaons_uncut);
  Double_t H_cal_etottracknorm_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("H_cal_etottracknorm", &H_cal_etottracknorm_kaons_uncut);
  Double_t H_cer_npeSum_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("H_cer_npeSum", &H_cer_npeSum_kaons_uncut);
  Double_t H_RFTime_Dist_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("H_RFTime_Dist", &H_RFTime_Dist_kaons_uncut);
  Double_t P_gtr_beta_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("P_gtr_beta", &P_gtr_beta_kaons_uncut);
  Double_t P_gtr_xp_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("P_gtr_xp", &P_gtr_xp_kaons_uncut);
  Double_t P_gtr_yp_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("P_gtr_yp", &P_gtr_yp_kaons_uncut);
  Double_t P_gtr_dp_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("P_gtr_dp", &P_gtr_dp_kaons_uncut);
  Double_t P_gtr_p_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("P_gtr_p", &P_gtr_p_kaons_uncut);
  Double_t P_hod_goodscinhit_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("P_hod_goodscinhit", &P_hod_goodscinhit_kaons_uncut);
  Double_t P_hod_goodstarttime_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("P_hod_goodstarttime", &P_hod_goodstarttime_kaons_uncut);
  Double_t P_cal_etotnorm_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("P_cal_etotnorm", &P_cal_etotnorm_kaons_uncut);
  Double_t P_cal_etottracknorm_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("P_cal_etottracknorm", &P_cal_etottracknorm_kaons_uncut);
  Double_t P_hgcer_npeSum_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("P_hgcer_npeSum", &P_hgcer_npeSum_kaons_uncut);
  Double_t P_hgcer_xAtCer_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("P_hgcer_xAtCer", &P_hgcer_xAtCer_kaons_uncut);
  Double_t P_hgcer_yAtCer_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("P_hgcer_yAtCer", &P_hgcer_yAtCer_kaons_uncut);
  Double_t P_aero_npeSum_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("P_aero_npeSum", &P_aero_npeSum_kaons_uncut);
  Double_t P_aero_xAtAero_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("P_aero_xAtAero", &P_aero_xAtAero_kaons_uncut);
  Double_t P_aero_yAtAero_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("P_aero_yAtAero", &P_aero_yAtAero_kaons_uncut);
  Double_t P_kin_MMK_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("MMK", &P_kin_MMK_kaons_uncut);
  Double_t P_RFTime_Dist_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("P_RF_Dist", &P_RFTime_Dist_kaons_uncut);
  Double_t CTime_eKCoinTime_ROC1_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("CTime_eKCoinTime_ROC1", &CTime_eKCoinTime_ROC1_kaons_uncut);
  Double_t Q2_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("Q2", &Q2_kaons_uncut);
  Double_t W_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("W", &W_kaons_uncut);
  Double_t epsilon_kaons_uncut; Uncut_Kaon_Events->SetBranchAddress("Q2", &epsilon_kaons_uncut);

  Double_t H_gtr_beta_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("H_gtr_beta", &H_gtr_beta_kaons_cut_all);
  Double_t H_gtr_xp_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("H_gtr_xp", &H_gtr_xp_kaons_cut_all);
  Double_t H_gtr_yp_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("H_gtr_yp", &H_gtr_yp_kaons_cut_all);
  Double_t H_gtr_dp_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("H_gtr_dp", &H_gtr_dp_kaons_cut_all);
  Double_t H_hod_goodscinhit_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("H_hod_goodscinhit", &H_hod_goodscinhit_kaons_cut_all);
  Double_t H_hod_goodstarttime_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("H_hod_goodstarttime", &H_hod_goodstarttime_kaons_cut_all);
  Double_t H_cal_etotnorm_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("H_cal_etotnorm", &H_cal_etotnorm_kaons_cut_all);
  Double_t H_cal_etottracknorm_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("H_cal_etottracknorm", &H_cal_etottracknorm_kaons_cut_all);
  Double_t H_cer_npeSum_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("H_cer_npeSum", &H_cer_npeSum_kaons_cut_all);
  Double_t H_RFTime_Dist_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("H_RFTime_Dist", &H_RFTime_Dist_kaons_cut_all);
  Double_t P_gtr_beta_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("P_gtr_beta", &P_gtr_beta_kaons_cut_all);
  Double_t P_gtr_xp_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("P_gtr_xp", &P_gtr_xp_kaons_cut_all);
  Double_t P_gtr_yp_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("P_gtr_yp", &P_gtr_yp_kaons_cut_all);
  Double_t P_gtr_dp_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("P_gtr_dp", &P_gtr_dp_kaons_cut_all);
  Double_t P_gtr_p_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("P_gtr_p", &P_gtr_p_kaons_cut_all);
  Double_t P_hod_goodscinhit_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("P_hod_goodscinhit", &P_hod_goodscinhit_kaons_cut_all);
  Double_t P_hod_goodstarttime_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("P_hod_goodstarttime", &P_hod_goodstarttime_kaons_cut_all);
  Double_t P_cal_etotnorm_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("P_cal_etotnorm", &P_cal_etotnorm_kaons_cut_all);
  Double_t P_cal_etottracknorm_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("P_cal_etottracknorm", &P_cal_etottracknorm_kaons_cut_all);
  Double_t P_hgcer_npeSum_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("P_hgcer_npeSum", &P_hgcer_npeSum_kaons_cut_all);
  Double_t P_hgcer_xAtCer_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("P_hgcer_xAtCer", &P_hgcer_xAtCer_kaons_cut_all);
  Double_t P_hgcer_yAtCer_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("P_hgcer_yAtCer", &P_hgcer_yAtCer_kaons_cut_all);
  Double_t P_aero_npeSum_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("P_aero_npeSum", &P_aero_npeSum_kaons_cut_all);
  Double_t P_aero_xAtAero_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("P_aero_xAtAero", &P_aero_xAtAero_kaons_cut_all);
  Double_t P_aero_yAtAero_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("P_aero_yAtAero", &P_aero_yAtAero_kaons_cut_all);
  Double_t P_kin_MMK_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("MMK", &P_kin_MMK_kaons_cut_all);
  Double_t P_RFTime_Dist_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("P_RF_Dist", &P_RFTime_Dist_kaons_cut_all);
  Double_t CTime_eKCoinTime_ROC1_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("CTime_eKCoinTime_ROC1", &CTime_eKCoinTime_ROC1_kaons_cut_all);
  Double_t Q2_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("Q2", &Q2_kaons_cut_all);
  Double_t W_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("W", &W_kaons_cut_all);
  Double_t epsilon_kaons_cut_all; Cut_Kaon_Events_All->SetBranchAddress("Q2", &epsilon_kaons_cut_all);

  Double_t P_gtr_beta_kaons_cut_pr; Cut_Kaon_Events_Prompt->SetBranchAddress("P_gtr_beta", &P_gtr_beta_kaons_cut_pr);
  Double_t P_RFTime_Dist_kaons_cut_pr; Cut_Kaon_Events_Prompt->SetBranchAddress("P_RF_Dist", &P_RFTime_Dist_kaons_cut_pr);
  Double_t CTime_eKCoinTime_ROC1_kaons_cut_pr; Cut_Kaon_Events_Prompt->SetBranchAddress("CTime_eKCoinTime_ROC1", &CTime_eKCoinTime_ROC1_kaons_cut_pr);
  Double_t P_kin_MMK_kaons_cut_pr; Cut_Kaon_Events_Prompt->SetBranchAddress("MMK", &P_kin_MMK_kaons_cut_pr);

  Double_t P_gtr_beta_kaons_cut_rn; Cut_Kaon_Events_Random->SetBranchAddress("P_gtr_beta", &P_gtr_beta_kaons_cut_rn);
  Double_t P_RFTime_Dist_kaons_cut_rn; Cut_Kaon_Events_Random->SetBranchAddress("P_RF_Dist", &P_RFTime_Dist_kaons_cut_rn);
  Double_t CTime_eKCoinTime_ROC1_kaons_cut_rn; Cut_Kaon_Events_Random->SetBranchAddress("CTime_eKCoinTime_ROC1", &CTime_eKCoinTime_ROC1_kaons_cut_rn);
  Double_t P_kin_MMK_kaons_cut_rn; Cut_Kaon_Events_Random->SetBranchAddress("MMK", &P_kin_MMK_kaons_cut_rn);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// For Protons
  Double_t H_gtr_beta_protons_uncut; Uncut_Proton_Events->SetBranchAddress("H_gtr_beta", &H_gtr_beta_protons_uncut);
  Double_t H_gtr_xp_protons_uncut; Uncut_Proton_Events->SetBranchAddress("H_gtr_xp", &H_gtr_xp_protons_uncut);
  Double_t H_gtr_yp_protons_uncut; Uncut_Proton_Events->SetBranchAddress("H_gtr_yp", &H_gtr_yp_protons_uncut);
  Double_t H_gtr_dp_protons_uncut; Uncut_Proton_Events->SetBranchAddress("H_gtr_dp", &H_gtr_dp_protons_uncut);
  Double_t H_hod_goodscinhit_protons_uncut; Uncut_Proton_Events->SetBranchAddress("H_hod_goodscinhit", &H_hod_goodscinhit_protons_uncut);
  Double_t H_hod_goodstarttime_protons_uncut; Uncut_Proton_Events->SetBranchAddress("H_hod_goodstarttime", &H_hod_goodstarttime_protons_uncut);
  Double_t H_cal_etotnorm_protons_uncut; Uncut_Proton_Events->SetBranchAddress("H_cal_etotnorm", &H_cal_etotnorm_protons_uncut);
  Double_t H_cal_etottracknorm_protons_uncut; Uncut_Proton_Events->SetBranchAddress("H_cal_etottracknorm", &H_cal_etottracknorm_protons_uncut);
  Double_t H_cer_npeSum_protons_uncut; Uncut_Proton_Events->SetBranchAddress("H_cer_npeSum", &H_cer_npeSum_protons_uncut);
  Double_t H_RFTime_Dist_protons_uncut; Uncut_Proton_Events->SetBranchAddress("H_RFTime_Dist", &H_RFTime_Dist_protons_uncut);
  Double_t P_gtr_beta_protons_uncut; Uncut_Proton_Events->SetBranchAddress("P_gtr_beta", &P_gtr_beta_protons_uncut);
  Double_t P_gtr_xp_protons_uncut; Uncut_Proton_Events->SetBranchAddress("P_gtr_xp", &P_gtr_xp_protons_uncut);
  Double_t P_gtr_yp_protons_uncut; Uncut_Proton_Events->SetBranchAddress("P_gtr_yp", &P_gtr_yp_protons_uncut);
  Double_t P_gtr_dp_protons_uncut; Uncut_Proton_Events->SetBranchAddress("P_gtr_dp", &P_gtr_dp_protons_uncut);
  Double_t P_gtr_p_protons_uncut; Uncut_Proton_Events->SetBranchAddress("P_gtr_p", &P_gtr_p_protons_uncut);
  Double_t P_hod_goodscinhit_protons_uncut; Uncut_Proton_Events->SetBranchAddress("P_hod_goodscinhit", &P_hod_goodscinhit_protons_uncut);
  Double_t P_hod_goodstarttime_protons_uncut; Uncut_Proton_Events->SetBranchAddress("P_hod_goodstarttime", &P_hod_goodstarttime_protons_uncut);
  Double_t P_cal_etotnorm_protons_uncut; Uncut_Proton_Events->SetBranchAddress("P_cal_etotnorm", &P_cal_etotnorm_protons_uncut);
  Double_t P_cal_etottracknorm_protons_uncut; Uncut_Proton_Events->SetBranchAddress("P_cal_etottracknorm", &P_cal_etottracknorm_protons_uncut);
  Double_t P_hgcer_npeSum_protons_uncut; Uncut_Proton_Events->SetBranchAddress("P_hgcer_npeSum", &P_hgcer_npeSum_protons_uncut);
  Double_t P_hgcer_xAtCer_protons_uncut; Uncut_Proton_Events->SetBranchAddress("P_hgcer_xAtCer", &P_hgcer_xAtCer_protons_uncut);
  Double_t P_hgcer_yAtCer_protons_uncut; Uncut_Proton_Events->SetBranchAddress("P_hgcer_yAtCer", &P_hgcer_yAtCer_protons_uncut);
  Double_t P_aero_npeSum_protons_uncut; Uncut_Proton_Events->SetBranchAddress("P_aero_npeSum", &P_aero_npeSum_protons_uncut);
  Double_t P_aero_xAtAero_protons_uncut; Uncut_Proton_Events->SetBranchAddress("P_aero_xAtAero", &P_aero_xAtAero_protons_uncut);
  Double_t P_aero_yAtAero_protons_uncut; Uncut_Proton_Events->SetBranchAddress("P_aero_yAtAero", &P_aero_yAtAero_protons_uncut);
  Double_t P_kin_MMp_protons_uncut; Uncut_Proton_Events->SetBranchAddress("MMp", &P_kin_MMp_protons_uncut);
  Double_t P_RFTime_Dist_protons_uncut; Uncut_Proton_Events->SetBranchAddress("P_RF_Dist", &P_RFTime_Dist_protons_uncut);
  Double_t CTime_epCoinTime_ROC1_protons_uncut; Uncut_Proton_Events->SetBranchAddress("CTime_epCoinTime_ROC1", &CTime_epCoinTime_ROC1_protons_uncut);
  Double_t Q2_protons_uncut; Uncut_Proton_Events->SetBranchAddress("Q2", &Q2_protons_uncut);
  Double_t W_protons_uncut; Uncut_Proton_Events->SetBranchAddress("W", &W_protons_uncut);
  Double_t epsilon_protons_uncut; Uncut_Proton_Events->SetBranchAddress("Q2", &epsilon_protons_uncut);

  Double_t H_gtr_beta_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("H_gtr_beta", &H_gtr_beta_protons_cut_all);
  Double_t H_gtr_xp_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("H_gtr_xp", &H_gtr_xp_protons_cut_all);
  Double_t H_gtr_yp_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("H_gtr_yp", &H_gtr_yp_protons_cut_all);
  Double_t H_gtr_dp_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("H_gtr_dp", &H_gtr_dp_protons_cut_all);
  Double_t H_hod_goodscinhit_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("H_hod_goodscinhit", &H_hod_goodscinhit_protons_cut_all);
  Double_t H_hod_goodstarttime_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("H_hod_goodstarttime", &H_hod_goodstarttime_protons_cut_all);
  Double_t H_cal_etotnorm_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("H_cal_etotnorm", &H_cal_etotnorm_protons_cut_all);
  Double_t H_cal_etottracknorm_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("H_cal_etottracknorm", &H_cal_etottracknorm_protons_cut_all);
  Double_t H_cer_npeSum_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("H_cer_npeSum", &H_cer_npeSum_protons_cut_all);
  Double_t H_RFTime_Dist_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("H_RFTime_Dist", &H_RFTime_Dist_protons_cut_all);
  Double_t P_gtr_beta_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("P_gtr_beta", &P_gtr_beta_protons_cut_all);
  Double_t P_gtr_xp_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("P_gtr_xp", &P_gtr_xp_protons_cut_all);
  Double_t P_gtr_yp_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("P_gtr_yp", &P_gtr_yp_protons_cut_all);
  Double_t P_gtr_dp_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("P_gtr_dp", &P_gtr_dp_protons_cut_all);
  Double_t P_gtr_p_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("P_gtr_p", &P_gtr_p_protons_cut_all);
  Double_t P_hod_goodscinhit_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("P_hod_goodscinhit", &P_hod_goodscinhit_protons_cut_all);
  Double_t P_hod_goodstarttime_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("P_hod_goodstarttime", &P_hod_goodstarttime_protons_cut_all);
  Double_t P_cal_etotnorm_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("P_cal_etotnorm", &P_cal_etotnorm_protons_cut_all);
  Double_t P_cal_etottracknorm_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("P_cal_etottracknorm", &P_cal_etottracknorm_protons_cut_all);
  Double_t P_hgcer_npeSum_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("P_hgcer_npeSum", &P_hgcer_npeSum_protons_cut_all);
  Double_t P_hgcer_xAtCer_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("P_hgcer_xAtCer", &P_hgcer_xAtCer_protons_cut_all);
  Double_t P_hgcer_yAtCer_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("P_hgcer_yAtCer", &P_hgcer_yAtCer_protons_cut_all);
  Double_t P_aero_npeSum_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("P_aero_npeSum", &P_aero_npeSum_protons_cut_all);
  Double_t P_aero_xAtAero_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("P_aero_xAtAero", &P_aero_xAtAero_protons_cut_all);
  Double_t P_aero_yAtAero_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("P_aero_yAtAero", &P_aero_yAtAero_protons_cut_all);
  Double_t P_kin_MMp_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("MMp", &P_kin_MMp_protons_cut_all);
  Double_t P_RFTime_Dist_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("P_RF_Dist", &P_RFTime_Dist_protons_cut_all);
  Double_t CTime_epCoinTime_ROC1_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("CTime_epCoinTime_ROC1", &CTime_epCoinTime_ROC1_protons_cut_all);
  Double_t Q2_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("Q2", &Q2_protons_cut_all);
  Double_t W_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("W", &W_protons_cut_all);
  Double_t epsilon_protons_cut_all; Cut_Proton_Events_All->SetBranchAddress("Q2", &epsilon_protons_cut_all);

  Double_t P_gtr_beta_protons_cut_pr; Cut_Proton_Events_Prompt->SetBranchAddress("P_gtr_beta", &P_gtr_beta_protons_cut_pr);
  Double_t P_RFTime_Dist_protons_cut_pr; Cut_Proton_Events_Prompt->SetBranchAddress("P_RF_Dist", &P_RFTime_Dist_protons_cut_pr);
  Double_t CTime_epCoinTime_ROC1_protons_cut_pr; Cut_Proton_Events_Prompt->SetBranchAddress("CTime_epCoinTime_ROC1", &CTime_epCoinTime_ROC1_protons_cut_pr);
  Double_t P_kin_MMp_protons_cut_pr; Cut_Proton_Events_Prompt->SetBranchAddress("MMp", &P_kin_MMp_protons_cut_pr);

  Double_t P_gtr_beta_protons_cut_rn; Cut_Proton_Events_Random->SetBranchAddress("P_gtr_beta", &P_gtr_beta_protons_cut_rn);
  Double_t P_RFTime_Dist_protons_cut_rn; Cut_Proton_Events_Random->SetBranchAddress("P_RF_Dist", &P_RFTime_Dist_protons_cut_rn);
  Double_t CTime_epCoinTime_ROC1_protons_cut_rn; Cut_Proton_Events_Random->SetBranchAddress("CTime_epCoinTime_ROC1", &CTime_epCoinTime_ROC1_protons_cut_rn);
  Double_t P_kin_MMp_protons_cut_rn; Cut_Proton_Events_Random->SetBranchAddress("MMp", &P_kin_MMp_protons_cut_rn);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Define Histograms for Pions
  TH1D *h1_H_gtr_beta_pions_Uncut = new TH1D("h1_H_gtr_beta_pions_Uncut", "H_gtr_beta", 220, 0.6, 1.4);
  TH1D *h1_H_gtr_xp_pions_Uncut = new TH1D("h1_H_gtr_xp_pions_Uncut", "H_gtr_xp", 220, -0.6, 0.6);
  TH1D *h1_H_gtr_yp_pions_Uncut = new TH1D("h1_H_gtr_yp_pions_Uncut", "H_gtr_yp", 220, -0.6, 0.6);
  TH1D *h1_H_gtr_dp_pions_Uncut = new TH1D("h1_H_gtr_dp_pions_Uncut", "H_gtr_dp", 220, -40, 40);
  TH1D *h1_H_hod_goodscinhit_pions_Uncut = new TH1D("h1_H_hod_goodscinhit_pions_Uncut", "H_hod_goodscinhit", 220, 0.6, 1.4);
  TH1D *h1_H_hod_goodstarttime_pions_Uncut = new TH1D("h1_H_hod_goodstarttime_pions_Uncut", "H_hod_goodstarttime", 220, 0.6, 1.4);
  TH1D *h1_H_cal_etotnorm_pions_Uncut = new TH1D("h1_H_cal_etotnorm_pions_Uncut", "H_cal_etotnorm", 220, 0, 1.5);
  TH1D *h1_H_cal_etottracknorm_pions_Uncut = new TH1D("h1_H_cal_etottracknorm_pions_Uncut", "H_cal_etottracknorm", 220, 0, 1.5); 
  TH1D *h1_H_cer_npeSum_pions_Uncut = new TH1D("h1_H_cer_npeSum_pions_Uncut", "H_cer_npeSum", 220, 0, 2);
  TH1D *h1_H_RFTime_Dist_pions_Uncut = new TH1D("h1_H_RFTime_Dist_pions_Uncut", "H_RFTime_Dist", 220, -1, 5);
  TH1D *h1_P_gtr_beta_pions_Uncut = new TH1D("h1_P_gtr_beta_pions_Uncut", "P_gtr_beta", 220, 0.6, 1.5);
  TH1D *h1_P_gtr_xp_pions_Uncut = new TH1D("h1_P_gtr_xp_pions_Uncut", "P_gtr_xp", 220, -0.6, 0.6);
  TH1D *h1_P_gtr_yp_pions_Uncut = new TH1D("h1_P_gtr_yp_pions_Uncut", "P_gtr_yp", 220, -0.8, 0.8);
  TH1D *h1_P_gtr_dp_pions_Uncut = new TH1D("h1_P_gtr_dp_pions_Uncut", "P_gtr_dp", 220, -40, 40);
  TH1D *h1_P_gtr_p_pions_Uncut = new TH1D("h1_P_gtr_p_pions_Uncut", "P_gtr_p", 220, -5, 15);
  TH1D *h1_P_hod_goodscinhit_pions_Uncut = new TH1D("h1_P_hod_goodscinhit_pions_Uncut", "P_hod_goodscinhit", 220, 0.6, 1.5);
  TH1D *h1_P_hod_goodstarttime_pions_Uncut = new TH1D("h1_P_hod_goodstarttime_pions_Uncut", "P_hod_goodstarttime", 220, 0.6, 1.5);
  TH1D *h1_P_cal_etotnorm_pions_Uncut = new TH1D("h1_P_cal_etotnorm_pions_Uncut", "P_cal_etotnorm", 220, 0, 1);
  TH1D *h1_P_cal_etottracknorm_pions_Uncut = new TH1D("h1_P_cal_etottracknorm_pions_Uncut", "P_cal_etottracknorm", 220, 0, 1);
  TH1D *h1_P_hgcer_npeSum_pions_Uncut = new TH1D("h1_P_hgcer_npeSum_pions_Uncut", "P_hgcer_npeSum", 220, 0, 0.5);
  TH1D *h1_P_hgcer_xAtCer_pions_Uncut = new TH1D("h1_P_hgcer_xAtCer_pions_Uncut", "P_hgcer_xAtCer", 220, -10, 10);
  TH1D *h1_P_hgcer_yAtCer_pions_Uncut = new TH1D("h1_P_hgcer_yAtCer_pions_Uncut", "P_hgcer_yAtCer", 220, -10, 10);
  TH1D *h1_P_aero_npeSum_pions_Uncut = new TH1D("h1_P_aero_npeSum_pions_Uncut", "P_aero_npeSum", 220, 0, 10);
  TH1D *h1_P_aero_xAtAero_pions_Uncut = new TH1D("h1_P_aero_xAtAero_pions_Uncut", "P_aero_xAtAero", 220, -10, 10);
  TH1D *h1_P_aero_yAtAero_pions_Uncut = new TH1D("h1_P_aero_yAtAero_pions_Uncut", "P_aero_yAtAero", 220, -10, 10);
  TH1D *h1_P_kin_MMpi_pions_Uncut = new TH1D("h1_P_kin_MMpi_pions_Uncut", "P_kin_MMpi", 220, 0, 2);
  TH1D *h1_P_RFTime_Dist_pions_Uncut = new TH1D("h1_P_RFTime_Dist_pions_Uncut", "P_RFTime_Dist", 220, -1, 5);
  TH1D *h1_CTime_ePiCoinTime_ROC1_pions_Uncut = new TH1D("h1_CTime_ePiCoinTime_ROC1_pions_Uncut", "CTime_ePiCoinTime_ROC1", 220, -100, 100);
  TH1D *h1_Q2_pions_Uncut = new TH1D("h1_Q2_pions_Uncut", "Q2", 220, -10, 10);
  TH1D *h1_W_pions_Uncut = new TH1D("h1_W_pions_Uncut", "W", 220, -10, 10);
  TH1D *h1_epsilon_pions_Uncut = new TH1D("h1_epsilon_pions_Uncut", "epsilon", 220, -10, 10);

  TH1D *h1_H_gtr_beta_pions_Cut_All = new TH1D("h1_H_gtr_beta_pions_Cut_All", "H_gtr_beta", 220, 0.6, 1.5);
  TH1D *h1_H_gtr_xp_pions_Cut_All = new TH1D("h1_H_gtr_xp_pions_Cut_All", "H_gtr_xp", 220, -0.6, 0.6);
  TH1D *h1_H_gtr_yp_pions_Cut_All = new TH1D("h1_H_gtr_yp_pions_Cut_All", "H_gtr_yp", 220, -0.6, 0.6);
  TH1D *h1_H_gtr_dp_pions_Cut_All = new TH1D("h1_H_gtr_dp_pions_Cut_All", "H_gtr_dp", 220, -15, 15);
  TH1D *h1_H_hod_goodscinhit_pions_Cut_All = new TH1D("h1_H_hod_goodscinhit_pions_Cut_All", "H_hod_goodscinhit", 220, 0.6, 1.5);
  TH1D *h1_H_hod_goodstarttime_pions_Cut_All = new TH1D("h1_H_hod_goodstarttime_pions_Cut_All", "H_hod_goodstarttime", 220, 0.6, 1.5);
  TH1D *h1_H_cal_etotnorm_pions_Cut_All = new TH1D("h1_H_cal_etotnorm_pions_Cut_All", "H_cal_etotnorm", 220, 0, 2);
  TH1D *h1_H_cal_etottracknorm_pions_Cut_All = new TH1D("h1_H_cal_etottracknorm_pions_Cut_All", "H_cal_etottracknorm", 220, 0, 2);
  TH1D *h1_H_cer_npeSum_pions_Cut_All = new TH1D("h1_H_cer_npeSum_pions_Cut_All", "H_cer_npeSum", 220, 0, 20);
  TH1D *h1_H_RFTime_Dist_pions_Cut_All = new TH1D("h1_H_RFTime_Dist_pions_Cut_All", "H_RFTime_Dist", 220, -1, 5);
  TH1D *h1_P_gtr_beta_pions_Cut_All = new TH1D("h1_P_gtr_beta_pions_Cut_All", "P_gtr_beta", 220, 0, 2);
  TH1D *h1_P_gtr_xp_pions_Cut_All = new TH1D("h1_P_gtr_xp_pions_Cut_All", "P_gtr_xp", 220, -0.6, 0.6);
  TH1D *h1_P_gtr_yp_pions_Cut_All = new TH1D("h1_P_gtr_yp_pions_Cut_All", "P_gtr_yp", 220, -0.6, 0.6);
  TH1D *h1_P_gtr_dp_pions_Cut_All = new TH1D("h1_P_gtr_dp_pions_Cut_All", "P_gtr_dp", 220, -20, 20);
  TH1D *h1_P_gtr_p_pions_Cut_All = new TH1D("h1_P_gtr_p_pions_Cut_All", "P_gtr_p", 220, 0, 12);
  TH1D *h1_P_hod_goodscinhit_pions_Cut_All = new TH1D("h1_P_hod_goodscinhit_pions_Cut_All", "P_hod_goodscinhit", 220, 0.6, 1.4);
  TH1D *h1_P_hod_goodstarttime_pions_Cut_All = new TH1D("h1_P_hod_goodstarttime_pions_Cut_All", "P_hod_goodstarttime", 220, 0.6, 1.4);
  TH1D *h1_P_cal_etotnorm_pions_Cut_All = new TH1D("h1_P_cal_etotnorm_pions_Cut_All", "P_cal_etotnorm", 220, 0, 1.5);
  TH1D *h1_P_cal_etottracknorm_pions_Cut_All = new TH1D("h1_P_cal_etottracknorm_pions_Cut_All", "P_cal_etottracknorm", 220, 0, 1.5);
  TH1D *h1_P_hgcer_npeSum_pions_Cut_All = new TH1D("h1_P_hgcer_npeSum_pions_Cut_All", "P_hgcer_npeSum", 220, -10, 50);
  TH1D *h1_P_hgcer_xAtCer_pions_Cut_All = new TH1D("h1_P_hgcer_xAtCer_pions_Cut_All", "P_hgcer_xAtCer", 220, -40, 40);
  TH1D *h1_P_hgcer_yAtCer_pions_Cut_All = new TH1D("h1_P_hgcer_yAtCer_pions_Cut_All", "P_hgcer_yAtCer", 220, -40, 40);
  TH1D *h1_P_aero_npeSum_pions_Cut_All = new TH1D("h1_P_aero_npeSum_pions_Cut_All", "P_aero_npeSum", 220, 1, 50);
  TH1D *h1_P_aero_xAtAero_pions_Cut_All = new TH1D("h1_P_aero_xAtAero_pions_Cut_All", "P_aero_xAtAero", 220, -40, 40);
  TH1D *h1_P_aero_yAtAero_pions_Cut_All = new TH1D("h1_P_aero_yAtAero_pions_Cut_All", "P_aero_yAtAero", 220, -40, 40);
  TH1D *h1_P_kin_MMpi_pions_Cut_All = new TH1D("h1_P_kin_MMpi_pions_Cut_All", "P_kin_MMpi", 220, 0.5, 1.8);
  TH1D *h1_P_RFTime_Dist_pions_Cut_All = new TH1D("h1_P_RFTime_Dist_pions_Cut_All", "P_RFTime_Dist", 220, -1, 5);
  TH1D *h1_CTime_ePiCoinTime_ROC1_pions_Cut_All = new TH1D("h1_CTime_ePiCoinTime_ROC1_pions_Cut_All", "CTime_ePiCoinTime_ROC1", 220, -50, 50);
  TH1D *h1_Q2_pions_Cut_All = new TH1D("h1_Q2_pions_Cut_All", "Q2", 220, -10, 10);
  TH1D *h1_W_pions_Cut_All = new TH1D("h1_W_pions_Cut_All", "W", 220, -2, 8);
  TH1D *h1_epsilon_pions_Cut_All = new TH1D("h1_epsilon_pions_Cut_All", "epsilon", 220, -2, 8);

  TH1D *h1_P_gtr_beta_pions_Cut_Prompt = new TH1D("h1_P_gtr_beta_pions_Cut_Prompt", "P_gtr_beta", 220, 0.6, 1.4);
  TH1D *h1_P_RFTime_Dist_pions_Cut_Prompt = new TH1D("h1_P_RFTime_Dist_pions_Cut_Prompt", "P_RFTime_Dist", 220, 0, 4);
  TH1D *h1_CTime_ePiCoinTime_ROC1_pions_Cut_Prompt = new TH1D("h1_CTime_ePiCoinTime_ROC1_pions_Cut_Prompt", "CTime_ePiCoinTime_ROC1", 220, -8, 8);
  TH1D *h1_P_kin_MMpi_pions_Cut_Prompt = new TH1D("h1_P_kin_MMpi_pions_Cut_Prompt", "P_kin_MMpi", 220, 0.5, 1.8);

  TH1D *h1_P_gtr_beta_pions_Cut_Random = new TH1D("h1_P_gtr_beta_pions_Cut_Random", "P_gtr_beta", 220, 0.6, 1.4);
  TH1D *h1_P_RFTime_Dist_pions_Cut_Random = new TH1D("h1_P_RFTime_Dist_pions_Cut_Random", "P_RFTime_Dist", 220, -1, 5);
  TH1D *h1_CTime_ePiCoinTime_ROC1_pions_Cut_Random = new TH1D("h1_CTime_ePiCoinTime_ROC1_pions_Cut_Random", "CTime_ePiCoinTime_ROC1", 220, -40, 40);
  TH1D *h1_P_kin_MMpi_pions_Cut_Random = new TH1D("h1_P_kin_MMpi_pions_Cut_Random", "P_kin_MMpi", 220, 0.5, 1.8);

  TH1D *h1_CTime_ePiCoinTime_ROC1_pions_Cut_Random_Scaled = new TH1D("h1_CTime_ePiCoinTime_ROC1_pions_Cut_Random_Scaled", "CTime_ePiCoinTime_ROC1", 220, -40, 40);
  TH1D *h1_CTime_ePiCoinTime_ROC1_pions_Rndm_Sub = new TH1D("h1_CTime_ePiCoinTime_ROC1_pions_Rndm_Sub", "CTime_ePiCoinTime_ROC1", 220, -8, 8);
  TH1D *h1_P_kin_MMpi_pions_Cut_Random_Scaled = new TH1D("h1_P_kin_MMpi_pions_Cut_Random_Scaled", "P_kin_MMpi", 220, 0.5, 1.8);
  TH1D *h1_P_kin_MMpi_pions_Rndm_Sub = new TH1D("h1_P_kin_MMpi_pions_Rndm_Sub", "P_kin_MMpi", 220, 0.5, 1.8);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Define Histograms for Kaons
  TH1D *h1_H_gtr_beta_kaons_Uncut = new TH1D("h1_H_gtr_beta_kaons_Uncut", "H_gtr_beta", 220, 0.6, 1.5);
  TH1D *h1_H_gtr_xp_kaons_Uncut = new TH1D("h1_H_gtr_xp_kaons_Uncut", "H_gtr_xp", 220, -0.6, 0.6);
  TH1D *h1_H_gtr_yp_kaons_Uncut = new TH1D("h1_H_gtr_yp_kaons_Uncut", "H_gtr_yp", 220, -0.6, 0.6);
  TH1D *h1_H_gtr_dp_kaons_Uncut = new TH1D("h1_H_gtr_dp_kaons_Uncut", "H_gtr_dp", 220, -40, 40);
  TH1D *h1_H_hod_goodscinhit_kaons_Uncut = new TH1D("h1_H_hod_goodscinhit_kaons_Uncut", "H_hod_goodscinhit", 220, 0.6, 1.4);
  TH1D *h1_H_hod_goodstarttime_kaons_Uncut = new TH1D("h1_H_hod_goodstarttime_kaons_Uncut", "H_hod_goodstarttime", 220, 0.6, 1.4);
  TH1D *h1_H_cal_etotnorm_kaons_Uncut = new TH1D("h1_H_cal_etotnorm_kaons_Uncut", "H_cal_etotnorm", 220, 0, 2);
  TH1D *h1_H_cal_etottracknorm_kaons_Uncut = new TH1D("h1_H_cal_etottracknorm_kaons_Uncut", "H_cal_etottracknorm", 220, 0, 2);
  TH1D *h1_H_cer_npeSum_kaons_Uncut = new TH1D("h1_H_cer_npeSum_kaons_Uncut", "H_cer_npeSum", 220, 0, 2);
  TH1D *h1_H_RFTime_Dist_kaons_Uncut = new TH1D("h1_H_RFTime_Dist_kaons_Uncut", "H_RFTime_Dist", 220, -1, 5);
  TH1D *h1_P_gtr_beta_kaons_Uncut = new TH1D("h1_P_gtr_beta_kaons_Uncut", "P_gtr_beta", 220, 0, 1.6);
  TH1D *h1_P_gtr_xp_kaons_Uncut = new TH1D("h1_P_gtr_xp_kaons_Uncut", "P_gtr_xp", 220, -0.6, 0.6);
  TH1D *h1_P_gtr_yp_kaons_Uncut = new TH1D("h1_P_gtr_yp_kaons_Uncut", "P_gtr_yp", 220, -0.8, 0.8);
  TH1D *h1_P_gtr_dp_kaons_Uncut = new TH1D("h1_P_gtr_dp_kaons_Uncut", "P_gtr_dp", 220, -40, 40);
  TH1D *h1_P_gtr_p_kaons_Uncut = new TH1D("h1_P_gtr_p_kaons_Uncut", "P_gtr_p", 220, -5, 15);
  TH1D *h1_P_hod_goodscinhit_kaons_Uncut = new TH1D("h1_P_hod_goodscinhit_kaons_Uncut", "P_hod_goodscinhit", 220, 0.6, 1.4);
  TH1D *h1_P_hod_goodstarttime_kaons_Uncut = new TH1D("h1_P_hod_goodstarttime_kaons_Uncut", "P_hod_goodstarttime", 220, 0.6, 1.4);
  TH1D *h1_P_cal_etotnorm_kaons_Uncut = new TH1D("h1_P_cal_etotnorm_kaons_Uncut", "P_cal_etotnorm", 220, 0, 1.6);
  TH1D *h1_P_cal_etottracknorm_kaons_Uncut = new TH1D("h1_P_cal_etottracknorm_kaons_Uncut", "P_cal_etottracknorm", 220, 0, 1.6);
  TH1D *h1_P_hgcer_npeSum_kaons_Uncut = new TH1D("h1_P_hgcer_npeSum_kaons_Uncut", "P_hgcer_npeSum", 220, 0, 0.5);
  TH1D *h1_P_hgcer_xAtCer_kaons_Uncut = new TH1D("h1_P_hgcer_xAtCer_kaons_Uncut", "P_hgcer_xAtCer", 220, -10, 10);
  TH1D *h1_P_hgcer_yAtCer_kaons_Uncut = new TH1D("h1_P_hgcer_yAtCer_kaons_Uncut", "P_hgcer_yAtCer", 220, -10, 10);
  TH1D *h1_P_aero_npeSum_kaons_Uncut = new TH1D("h1_P_aero_npeSum_kaons_Uncut", "P_aero_npeSum", 220, 0, 10);
  TH1D *h1_P_aero_xAtAero_kaons_Uncut = new TH1D("h1_P_aero_xAtAero_kaons_Uncut", "P_aero_xAtAero", 220, -10, 10);
  TH1D *h1_P_aero_yAtAero_kaons_Uncut = new TH1D("h1_P_aero_yAtAero_kaons_Uncut", "P_aero_yAtAero", 220, -10, 10);
  TH1D *h1_P_kin_MMK_kaons_Uncut = new TH1D("h1_P_kin_MMK_kaons_Uncut", "P_kin_MMK", 220, 0, 2);
  TH1D *h1_P_RFTime_Dist_kaons_Uncut = new TH1D("h1_P_RFTime_Dist_kaons_Uncut", "P_RFTime_Dist", 220, -1, 5);
  TH1D *h1_CTime_eKCoinTime_ROC1_kaons_Uncut = new TH1D("h1_CTime_eKCoinTime_ROC1_kaons_Uncut", "CTime_eKCoinTime_ROC1", 220, -100, 100);
  TH1D *h1_Q2_kaons_Uncut = new TH1D("h1_Q2_kaons_Uncut", "Q2", 220, -10, 10);
  TH1D *h1_W_kaons_Uncut = new TH1D("h1_W_kaons_Uncut", "W", 220, -10, 10);
  TH1D *h1_epsilon_kaons_Uncut = new TH1D("h1_epsilon_kaons_Uncut", "epsilon", 220, -10, 10);

  TH1D *h1_H_gtr_beta_kaons_Cut_All = new TH1D("h1_H_gtr_beta_kaons_Cut_All", "H_gtr_beta", 220, 0.6, 1.4);
  TH1D *h1_H_gtr_xp_kaons_Cut_All = new TH1D("h1_H_gtr_xp_kaons_Cut_All", "H_gtr_xp", 220, -0.6, 0.6);
  TH1D *h1_H_gtr_yp_kaons_Cut_All = new TH1D("h1_H_gtr_yp_kaons_Cut_All", "H_gtr_yp", 220, -0.6, 0.6);
  TH1D *h1_H_gtr_dp_kaons_Cut_All = new TH1D("h1_H_gtr_dp_kaons_Cut_All", "H_gtr_dp", 220, -15, 15);
  TH1D *h1_H_hod_goodscinhit_kaons_Cut_All = new TH1D("h1_H_hod_goodscinhit_kaons_Cut_All", "H_hod_goodscinhit", 220, 0.6, 1.4);
  TH1D *h1_H_hod_goodstarttime_kaons_Cut_All = new TH1D("h1_H_hod_goodstarttime_kaons_Cut_All", "H_hod_goodstarttime", 220, 0.6, 1.4);
  TH1D *h1_H_cal_etotnorm_kaons_Cut_All = new TH1D("h1_H_cal_etotnorm_kaons_Cut_All", "H_cal_etotnorm", 220, 0, 2);
  TH1D *h1_H_cal_etottracknorm_kaons_Cut_All = new TH1D("h1_H_cal_etottracknorm_kaons_Cut_All", "H_cal_etottracknorm", 220, 0, 2);
  TH1D *h1_H_cer_npeSum_kaons_Cut_All = new TH1D("h1_H_cer_npeSum_kaons_Cut_All", "H_cer_npeSum", 220, 0, 20);
  TH1D *h1_H_RFTime_Dist_kaons_Cut_All = new TH1D("h1_H_RFTime_Dist_kaons_Cut_All", "H_RFTime_Dist", 220, -1, 5);
  TH1D *h1_P_gtr_beta_kaons_Cut_All = new TH1D("h1_P_gtr_beta_kaons_Cut_All", "P_gtr_beta", 220, 0.6, 1.4);
  TH1D *h1_P_gtr_xp_kaons_Cut_All = new TH1D("h1_P_gtr_xp_kaons_Cut_All", "P_gtr_xp", 220, -0.6, 0.6);
  TH1D *h1_P_gtr_yp_kaons_Cut_All = new TH1D("h1_P_gtr_yp_kaons_Cut_All", "P_gtr_yp", 220, -0.6, 0.6);
  TH1D *h1_P_gtr_dp_kaons_Cut_All = new TH1D("h1_P_gtr_dp_kaons_Cut_All", "P_gtr_dp", 220, -20, 20);
  TH1D *h1_P_gtr_p_kaons_Cut_All = new TH1D("h1_P_gtr_p_kaons_Cut_All", "P_gtr_p", 220, 4, 8);
  TH1D *h1_P_hod_goodscinhit_kaons_Cut_All = new TH1D("h1_P_hod_goodscinhit_kaons_Cut_All", "P_hod_goodscinhit", 220, 0.6, 1.4);
  TH1D *h1_P_hod_goodstarttime_kaons_Cut_All = new TH1D("h1_P_hod_goodstarttime_kaons_Cut_All", "P_hod_goodstarttime", 220, 0.6, 1.4);
  TH1D *h1_P_cal_etotnorm_kaons_Cut_All = new TH1D("h1_P_cal_etotnorm_kaons_Cut_All", "P_cal_etotnorm", 220, 0, 1.4);
  TH1D *h1_P_cal_etottracknorm_kaons_Cut_All = new TH1D("h1_P_cal_etottracknorm_kaons_Cut_All", "P_cal_etottracknorm", 220, 0, 1.4);
  TH1D *h1_P_hgcer_npeSum_kaons_Cut_All = new TH1D("h1_P_hgcer_npeSum_kaons_Cut_All", "P_hgcer_npeSum", 220, 0, 0.5);
  TH1D *h1_P_hgcer_xAtCer_kaons_Cut_All = new TH1D("h1_P_hgcer_xAtCer_kaons_Cut_All", "P_hgcer_xAtCer", 220, -40, 40);
  TH1D *h1_P_hgcer_yAtCer_kaons_Cut_All = new TH1D("h1_P_hgcer_yAtCer_kaons_Cut_All", "P_hgcer_yAtCer", 220, -40, 40);
  TH1D *h1_P_aero_npeSum_kaons_Cut_All = new TH1D("h1_P_aero_npeSum_kaons_Cut_All", "P_aero_npeSum", 220, 1, 50);
  TH1D *h1_P_aero_xAtAero_kaons_Cut_All = new TH1D("h1_P_aero_xAtAero_kaons_Cut_All", "P_aero_xAtAero", 220, -40, 40);
  TH1D *h1_P_aero_yAtAero_kaons_Cut_All = new TH1D("h1_P_aero_yAtAero_kaons_Cut_All", "P_aero_yAtAero", 220, -40, 40);
  TH1D *h1_P_kin_MMK_kaons_Cut_All = new TH1D("h1_P_kin_MMK_kaons_Cut_All", "P_kin_MMK", 220, 0.5, 1.8);
  TH1D *h1_P_RFTime_Dist_kaons_Cut_All = new TH1D("h1_P_RFTime_Dist_kaons_Cut_All", "P_RFTime_Dist", 220, -1, 5);
  TH1D *h1_CTime_eKCoinTime_ROC1_kaons_Cut_All = new TH1D("h1_CTime_eKCoinTime_ROC1_kaons_Cut_All", "CTime_eKCoinTime_ROC1", 220, -50, 50);
  TH1D *h1_Q2_kaons_Cut_All = new TH1D("h1_Q2_kaons_Cut_All", "Q2", 220, -10, 10);
  TH1D *h1_W_kaons_Cut_All = new TH1D("h1_W_kaons_Cut_All", "W", 220, -4, 10);
  TH1D *h1_epsilon_kaons_Cut_All = new TH1D("h1_epsilon_kaons_Cut_All", "epsilon", 220, -4, 10);

  TH1D *h1_P_gtr_beta_kaons_Cut_Prompt = new TH1D("h1_P_gtr_beta_kaons_Cut_Prompt", "P_gtr_beta", 220, 0.6, 1.4);
  TH1D *h1_P_RFTime_Dist_kaons_Cut_Prompt = new TH1D("h1_P_RFTime_Dist_kaons_Cut_Prompt", "P_RFTime_Dist", 220, -1, 5);
  TH1D *h1_CTime_eKCoinTime_ROC1_kaons_Cut_Prompt = new TH1D("h1_CTime_eKCoinTime_ROC1_kaons_Cut_Prompt", "CTime_eKCoinTime_ROC1", 220, -8, 8);
  TH1D *h1_P_kin_MMK_kaons_Cut_Prompt = new TH1D("h1_P_kin_MMK_kaons_Cut_Prompt", "P_kin_MMK", 220, 0.5, 1.8);

  TH1D *h1_P_gtr_beta_kaons_Cut_Random = new TH1D("h1_P_gtr_beta_kaons_Cut_Random", "P_gtr_beta", 220, 0.6, 1.4);
  TH1D *h1_P_RFTime_Dist_kaons_Cut_Random = new TH1D("h1_P_RFTime_Dist_kaons_Cut_Random", "P_RFTime_Dist", 220, -1, 5);
  TH1D *h1_CTime_eKCoinTime_ROC1_kaons_Cut_Random = new TH1D("h1_CTime_eKCoinTime_ROC1_kaons_Cut_Random", "CTime_eKCoinTime_ROC1", 220, -40, 40);
  TH1D *h1_P_kin_MMK_kaons_Cut_Random = new TH1D("h1_P_kin_MMK_kaons_Cut_Random", "P_kin_MMK", 220, 0.5, 1.8);

  TH1D *h1_CTime_eKCoinTime_ROC1_kaons_Cut_Random_Scaled = new TH1D("h1_CTime_eKCoinTime_ROC1_kaons_Cut_Random_Scaled", "CTime_eKCoinTime_ROC1", 220, -40, 40);
  TH1D *h1_CTime_eKCoinTime_ROC1_kaons_Rndm_Sub = new TH1D("h1_CTime_eKCoinTime_ROC1_kaons_Rndm_Sub", "CTime_eKCoinTime_ROC1", 220, -8, 8);
  TH1D *h1_P_kin_MMK_kaons_Cut_Random_Scaled = new TH1D("h1_P_kin_MMK_kaons_Cut_Random_Scaled", "P_kin_MMK", 220, 0.5, 1.8);
  TH1D *h1_P_kin_MMK_kaons_Rndm_Sub = new TH1D("h1_P_kin_MMK_kaons_Rndm_Sub", "P_kin_MMK", 220, 0.5, 1.8);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Define Histograms for Protons
  TH1D *h1_H_gtr_beta_protons_Uncut = new TH1D("h1_H_gtr_beta_protons_Uncut", "H_gtr_beta", 220, 0.6, 1.4);
  TH1D *h1_H_gtr_xp_protons_Uncut = new TH1D("h1_H_gtr_xp_protons_Uncut", "H_gtr_xp", 220, -0.6, 0.6);
  TH1D *h1_H_gtr_yp_protons_Uncut = new TH1D("h1_H_gtr_yp_protons_Uncut", "H_gtr_yp", 220, -0.6, 0.6);
  TH1D *h1_H_gtr_dp_protons_Uncut = new TH1D("h1_H_gtr_dp_protons_Uncut", "H_gtr_dp", 220, -40, 40);
  TH1D *h1_H_hod_goodscinhit_protons_Uncut = new TH1D("h1_H_hod_goodscinhit_protons_Uncut", "H_hod_goodscinhit", 220, 0.6, 1.4);
  TH1D *h1_H_hod_goodstarttime_protons_Uncut = new TH1D("h1_H_hod_goodstarttime_protons_Uncut", "H_hod_goodstarttime", 220, 0.6, 1.4);
  TH1D *h1_H_cal_etotnorm_protons_Uncut = new TH1D("h1_H_cal_etotnorm_protons_Uncut", "H_cal_etotnorm", 220, 0, 2);
  TH1D *h1_H_cal_etottracknorm_protons_Uncut = new TH1D("h1_H_cal_etottracknorm_protons_Uncut", "H_cal_etottracknorm", 220, 0, 2);
  TH1D *h1_H_cer_npeSum_protons_Uncut = new TH1D("h1_H_cer_npeSum_protons_Uncut", "H_cer_npeSum", 220, 0, 2);
  TH1D *h1_H_RFTime_Dist_protons_Uncut = new TH1D("h1_H_RFTime_Dist_protons_Uncut", "H_RFTime_Dist", 220, -1, 5);
  TH1D *h1_P_gtr_beta_protons_Uncut = new TH1D("h1_P_gtr_beta_protons_Uncut", "P_gtr_beta", 220, 0.6, 1.4);
  TH1D *h1_P_gtr_xp_protons_Uncut = new TH1D("h1_P_gtr_xp_protons_Uncut", "P_gtr_xp", 220, -0.6, 0.6);
  TH1D *h1_P_gtr_yp_protons_Uncut = new TH1D("h1_P_gtr_yp_protons_Uncut", "P_gtr_yp", 220, -0.8, 0.8);
  TH1D *h1_P_gtr_dp_protons_Uncut = new TH1D("h1_P_gtr_dp_protons_Uncut", "P_gtr_dp", 220, -50, 50);
  TH1D *h1_P_gtr_p_protons_Uncut = new TH1D("h1_P_gtr_p_protons_Uncut", "P_gtr_p", 220, -10, 20);
  TH1D *h1_P_hod_goodscinhit_protons_Uncut = new TH1D("h1_P_hod_goodscinhit_protons_Uncut", "P_hod_goodscinhit", 220, 0.6, 1.4);
  TH1D *h1_P_hod_goodstarttime_protons_Uncut = new TH1D("h1_P_hod_goodstarttime_protons_Uncut", "P_hod_goodstarttime", 220, 0.6, 1.4);
  TH1D *h1_P_cal_etotnorm_protons_Uncut = new TH1D("h1_P_cal_etotnorm_protons_Uncut", "P_cal_etotnorm", 220, 0, 1.6);
  TH1D *h1_P_cal_etottracknorm_protons_Uncut = new TH1D("h1_P_cal_etottracknorm_protons_Uncut", "P_cal_etottracknorm", 220, 0, 1.6);
  TH1D *h1_P_hgcer_npeSum_protons_Uncut = new TH1D("h1_P_hgcer_npeSum_protons_Uncut", "P_hgcer_npeSum", 220, 0, 0.5);
  TH1D *h1_P_hgcer_xAtCer_protons_Uncut = new TH1D("h1_P_hgcer_xAtCer_protons_Uncut", "P_hgcer_xAtCer", 220, -10, 10);
  TH1D *h1_P_hgcer_yAtCer_protons_Uncut = new TH1D("h1_P_hgcer_yAtCer_protons_Uncut", "P_hgcer_yAtCer", 220, -10, 10);
  TH1D *h1_P_aero_npeSum_protons_Uncut = new TH1D("h1_P_aero_npeSum_protons_Uncut", "P_aero_npeSum", 220, 0, 10);
  TH1D *h1_P_aero_xAtAero_protons_Uncut = new TH1D("h1_P_aero_xAtAero_protons_Uncut", "P_aero_xAtAero", 220, -10, 10);
  TH1D *h1_P_aero_yAtAero_protons_Uncut = new TH1D("h1_P_aero_yAtAero_protons_Uncut", "P_aero_yAtAero", 220, -10, 10);
  TH1D *h1_P_kin_MMp_protons_Uncut = new TH1D("h1_P_kin_MMp_protons_Uncut", "P_kin_MMp", 220, 0, 2);
  TH1D *h1_P_RFTime_Dist_protons_Uncut = new TH1D("h1_P_RFTime_Dist_protons_Uncut", "P_RFTime_Dist", 220, -1, 5);
  TH1D *h1_CTime_epCoinTime_ROC1_protons_Uncut = new TH1D("h1_CTime_epCoinTime_ROC1_protons_Uncut", "CTime_epCoinTime_ROC1", 220, -100, 100);
  TH1D *h1_Q2_protons_Uncut = new TH1D("h1_Q2_protons_Uncut", "Q2", 220, -10, 10);
  TH1D *h1_W_protons_Uncut = new TH1D("h1_W_protons_Uncut", "W", 220, -10, 10);
  TH1D *h1_epsilon_protons_Uncut = new TH1D("h1_epsilon_protons_Uncut", "epsilon", 220, -10, 10);

  TH1D *h1_H_gtr_beta_protons_Cut_All = new TH1D("h1_H_gtr_beta_protons_Cut_All", "H_gtr_beta", 220, 0.6, 1.4);
  TH1D *h1_H_gtr_xp_protons_Cut_All = new TH1D("h1_H_gtr_xp_protons_Cut_All", "H_gtr_xp", 220, -0.6, 0.6);
  TH1D *h1_H_gtr_yp_protons_Cut_All = new TH1D("h1_H_gtr_yp_protons_Cut_All", "H_gtr_yp", 220, -0.6, 0.6);
  TH1D *h1_H_gtr_dp_protons_Cut_All = new TH1D("h1_H_gtr_dp_protons_Cut_All", "H_gtr_dp", 220, -15, 15);
  TH1D *h1_H_hod_goodscinhit_protons_Cut_All = new TH1D("h1_H_hod_goodscinhit_protons_Cut_All", "H_hod_goodscinhit", 220, 0.6, 1.4);
  TH1D *h1_H_hod_goodstarttime_protons_Cut_All = new TH1D("h1_H_hod_goodstarttime_protons_Cut_All", "H_hod_goodstarttime", 220, 0.6, 1.4);
  TH1D *h1_H_cal_etotnorm_protons_Cut_All = new TH1D("h1_H_cal_etotnorm_protons_Cut_All", "H_cal_etotnorm", 220, 0, 2);
  TH1D *h1_H_cal_etottracknorm_protons_Cut_All = new TH1D("h1_H_cal_etottracknorm_protons_Cut_All", "H_cal_etottracknorm", 220, 0, 2);
  TH1D *h1_H_cer_npeSum_protons_Cut_All = new TH1D("h1_H_cer_npeSum_protons_Cut_All", "H_cer_npeSum", 220, 0, 20);
  TH1D *h1_H_RFTime_Dist_protons_Cut_All = new TH1D("h1_H_RFTime_Dist_protons_Cut_All", "H_RFTime_Dist", 220, -1, 5);
  TH1D *h1_P_gtr_beta_protons_Cut_All = new TH1D("h1_P_gtr_beta_protons_Cut_All", "P_gtr_beta", 220, 0.6, 1.4);
  TH1D *h1_P_gtr_xp_protons_Cut_All = new TH1D("h1_P_gtr_xp_protons_Cut_All", "P_gtr_xp", 220, -0.6, 0.6);
  TH1D *h1_P_gtr_yp_protons_Cut_All = new TH1D("h1_P_gtr_yp_protons_Cut_All", "P_gtr_yp", 220, -0.6, 0.6);
  TH1D *h1_P_gtr_dp_protons_Cut_All = new TH1D("h1_P_gtr_dp_protons_Cut_All", "P_gtr_dp", 220, -20, 20);
  TH1D *h1_P_gtr_p_protons_Cut_All = new TH1D("h1_P_gtr_p_protons_Cut_All", "P_gtr_p", 220, 0, 12);
  TH1D *h1_P_hod_goodscinhit_protons_Cut_All = new TH1D("h1_P_hod_goodscinhit_protons_Cut_All", "P_hod_goodscinhit", 220, 0.6, 1.4);
  TH1D *h1_P_hod_goodstarttime_protons_Cut_All = new TH1D("h1_P_hod_goodstarttime_protons_Cut_All", "P_hod_goodstarttime", 220, 0.6, 1.4);
  TH1D *h1_P_cal_etotnorm_protons_Cut_All = new TH1D("h1_P_cal_etotnorm_protons_Cut_All", "P_cal_etotnorm", 220, 0, 1.5);
  TH1D *h1_P_cal_etottracknorm_protons_Cut_All = new TH1D("h1_P_cal_etottracknorm_protons_Cut_All", "P_cal_etottracknorm", 220, 0, 1.5);
  TH1D *h1_P_hgcer_npeSum_protons_Cut_All = new TH1D("h1_P_hgcer_npeSum_protons_Cut_All", "P_hgcer_npeSum", 220, 0, 0.5);
  TH1D *h1_P_hgcer_xAtCer_protons_Cut_All = new TH1D("h1_P_hgcer_xAtCer_protons_Cut_All", "P_hgcer_xAtCer", 220, -40, 40);
  TH1D *h1_P_hgcer_yAtCer_protons_Cut_All = new TH1D("h1_P_hgcer_yAtCer_protons_Cut_All", "P_hgcer_yAtCer", 220, -40, 40);
  TH1D *h1_P_aero_npeSum_protons_Cut_All = new TH1D("h1_P_aero_npeSum_protons_Cut_All", "P_aero_npeSum", 220, 1, 50);
  TH1D *h1_P_aero_xAtAero_protons_Cut_All = new TH1D("h1_P_aero_xAtAero_protons_Cut_All", "P_aero_xAtAero", 220, -40, 40);
  TH1D *h1_P_aero_yAtAero_protons_Cut_All = new TH1D("h1_P_aero_yAtAero_protons_Cut_All", "P_aero_yAtAero", 220, -40, 40);
  TH1D *h1_P_kin_MMp_protons_Cut_All = new TH1D("h1_P_kin_MMp_protons_Cut_All", "P_kin_MMp", 220, 0.5, 1.8);
  TH1D *h1_P_RFTime_Dist_protons_Cut_All = new TH1D("h1_P_RFTime_Dist_protons_Cut_All", "P_RFTime_Dist", 220, -1, 5);
  TH1D *h1_CTime_epCoinTime_ROC1_protons_Cut_All = new TH1D("h1_CTime_epCoinTime_ROC1_protons_Cut_All", "CTime_epCoinTime_ROC1", 220, -50, 50);
  TH1D *h1_Q2_protons_Cut_All = new TH1D("h1_Q2_protons_Cut_All", "Q2", 220, -10, 10);
  TH1D *h1_W_protons_Cut_All = new TH1D("h1_W_protons_Cut_All", "W", 220, -4, 10);
  TH1D *h1_epsilon_protons_Cut_All = new TH1D("h1_epsilon_protons_Cut_All", "epsilon", 220, -4, 10);

  TH1D *h1_P_gtr_beta_protons_Cut_Prompt = new TH1D("h1_P_gtr_beta_protons_Cut_Prompt", "P_gtr_beta", 220, 0.6, 1.4);
  TH1D *h1_P_RFTime_Dist_protons_Cut_Prompt = new TH1D("h1_P_RFTime_Dist_protons_Cut_Prompt", "P_RFTime_Dist", 220, 0, 4);
  TH1D *h1_CTime_epCoinTime_ROC1_protons_Cut_Prompt = new TH1D("h1_CTime_epCoinTime_ROC1_protons_Cut_Prompt", "CTime_epCoinTime_ROC1", 220, -8, 8);
  TH1D *h1_P_kin_MMp_protons_Cut_Prompt = new TH1D("h1_P_kin_MMp_protons_Cut_Prompt", "P_kin_MMp", 220, 0.5, 1.8);

  TH1D *h1_P_gtr_beta_protons_Cut_Random = new TH1D("h1_P_gtr_beta_protons_Cut_Random", "P_gtr_beta", 220, 0.6, 1.4);
  TH1D *h1_P_RFTime_Dist_protons_Cut_Random = new TH1D("h1_P_RFTime_Dist_protons_Cut_Random", "P_RFTime_Dist", 220, -1, 5);
  TH1D *h1_CTime_epCoinTime_ROC1_protons_Cut_Random = new TH1D("h1_CTime_epCoinTime_ROC1_protons_Cut_Random", "CTime_epCoinTime_ROC1", 220, -40, 40);
  TH1D *h1_P_kin_MMp_protons_Cut_Random = new TH1D("h1_P_kin_MMp_protons_Cut_Random", "P_kin_MMp", 220, 0.5, 1.8);

  TH1D *h1_CTime_epCoinTime_ROC1_protons_Cut_Random_Scaled = new TH1D("h1_CTime_epCoinTime_ROC1_protons_Cut_Random_Scaled", "CTime_epCoinTime_ROC1", 220, -40, 40);
  TH1D *h1_CTime_epCoinTime_ROC1_protons_Rndm_Sub = new TH1D("h1_CTime_epCoinTime_ROC1_protons_Rndm_Sub", "CTime_epCoinTime_ROC1", 220, -8, 8);
  TH1D *h1_P_kin_MMp_protons_Cut_Random_Scaled = new TH1D("h1_P_kin_MMp_protons_Cut_Random_Scaled", "P_kin_MMp", 220, 0.5, 1.8);
  TH1D *h1_P_kin_MMp_protons_Rndm_Sub = new TH1D("h1_P_kin_MMp_protons_Rndm_Sub", "P_kin_MMp", 220, 0.5, 1.8);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TH2D *h2_P_hgcer_npeSum_vs_aero_npeSum_pions_Uncut = new TH2D("h2_P_hgcer_npeSum_vs_aero_npeSum_pions_Uncut","P_hgcer_npeSum vs P_aero_npeSum",100,0,50,100,0,50);
  TH2D *h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Uncut = new TH2D("h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Uncut","CTime_ePiCoinTime_ROC1 vs P_kin_MMpi",100,-2,2,100,0,3);
  TH2D *h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_pions_Uncut = new TH2D("h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_pions_Uncut","P_hgcer_xAtCer vs P_hgcer_yAtCer",100,-10,10,100,-10,10);

  TH2D *h2_P_hgcer_npeSum_vs_aero_npeSum_pions_Cut_All = new TH2D("h2_P_hgcer_npeSum_vs_aero_npeSum_pions_Cut_All","P_hgcer_npeSum vs P_aero_npeSum",100,0,50,100,0,50);
  TH2D *h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Cut_All = new TH2D("h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Cut_All","CTime_ePiCoinTime_ROC1 vs P_kin_MMpi",100,-1,1,100,0.8,1.6);
  TH2D *h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_pions_Cut_All = new TH2D("h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_pions_Cut_All","P_hgcer_xAtCer vs P_hgcer_yAtCer",100,-32,32,100,-30,30);

//  TH2D *h2_P_hgcer_npeSum_vs_aero_npeSum_pions_Cut_Prompt = new TH2D("h2_P_hgcer_npeSum_vs_aero_npeSum_pions_Cut_Prompt","P_hgcer_npeSum vs P_aero_npeSum",100,0,50,100,0,50);
  TH2D *h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Cut_Prompt = new TH2D("h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Cut_Prompt","CTime_ePiCoinTime_ROC1 vs P_kin_MMpi",100,-1,1,100,0.8,1.6);

  TH2D *h2_P_hgcer_npeSum_vs_aero_npeSum_kaons_Uncut = new TH2D("h2_P_hgcer_npeSum_vs_aero_npeSum_kaons_Uncut","P_hgcer_npeSum vs P_aero_npeSum",100,0,50,100,0,50);
  TH2D *h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Uncut = new TH2D("h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Uncut","CTime_eKCoinTime_ROC1 vs P_kin_MMK",100,-2,2,100,0,3);
  TH2D *h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_kaons_Uncut = new TH2D("h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_kaons_Uncut","P_hgcer_xAtCer vs P_hgcer_yAtCer",100,-10,10,100,-10,10);

  TH2D *h2_P_hgcer_npeSum_vs_aero_npeSum_kaons_Cut_All = new TH2D("h2_P_hgcer_npeSum_vs_aero_npeSum_kaons_Cut_All","P_hgcer_npeSum vs P_aero_npeSum",100,0,50,100,0,50);
  TH2D *h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Cut_All = new TH2D("h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Cut_All","CTime_eKCoinTime_ROC1 vs P_kin_MMK",100,-1,1,100,0.6,1.6);
  TH2D *h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_kaons_Cut_All = new TH2D("h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_kaons_Cut_All","P_hgcer_xAtCer vs P_hgcer_yAtCer",100,-32,32,100,-30,30);

//  TH2D *h2_P_hgcer_npeSum_vs_aero_npeSum_kaons_Cut_Prompt = new TH2D("h2_P_hgcer_npeSum_vs_aero_npeSum_kaons_Cut_Prompt","P_hgcer_npeSum vs P_aero_npeSum",100,0,50,100,0,50);
  TH2D *h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Cut_Prompt = new TH2D("h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Cut_Prompt","CTime_eKCoinTime_ROC1 vs P_kin_MMK",100,-1,1,100,0.6,1.6);

  TH2D *h2_P_hgcer_npeSum_vs_aero_npeSum_protons_Uncut = new TH2D("h2_P_hgcer_npeSum_vs_aero_npeSum_protons_Uncut","P_hgcer_npeSum vs P_aero_npeSum",100,0,50,100,0,50);
  TH2D *h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Uncut = new TH2D("h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Uncut","CTime_epCoinTime_ROC1 vs P_kin_MMp",100,-2,2,100,0,3);
  TH2D *h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_protons_Uncut = new TH2D("h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_protons_Uncut","P_hgcer_xAtCer vs P_hgcer_yAtCer",100,-10,10,100,-10,10);

  TH2D *h2_P_hgcer_npeSum_vs_aero_npeSum_protons_Cut_All = new TH2D("h2_P_hgcer_npeSum_vs_aero_npeSum_protons_Cut_All","P_hgcer_npeSum vs P_aero_npeSum",100,0,50,100,0,50);
  TH2D *h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Cut_All = new TH2D("h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Cut_All","CTime_epCoinTime_ROC1 vs P_kin_MMp",100,-1,1,100,0.4,1.6);
  TH2D *h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_protons_Cut_All = new TH2D("h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_protons_Cut_All","P_hgcer_xAtCer vs P_hgcer_yAtCer",100,-32,32,100,-30,30);

//  TH2D *h2_P_hgcer_npeSum_vs_aero_npeSum_protons_Cut_Prompt = new TH2D("h2_P_hgcer_npeSum_vs_aero_npeSum_protons_Cut_Prompt","P_hgcer_npeSum vs P_aero_npeSum",100,0,50,100,0,50);
  TH2D *h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Cut_Prompt = new TH2D("h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Cut_Prompt","CTime_epCoinTime_ROC1 vs P_kin_MMp",100,-1,1,100,0.4,1.6);

/////////////////////////////////////////////////////////////////////////////////////////////

  // Loop over pion events in tree 
  for(Long64_t i = 0; i < nEntries_Uncut_Pion_Events; i++){
    Uncut_Pion_Events->GetEntry(i);

    h1_H_gtr_beta_pions_Uncut->Fill(H_gtr_beta_pions_uncut);
    h1_H_gtr_beta_pions_Uncut->GetXaxis()->SetTitle("H_gtr_beta");
    h1_H_gtr_beta_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_beta_pions_Uncut->GetXaxis()->CenterTitle();
    h1_H_gtr_beta_pions_Uncut->GetYaxis()->CenterTitle();
    h1_H_gtr_xp_pions_Uncut->Fill(H_gtr_xp_pions_uncut);
    h1_H_gtr_xp_pions_Uncut->GetXaxis()->SetTitle("H_gtr_xp");
    h1_H_gtr_xp_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_xp_pions_Uncut->GetXaxis()->CenterTitle();
    h1_H_gtr_xp_pions_Uncut->GetYaxis()->CenterTitle();
    h1_H_gtr_yp_pions_Uncut->Fill(H_gtr_yp_pions_uncut);
    h1_H_gtr_yp_pions_Uncut->GetXaxis()->SetTitle("H_gtr_yp");
    h1_H_gtr_yp_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_yp_pions_Uncut->GetXaxis()->CenterTitle();
    h1_H_gtr_yp_pions_Uncut->GetYaxis()->CenterTitle();
    h1_H_gtr_dp_pions_Uncut->Fill(H_gtr_dp_pions_uncut);
    h1_H_gtr_dp_pions_Uncut->GetXaxis()->SetTitle("H_gtr_dp");
    h1_H_gtr_dp_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_dp_pions_Uncut->GetXaxis()->CenterTitle();
    h1_H_gtr_dp_pions_Uncut->GetYaxis()->CenterTitle();
    h1_H_hod_goodscinhit_pions_Uncut->Fill(H_hod_goodscinhit_pions_uncut);
    h1_H_hod_goodscinhit_pions_Uncut->GetXaxis()->SetTitle("H_hod_goodscinhit");
    h1_H_hod_goodscinhit_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_hod_goodscinhit_pions_Uncut->GetXaxis()->CenterTitle();
    h1_H_hod_goodscinhit_pions_Uncut->GetYaxis()->CenterTitle();
    h1_H_hod_goodstarttime_pions_Uncut->Fill(H_hod_goodstarttime_pions_uncut);
    h1_H_hod_goodstarttime_pions_Uncut->GetXaxis()->SetTitle("H_hod_goodstarttime");
    h1_H_hod_goodstarttime_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_hod_goodstarttime_pions_Uncut->GetXaxis()->CenterTitle();
    h1_H_hod_goodstarttime_pions_Uncut->GetYaxis()->CenterTitle();
    h1_H_cal_etotnorm_pions_Uncut->Fill(H_cal_etotnorm_pions_uncut);
    h1_H_cal_etotnorm_pions_Uncut->GetXaxis()->SetTitle("H_cal_etotnorm");
    h1_H_cal_etotnorm_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_cal_etotnorm_pions_Uncut->GetXaxis()->CenterTitle();
    h1_H_cal_etotnorm_pions_Uncut->GetYaxis()->CenterTitle();
    h1_H_cal_etottracknorm_pions_Uncut->Fill(H_cal_etottracknorm_pions_uncut);
    h1_H_cal_etottracknorm_pions_Uncut->GetXaxis()->SetTitle("H_cal_etottracknorm");
    h1_H_cal_etottracknorm_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_cal_etottracknorm_pions_Uncut->GetXaxis()->CenterTitle();
    h1_H_cal_etottracknorm_pions_Uncut->GetYaxis()->CenterTitle();
    h1_H_cer_npeSum_pions_Uncut->Fill(H_cer_npeSum_pions_uncut);
    h1_H_cer_npeSum_pions_Uncut->GetXaxis()->SetTitle("H_cer_npeSum");
    h1_H_cer_npeSum_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_cer_npeSum_pions_Uncut->GetXaxis()->CenterTitle();
    h1_H_cer_npeSum_pions_Uncut->GetYaxis()->CenterTitle();
    h1_H_RFTime_Dist_pions_Uncut->Fill(H_RFTime_Dist_pions_uncut);
    h1_H_RFTime_Dist_pions_Uncut->GetXaxis()->SetTitle("H_RFTime_Dist");
    h1_H_RFTime_Dist_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_RFTime_Dist_pions_Uncut->GetXaxis()->CenterTitle();
    h1_H_RFTime_Dist_pions_Uncut->GetYaxis()->CenterTitle();

    h1_P_gtr_beta_pions_Uncut->Fill(P_gtr_beta_pions_uncut);
    h1_P_gtr_beta_pions_Uncut->GetXaxis()->SetTitle("P_gtr_beta");
    h1_P_gtr_beta_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_beta_pions_Uncut->GetXaxis()->CenterTitle();
    h1_P_gtr_beta_pions_Uncut->GetYaxis()->CenterTitle();
    h1_P_gtr_xp_pions_Uncut->Fill(P_gtr_xp_pions_uncut);
    h1_P_gtr_xp_pions_Uncut->GetXaxis()->SetTitle("P_gtr_xp");
    h1_P_gtr_xp_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_xp_pions_Uncut->GetXaxis()->CenterTitle();
    h1_P_gtr_xp_pions_Uncut->GetYaxis()->CenterTitle();
    h1_P_gtr_yp_pions_Uncut->Fill(P_gtr_yp_pions_uncut);
    h1_P_gtr_yp_pions_Uncut->GetXaxis()->SetTitle("P_gtr_yp");
    h1_P_gtr_yp_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_yp_pions_Uncut->GetXaxis()->CenterTitle();
    h1_P_gtr_yp_pions_Uncut->GetYaxis()->CenterTitle();
    h1_P_gtr_dp_pions_Uncut->Fill(P_gtr_dp_pions_uncut);
    h1_P_gtr_dp_pions_Uncut->GetXaxis()->SetTitle("P_gtr_dp");
    h1_P_gtr_dp_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_dp_pions_Uncut->GetXaxis()->CenterTitle();
    h1_P_gtr_dp_pions_Uncut->GetYaxis()->CenterTitle();
    h1_P_gtr_p_pions_Uncut->Fill(P_gtr_p_pions_uncut);
    h1_P_gtr_p_pions_Uncut->GetXaxis()->SetTitle("P_gtr_p");
    h1_P_gtr_p_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_p_pions_Uncut->GetXaxis()->CenterTitle();
    h1_P_gtr_p_pions_Uncut->GetYaxis()->CenterTitle();
    h1_P_hod_goodscinhit_pions_Uncut->Fill(P_hod_goodscinhit_pions_uncut);
    h1_P_hod_goodscinhit_pions_Uncut->GetXaxis()->SetTitle("P_hod_goodscinhit");
    h1_P_hod_goodscinhit_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_hod_goodscinhit_pions_Uncut->GetXaxis()->CenterTitle();
    h1_P_hod_goodscinhit_pions_Uncut->GetYaxis()->CenterTitle();
    h1_P_hod_goodstarttime_pions_Uncut->Fill(P_hod_goodstarttime_pions_uncut);
    h1_P_hod_goodstarttime_pions_Uncut->GetXaxis()->SetTitle("P_hod_goodstarttime");
    h1_P_hod_goodstarttime_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_hod_goodstarttime_pions_Uncut->GetXaxis()->CenterTitle();
    h1_P_hod_goodstarttime_pions_Uncut->GetYaxis()->CenterTitle();
    h1_P_cal_etotnorm_pions_Uncut->Fill(P_cal_etotnorm_pions_uncut);
    h1_P_cal_etotnorm_pions_Uncut->GetXaxis()->SetTitle("P_cal_etotnorm");
    h1_P_cal_etotnorm_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_cal_etotnorm_pions_Uncut->GetXaxis()->CenterTitle();
    h1_P_cal_etotnorm_pions_Uncut->GetYaxis()->CenterTitle();
    h1_P_cal_etottracknorm_pions_Uncut->Fill(P_cal_etottracknorm_pions_uncut);
    h1_P_cal_etottracknorm_pions_Uncut->GetXaxis()->SetTitle("P_cal_etottracknorm");
    h1_P_cal_etottracknorm_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_cal_etottracknorm_pions_Uncut->GetXaxis()->CenterTitle();
    h1_P_cal_etottracknorm_pions_Uncut->GetYaxis()->CenterTitle();
    h1_P_hgcer_npeSum_pions_Uncut->Fill(P_hgcer_npeSum_pions_uncut);
    h1_P_hgcer_npeSum_pions_Uncut->GetXaxis()->SetTitle("P_hgcer_npeSum");
    h1_P_hgcer_npeSum_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_hgcer_npeSum_pions_Uncut->GetXaxis()->CenterTitle();
    h1_P_hgcer_npeSum_pions_Uncut->GetYaxis()->CenterTitle();
    h1_P_hgcer_xAtCer_pions_Uncut->Fill(P_hgcer_xAtCer_pions_uncut);
    h1_P_hgcer_xAtCer_pions_Uncut->GetXaxis()->SetTitle("P_hgcer_xAtCer");
    h1_P_hgcer_xAtCer_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_hgcer_xAtCer_pions_Uncut->GetXaxis()->CenterTitle();
    h1_P_hgcer_xAtCer_pions_Uncut->GetYaxis()->CenterTitle();
    h1_P_hgcer_yAtCer_pions_Uncut->Fill(P_hgcer_yAtCer_pions_uncut);
    h1_P_hgcer_yAtCer_pions_Uncut->GetXaxis()->SetTitle("P_hgcer_yAtCer");
    h1_P_hgcer_yAtCer_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_hgcer_yAtCer_pions_Uncut->GetXaxis()->CenterTitle();
    h1_P_hgcer_yAtCer_pions_Uncut->GetYaxis()->CenterTitle();
    h1_P_aero_npeSum_pions_Uncut->Fill(P_aero_npeSum_pions_uncut);
    h1_P_aero_npeSum_pions_Uncut->GetXaxis()->SetTitle("P_aero_npeSum");
    h1_P_aero_npeSum_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_aero_npeSum_pions_Uncut->GetXaxis()->CenterTitle();
    h1_P_aero_npeSum_pions_Uncut->GetYaxis()->CenterTitle();
    h1_P_aero_xAtAero_pions_Uncut->Fill(P_aero_xAtAero_pions_uncut);
    h1_P_aero_xAtAero_pions_Uncut->GetXaxis()->SetTitle("P_aero_xAtAero");
    h1_P_aero_xAtAero_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_aero_xAtAero_pions_Uncut->GetXaxis()->CenterTitle();
    h1_P_aero_xAtAero_pions_Uncut->GetYaxis()->CenterTitle();
    h1_P_aero_yAtAero_pions_Uncut->Fill(P_aero_yAtAero_pions_uncut);
    h1_P_aero_yAtAero_pions_Uncut->GetXaxis()->SetTitle("P_aero_yAtAero");
    h1_P_aero_yAtAero_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_aero_yAtAero_pions_Uncut->GetXaxis()->CenterTitle();
    h1_P_aero_yAtAero_pions_Uncut->GetYaxis()->CenterTitle();
    h1_P_kin_MMpi_pions_Uncut->Fill(P_kin_MMpi_pions_uncut);
    h1_P_kin_MMpi_pions_Uncut->GetXaxis()->SetTitle("Missing_Mass_Pions_uncut");
    h1_P_kin_MMpi_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_kin_MMpi_pions_Uncut->GetXaxis()->CenterTitle();
    h1_P_kin_MMpi_pions_Uncut->GetYaxis()->CenterTitle();
    h1_P_RFTime_Dist_pions_Uncut->Fill(P_RFTime_Dist_pions_uncut);
    h1_P_RFTime_Dist_pions_Uncut->GetXaxis()->SetTitle("P_RFTime_Dist");
    h1_P_RFTime_Dist_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_RFTime_Dist_pions_Uncut->GetXaxis()->CenterTitle();
    h1_P_RFTime_Dist_pions_Uncut->GetYaxis()->CenterTitle();
    h1_CTime_ePiCoinTime_ROC1_pions_Uncut->Fill(CTime_ePiCoinTime_ROC1_pions_uncut);
    h1_CTime_ePiCoinTime_ROC1_pions_Uncut->GetXaxis()->SetTitle("CTime_ePiCoinTime_ROC1");
    h1_CTime_ePiCoinTime_ROC1_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_CTime_ePiCoinTime_ROC1_pions_Uncut->GetXaxis()->CenterTitle();
    h1_CTime_ePiCoinTime_ROC1_pions_Uncut->GetYaxis()->CenterTitle();
    h1_Q2_pions_Uncut->Fill(Q2_pions_uncut);
    h1_Q2_pions_Uncut->GetXaxis()->SetTitle("Q2");
    h1_Q2_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_Q2_pions_Uncut->GetXaxis()->CenterTitle();
    h1_Q2_pions_Uncut->GetYaxis()->CenterTitle();
    h1_W_pions_Uncut->Fill(W_pions_uncut);
    h1_W_pions_Uncut->GetXaxis()->SetTitle("W");
    h1_W_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_W_pions_Uncut->GetXaxis()->CenterTitle();
    h1_W_pions_Uncut->GetYaxis()->CenterTitle();
    h1_epsilon_pions_Uncut->Fill(epsilon_pions_uncut);
    h1_epsilon_pions_Uncut->GetXaxis()->SetTitle("epsilon");
    h1_epsilon_pions_Uncut->GetYaxis()->SetTitle("Entries");
    h1_epsilon_pions_Uncut->GetXaxis()->CenterTitle();
    h1_epsilon_pions_Uncut->GetYaxis()->CenterTitle();
  }

  for(Long64_t i = 0; i < nEntries_Cut_Pion_Events_All; i++){
    Cut_Pion_Events_All->GetEntry(i);

    h1_H_gtr_beta_pions_Cut_All->Fill(H_gtr_beta_pions_cut_all);
    h1_H_gtr_beta_pions_Cut_All->GetXaxis()->SetTitle("H_gtr_beta");
    h1_H_gtr_beta_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_beta_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_H_gtr_beta_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_H_gtr_xp_pions_Cut_All->Fill(H_gtr_xp_pions_cut_all);
    h1_H_gtr_xp_pions_Cut_All->GetXaxis()->SetTitle("H_gtr_xp");
    h1_H_gtr_xp_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_xp_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_H_gtr_xp_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_H_gtr_yp_pions_Cut_All->Fill(H_gtr_yp_pions_cut_all);
    h1_H_gtr_yp_pions_Cut_All->GetXaxis()->SetTitle("H_gtr_yp");
    h1_H_gtr_yp_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_yp_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_H_gtr_yp_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_H_gtr_dp_pions_Cut_All->Fill(H_gtr_dp_pions_cut_all);
    h1_H_gtr_dp_pions_Cut_All->GetXaxis()->SetTitle("H_gtr_dp");
    h1_H_gtr_dp_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_dp_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_H_gtr_dp_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_H_hod_goodscinhit_pions_Cut_All->Fill(H_hod_goodscinhit_pions_cut_all);
    h1_H_hod_goodscinhit_pions_Cut_All->GetXaxis()->SetTitle("H_hod_goodscinhit");
    h1_H_hod_goodscinhit_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_hod_goodscinhit_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_H_hod_goodscinhit_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_H_hod_goodstarttime_pions_Cut_All->Fill(H_hod_goodstarttime_pions_cut_all);
    h1_H_hod_goodstarttime_pions_Cut_All->GetXaxis()->SetTitle("H_hod_goodstarttime");
    h1_H_hod_goodstarttime_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_hod_goodstarttime_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_H_hod_goodstarttime_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_H_cal_etotnorm_pions_Cut_All->Fill(H_cal_etotnorm_pions_cut_all);
    h1_H_cal_etotnorm_pions_Cut_All->GetXaxis()->SetTitle("H_cal_etotnorm");
    h1_H_cal_etotnorm_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_cal_etotnorm_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_H_cal_etotnorm_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_H_cal_etottracknorm_pions_Cut_All->Fill(H_cal_etottracknorm_pions_cut_all);
    h1_H_cal_etottracknorm_pions_Cut_All->GetXaxis()->SetTitle("H_cal_etottracknorm");
    h1_H_cal_etottracknorm_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_cal_etottracknorm_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_H_cal_etottracknorm_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_H_cer_npeSum_pions_Cut_All->Fill(H_cer_npeSum_pions_cut_all);
    h1_H_cer_npeSum_pions_Cut_All->GetXaxis()->SetTitle("H_cer_npeSum");
    h1_H_cer_npeSum_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_cer_npeSum_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_H_cer_npeSum_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_H_RFTime_Dist_pions_Cut_All->Fill(H_RFTime_Dist_pions_cut_all);
    h1_H_RFTime_Dist_pions_Cut_All->GetXaxis()->SetTitle("H_RFTime_Dist");
    h1_H_RFTime_Dist_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_RFTime_Dist_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_H_RFTime_Dist_pions_Cut_All->GetYaxis()->CenterTitle();

    h1_P_gtr_beta_pions_Cut_All->Fill(P_gtr_beta_pions_cut_all);
    h1_P_gtr_beta_pions_Cut_All->GetXaxis()->SetTitle("P_gtr_beta");
    h1_P_gtr_beta_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_beta_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_P_gtr_beta_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_P_gtr_xp_pions_Cut_All->Fill(P_gtr_xp_pions_cut_all);
    h1_P_gtr_xp_pions_Cut_All->GetXaxis()->SetTitle("P_gtr_xp");
    h1_P_gtr_xp_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_xp_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_P_gtr_xp_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_P_gtr_yp_pions_Cut_All->Fill(P_gtr_yp_pions_cut_all);
    h1_P_gtr_yp_pions_Cut_All->GetXaxis()->SetTitle("P_gtr_yp");
    h1_P_gtr_yp_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_yp_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_P_gtr_yp_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_P_gtr_dp_pions_Cut_All->Fill(P_gtr_dp_pions_cut_all);
    h1_P_gtr_dp_pions_Cut_All->GetXaxis()->SetTitle("P_gtr_dp");
    h1_P_gtr_dp_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_dp_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_P_gtr_dp_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_P_gtr_p_pions_Cut_All->Fill(P_gtr_p_pions_cut_all);
    h1_P_gtr_p_pions_Cut_All->GetXaxis()->SetTitle("P_gtr_p");
    h1_P_gtr_p_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_p_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_P_gtr_p_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_P_hod_goodscinhit_pions_Cut_All->Fill(P_hod_goodscinhit_pions_cut_all);
    h1_P_hod_goodscinhit_pions_Cut_All->GetXaxis()->SetTitle("P_hod_goodscinhit");
    h1_P_hod_goodscinhit_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_hod_goodscinhit_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_P_hod_goodscinhit_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_P_hod_goodstarttime_pions_Cut_All->Fill(P_hod_goodstarttime_pions_cut_all);
    h1_P_hod_goodstarttime_pions_Cut_All->GetXaxis()->SetTitle("P_hod_goodstarttime");
    h1_P_hod_goodstarttime_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_hod_goodstarttime_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_P_hod_goodstarttime_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_P_cal_etotnorm_pions_Cut_All->Fill(P_cal_etotnorm_pions_cut_all);
    h1_P_cal_etotnorm_pions_Cut_All->GetXaxis()->SetTitle("P_cal_etotnorm");
    h1_P_cal_etotnorm_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_cal_etotnorm_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_P_cal_etotnorm_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_P_cal_etottracknorm_pions_Cut_All->Fill(P_cal_etottracknorm_pions_cut_all);
    h1_P_cal_etottracknorm_pions_Cut_All->GetXaxis()->SetTitle("P_cal_etottracknorm");
    h1_P_cal_etottracknorm_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_cal_etottracknorm_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_P_cal_etottracknorm_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_P_hgcer_npeSum_pions_Cut_All->Fill(P_hgcer_npeSum_pions_cut_all);
    h1_P_hgcer_npeSum_pions_Cut_All->GetXaxis()->SetTitle("P_hgcer_npeSum");
    h1_P_hgcer_npeSum_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_hgcer_npeSum_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_P_hgcer_npeSum_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_P_hgcer_xAtCer_pions_Cut_All->Fill(P_hgcer_xAtCer_pions_cut_all);
    h1_P_hgcer_xAtCer_pions_Cut_All->GetXaxis()->SetTitle("P_hgcer_xAtCer");
    h1_P_hgcer_xAtCer_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_hgcer_xAtCer_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_P_hgcer_xAtCer_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_P_hgcer_yAtCer_pions_Cut_All->Fill(P_hgcer_yAtCer_pions_cut_all);
    h1_P_hgcer_yAtCer_pions_Cut_All->GetXaxis()->SetTitle("P_hgcer_yAtCer");
    h1_P_hgcer_yAtCer_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_hgcer_yAtCer_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_P_hgcer_yAtCer_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_P_aero_npeSum_pions_Cut_All->Fill(P_aero_npeSum_pions_cut_all);
    h1_P_aero_npeSum_pions_Cut_All->GetXaxis()->SetTitle("P_aero_npeSum");
    h1_P_aero_npeSum_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_aero_npeSum_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_P_aero_npeSum_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_P_aero_xAtAero_pions_Cut_All->Fill(P_aero_xAtAero_pions_cut_all);
    h1_P_aero_xAtAero_pions_Cut_All->GetXaxis()->SetTitle("P_aero_xAtAero");
    h1_P_aero_xAtAero_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_aero_xAtAero_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_P_aero_xAtAero_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_P_aero_yAtAero_pions_Cut_All->Fill(P_aero_yAtAero_pions_cut_all);
    h1_P_aero_yAtAero_pions_Cut_All->GetXaxis()->SetTitle("P_aero_yAtAero");
    h1_P_aero_yAtAero_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_aero_yAtAero_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_P_aero_yAtAero_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_P_kin_MMpi_pions_Cut_All->Fill(P_kin_MMpi_pions_cut_all);
    h1_P_kin_MMpi_pions_Cut_All->GetXaxis()->SetTitle("Missing_Mass_Pions_cut_all");
    h1_P_kin_MMpi_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_kin_MMpi_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_P_kin_MMpi_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_P_RFTime_Dist_pions_Cut_All->Fill(P_RFTime_Dist_pions_cut_all);
    h1_P_RFTime_Dist_pions_Cut_All->GetXaxis()->SetTitle("P_RFTime_Dist");
    h1_P_RFTime_Dist_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_RFTime_Dist_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_P_RFTime_Dist_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_CTime_ePiCoinTime_ROC1_pions_Cut_All->Fill(CTime_ePiCoinTime_ROC1_pions_cut_all);
    h1_CTime_ePiCoinTime_ROC1_pions_Cut_All->GetXaxis()->SetTitle("CTime_ePiCoinTime_ROC1");
    h1_CTime_ePiCoinTime_ROC1_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_CTime_ePiCoinTime_ROC1_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_CTime_ePiCoinTime_ROC1_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_Q2_pions_Cut_All->Fill(Q2_pions_cut_all);
    h1_Q2_pions_Cut_All->GetXaxis()->SetTitle("Q2");
    h1_Q2_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_Q2_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_Q2_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_W_pions_Cut_All->Fill(W_pions_cut_all);
    h1_W_pions_Cut_All->GetXaxis()->SetTitle("W");
    h1_W_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_W_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_W_pions_Cut_All->GetYaxis()->CenterTitle();
    h1_epsilon_pions_Cut_All->Fill(epsilon_pions_cut_all);
    h1_epsilon_pions_Cut_All->GetXaxis()->SetTitle("epsilon");
    h1_epsilon_pions_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_epsilon_pions_Cut_All->GetXaxis()->CenterTitle();
    h1_epsilon_pions_Cut_All->GetYaxis()->CenterTitle();
  }

  for(Long64_t i = 0; i < nEntries_Cut_Pion_Events_Prompt; i++){
    Cut_Pion_Events_Prompt->GetEntry(i);

    h1_P_gtr_beta_pions_Cut_Prompt->Fill(P_gtr_beta_pions_cut_pr);
    h1_P_gtr_beta_pions_Cut_Prompt->GetXaxis()->SetTitle("P_gtr_beta");
    h1_P_gtr_beta_pions_Cut_Prompt->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_beta_pions_Cut_Prompt->GetXaxis()->CenterTitle();
    h1_P_gtr_beta_pions_Cut_Prompt->GetYaxis()->CenterTitle();
    h1_P_RFTime_Dist_pions_Cut_Prompt->Fill(P_RFTime_Dist_pions_cut_pr);
    h1_P_RFTime_Dist_pions_Cut_Prompt->GetXaxis()->SetTitle("P_RFTime_Dist");
    h1_P_RFTime_Dist_pions_Cut_Prompt->GetYaxis()->SetTitle("Entries");
    h1_P_RFTime_Dist_pions_Cut_Prompt->GetXaxis()->CenterTitle();
    h1_P_RFTime_Dist_pions_Cut_Prompt->GetYaxis()->CenterTitle();
    h1_CTime_ePiCoinTime_ROC1_pions_Cut_Prompt->Fill(CTime_ePiCoinTime_ROC1_pions_cut_pr);
    h1_CTime_ePiCoinTime_ROC1_pions_Cut_Prompt->GetXaxis()->SetTitle("CTime_ePiCoinTime_ROC1");
    h1_CTime_ePiCoinTime_ROC1_pions_Cut_Prompt->GetYaxis()->SetTitle("Entries");
    h1_CTime_ePiCoinTime_ROC1_pions_Cut_Prompt->GetXaxis()->CenterTitle();
    h1_CTime_ePiCoinTime_ROC1_pions_Cut_Prompt->GetYaxis()->CenterTitle();
    h1_P_kin_MMpi_pions_Cut_Prompt->Fill(P_kin_MMpi_pions_cut_pr);
    h1_P_kin_MMpi_pions_Cut_Prompt->GetXaxis()->SetTitle("Missing_Mass_Pions_prompt");
    h1_P_kin_MMpi_pions_Cut_Prompt->GetYaxis()->SetTitle("Entries");
    h1_P_kin_MMpi_pions_Cut_Prompt->GetXaxis()->CenterTitle();
    h1_P_kin_MMpi_pions_Cut_Prompt->GetYaxis()->CenterTitle();
  }

  for(Long64_t i = 0; i < nEntries_Cut_Pion_Events_Random; i++){
    Cut_Pion_Events_Random->GetEntry(i);

    h1_P_gtr_beta_pions_Cut_Random->Fill(P_gtr_beta_pions_cut_rn);
    h1_P_gtr_beta_pions_Cut_Random->GetXaxis()->SetTitle("P_gtr_beta");
    h1_P_gtr_beta_pions_Cut_Random->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_beta_pions_Cut_Random->GetXaxis()->CenterTitle();
    h1_P_gtr_beta_pions_Cut_Random->GetYaxis()->CenterTitle();
    h1_P_RFTime_Dist_pions_Cut_Random->Fill(P_RFTime_Dist_pions_cut_rn);
    h1_P_RFTime_Dist_pions_Cut_Random->GetXaxis()->SetTitle("P_RFTime_Dist");
    h1_P_RFTime_Dist_pions_Cut_Random->GetYaxis()->SetTitle("Entries");
    h1_P_RFTime_Dist_pions_Cut_Random->GetXaxis()->CenterTitle();
    h1_P_RFTime_Dist_pions_Cut_Random->GetYaxis()->CenterTitle();
    h1_CTime_ePiCoinTime_ROC1_pions_Cut_Random->Fill(CTime_ePiCoinTime_ROC1_pions_cut_rn);
    h1_CTime_ePiCoinTime_ROC1_pions_Cut_Random->GetXaxis()->SetTitle("CTime_ePiCoinTime_ROC1");
    h1_CTime_ePiCoinTime_ROC1_pions_Cut_Random->GetYaxis()->SetTitle("Entries");
    h1_CTime_ePiCoinTime_ROC1_pions_Cut_Random->GetXaxis()->CenterTitle();
    h1_CTime_ePiCoinTime_ROC1_pions_Cut_Random->GetYaxis()->CenterTitle();
    h1_P_kin_MMpi_pions_Cut_Random->Fill(P_kin_MMpi_pions_cut_rn);
    h1_P_kin_MMpi_pions_Cut_Random->GetXaxis()->SetTitle("Missing_Mass_Pions_random");
    h1_P_kin_MMpi_pions_Cut_Random->GetYaxis()->SetTitle("Entries");
    h1_P_kin_MMpi_pions_Cut_Random->GetXaxis()->CenterTitle();
    h1_P_kin_MMpi_pions_Cut_Random->GetYaxis()->CenterTitle();
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Loop over kaon events in tree
  for(Long64_t i = 0; i < nEntries_Uncut_Kaon_Events; i++){
    Uncut_Kaon_Events->GetEntry(i);

    h1_H_gtr_beta_kaons_Uncut->Fill(H_gtr_beta_kaons_uncut);
    h1_H_gtr_beta_kaons_Uncut->GetXaxis()->SetTitle("H_gtr_beta");
    h1_H_gtr_beta_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_beta_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_H_gtr_beta_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_H_gtr_xp_kaons_Uncut->Fill(H_gtr_xp_kaons_uncut);
    h1_H_gtr_xp_kaons_Uncut->GetXaxis()->SetTitle("H_gtr_xp");
    h1_H_gtr_xp_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_xp_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_H_gtr_xp_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_H_gtr_yp_kaons_Uncut->Fill(H_gtr_yp_kaons_uncut);
    h1_H_gtr_yp_kaons_Uncut->GetXaxis()->SetTitle("H_gtr_yp");
    h1_H_gtr_yp_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_yp_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_H_gtr_yp_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_H_gtr_dp_kaons_Uncut->Fill(H_gtr_dp_kaons_uncut);
    h1_H_gtr_dp_kaons_Uncut->GetXaxis()->SetTitle("H_gtr_dp");
    h1_H_gtr_dp_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_dp_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_H_gtr_dp_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_H_hod_goodscinhit_kaons_Uncut->Fill(H_hod_goodscinhit_kaons_uncut);
    h1_H_hod_goodscinhit_kaons_Uncut->GetXaxis()->SetTitle("H_hod_goodscinhit");
    h1_H_hod_goodscinhit_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_hod_goodscinhit_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_H_hod_goodscinhit_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_H_hod_goodstarttime_kaons_Uncut->Fill(H_hod_goodstarttime_kaons_uncut);
    h1_H_hod_goodstarttime_kaons_Uncut->GetXaxis()->SetTitle("H_hod_goodstarttime");
    h1_H_hod_goodstarttime_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_hod_goodstarttime_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_H_hod_goodstarttime_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_H_cal_etotnorm_kaons_Uncut->Fill(H_cal_etotnorm_kaons_uncut);
    h1_H_cal_etotnorm_kaons_Uncut->GetXaxis()->SetTitle("H_cal_etotnorm");
    h1_H_cal_etotnorm_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_cal_etotnorm_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_H_cal_etotnorm_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_H_cal_etottracknorm_kaons_Uncut->Fill(H_cal_etottracknorm_kaons_uncut);
    h1_H_cal_etottracknorm_kaons_Uncut->GetXaxis()->SetTitle("H_cal_etottracknorm");
    h1_H_cal_etottracknorm_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_cal_etottracknorm_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_H_cal_etottracknorm_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_H_cer_npeSum_kaons_Uncut->Fill(H_cer_npeSum_kaons_uncut);
    h1_H_cer_npeSum_kaons_Uncut->GetXaxis()->SetTitle("H_cer_npeSum");
    h1_H_cer_npeSum_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_cer_npeSum_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_H_cer_npeSum_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_H_RFTime_Dist_kaons_Uncut->Fill(H_RFTime_Dist_kaons_uncut);
    h1_H_RFTime_Dist_kaons_Uncut->GetXaxis()->SetTitle("H_RFTime_Dist");
    h1_H_RFTime_Dist_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_RFTime_Dist_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_H_RFTime_Dist_kaons_Uncut->GetYaxis()->CenterTitle();

    h1_P_gtr_beta_kaons_Uncut->Fill(P_gtr_beta_kaons_uncut);
    h1_P_gtr_beta_kaons_Uncut->GetXaxis()->SetTitle("P_gtr_beta");
    h1_P_gtr_beta_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_beta_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_P_gtr_beta_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_P_gtr_xp_kaons_Uncut->Fill(P_gtr_xp_kaons_uncut);
    h1_P_gtr_xp_kaons_Uncut->GetXaxis()->SetTitle("P_gtr_xp");
    h1_P_gtr_xp_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_xp_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_P_gtr_xp_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_P_gtr_yp_kaons_Uncut->Fill(P_gtr_yp_kaons_uncut);
    h1_P_gtr_yp_kaons_Uncut->GetXaxis()->SetTitle("P_gtr_yp");
    h1_P_gtr_yp_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_yp_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_P_gtr_yp_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_P_gtr_dp_kaons_Uncut->Fill(P_gtr_dp_kaons_uncut);
    h1_P_gtr_dp_kaons_Uncut->GetXaxis()->SetTitle("P_gtr_dp");
    h1_P_gtr_dp_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_dp_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_P_gtr_dp_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_P_gtr_p_kaons_Uncut->Fill(P_gtr_p_kaons_uncut);
    h1_P_gtr_p_kaons_Uncut->GetXaxis()->SetTitle("P_gtr_p");
    h1_P_gtr_p_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_p_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_P_gtr_p_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_P_hod_goodscinhit_kaons_Uncut->Fill(P_hod_goodscinhit_kaons_uncut);
    h1_P_hod_goodscinhit_kaons_Uncut->GetXaxis()->SetTitle("P_hod_goodscinhit");
    h1_P_hod_goodscinhit_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_hod_goodscinhit_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_P_hod_goodscinhit_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_P_hod_goodstarttime_kaons_Uncut->Fill(P_hod_goodstarttime_kaons_uncut);
    h1_P_hod_goodstarttime_kaons_Uncut->GetXaxis()->SetTitle("P_hod_goodstarttime");
    h1_P_hod_goodstarttime_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_hod_goodstarttime_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_P_hod_goodstarttime_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_P_cal_etotnorm_kaons_Uncut->Fill(P_cal_etotnorm_kaons_uncut);
    h1_P_cal_etotnorm_kaons_Uncut->GetXaxis()->SetTitle("P_cal_etotnorm");
    h1_P_cal_etotnorm_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_cal_etotnorm_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_P_cal_etotnorm_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_P_cal_etottracknorm_kaons_Uncut->Fill(P_cal_etottracknorm_kaons_uncut);
    h1_P_cal_etottracknorm_kaons_Uncut->GetXaxis()->SetTitle("P_cal_etottracknorm");
    h1_P_cal_etottracknorm_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_cal_etottracknorm_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_P_cal_etottracknorm_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_P_hgcer_npeSum_kaons_Uncut->Fill(P_hgcer_npeSum_kaons_uncut);
    h1_P_hgcer_npeSum_kaons_Uncut->GetXaxis()->SetTitle("P_hgcer_npeSum");
    h1_P_hgcer_npeSum_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_hgcer_npeSum_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_P_hgcer_npeSum_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_P_hgcer_xAtCer_kaons_Uncut->Fill(P_hgcer_xAtCer_kaons_uncut);
    h1_P_hgcer_xAtCer_kaons_Uncut->GetXaxis()->SetTitle("P_hgcer_xAtCer");
    h1_P_hgcer_xAtCer_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_hgcer_xAtCer_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_P_hgcer_xAtCer_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_P_hgcer_yAtCer_kaons_Uncut->Fill(P_hgcer_yAtCer_kaons_uncut);
    h1_P_hgcer_yAtCer_kaons_Uncut->GetXaxis()->SetTitle("P_hgcer_yAtCer");
    h1_P_hgcer_yAtCer_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_hgcer_yAtCer_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_P_hgcer_yAtCer_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_P_aero_npeSum_kaons_Uncut->Fill(P_aero_npeSum_kaons_uncut);
    h1_P_aero_npeSum_kaons_Uncut->GetXaxis()->SetTitle("P_aero_npeSum");
    h1_P_aero_npeSum_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_aero_npeSum_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_P_aero_npeSum_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_P_aero_xAtAero_kaons_Uncut->Fill(P_aero_xAtAero_kaons_uncut);
    h1_P_aero_xAtAero_kaons_Uncut->GetXaxis()->SetTitle("P_aero_xAtAero");
    h1_P_aero_xAtAero_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_aero_xAtAero_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_P_aero_xAtAero_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_P_aero_yAtAero_kaons_Uncut->Fill(P_aero_yAtAero_kaons_uncut);
    h1_P_aero_yAtAero_kaons_Uncut->GetXaxis()->SetTitle("P_aero_yAtAero");
    h1_P_aero_yAtAero_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_aero_yAtAero_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_P_aero_yAtAero_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_P_kin_MMK_kaons_Uncut->Fill(P_kin_MMK_kaons_uncut);
    h1_P_kin_MMK_kaons_Uncut->GetXaxis()->SetTitle("Missing_Mass_Kaons_uncut");
    h1_P_kin_MMK_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_kin_MMK_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_P_kin_MMK_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_P_RFTime_Dist_kaons_Uncut->Fill(P_RFTime_Dist_kaons_uncut);
    h1_P_RFTime_Dist_kaons_Uncut->GetXaxis()->SetTitle("P_RFTime_Dist");
    h1_P_RFTime_Dist_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_RFTime_Dist_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_P_RFTime_Dist_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_CTime_eKCoinTime_ROC1_kaons_Uncut->Fill(CTime_eKCoinTime_ROC1_kaons_uncut);
    h1_CTime_eKCoinTime_ROC1_kaons_Uncut->GetXaxis()->SetTitle("CTime_eKCoinTime_ROC1");
    h1_CTime_eKCoinTime_ROC1_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_CTime_eKCoinTime_ROC1_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_CTime_eKCoinTime_ROC1_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_Q2_kaons_Uncut->Fill(Q2_kaons_uncut);
    h1_Q2_kaons_Uncut->GetXaxis()->SetTitle("Q2");
    h1_Q2_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_Q2_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_Q2_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_W_kaons_Uncut->Fill(W_kaons_uncut);
    h1_W_kaons_Uncut->GetXaxis()->SetTitle("W");
    h1_W_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_W_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_W_kaons_Uncut->GetYaxis()->CenterTitle();
    h1_epsilon_kaons_Uncut->Fill(epsilon_kaons_uncut);
    h1_epsilon_kaons_Uncut->GetXaxis()->SetTitle("epsilon");
    h1_epsilon_kaons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_epsilon_kaons_Uncut->GetXaxis()->CenterTitle();
    h1_epsilon_kaons_Uncut->GetYaxis()->CenterTitle();
  }

  for(Long64_t i = 0; i < nEntries_Cut_Kaon_Events_All; i++){
    Cut_Kaon_Events_All->GetEntry(i);

    h1_H_gtr_beta_kaons_Cut_All->Fill(H_gtr_beta_kaons_cut_all);
    h1_H_gtr_beta_kaons_Cut_All->GetXaxis()->SetTitle("H_gtr_beta");
    h1_H_gtr_beta_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_beta_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_H_gtr_beta_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_H_gtr_xp_kaons_Cut_All->Fill(H_gtr_xp_kaons_cut_all);
    h1_H_gtr_xp_kaons_Cut_All->GetXaxis()->SetTitle("H_gtr_xp");
    h1_H_gtr_xp_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_xp_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_H_gtr_xp_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_H_gtr_yp_kaons_Cut_All->Fill(H_gtr_yp_kaons_cut_all);
    h1_H_gtr_yp_kaons_Cut_All->GetXaxis()->SetTitle("H_gtr_yp");
    h1_H_gtr_yp_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_yp_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_H_gtr_yp_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_H_gtr_dp_kaons_Cut_All->Fill(H_gtr_dp_kaons_cut_all);
    h1_H_gtr_dp_kaons_Cut_All->GetXaxis()->SetTitle("H_gtr_dp");
    h1_H_gtr_dp_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_dp_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_H_gtr_dp_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_H_hod_goodscinhit_kaons_Cut_All->Fill(H_hod_goodscinhit_kaons_cut_all);
    h1_H_hod_goodscinhit_kaons_Cut_All->GetXaxis()->SetTitle("H_hod_goodscinhit");
    h1_H_hod_goodscinhit_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_hod_goodscinhit_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_H_hod_goodscinhit_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_H_hod_goodstarttime_kaons_Cut_All->Fill(H_hod_goodstarttime_kaons_cut_all);
    h1_H_hod_goodstarttime_kaons_Cut_All->GetXaxis()->SetTitle("H_hod_goodstarttime");
    h1_H_hod_goodstarttime_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_hod_goodstarttime_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_H_hod_goodstarttime_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_H_cal_etotnorm_kaons_Cut_All->Fill(H_cal_etotnorm_kaons_cut_all);
    h1_H_cal_etotnorm_kaons_Cut_All->GetXaxis()->SetTitle("H_cal_etotnorm");
    h1_H_cal_etotnorm_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_cal_etotnorm_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_H_cal_etotnorm_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_H_cal_etottracknorm_kaons_Cut_All->Fill(H_cal_etottracknorm_kaons_cut_all);
    h1_H_cal_etottracknorm_kaons_Cut_All->GetXaxis()->SetTitle("H_cal_etottracknorm");
    h1_H_cal_etottracknorm_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_cal_etottracknorm_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_H_cal_etottracknorm_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_H_cer_npeSum_kaons_Cut_All->Fill(H_cer_npeSum_kaons_cut_all);
    h1_H_cer_npeSum_kaons_Cut_All->GetXaxis()->SetTitle("H_cer_npeSum");
    h1_H_cer_npeSum_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_cer_npeSum_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_H_cer_npeSum_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_H_RFTime_Dist_kaons_Cut_All->Fill(H_RFTime_Dist_kaons_cut_all);
    h1_H_RFTime_Dist_kaons_Cut_All->GetXaxis()->SetTitle("H_RFTime_Dist");
    h1_H_RFTime_Dist_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_RFTime_Dist_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_H_RFTime_Dist_kaons_Cut_All->GetYaxis()->CenterTitle();

    h1_P_gtr_beta_kaons_Cut_All->Fill(P_gtr_beta_kaons_cut_all);
    h1_P_gtr_beta_kaons_Cut_All->GetXaxis()->SetTitle("P_gtr_beta");
    h1_P_gtr_beta_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_beta_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_gtr_beta_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_gtr_xp_kaons_Cut_All->Fill(P_gtr_xp_kaons_cut_all);
    h1_P_gtr_xp_kaons_Cut_All->GetXaxis()->SetTitle("P_gtr_xp");
    h1_P_gtr_xp_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_xp_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_gtr_xp_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_gtr_yp_kaons_Cut_All->Fill(P_gtr_yp_kaons_cut_all);
    h1_P_gtr_yp_kaons_Cut_All->GetXaxis()->SetTitle("P_gtr_yp");
    h1_P_gtr_yp_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_yp_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_gtr_yp_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_gtr_dp_kaons_Cut_All->Fill(P_gtr_dp_kaons_cut_all);
    h1_P_gtr_dp_kaons_Cut_All->GetXaxis()->SetTitle("P_gtr_dp");
    h1_P_gtr_dp_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_dp_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_gtr_dp_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_gtr_p_kaons_Cut_All->Fill(P_gtr_p_kaons_cut_all);
    h1_P_gtr_p_kaons_Cut_All->GetXaxis()->SetTitle("P_gtr_p");
    h1_P_gtr_p_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_p_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_gtr_p_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_hod_goodscinhit_kaons_Cut_All->Fill(P_hod_goodscinhit_kaons_cut_all);
    h1_P_hod_goodscinhit_kaons_Cut_All->GetXaxis()->SetTitle("P_hod_goodscinhit");
    h1_P_hod_goodscinhit_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_hod_goodscinhit_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_hod_goodscinhit_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_hod_goodstarttime_kaons_Cut_All->Fill(P_hod_goodstarttime_kaons_cut_all);
    h1_P_hod_goodstarttime_kaons_Cut_All->GetXaxis()->SetTitle("P_hod_goodstarttime");
    h1_P_hod_goodstarttime_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_hod_goodstarttime_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_hod_goodstarttime_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_cal_etotnorm_kaons_Cut_All->Fill(P_cal_etotnorm_kaons_cut_all);
    h1_P_cal_etotnorm_kaons_Cut_All->GetXaxis()->SetTitle("P_cal_etotnorm");
    h1_P_cal_etotnorm_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_cal_etotnorm_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_cal_etotnorm_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_cal_etottracknorm_kaons_Cut_All->Fill(P_cal_etottracknorm_kaons_cut_all);
    h1_P_cal_etottracknorm_kaons_Cut_All->GetXaxis()->SetTitle("P_cal_etottracknorm");
    h1_P_cal_etottracknorm_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_cal_etottracknorm_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_cal_etottracknorm_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_hgcer_npeSum_kaons_Cut_All->Fill(P_hgcer_npeSum_kaons_cut_all);
    h1_P_hgcer_npeSum_kaons_Cut_All->GetXaxis()->SetTitle("P_hgcer_npeSum");
    h1_P_hgcer_npeSum_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_hgcer_npeSum_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_hgcer_npeSum_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_hgcer_xAtCer_kaons_Cut_All->Fill(P_hgcer_xAtCer_kaons_cut_all);
    h1_P_hgcer_xAtCer_kaons_Cut_All->GetXaxis()->SetTitle("P_hgcer_xAtCer");
    h1_P_hgcer_xAtCer_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_hgcer_xAtCer_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_hgcer_xAtCer_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_hgcer_yAtCer_kaons_Cut_All->Fill(P_hgcer_yAtCer_kaons_cut_all);
    h1_P_hgcer_yAtCer_kaons_Cut_All->GetXaxis()->SetTitle("P_hgcer_yAtCer");
    h1_P_hgcer_yAtCer_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_hgcer_yAtCer_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_hgcer_yAtCer_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_aero_npeSum_kaons_Cut_All->Fill(P_aero_npeSum_kaons_cut_all);
    h1_P_aero_npeSum_kaons_Cut_All->GetXaxis()->SetTitle("P_aero_npeSum");
    h1_P_aero_npeSum_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_aero_npeSum_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_aero_npeSum_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_aero_xAtAero_kaons_Cut_All->Fill(P_aero_xAtAero_kaons_cut_all);
    h1_P_aero_xAtAero_kaons_Cut_All->GetXaxis()->SetTitle("P_aero_xAtAero");
    h1_P_aero_xAtAero_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_aero_xAtAero_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_aero_xAtAero_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_aero_yAtAero_kaons_Cut_All->Fill(P_aero_yAtAero_kaons_cut_all);
    h1_P_aero_yAtAero_kaons_Cut_All->GetXaxis()->SetTitle("P_aero_yAtAero");
    h1_P_aero_yAtAero_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_aero_yAtAero_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_aero_yAtAero_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_kin_MMK_kaons_Cut_All->Fill(P_kin_MMK_kaons_cut_all);
    h1_P_kin_MMK_kaons_Cut_All->GetXaxis()->SetTitle("Missing_Mass_Kaons_cut_all");
    h1_P_kin_MMK_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_kin_MMK_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_kin_MMK_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_RFTime_Dist_kaons_Cut_All->Fill(P_RFTime_Dist_kaons_cut_all);
    h1_P_RFTime_Dist_kaons_Cut_All->GetXaxis()->SetTitle("P_RFTime_Dist");
    h1_P_RFTime_Dist_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_RFTime_Dist_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_RFTime_Dist_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_CTime_eKCoinTime_ROC1_kaons_Cut_All->Fill(CTime_eKCoinTime_ROC1_kaons_cut_all);
    h1_CTime_eKCoinTime_ROC1_kaons_Cut_All->GetXaxis()->SetTitle("CTime_eKCoinTime_ROC1");
    h1_CTime_eKCoinTime_ROC1_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_CTime_eKCoinTime_ROC1_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_CTime_eKCoinTime_ROC1_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_Q2_kaons_Cut_All->Fill(Q2_kaons_cut_all);
    h1_Q2_kaons_Cut_All->GetXaxis()->SetTitle("Q2");
    h1_Q2_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_Q2_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_Q2_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_W_kaons_Cut_All->Fill(W_kaons_cut_all);
    h1_W_kaons_Cut_All->GetXaxis()->SetTitle("W");
    h1_W_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_W_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_W_kaons_Cut_All->GetYaxis()->CenterTitle();
    h1_epsilon_kaons_Cut_All->Fill(epsilon_kaons_cut_all);
    h1_epsilon_kaons_Cut_All->GetXaxis()->SetTitle("epsilon");
    h1_epsilon_kaons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_epsilon_kaons_Cut_All->GetXaxis()->CenterTitle();
    h1_epsilon_kaons_Cut_All->GetYaxis()->CenterTitle();
  }

  for(Long64_t i = 0; i < nEntries_Cut_Kaon_Events_Prompt; i++){
    Cut_Kaon_Events_Prompt->GetEntry(i);

    h1_P_gtr_beta_kaons_Cut_Prompt->Fill(P_gtr_beta_kaons_cut_pr);
    h1_P_gtr_beta_kaons_Cut_Prompt->GetXaxis()->SetTitle("P_gtr_beta");
    h1_P_gtr_beta_kaons_Cut_Prompt->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_beta_kaons_Cut_Prompt->GetXaxis()->CenterTitle();
    h1_P_gtr_beta_kaons_Cut_Prompt->GetYaxis()->CenterTitle();
    h1_P_RFTime_Dist_kaons_Cut_Prompt->Fill(P_RFTime_Dist_kaons_cut_pr);
    h1_P_RFTime_Dist_kaons_Cut_Prompt->GetXaxis()->SetTitle("P_RFTime_Dist");
    h1_P_RFTime_Dist_kaons_Cut_Prompt->GetYaxis()->SetTitle("Entries");
    h1_P_RFTime_Dist_kaons_Cut_Prompt->GetXaxis()->CenterTitle();
    h1_P_RFTime_Dist_kaons_Cut_Prompt->GetYaxis()->CenterTitle();
    h1_CTime_eKCoinTime_ROC1_kaons_Cut_Prompt->Fill(CTime_eKCoinTime_ROC1_kaons_cut_pr);
    h1_CTime_eKCoinTime_ROC1_kaons_Cut_Prompt->GetXaxis()->SetTitle("CTime_eKCoinTime_ROC1");
    h1_CTime_eKCoinTime_ROC1_kaons_Cut_Prompt->GetYaxis()->SetTitle("Entries");
    h1_CTime_eKCoinTime_ROC1_kaons_Cut_Prompt->GetXaxis()->CenterTitle();
    h1_CTime_eKCoinTime_ROC1_kaons_Cut_Prompt->GetYaxis()->CenterTitle();
    h1_P_kin_MMK_kaons_Cut_Prompt->Fill(P_kin_MMK_kaons_cut_pr);
    h1_P_kin_MMK_kaons_Cut_Prompt->GetXaxis()->SetTitle("Missing_Mass_Kaons_prompt");
    h1_P_kin_MMK_kaons_Cut_Prompt->GetYaxis()->SetTitle("Entries");
    h1_P_kin_MMK_kaons_Cut_Prompt->GetXaxis()->CenterTitle();
    h1_P_kin_MMK_kaons_Cut_Prompt->GetYaxis()->CenterTitle();
  }

  for(Long64_t i = 0; i < nEntries_Cut_Kaon_Events_Random; i++){
    Cut_Kaon_Events_Random->GetEntry(i);

    h1_P_gtr_beta_kaons_Cut_Random->Fill(P_gtr_beta_kaons_cut_rn);
    h1_P_gtr_beta_kaons_Cut_Random->GetXaxis()->SetTitle("P_gtr_beta");
    h1_P_gtr_beta_kaons_Cut_Random->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_beta_kaons_Cut_Random->GetXaxis()->CenterTitle();
    h1_P_gtr_beta_kaons_Cut_Random->GetYaxis()->CenterTitle();
    h1_P_RFTime_Dist_kaons_Cut_Random->Fill(P_RFTime_Dist_kaons_cut_rn);
    h1_P_RFTime_Dist_kaons_Cut_Random->GetXaxis()->SetTitle("P_RFTime_Dist");
    h1_P_RFTime_Dist_kaons_Cut_Random->GetYaxis()->SetTitle("Entries");
    h1_P_RFTime_Dist_kaons_Cut_Random->GetXaxis()->CenterTitle();
    h1_P_RFTime_Dist_kaons_Cut_Random->GetYaxis()->CenterTitle();
    h1_CTime_eKCoinTime_ROC1_kaons_Cut_Random->Fill(CTime_eKCoinTime_ROC1_kaons_cut_rn);
    h1_CTime_eKCoinTime_ROC1_kaons_Cut_Random->GetXaxis()->SetTitle("CTime_eKCoinTime_ROC1");
    h1_CTime_eKCoinTime_ROC1_kaons_Cut_Random->GetYaxis()->SetTitle("Entries");
    h1_CTime_eKCoinTime_ROC1_kaons_Cut_Random->GetXaxis()->CenterTitle();
    h1_CTime_eKCoinTime_ROC1_kaons_Cut_Random->GetYaxis()->CenterTitle();
    h1_P_kin_MMK_kaons_Cut_Random->Fill(P_kin_MMK_kaons_cut_rn);
    h1_P_kin_MMK_kaons_Cut_Random->GetXaxis()->SetTitle("Missing_Mass_Kaons_random");
    h1_P_kin_MMK_kaons_Cut_Random->GetYaxis()->SetTitle("Entries");
    h1_P_kin_MMK_kaons_Cut_Random->GetXaxis()->CenterTitle();
    h1_P_kin_MMK_kaons_Cut_Random->GetYaxis()->CenterTitle();
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Loop over proton events in tree
  for(Long64_t i = 0; i < nEntries_Uncut_Proton_Events; i++){
    Uncut_Proton_Events->GetEntry(i);

    h1_H_gtr_beta_protons_Uncut->Fill(H_gtr_beta_protons_uncut);
    h1_H_gtr_beta_protons_Uncut->GetXaxis()->SetTitle("H_gtr_beta");
    h1_H_gtr_beta_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_beta_protons_Uncut->GetXaxis()->CenterTitle();
    h1_H_gtr_beta_protons_Uncut->GetYaxis()->CenterTitle();
    h1_H_gtr_xp_protons_Uncut->Fill(H_gtr_xp_protons_uncut);
    h1_H_gtr_xp_protons_Uncut->GetXaxis()->SetTitle("H_gtr_xp");
    h1_H_gtr_xp_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_xp_protons_Uncut->GetXaxis()->CenterTitle();
    h1_H_gtr_xp_protons_Uncut->GetYaxis()->CenterTitle();
    h1_H_gtr_yp_protons_Uncut->Fill(H_gtr_yp_protons_uncut);
    h1_H_gtr_yp_protons_Uncut->GetXaxis()->SetTitle("H_gtr_yp");
    h1_H_gtr_yp_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_yp_protons_Uncut->GetXaxis()->CenterTitle();
    h1_H_gtr_yp_protons_Uncut->GetYaxis()->CenterTitle();
    h1_H_gtr_dp_protons_Uncut->Fill(H_gtr_dp_protons_uncut);
    h1_H_gtr_dp_protons_Uncut->GetXaxis()->SetTitle("H_gtr_dp");
    h1_H_gtr_dp_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_dp_protons_Uncut->GetXaxis()->CenterTitle();
    h1_H_gtr_dp_protons_Uncut->GetYaxis()->CenterTitle();
    h1_H_hod_goodscinhit_protons_Uncut->Fill(H_hod_goodscinhit_protons_uncut);
    h1_H_hod_goodscinhit_protons_Uncut->GetXaxis()->SetTitle("H_hod_goodscinhit");
    h1_H_hod_goodscinhit_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_hod_goodscinhit_protons_Uncut->GetXaxis()->CenterTitle();
    h1_H_hod_goodscinhit_protons_Uncut->GetYaxis()->CenterTitle();
    h1_H_hod_goodstarttime_protons_Uncut->Fill(H_hod_goodstarttime_protons_uncut);
    h1_H_hod_goodstarttime_protons_Uncut->GetXaxis()->SetTitle("H_hod_goodstarttime");
    h1_H_hod_goodstarttime_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_hod_goodstarttime_protons_Uncut->GetXaxis()->CenterTitle();
    h1_H_hod_goodstarttime_protons_Uncut->GetYaxis()->CenterTitle();
    h1_H_cal_etotnorm_protons_Uncut->Fill(H_cal_etotnorm_protons_uncut);
    h1_H_cal_etotnorm_protons_Uncut->GetXaxis()->SetTitle("H_cal_etotnorm");
    h1_H_cal_etotnorm_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_cal_etotnorm_protons_Uncut->GetXaxis()->CenterTitle();
    h1_H_cal_etotnorm_protons_Uncut->GetYaxis()->CenterTitle();
    h1_H_cal_etottracknorm_protons_Uncut->Fill(H_cal_etottracknorm_protons_uncut);
    h1_H_cal_etottracknorm_protons_Uncut->GetXaxis()->SetTitle("H_cal_etottracknorm");
    h1_H_cal_etottracknorm_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_cal_etottracknorm_protons_Uncut->GetXaxis()->CenterTitle();
    h1_H_cal_etottracknorm_protons_Uncut->GetYaxis()->CenterTitle();
    h1_H_cer_npeSum_protons_Uncut->Fill(H_cer_npeSum_protons_uncut);
    h1_H_cer_npeSum_protons_Uncut->GetXaxis()->SetTitle("H_cer_npeSum");
    h1_H_cer_npeSum_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_cer_npeSum_protons_Uncut->GetXaxis()->CenterTitle();
    h1_H_cer_npeSum_protons_Uncut->GetYaxis()->CenterTitle();
    h1_H_RFTime_Dist_protons_Uncut->Fill(H_RFTime_Dist_protons_uncut);
    h1_H_RFTime_Dist_protons_Uncut->GetXaxis()->SetTitle("H_RFTime_Dist");
    h1_H_RFTime_Dist_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_H_RFTime_Dist_protons_Uncut->GetXaxis()->CenterTitle();
    h1_H_RFTime_Dist_protons_Uncut->GetYaxis()->CenterTitle();

    h1_P_gtr_beta_protons_Uncut->Fill(P_gtr_beta_protons_uncut);
    h1_P_gtr_beta_protons_Uncut->GetXaxis()->SetTitle("P_gtr_beta");
    h1_P_gtr_beta_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_beta_protons_Uncut->GetXaxis()->CenterTitle();
    h1_P_gtr_beta_protons_Uncut->GetYaxis()->CenterTitle();
    h1_P_gtr_xp_protons_Uncut->Fill(P_gtr_xp_protons_uncut);
    h1_P_gtr_xp_protons_Uncut->GetXaxis()->SetTitle("P_gtr_xp");
    h1_P_gtr_xp_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_xp_protons_Uncut->GetXaxis()->CenterTitle();
    h1_P_gtr_xp_protons_Uncut->GetYaxis()->CenterTitle();
    h1_P_gtr_yp_protons_Uncut->Fill(P_gtr_yp_protons_uncut);
    h1_P_gtr_yp_protons_Uncut->GetXaxis()->SetTitle("P_gtr_yp");
    h1_P_gtr_yp_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_yp_protons_Uncut->GetXaxis()->CenterTitle();
    h1_P_gtr_yp_protons_Uncut->GetYaxis()->CenterTitle();
    h1_P_gtr_dp_protons_Uncut->Fill(P_gtr_dp_protons_uncut);
    h1_P_gtr_dp_protons_Uncut->GetXaxis()->SetTitle("P_gtr_dp");
    h1_P_gtr_dp_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_dp_protons_Uncut->GetXaxis()->CenterTitle();
    h1_P_gtr_dp_protons_Uncut->GetYaxis()->CenterTitle();
    h1_P_gtr_p_protons_Uncut->Fill(P_gtr_p_protons_uncut);
    h1_P_gtr_p_protons_Uncut->GetXaxis()->SetTitle("P_gtr_p");
    h1_P_gtr_p_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_p_protons_Uncut->GetXaxis()->CenterTitle();
    h1_P_gtr_p_protons_Uncut->GetYaxis()->CenterTitle();
    h1_P_hod_goodscinhit_protons_Uncut->Fill(P_hod_goodscinhit_protons_uncut);
    h1_P_hod_goodscinhit_protons_Uncut->GetXaxis()->SetTitle("P_hod_goodscinhit");
    h1_P_hod_goodscinhit_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_hod_goodscinhit_protons_Uncut->GetXaxis()->CenterTitle();
    h1_P_hod_goodscinhit_protons_Uncut->GetYaxis()->CenterTitle();
    h1_P_hod_goodstarttime_protons_Uncut->Fill(P_hod_goodstarttime_protons_uncut);
    h1_P_hod_goodstarttime_protons_Uncut->GetXaxis()->SetTitle("P_hod_goodstarttime");
    h1_P_hod_goodstarttime_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_hod_goodstarttime_protons_Uncut->GetXaxis()->CenterTitle();
    h1_P_hod_goodstarttime_protons_Uncut->GetYaxis()->CenterTitle();
    h1_P_cal_etotnorm_protons_Uncut->Fill(P_cal_etotnorm_protons_uncut);
    h1_P_cal_etotnorm_protons_Uncut->GetXaxis()->SetTitle("P_cal_etotnorm");
    h1_P_cal_etotnorm_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_cal_etotnorm_protons_Uncut->GetXaxis()->CenterTitle();
    h1_P_cal_etotnorm_protons_Uncut->GetYaxis()->CenterTitle();
    h1_P_cal_etottracknorm_protons_Uncut->Fill(P_cal_etottracknorm_protons_uncut);
    h1_P_cal_etottracknorm_protons_Uncut->GetXaxis()->SetTitle("P_cal_etottracknorm");
    h1_P_cal_etottracknorm_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_cal_etottracknorm_protons_Uncut->GetXaxis()->CenterTitle();
    h1_P_cal_etottracknorm_protons_Uncut->GetYaxis()->CenterTitle();
    h1_P_hgcer_npeSum_protons_Uncut->Fill(P_hgcer_npeSum_protons_uncut);
    h1_P_hgcer_npeSum_protons_Uncut->GetXaxis()->SetTitle("P_hgcer_npeSum");
    h1_P_hgcer_npeSum_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_hgcer_npeSum_protons_Uncut->GetXaxis()->CenterTitle();
    h1_P_hgcer_npeSum_protons_Uncut->GetYaxis()->CenterTitle();
    h1_P_hgcer_xAtCer_protons_Uncut->Fill(P_hgcer_xAtCer_protons_uncut);
    h1_P_hgcer_xAtCer_protons_Uncut->GetXaxis()->SetTitle("P_hgcer_xAtCer");
    h1_P_hgcer_xAtCer_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_hgcer_xAtCer_protons_Uncut->GetXaxis()->CenterTitle();
    h1_P_hgcer_xAtCer_protons_Uncut->GetYaxis()->CenterTitle();
    h1_P_hgcer_yAtCer_protons_Uncut->Fill(P_hgcer_yAtCer_protons_uncut);
    h1_P_hgcer_yAtCer_protons_Uncut->GetXaxis()->SetTitle("P_hgcer_yAtCer");
    h1_P_hgcer_yAtCer_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_hgcer_yAtCer_protons_Uncut->GetXaxis()->CenterTitle();
    h1_P_hgcer_yAtCer_protons_Uncut->GetYaxis()->CenterTitle();
    h1_P_aero_npeSum_protons_Uncut->Fill(P_aero_npeSum_protons_uncut);
    h1_P_aero_npeSum_protons_Uncut->GetXaxis()->SetTitle("P_aero_npeSum");
    h1_P_aero_npeSum_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_aero_npeSum_protons_Uncut->GetXaxis()->CenterTitle();
    h1_P_aero_npeSum_protons_Uncut->GetYaxis()->CenterTitle();
    h1_P_aero_xAtAero_protons_Uncut->Fill(P_aero_xAtAero_protons_uncut);
    h1_P_aero_xAtAero_protons_Uncut->GetXaxis()->SetTitle("P_aero_xAtAero");
    h1_P_aero_xAtAero_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_aero_xAtAero_protons_Uncut->GetXaxis()->CenterTitle();
    h1_P_aero_xAtAero_protons_Uncut->GetYaxis()->CenterTitle();
    h1_P_aero_yAtAero_protons_Uncut->Fill(P_aero_yAtAero_protons_uncut);
    h1_P_aero_yAtAero_protons_Uncut->GetXaxis()->SetTitle("P_aero_yAtAero");
    h1_P_aero_yAtAero_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_aero_yAtAero_protons_Uncut->GetXaxis()->CenterTitle();
    h1_P_aero_yAtAero_protons_Uncut->GetYaxis()->CenterTitle();
    h1_P_kin_MMp_protons_Uncut->Fill(P_kin_MMp_protons_uncut);
    h1_P_kin_MMp_protons_Uncut->GetXaxis()->SetTitle("Missing_Mass_Protons_uncut");
    h1_P_kin_MMp_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_kin_MMp_protons_Uncut->GetXaxis()->CenterTitle();
    h1_P_kin_MMp_protons_Uncut->GetYaxis()->CenterTitle();
    h1_P_RFTime_Dist_protons_Uncut->Fill(P_RFTime_Dist_protons_uncut);
    h1_P_RFTime_Dist_protons_Uncut->GetXaxis()->SetTitle("P_RFTime_Dist");
    h1_P_RFTime_Dist_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_P_RFTime_Dist_protons_Uncut->GetXaxis()->CenterTitle();
    h1_P_RFTime_Dist_protons_Uncut->GetYaxis()->CenterTitle();
    h1_CTime_epCoinTime_ROC1_protons_Uncut->Fill(CTime_epCoinTime_ROC1_protons_uncut);
    h1_CTime_epCoinTime_ROC1_protons_Uncut->GetXaxis()->SetTitle("CTime_epCoinTime_ROC1");
    h1_CTime_epCoinTime_ROC1_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_CTime_epCoinTime_ROC1_protons_Uncut->GetXaxis()->CenterTitle();
    h1_CTime_epCoinTime_ROC1_protons_Uncut->GetYaxis()->CenterTitle();
    h1_Q2_protons_Uncut->Fill(Q2_protons_uncut);
    h1_Q2_protons_Uncut->GetXaxis()->SetTitle("Q2");
    h1_Q2_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_Q2_protons_Uncut->GetXaxis()->CenterTitle();
    h1_Q2_protons_Uncut->GetYaxis()->CenterTitle();
    h1_W_protons_Uncut->Fill(W_protons_uncut);
    h1_W_protons_Uncut->GetXaxis()->SetTitle("W");
    h1_W_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_W_protons_Uncut->GetXaxis()->CenterTitle();
    h1_W_protons_Uncut->GetYaxis()->CenterTitle();
    h1_epsilon_protons_Uncut->Fill(epsilon_protons_uncut);
    h1_epsilon_protons_Uncut->GetXaxis()->SetTitle("epsilon");
    h1_epsilon_protons_Uncut->GetYaxis()->SetTitle("Entries");
    h1_epsilon_protons_Uncut->GetXaxis()->CenterTitle();
    h1_epsilon_protons_Uncut->GetYaxis()->CenterTitle();
  }

  for(Long64_t i = 0; i < nEntries_Cut_Proton_Events_All; i++){
    Cut_Proton_Events_All->GetEntry(i);

    h1_H_gtr_beta_protons_Cut_All->Fill(H_gtr_beta_protons_cut_all);
    h1_H_gtr_beta_protons_Cut_All->GetXaxis()->SetTitle("H_gtr_beta");
    h1_H_gtr_beta_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_beta_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_H_gtr_beta_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_H_gtr_xp_protons_Cut_All->Fill(H_gtr_xp_protons_cut_all);
    h1_H_gtr_xp_protons_Cut_All->GetXaxis()->SetTitle("H_gtr_xp");
    h1_H_gtr_xp_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_xp_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_H_gtr_xp_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_H_gtr_yp_protons_Cut_All->Fill(H_gtr_yp_protons_cut_all);
    h1_H_gtr_yp_protons_Cut_All->GetXaxis()->SetTitle("H_gtr_yp");
    h1_H_gtr_yp_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_yp_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_H_gtr_yp_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_H_gtr_dp_protons_Cut_All->Fill(H_gtr_dp_protons_cut_all);
    h1_H_gtr_dp_protons_Cut_All->GetXaxis()->SetTitle("H_gtr_dp");
    h1_H_gtr_dp_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_gtr_dp_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_H_gtr_dp_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_H_hod_goodscinhit_protons_Cut_All->Fill(H_hod_goodscinhit_protons_cut_all);
    h1_H_hod_goodscinhit_protons_Cut_All->GetXaxis()->SetTitle("H_hod_goodscinhit");
    h1_H_hod_goodscinhit_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_hod_goodscinhit_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_H_hod_goodscinhit_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_H_hod_goodstarttime_protons_Cut_All->Fill(H_hod_goodstarttime_protons_cut_all);
    h1_H_hod_goodstarttime_protons_Cut_All->GetXaxis()->SetTitle("H_hod_goodstarttime");
    h1_H_hod_goodstarttime_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_hod_goodstarttime_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_H_hod_goodstarttime_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_H_cal_etotnorm_protons_Cut_All->Fill(H_cal_etotnorm_protons_cut_all);
    h1_H_cal_etotnorm_protons_Cut_All->GetXaxis()->SetTitle("H_cal_etotnorm");
    h1_H_cal_etotnorm_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_cal_etotnorm_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_H_cal_etotnorm_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_H_cal_etottracknorm_protons_Cut_All->Fill(H_cal_etottracknorm_protons_cut_all);
    h1_H_cal_etottracknorm_protons_Cut_All->GetXaxis()->SetTitle("H_cal_etottracknorm");
    h1_H_cal_etottracknorm_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_cal_etottracknorm_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_H_cal_etottracknorm_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_H_cer_npeSum_protons_Cut_All->Fill(H_cer_npeSum_protons_cut_all);
    h1_H_cer_npeSum_protons_Cut_All->GetXaxis()->SetTitle("H_cer_npeSum");
    h1_H_cer_npeSum_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_cer_npeSum_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_H_cer_npeSum_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_H_RFTime_Dist_protons_Cut_All->Fill(H_RFTime_Dist_protons_cut_all);
    h1_H_RFTime_Dist_protons_Cut_All->GetXaxis()->SetTitle("H_RFTime_Dist");
    h1_H_RFTime_Dist_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_H_RFTime_Dist_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_H_RFTime_Dist_protons_Cut_All->GetYaxis()->CenterTitle();

    h1_P_gtr_beta_protons_Cut_All->Fill(P_gtr_beta_protons_cut_all);
    h1_P_gtr_beta_protons_Cut_All->GetXaxis()->SetTitle("P_gtr_beta");
    h1_P_gtr_beta_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_beta_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_gtr_beta_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_gtr_xp_protons_Cut_All->Fill(P_gtr_xp_protons_cut_all);
    h1_P_gtr_xp_protons_Cut_All->GetXaxis()->SetTitle("P_gtr_xp");
    h1_P_gtr_xp_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_xp_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_gtr_xp_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_gtr_yp_protons_Cut_All->Fill(P_gtr_yp_protons_cut_all);
    h1_P_gtr_yp_protons_Cut_All->GetXaxis()->SetTitle("P_gtr_yp");
    h1_P_gtr_yp_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_yp_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_gtr_yp_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_gtr_dp_protons_Cut_All->Fill(P_gtr_dp_protons_cut_all);
    h1_P_gtr_dp_protons_Cut_All->GetXaxis()->SetTitle("P_gtr_dp");
    h1_P_gtr_dp_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_dp_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_gtr_dp_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_gtr_p_protons_Cut_All->Fill(P_gtr_p_protons_cut_all);
    h1_P_gtr_p_protons_Cut_All->GetXaxis()->SetTitle("P_gtr_p");
    h1_P_gtr_p_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_p_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_gtr_p_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_hod_goodscinhit_protons_Cut_All->Fill(P_hod_goodscinhit_protons_cut_all);
    h1_P_hod_goodscinhit_protons_Cut_All->GetXaxis()->SetTitle("P_hod_goodscinhit");
    h1_P_hod_goodscinhit_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_hod_goodscinhit_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_hod_goodscinhit_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_hod_goodstarttime_protons_Cut_All->Fill(P_hod_goodstarttime_protons_cut_all);
    h1_P_hod_goodstarttime_protons_Cut_All->GetXaxis()->SetTitle("P_hod_goodstarttime");
    h1_P_hod_goodstarttime_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_hod_goodstarttime_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_hod_goodstarttime_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_cal_etotnorm_protons_Cut_All->Fill(P_cal_etotnorm_protons_cut_all);
    h1_P_cal_etotnorm_protons_Cut_All->GetXaxis()->SetTitle("P_cal_etotnorm");
    h1_P_cal_etotnorm_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_cal_etotnorm_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_cal_etotnorm_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_cal_etottracknorm_protons_Cut_All->Fill(P_cal_etottracknorm_protons_cut_all);
    h1_P_cal_etottracknorm_protons_Cut_All->GetXaxis()->SetTitle("P_cal_etottracknorm");
    h1_P_cal_etottracknorm_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_cal_etottracknorm_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_cal_etottracknorm_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_hgcer_npeSum_protons_Cut_All->Fill(P_hgcer_npeSum_protons_cut_all);
    h1_P_hgcer_npeSum_protons_Cut_All->GetXaxis()->SetTitle("P_hgcer_npeSum");
    h1_P_hgcer_npeSum_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_hgcer_npeSum_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_hgcer_npeSum_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_hgcer_xAtCer_protons_Cut_All->Fill(P_hgcer_xAtCer_protons_cut_all);
    h1_P_hgcer_xAtCer_protons_Cut_All->GetXaxis()->SetTitle("P_hgcer_xAtCer");
    h1_P_hgcer_xAtCer_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_hgcer_xAtCer_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_hgcer_xAtCer_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_hgcer_yAtCer_protons_Cut_All->Fill(P_hgcer_yAtCer_protons_cut_all);
    h1_P_hgcer_yAtCer_protons_Cut_All->GetXaxis()->SetTitle("P_hgcer_yAtCer");
    h1_P_hgcer_yAtCer_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_hgcer_yAtCer_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_hgcer_yAtCer_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_aero_npeSum_protons_Cut_All->Fill(P_aero_npeSum_protons_cut_all);
    h1_P_aero_npeSum_protons_Cut_All->GetXaxis()->SetTitle("P_aero_npeSum");
    h1_P_aero_npeSum_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_aero_npeSum_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_aero_npeSum_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_aero_xAtAero_protons_Cut_All->Fill(P_aero_xAtAero_protons_cut_all);
    h1_P_aero_xAtAero_protons_Cut_All->GetXaxis()->SetTitle("P_aero_xAtAero");
    h1_P_aero_xAtAero_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_aero_xAtAero_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_aero_xAtAero_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_aero_yAtAero_protons_Cut_All->Fill(P_aero_yAtAero_protons_cut_all);
    h1_P_aero_yAtAero_protons_Cut_All->GetXaxis()->SetTitle("P_aero_yAtAero");
    h1_P_aero_yAtAero_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_aero_yAtAero_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_aero_yAtAero_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_kin_MMp_protons_Cut_All->Fill(P_kin_MMp_protons_cut_all);
    h1_P_kin_MMp_protons_Cut_All->GetXaxis()->SetTitle("Missing_Mass_Protons_cut_all");
    h1_P_kin_MMp_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_kin_MMp_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_kin_MMp_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_P_RFTime_Dist_protons_Cut_All->Fill(P_RFTime_Dist_protons_cut_all);
    h1_P_RFTime_Dist_protons_Cut_All->GetXaxis()->SetTitle("P_RFTime_Dist");
    h1_P_RFTime_Dist_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_P_RFTime_Dist_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_P_RFTime_Dist_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_CTime_epCoinTime_ROC1_protons_Cut_All->Fill(CTime_epCoinTime_ROC1_protons_cut_all);
    h1_CTime_epCoinTime_ROC1_protons_Cut_All->GetXaxis()->SetTitle("CTime_epCoinTime_ROC1");
    h1_CTime_epCoinTime_ROC1_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_CTime_epCoinTime_ROC1_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_CTime_epCoinTime_ROC1_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_Q2_protons_Cut_All->Fill(Q2_protons_cut_all);
    h1_Q2_protons_Cut_All->GetXaxis()->SetTitle("Q2");
    h1_Q2_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_Q2_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_Q2_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_W_protons_Cut_All->Fill(W_protons_cut_all);
    h1_W_protons_Cut_All->GetXaxis()->SetTitle("W");
    h1_W_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_W_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_W_protons_Cut_All->GetYaxis()->CenterTitle();
    h1_epsilon_protons_Cut_All->Fill(epsilon_protons_cut_all);
    h1_epsilon_protons_Cut_All->GetXaxis()->SetTitle("epsilon");
    h1_epsilon_protons_Cut_All->GetYaxis()->SetTitle("Entries");
    h1_epsilon_protons_Cut_All->GetXaxis()->CenterTitle();
    h1_epsilon_protons_Cut_All->GetYaxis()->CenterTitle();
  }

  for(Long64_t i = 0; i < nEntries_Cut_Proton_Events_Prompt; i++){
    Cut_Proton_Events_Prompt->GetEntry(i);

    h1_P_gtr_beta_protons_Cut_Prompt->Fill(P_gtr_beta_protons_cut_pr);
    h1_P_gtr_beta_protons_Cut_Prompt->GetXaxis()->SetTitle("P_gtr_beta");
    h1_P_gtr_beta_protons_Cut_Prompt->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_beta_protons_Cut_Prompt->GetXaxis()->CenterTitle();
    h1_P_gtr_beta_protons_Cut_Prompt->GetYaxis()->CenterTitle();
    h1_P_RFTime_Dist_protons_Cut_Prompt->Fill(P_RFTime_Dist_protons_cut_pr);
    h1_P_RFTime_Dist_protons_Cut_Prompt->GetXaxis()->SetTitle("P_RFTime_Dist");
    h1_P_RFTime_Dist_protons_Cut_Prompt->GetYaxis()->SetTitle("Entries");
    h1_P_RFTime_Dist_protons_Cut_Prompt->GetXaxis()->CenterTitle();
    h1_P_RFTime_Dist_protons_Cut_Prompt->GetYaxis()->CenterTitle();
    h1_CTime_epCoinTime_ROC1_protons_Cut_Prompt->Fill(CTime_epCoinTime_ROC1_protons_cut_pr);
    h1_CTime_epCoinTime_ROC1_protons_Cut_Prompt->GetXaxis()->SetTitle("CTime_epCoinTime_ROC1");
    h1_CTime_epCoinTime_ROC1_protons_Cut_Prompt->GetYaxis()->SetTitle("Entries");
    h1_CTime_epCoinTime_ROC1_protons_Cut_Prompt->GetXaxis()->CenterTitle();
    h1_CTime_epCoinTime_ROC1_protons_Cut_Prompt->GetYaxis()->CenterTitle();
    h1_P_kin_MMp_protons_Cut_Prompt->Fill(P_kin_MMp_protons_cut_pr);
    h1_P_kin_MMp_protons_Cut_Prompt->GetXaxis()->SetTitle("Missing_Mass_Protons_prompt");
    h1_P_kin_MMp_protons_Cut_Prompt->GetYaxis()->SetTitle("Entries");
    h1_P_kin_MMp_protons_Cut_Prompt->GetXaxis()->CenterTitle();
    h1_P_kin_MMp_protons_Cut_Prompt->GetYaxis()->CenterTitle();
  }

  for(Long64_t i = 0; i < nEntries_Cut_Proton_Events_Random; i++){
    Cut_Proton_Events_Random->GetEntry(i);

    h1_P_gtr_beta_protons_Cut_Random->Fill(P_gtr_beta_protons_cut_rn);
    h1_P_gtr_beta_protons_Cut_Random->GetXaxis()->SetTitle("P_gtr_beta");
    h1_P_gtr_beta_protons_Cut_Random->GetYaxis()->SetTitle("Entries");
    h1_P_gtr_beta_protons_Cut_Random->GetXaxis()->CenterTitle();
    h1_P_gtr_beta_protons_Cut_Random->GetYaxis()->CenterTitle();
    h1_P_RFTime_Dist_protons_Cut_Random->Fill(P_RFTime_Dist_protons_cut_rn);
    h1_P_RFTime_Dist_protons_Cut_Random->GetXaxis()->SetTitle("P_RFTime_Dist");
    h1_P_RFTime_Dist_protons_Cut_Random->GetYaxis()->SetTitle("Entries");
    h1_P_RFTime_Dist_protons_Cut_Random->GetXaxis()->CenterTitle();
    h1_P_RFTime_Dist_protons_Cut_Random->GetYaxis()->CenterTitle();
    h1_CTime_epCoinTime_ROC1_protons_Cut_Random->Fill(CTime_epCoinTime_ROC1_protons_cut_rn);
    h1_CTime_epCoinTime_ROC1_protons_Cut_Random->GetXaxis()->SetTitle("CTime_epCoinTime_ROC1");
    h1_CTime_epCoinTime_ROC1_protons_Cut_Random->GetYaxis()->SetTitle("Entries");
    h1_CTime_epCoinTime_ROC1_protons_Cut_Random->GetXaxis()->CenterTitle();
    h1_CTime_epCoinTime_ROC1_protons_Cut_Random->GetYaxis()->CenterTitle();
    h1_P_kin_MMp_protons_Cut_Random->Fill(P_kin_MMp_protons_cut_rn);
    h1_P_kin_MMp_protons_Cut_Random->GetXaxis()->SetTitle("Missing_Mass_Protons_random");
    h1_P_kin_MMp_protons_Cut_Random->GetYaxis()->SetTitle("Entries");
    h1_P_kin_MMp_protons_Cut_Random->GetXaxis()->CenterTitle();
    h1_P_kin_MMp_protons_Cut_Random->GetYaxis()->CenterTitle();
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for(Long64_t i = 0; i < nEntries_Uncut_Pion_Events; i++){
    Uncut_Pion_Events->GetEntry(i);

  h2_P_hgcer_npeSum_vs_aero_npeSum_pions_Uncut->Fill(P_hgcer_npeSum_pions_uncut, P_aero_npeSum_pions_uncut);
  h2_P_hgcer_npeSum_vs_aero_npeSum_pions_Uncut->GetXaxis()->SetTitle("P_hgcer_npeSum");
  h2_P_hgcer_npeSum_vs_aero_npeSum_pions_Uncut->GetYaxis()->SetTitle("P_aero_npeSum");
  h2_P_hgcer_npeSum_vs_aero_npeSum_pions_Uncut->GetXaxis()->CenterTitle();
  h2_P_hgcer_npeSum_vs_aero_npeSum_pions_Uncut->GetYaxis()->CenterTitle();
  
  h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Uncut->Fill(CTime_ePiCoinTime_ROC1_pions_uncut, P_kin_MMpi_pions_uncut);
  h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Uncut->GetXaxis()->SetTitle("CTime_ePiCoinTime_ROC1");
  h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Uncut->GetYaxis()->SetTitle("P_kin_MMpi");
  h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Uncut->GetXaxis()->CenterTitle();
  h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Uncut->GetYaxis()->CenterTitle();

  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_pions_Uncut->Fill(P_hgcer_xAtCer_pions_uncut, P_hgcer_yAtCer_pions_uncut);
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_pions_Uncut->GetXaxis()->SetTitle("P_hgcer_xAtCer");
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_pions_Uncut->GetYaxis()->SetTitle("P_hgcer_yAtCer");
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_pions_Uncut->GetXaxis()->CenterTitle();
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_pions_Uncut->GetYaxis()->CenterTitle();
  }

  for(Long64_t i = 0; i < nEntries_Cut_Pion_Events_All; i++){
    Cut_Pion_Events_All->GetEntry(i);

  h2_P_hgcer_npeSum_vs_aero_npeSum_pions_Cut_All->Fill(P_hgcer_npeSum_pions_cut_all, P_aero_npeSum_pions_cut_all);
  h2_P_hgcer_npeSum_vs_aero_npeSum_pions_Cut_All->GetXaxis()->SetTitle("P_hgcer_npeSum");
  h2_P_hgcer_npeSum_vs_aero_npeSum_pions_Cut_All->GetYaxis()->SetTitle("P_aero_npeSum");
  h2_P_hgcer_npeSum_vs_aero_npeSum_pions_Cut_All->GetXaxis()->CenterTitle();
  h2_P_hgcer_npeSum_vs_aero_npeSum_pions_Cut_All->GetYaxis()->CenterTitle();

  h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Cut_All->Fill(CTime_ePiCoinTime_ROC1_pions_cut_all, P_kin_MMpi_pions_cut_all);
  h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Cut_All->GetXaxis()->SetTitle("CTime_ePiCoinTime_ROC1");
  h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Cut_All->GetYaxis()->SetTitle("P_kin_MMpi");
  h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Cut_All->GetXaxis()->CenterTitle();
  h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Cut_All->GetYaxis()->CenterTitle();

  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_pions_Cut_All->Fill(P_hgcer_xAtCer_pions_cut_all, P_hgcer_yAtCer_pions_cut_all);
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_pions_Cut_All->GetXaxis()->SetTitle("P_hgcer_xAtCer");
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_pions_Cut_All->GetYaxis()->SetTitle("P_hgcer_yAtCer");
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_pions_Cut_All->GetXaxis()->CenterTitle();
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_pions_Cut_All->GetYaxis()->CenterTitle();
  }

  for(Long64_t i = 0; i < nEntries_Cut_Pion_Events_Prompt; i++){
    Cut_Pion_Events_Prompt->GetEntry(i);

  h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Cut_Prompt->Fill(CTime_ePiCoinTime_ROC1_pions_cut_pr, P_kin_MMpi_pions_cut_pr);
  h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Cut_Prompt->GetXaxis()->SetTitle("CTime_ePiCoinTime_ROC1");
  h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Cut_Prompt->GetYaxis()->SetTitle("P_kin_MMpi");
  h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Cut_Prompt->GetXaxis()->CenterTitle();
  h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Cut_Prompt->GetYaxis()->CenterTitle();
  }

  for(Long64_t i = 0; i < nEntries_Uncut_Kaon_Events; i++){
    Uncut_Kaon_Events->GetEntry(i);

  h2_P_hgcer_npeSum_vs_aero_npeSum_kaons_Uncut->Fill(P_hgcer_npeSum_kaons_uncut, P_aero_npeSum_kaons_uncut);
  h2_P_hgcer_npeSum_vs_aero_npeSum_kaons_Uncut->GetXaxis()->SetTitle("P_hgcer_npeSum");
  h2_P_hgcer_npeSum_vs_aero_npeSum_kaons_Uncut->GetYaxis()->SetTitle("P_aero_npeSum");
  h2_P_hgcer_npeSum_vs_aero_npeSum_kaons_Uncut->GetXaxis()->CenterTitle();
  h2_P_hgcer_npeSum_vs_aero_npeSum_kaons_Uncut->GetYaxis()->CenterTitle();

  h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Uncut->Fill(CTime_eKCoinTime_ROC1_kaons_uncut, P_kin_MMK_kaons_uncut);
  h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Uncut->GetXaxis()->SetTitle("CTime_eKCoinTime_ROC1");
  h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Uncut->GetYaxis()->SetTitle("P_kin_MMK");
  h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Uncut->GetXaxis()->CenterTitle();
  h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Uncut->GetYaxis()->CenterTitle();

  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_kaons_Uncut->Fill(P_hgcer_xAtCer_kaons_uncut, P_hgcer_yAtCer_kaons_uncut);
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_kaons_Uncut->GetXaxis()->SetTitle("P_hgcer_xAtCer");
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_kaons_Uncut->GetYaxis()->SetTitle("P_hgcer_yAtCer");
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_kaons_Uncut->GetXaxis()->CenterTitle();
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_kaons_Uncut->GetYaxis()->CenterTitle();
  }

  for(Long64_t i = 0; i < nEntries_Cut_Kaon_Events_All; i++){
    Cut_Kaon_Events_All->GetEntry(i);

  h2_P_hgcer_npeSum_vs_aero_npeSum_kaons_Cut_All->Fill(P_hgcer_npeSum_kaons_cut_all, P_aero_npeSum_kaons_cut_all);
  h2_P_hgcer_npeSum_vs_aero_npeSum_kaons_Cut_All->GetXaxis()->SetTitle("P_hgcer_npeSum");
  h2_P_hgcer_npeSum_vs_aero_npeSum_kaons_Cut_All->GetYaxis()->SetTitle("P_aero_npeSum");
  h2_P_hgcer_npeSum_vs_aero_npeSum_kaons_Cut_All->GetXaxis()->CenterTitle();
  h2_P_hgcer_npeSum_vs_aero_npeSum_kaons_Cut_All->GetYaxis()->CenterTitle();

  h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Cut_All->Fill(CTime_eKCoinTime_ROC1_kaons_cut_all, P_kin_MMK_kaons_cut_all);
  h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Cut_All->GetXaxis()->SetTitle("CTime_eKCoinTime_ROC1");
  h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Cut_All->GetYaxis()->SetTitle("P_kin_MMK");
  h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Cut_All->GetXaxis()->CenterTitle();
  h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Cut_All->GetYaxis()->CenterTitle();

  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_kaons_Cut_All->Fill(P_hgcer_xAtCer_kaons_cut_all, P_hgcer_yAtCer_kaons_cut_all);
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_kaons_Cut_All->GetXaxis()->SetTitle("P_hgcer_xAtCer");
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_kaons_Cut_All->GetYaxis()->SetTitle("P_hgcer_yAtCer");
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_kaons_Cut_All->GetXaxis()->CenterTitle();
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_kaons_Cut_All->GetYaxis()->CenterTitle();
  }

  for(Long64_t i = 0; i < nEntries_Cut_Kaon_Events_Prompt; i++){
    Cut_Kaon_Events_Prompt->GetEntry(i);

  h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Cut_Prompt->Fill(CTime_eKCoinTime_ROC1_kaons_cut_pr, P_kin_MMK_kaons_cut_pr);
  h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Cut_Prompt->GetXaxis()->SetTitle("CTime_eKCoinTime_ROC1");
  h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Cut_Prompt->GetYaxis()->SetTitle("P_kin_MMK");
  h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Cut_Prompt->GetXaxis()->CenterTitle();
  h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Cut_Prompt->GetYaxis()->CenterTitle();
  }

  for(Long64_t i = 0; i < nEntries_Uncut_Proton_Events; i++){
    Uncut_Proton_Events->GetEntry(i);

  h2_P_hgcer_npeSum_vs_aero_npeSum_protons_Uncut->Fill(P_hgcer_npeSum_protons_uncut, P_aero_npeSum_protons_uncut);
  h2_P_hgcer_npeSum_vs_aero_npeSum_protons_Uncut->GetXaxis()->SetTitle("P_hgcer_npeSum");
  h2_P_hgcer_npeSum_vs_aero_npeSum_protons_Uncut->GetYaxis()->SetTitle("P_aero_npeSum");
  h2_P_hgcer_npeSum_vs_aero_npeSum_protons_Uncut->GetXaxis()->CenterTitle();
  h2_P_hgcer_npeSum_vs_aero_npeSum_protons_Uncut->GetYaxis()->CenterTitle();

  h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Uncut->Fill(CTime_epCoinTime_ROC1_protons_uncut, P_kin_MMp_protons_uncut);
  h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Uncut->GetXaxis()->SetTitle("CTime_epCoinTime_ROC1");
  h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Uncut->GetYaxis()->SetTitle("P_kin_MMp");
  h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Uncut->GetXaxis()->CenterTitle();
  h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Uncut->GetYaxis()->CenterTitle();

  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_protons_Uncut->Fill(P_hgcer_xAtCer_protons_uncut, P_hgcer_yAtCer_protons_uncut);
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_protons_Uncut->GetXaxis()->SetTitle("P_hgcer_xAtCer");
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_protons_Uncut->GetYaxis()->SetTitle("P_hgcer_yAtCer");
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_protons_Uncut->GetXaxis()->CenterTitle();
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_protons_Uncut->GetYaxis()->CenterTitle();
  }

  for(Long64_t i = 0; i < nEntries_Cut_Proton_Events_All; i++){
    Cut_Proton_Events_All->GetEntry(i);

  h2_P_hgcer_npeSum_vs_aero_npeSum_protons_Cut_All->Fill(P_hgcer_npeSum_protons_cut_all, P_aero_npeSum_protons_cut_all);
  h2_P_hgcer_npeSum_vs_aero_npeSum_protons_Cut_All->GetXaxis()->SetTitle("P_hgcer_npeSum");
  h2_P_hgcer_npeSum_vs_aero_npeSum_protons_Cut_All->GetYaxis()->SetTitle("P_aero_npeSum");
  h2_P_hgcer_npeSum_vs_aero_npeSum_protons_Cut_All->GetXaxis()->CenterTitle();
  h2_P_hgcer_npeSum_vs_aero_npeSum_protons_Cut_All->GetYaxis()->CenterTitle();

  h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Cut_All->Fill(CTime_epCoinTime_ROC1_protons_cut_all, P_kin_MMp_protons_cut_all);
  h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Cut_All->GetXaxis()->SetTitle("CTime_epCoinTime_ROC1");
  h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Cut_All->GetYaxis()->SetTitle("P_kin_MMp");
  h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Cut_All->GetXaxis()->CenterTitle();
  h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Cut_All->GetYaxis()->CenterTitle();

  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_protons_Cut_All->Fill(P_hgcer_xAtCer_protons_cut_all, P_hgcer_yAtCer_protons_cut_all);
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_protons_Cut_All->GetXaxis()->SetTitle("P_hgcer_xAtCer");
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_protons_Cut_All->GetYaxis()->SetTitle("P_hgcer_yAtCer");
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_protons_Cut_All->GetXaxis()->CenterTitle();
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_protons_Cut_All->GetYaxis()->CenterTitle();
  }

  for(Long64_t i = 0; i < nEntries_Cut_Proton_Events_Prompt; i++){
    Cut_Proton_Events_Prompt->GetEntry(i);

  h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Cut_Prompt->Fill(CTime_epCoinTime_ROC1_protons_cut_pr, P_kin_MMp_protons_cut_pr);
  h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Cut_Prompt->GetXaxis()->SetTitle("CTime_epCoinTime_ROC1");
  h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Cut_Prompt->GetYaxis()->SetTitle("P_kin_MMp");
  h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Cut_Prompt->GetXaxis()->CenterTitle();
  h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Cut_Prompt->GetYaxis()->CenterTitle();
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for(Long64_t i = 0; i < nEntries_Cut_Pion_Events_Random; i++){
    Cut_Pion_Events_Random->GetEntry(i);
  h1_CTime_ePiCoinTime_ROC1_pions_Cut_Random_Scaled->Fill(CTime_ePiCoinTime_ROC1_pions_cut_rn);
  h1_CTime_ePiCoinTime_ROC1_pions_Cut_Random_Scaled->Scale(1.0/nWindows);
  }
  h1_CTime_ePiCoinTime_ROC1_pions_Rndm_Sub->Add(h1_CTime_ePiCoinTime_ROC1_pions_Cut_Prompt, h1_CTime_ePiCoinTime_ROC1_pions_Cut_Random_Scaled, 1, -1);
  h1_CTime_ePiCoinTime_ROC1_pions_Rndm_Sub->GetXaxis()->SetTitle("CTime_ePiCoinTime_ROC1");
  h1_CTime_ePiCoinTime_ROC1_pions_Rndm_Sub->GetYaxis()->SetTitle("Entries");
  h1_CTime_ePiCoinTime_ROC1_pions_Rndm_Sub->GetXaxis()->CenterTitle();
  h1_CTime_ePiCoinTime_ROC1_pions_Rndm_Sub->GetYaxis()->CenterTitle();

  for(Long64_t i = 0; i < nEntries_Cut_Pion_Events_Random; i++){
    Cut_Pion_Events_Random->GetEntry(i);
  h1_P_kin_MMpi_pions_Cut_Random_Scaled->Fill(P_kin_MMpi_pions_cut_rn);
  h1_P_kin_MMpi_pions_Cut_Random_Scaled->Scale(1.0/nWindows);
  }
  h1_P_kin_MMpi_pions_Rndm_Sub->Add(h1_P_kin_MMpi_pions_Cut_Prompt, h1_P_kin_MMpi_pions_Cut_Random_Scaled, 1, -1);
  h1_P_kin_MMpi_pions_Rndm_Sub->GetXaxis()->SetTitle("P_kin_MMpi");
  h1_P_kin_MMpi_pions_Rndm_Sub->GetYaxis()->SetTitle("Entries");
  h1_P_kin_MMpi_pions_Rndm_Sub->GetXaxis()->CenterTitle();
  h1_P_kin_MMpi_pions_Rndm_Sub->GetYaxis()->CenterTitle();

  for(Long64_t i = 0; i < nEntries_Cut_Kaon_Events_Random; i++){
    Cut_Kaon_Events_Random->GetEntry(i);
  h1_CTime_eKCoinTime_ROC1_kaons_Cut_Random_Scaled->Fill(CTime_eKCoinTime_ROC1_kaons_cut_rn);
  h1_CTime_eKCoinTime_ROC1_kaons_Cut_Random_Scaled->Scale(1.0/nWindows);
  }
  h1_CTime_eKCoinTime_ROC1_kaons_Rndm_Sub->Add(h1_CTime_eKCoinTime_ROC1_kaons_Cut_Prompt, h1_CTime_eKCoinTime_ROC1_kaons_Cut_Random_Scaled, 1, -1);
  h1_CTime_eKCoinTime_ROC1_kaons_Rndm_Sub->GetXaxis()->SetTitle("CTime_eKCoinTime_ROC1");
  h1_CTime_eKCoinTime_ROC1_kaons_Rndm_Sub->GetYaxis()->SetTitle("Entries");
  h1_CTime_eKCoinTime_ROC1_kaons_Rndm_Sub->GetXaxis()->CenterTitle();
  h1_CTime_eKCoinTime_ROC1_kaons_Rndm_Sub->GetYaxis()->CenterTitle();

  for(Long64_t i = 0; i < nEntries_Cut_Kaon_Events_Random; i++){
    Cut_Kaon_Events_Random->GetEntry(i);
  h1_P_kin_MMK_kaons_Cut_Random_Scaled->Fill(P_kin_MMK_kaons_cut_rn);
  h1_P_kin_MMK_kaons_Cut_Random_Scaled->Scale(1.0/nWindows);
  }
  h1_P_kin_MMK_kaons_Rndm_Sub->Add(h1_P_kin_MMK_kaons_Cut_Prompt, h1_P_kin_MMK_kaons_Cut_Random_Scaled, 1, -1);
  h1_P_kin_MMK_kaons_Rndm_Sub->GetXaxis()->SetTitle("P_kin_MMK");
  h1_P_kin_MMK_kaons_Rndm_Sub->GetYaxis()->SetTitle("Entries");
  h1_P_kin_MMK_kaons_Rndm_Sub->GetXaxis()->CenterTitle();
  h1_P_kin_MMK_kaons_Rndm_Sub->GetYaxis()->CenterTitle();

  for(Long64_t i = 0; i < nEntries_Cut_Proton_Events_Random; i++){
    Cut_Proton_Events_Random->GetEntry(i);
  h1_CTime_epCoinTime_ROC1_protons_Cut_Random_Scaled->Fill(CTime_epCoinTime_ROC1_protons_cut_rn);
  h1_CTime_epCoinTime_ROC1_protons_Cut_Random_Scaled->Scale(1.0/nWindows);
  }
  h1_CTime_epCoinTime_ROC1_protons_Rndm_Sub->Add(h1_CTime_epCoinTime_ROC1_protons_Cut_Prompt, h1_CTime_epCoinTime_ROC1_protons_Cut_Random_Scaled, 1, -1);
  h1_CTime_epCoinTime_ROC1_protons_Rndm_Sub->GetXaxis()->SetTitle("CTime_epCoinTime_ROC1");
  h1_CTime_epCoinTime_ROC1_protons_Rndm_Sub->GetYaxis()->SetTitle("Entries");
  h1_CTime_epCoinTime_ROC1_protons_Rndm_Sub->GetXaxis()->CenterTitle();
  h1_CTime_epCoinTime_ROC1_protons_Rndm_Sub->GetYaxis()->CenterTitle();

  for(Long64_t i = 0; i < nEntries_Cut_Proton_Events_Random; i++){
    Cut_Proton_Events_Random->GetEntry(i);
  h1_P_kin_MMp_protons_Cut_Random_Scaled->Fill(P_kin_MMp_protons_cut_rn);
  h1_P_kin_MMp_protons_Cut_Random_Scaled->Scale(1.0/nWindows);
  }
  h1_P_kin_MMp_protons_Rndm_Sub->Add(h1_P_kin_MMp_protons_Cut_Prompt, h1_P_kin_MMp_protons_Cut_Random_Scaled, 1, -1);
  h1_P_kin_MMp_protons_Rndm_Sub->GetXaxis()->SetTitle("P_kin_MMp");
  h1_P_kin_MMp_protons_Rndm_Sub->GetYaxis()->SetTitle("Entries");
  h1_P_kin_MMp_protons_Rndm_Sub->GetXaxis()->CenterTitle();
  h1_P_kin_MMp_protons_Rndm_Sub->GetYaxis()->CenterTitle();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 // Horrible block of stuff for polar plotting
  gPad->SetTheta(90); gPad->SetPhi(180);
  TPaveText *tvsphi_title = new TPaveText(0.0277092,0.89779,0.096428,0.991854,"NDC");
  tvsphi_title->AddText("-t vs #phi"); tvsphi_title->Draw();
  TPaveText *ptphizero = new TPaveText(0.923951,0.513932,0.993778,0.574551,"NDC");
  ptphizero->AddText("#phi = 0"); ptphizero->Draw();
  TLine *phihalfpi = new TLine(0,0,0,0.6); 
  phihalfpi->SetLineColor(kBlack); phihalfpi->SetLineWidth(2); phihalfpi->Draw();  
  TPaveText *ptphihalfpi = new TPaveText(0.417855,0.901876,0.486574,0.996358,"NDC");
  ptphihalfpi->AddText("#phi = #frac{#pi}{2}"); ptphihalfpi->Draw();
  TLine *phipi = new TLine(0,0,-0.6,0); 
  phipi->SetLineColor(kBlack); phipi->SetLineWidth(2); phipi->Draw();  
  TPaveText *ptphipi = new TPaveText(0.0277092,0.514217,0.096428,0.572746,"NDC");
  ptphipi->AddText("#phi = #pi"); ptphipi->Draw();
  TLine *phithreepi = new TLine(0,0,0,-0.6); 
  phithreepi->SetLineColor(kBlack); phithreepi->SetLineWidth(2); phithreepi->Draw();  
  TPaveText *ptphithreepi = new TPaveText(0.419517,0.00514928,0.487128,0.0996315,"NDC");
  ptphithreepi->AddText("#phi = #frac{3#pi}{2}"); ptphithreepi->Draw();
  TArc *Arc[10];

  for (Int_t k = 0; k < 10; k++){
    Arc[k] = new TArc(); 
    Arc[k]->SetFillStyle(0);
    Arc[k]->SetLineWidth(2);
    Arc[k]->DrawArc(0,0,0.575*(k+1)/(10),0.,360.,"same"); 
  }
  TGaxis *tradius = new TGaxis(0,0,0.575,0,0,2.0,10,"-+"); 
  tradius->SetLineColor(2);tradius->SetLabelColor(2);tradius->Draw();
  TLine *phizero = new TLine(0,0,0.6,0); 
:w  phizero->SetLineColor(kBlack); phizero->SetLineWidth(2); phizero->Draw();
  c_Kine->cd(4);
  h1_MMpi_BGSub->Draw("HIST");
  c_Kine->Print(foutpdf + ')');

*/

////////////////////////////////////////////////////////////////////////////////////////////////

  TFile *OutHisto_file = new TFile(foutname,"RECREATE");
  TDirectory *d_Uncut_Pion_Events = OutHisto_file->mkdir("Uncut_Pion_Events");
  TDirectory *d_Cut_Pion_Events_All = OutHisto_file->mkdir("Cut_Pion_Events_All");
  TDirectory *d_Cut_Pion_Events_Prompt = OutHisto_file->mkdir("Cut_Pion_Events_Prompt");
  TDirectory *d_Cut_Pion_Events_Random = OutHisto_file->mkdir("Cut_Pion_Events_Random");
  TDirectory *d_Uncut_Kaon_Events = OutHisto_file->mkdir("Uncut_Kaon_Events");
  TDirectory *d_Cut_Kaon_Events_All = OutHisto_file->mkdir("Cut_Kaon_Events_All");
  TDirectory *d_Cut_Kaon_Events_Prompt = OutHisto_file->mkdir("Cut_Kaon_Events_Prompt");
  TDirectory *d_Cut_Kaon_Events_Random = OutHisto_file->mkdir("Cut_Kaon_Events_Random");
  TDirectory *d_Uncut_Proton_Events = OutHisto_file->mkdir("Uncut_Proton_Events");
  TDirectory *d_Cut_Proton_Events_All = OutHisto_file->mkdir("Cut_Proton_Events_All");
  TDirectory *d_Cut_Proton_Events_Prompt = OutHisto_file->mkdir("Cut_Proton_Events_Prompt");
  TDirectory *d_Cut_Proton_Events_Random = OutHisto_file->mkdir("Cut_Proton_Events_Random");

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  d_Uncut_Pion_Events->cd();
  h1_H_gtr_beta_pions_Uncut->Write();
  h1_H_gtr_xp_pions_Uncut->Write();
  h1_H_gtr_yp_pions_Uncut->Write();
  h1_H_gtr_dp_pions_Uncut->Write();
  h1_H_hod_goodscinhit_pions_Uncut->Write();
  h1_H_hod_goodstarttime_pions_Uncut->Write();
  h1_H_cal_etotnorm_pions_Uncut->Write();
  h1_H_cal_etottracknorm_pions_Uncut->Write();
  h1_H_cer_npeSum_pions_Uncut->Write();
  h1_H_RFTime_Dist_pions_Uncut->Write();
  h1_P_gtr_beta_pions_Uncut->Write(); 
  h1_P_gtr_xp_pions_Uncut->Write();
  h1_P_gtr_yp_pions_Uncut->Write();
  h1_P_gtr_dp_pions_Uncut->Write();
  h1_P_gtr_p_pions_Uncut->Write();
  h1_P_hod_goodscinhit_pions_Uncut->Write();
  h1_P_hod_goodstarttime_pions_Uncut->Write();
  h1_P_cal_etotnorm_pions_Uncut->Write();
  h1_P_cal_etottracknorm_pions_Uncut->Write();
  h1_P_hgcer_npeSum_pions_Uncut->Write();
  h1_P_hgcer_xAtCer_pions_Uncut->Write();
  h1_P_hgcer_yAtCer_pions_Uncut->Write();
  h1_P_aero_npeSum_pions_Uncut->Write();
  h1_P_aero_xAtAero_pions_Uncut->Write();
  h1_P_aero_yAtAero_pions_Uncut->Write();
  h1_P_kin_MMpi_pions_Uncut->Write();
  h1_P_RFTime_Dist_pions_Uncut->Write();      
  h1_CTime_ePiCoinTime_ROC1_pions_Uncut->Write();
  h1_Q2_pions_Uncut->Write();
  h1_W_pions_Uncut->Write();
  h1_epsilon_pions_Uncut->Write();

  h2_P_hgcer_npeSum_vs_aero_npeSum_pions_Uncut->Write();
  h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Uncut->Write();
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_pions_Uncut->Write();

  d_Cut_Pion_Events_All->cd();
  h1_H_gtr_beta_pions_Cut_All->Write();
  h1_H_gtr_xp_pions_Cut_All->Write();
  h1_H_gtr_yp_pions_Cut_All->Write();
  h1_H_gtr_dp_pions_Cut_All->Write();
  h1_H_hod_goodscinhit_pions_Cut_All->Write();
  h1_H_hod_goodstarttime_pions_Cut_All->Write();
  h1_H_cal_etotnorm_pions_Cut_All->Write();
  h1_H_cal_etottracknorm_pions_Cut_All->Write();
  h1_H_cer_npeSum_pions_Cut_All->Write();
  h1_H_RFTime_Dist_pions_Cut_All->Write();
  h1_P_gtr_beta_pions_Cut_All->Write();
  h1_P_gtr_xp_pions_Cut_All->Write();
  h1_P_gtr_yp_pions_Cut_All->Write();
  h1_P_gtr_dp_pions_Cut_All->Write();
  h1_P_gtr_p_pions_Cut_All->Write();
  h1_P_hod_goodscinhit_pions_Cut_All->Write();
  h1_P_hod_goodstarttime_pions_Cut_All->Write();
  h1_P_cal_etotnorm_pions_Cut_All->Write();
  h1_P_cal_etottracknorm_pions_Cut_All->Write();
  h1_P_hgcer_npeSum_pions_Cut_All->Write();
  h1_P_hgcer_xAtCer_pions_Cut_All->Write();
  h1_P_hgcer_yAtCer_pions_Cut_All->Write();
  h1_P_aero_npeSum_pions_Cut_All->Write();
  h1_P_aero_xAtAero_pions_Cut_All->Write();
  h1_P_aero_yAtAero_pions_Cut_All->Write();
  h1_P_kin_MMpi_pions_Cut_All->Write();
  h1_P_RFTime_Dist_pions_Cut_All->Write();         
  h1_CTime_ePiCoinTime_ROC1_pions_Cut_All->Write();
  h1_Q2_pions_Cut_All->Write();
  h1_W_pions_Cut_All->Write();
  h1_epsilon_pions_Cut_All->Write();

  h1_CTime_ePiCoinTime_ROC1_pions_Rndm_Sub->Write();
  h1_P_kin_MMpi_pions_Rndm_Sub->Write();

  h2_P_hgcer_npeSum_vs_aero_npeSum_pions_Cut_All->Write();
  h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Cut_All->Write();
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_pions_Cut_All->Write();

  d_Cut_Pion_Events_Prompt->cd();
  h1_P_gtr_beta_pions_Cut_Prompt->Write();
  h1_P_RFTime_Dist_pions_Cut_Prompt->Write();
  h1_CTime_ePiCoinTime_ROC1_pions_Cut_Prompt->Write();
  h1_P_kin_MMpi_pions_Cut_Prompt->Write();

//  h2_P_hgcer_npeSum_vs_aero_npeSum_pions_Cut_Prompt->Write();
  h2_CTime_ePiCoinTime_ROC1_vs_P_kin_MMpi_pions_Cut_Prompt->Write();

  d_Cut_Pion_Events_Random->cd();
  h1_P_gtr_beta_pions_Cut_Random->Write();
  h1_P_RFTime_Dist_pions_Cut_Random->Write();
  h1_CTime_ePiCoinTime_ROC1_pions_Cut_Random->Write();
  h1_P_kin_MMpi_pions_Cut_Random->Write();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  d_Uncut_Kaon_Events->cd();
  h1_H_gtr_beta_kaons_Uncut->Write();
  h1_H_gtr_xp_kaons_Uncut->Write();
  h1_H_gtr_yp_kaons_Uncut->Write();
  h1_H_gtr_dp_kaons_Uncut->Write();
  h1_H_hod_goodscinhit_kaons_Uncut->Write();
  h1_H_hod_goodstarttime_kaons_Uncut->Write();
  h1_H_cal_etotnorm_kaons_Uncut->Write();
  h1_H_cal_etottracknorm_kaons_Uncut->Write();
  h1_H_cer_npeSum_kaons_Uncut->Write();
  h1_H_RFTime_Dist_kaons_Uncut->Write();
  h1_P_gtr_beta_kaons_Uncut->Write();
  h1_P_gtr_xp_kaons_Uncut->Write();
  h1_P_gtr_yp_kaons_Uncut->Write();
  h1_P_gtr_dp_kaons_Uncut->Write();
  h1_P_gtr_p_kaons_Uncut->Write();
  h1_P_hod_goodscinhit_kaons_Uncut->Write();
  h1_P_hod_goodstarttime_kaons_Uncut->Write();
  h1_P_cal_etotnorm_kaons_Uncut->Write();
  h1_P_cal_etottracknorm_kaons_Uncut->Write();
  h1_P_hgcer_npeSum_kaons_Uncut->Write();
  h1_P_hgcer_xAtCer_kaons_Uncut->Write();
  h1_P_hgcer_yAtCer_kaons_Uncut->Write();
  h1_P_aero_npeSum_kaons_Uncut->Write();
  h1_P_aero_xAtAero_kaons_Uncut->Write();
  h1_P_aero_yAtAero_kaons_Uncut->Write();
  h1_P_kin_MMK_kaons_Uncut->Write();
  h1_P_RFTime_Dist_kaons_Uncut->Write();
  h1_CTime_eKCoinTime_ROC1_kaons_Uncut->Write();
  h1_Q2_kaons_Uncut->Write();
  h1_W_kaons_Uncut->Write();
  h1_epsilon_kaons_Uncut->Write();

  h2_P_hgcer_npeSum_vs_aero_npeSum_kaons_Uncut->Write();
  h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Uncut->Write();
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_kaons_Uncut->Write();

  d_Cut_Kaon_Events_All->cd();
  h1_H_gtr_beta_kaons_Cut_All->Write();
  h1_H_gtr_xp_kaons_Cut_All->Write();
  h1_H_gtr_yp_kaons_Cut_All->Write();
  h1_H_gtr_dp_kaons_Cut_All->Write();
  h1_H_hod_goodscinhit_kaons_Cut_All->Write();
  h1_H_hod_goodstarttime_kaons_Cut_All->Write();
  h1_H_cal_etotnorm_kaons_Cut_All->Write();
  h1_H_cal_etottracknorm_kaons_Cut_All->Write();
  h1_H_cer_npeSum_kaons_Cut_All->Write();
  h1_H_RFTime_Dist_kaons_Cut_All->Write();
  h1_P_gtr_beta_kaons_Cut_All->Write();
  h1_P_gtr_xp_kaons_Cut_All->Write();
  h1_P_gtr_yp_kaons_Cut_All->Write();
  h1_P_gtr_dp_kaons_Cut_All->Write();
  h1_P_gtr_p_kaons_Cut_All->Write();
  h1_P_hod_goodscinhit_kaons_Cut_All->Write();
  h1_P_hod_goodstarttime_kaons_Cut_All->Write();
  h1_P_cal_etotnorm_kaons_Cut_All->Write();
  h1_P_cal_etottracknorm_kaons_Cut_All->Write();
  h1_P_hgcer_npeSum_kaons_Cut_All->Write();
  h1_P_hgcer_xAtCer_kaons_Cut_All->Write();
  h1_P_hgcer_yAtCer_kaons_Cut_All->Write();
  h1_P_aero_npeSum_kaons_Cut_All->Write();
  h1_P_aero_xAtAero_kaons_Cut_All->Write();
  h1_P_aero_yAtAero_kaons_Cut_All->Write();
  h1_P_kin_MMK_kaons_Cut_All->Write();
  h1_P_RFTime_Dist_kaons_Cut_All->Write();
  h1_CTime_eKCoinTime_ROC1_kaons_Cut_All->Write();
  h1_Q2_kaons_Cut_All->Write();
  h1_W_kaons_Cut_All->Write();
  h1_epsilon_kaons_Cut_All->Write();

  h1_CTime_eKCoinTime_ROC1_kaons_Rndm_Sub->Write();
  h1_P_kin_MMK_kaons_Rndm_Sub->Write();

  h2_P_hgcer_npeSum_vs_aero_npeSum_kaons_Cut_All->Write();
  h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Cut_All->Write();
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_kaons_Cut_All->Write();

  d_Cut_Kaon_Events_Prompt->cd();
  h1_P_gtr_beta_kaons_Cut_Prompt->Write();
  h1_P_RFTime_Dist_kaons_Cut_Prompt->Write();
  h1_CTime_eKCoinTime_ROC1_kaons_Cut_Prompt->Write();
  h1_P_kin_MMK_kaons_Cut_Prompt->Write();

//  h2_P_hgcer_npeSum_vs_aero_npeSum_kaons_Cut_Prompt->Write();
  h2_CTime_eKCoinTime_ROC1_vs_P_kin_MMK_kaons_Cut_Prompt->Write();

  d_Cut_Kaon_Events_Random->cd();
  h1_P_gtr_beta_kaons_Cut_Random->Write();
  h1_P_RFTime_Dist_kaons_Cut_Random->Write();
  h1_CTime_eKCoinTime_ROC1_kaons_Cut_Random->Write();
  h1_P_kin_MMK_kaons_Cut_Random->Write();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  d_Uncut_Proton_Events->cd();
  h1_H_gtr_beta_protons_Uncut->Write();
  h1_H_gtr_xp_protons_Uncut->Write();
  h1_H_gtr_yp_protons_Uncut->Write();
  h1_H_gtr_dp_protons_Uncut->Write();
  h1_H_hod_goodscinhit_protons_Uncut->Write();
  h1_H_hod_goodstarttime_protons_Uncut->Write();
  h1_H_cal_etotnorm_protons_Uncut->Write();
  h1_H_cal_etottracknorm_protons_Uncut->Write();
  h1_H_cer_npeSum_protons_Uncut->Write();
  h1_H_RFTime_Dist_protons_Uncut->Write();
  h1_P_gtr_beta_protons_Uncut->Write();
  h1_P_gtr_xp_protons_Uncut->Write();
  h1_P_gtr_yp_protons_Uncut->Write();
  h1_P_gtr_dp_protons_Uncut->Write();
  h1_P_gtr_p_protons_Uncut->Write();
  h1_P_hod_goodscinhit_protons_Uncut->Write();
  h1_P_hod_goodstarttime_protons_Uncut->Write();
  h1_P_cal_etotnorm_protons_Uncut->Write();
  h1_P_cal_etottracknorm_protons_Uncut->Write();
  h1_P_hgcer_npeSum_protons_Uncut->Write();
  h1_P_hgcer_xAtCer_protons_Uncut->Write();
  h1_P_hgcer_yAtCer_protons_Uncut->Write();
  h1_P_aero_npeSum_protons_Uncut->Write();
  h1_P_aero_xAtAero_protons_Uncut->Write();
  h1_P_aero_yAtAero_protons_Uncut->Write();
  h1_P_kin_MMp_protons_Uncut->Write();
  h1_P_RFTime_Dist_protons_Uncut->Write();
  h1_CTime_epCoinTime_ROC1_protons_Uncut->Write();
  h1_Q2_protons_Uncut->Write();
  h1_W_protons_Uncut->Write();
  h1_epsilon_protons_Uncut->Write();

  h2_P_hgcer_npeSum_vs_aero_npeSum_protons_Uncut->Write();
  h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Uncut->Write();
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_protons_Uncut->Write();

  d_Cut_Proton_Events_All->cd();
  h1_H_gtr_beta_protons_Cut_All->Write();
  h1_H_gtr_xp_protons_Cut_All->Write();
  h1_H_gtr_yp_protons_Cut_All->Write();
  h1_H_gtr_dp_protons_Cut_All->Write();
  h1_H_hod_goodscinhit_protons_Cut_All->Write();
  h1_H_hod_goodstarttime_protons_Cut_All->Write();
  h1_H_cal_etotnorm_protons_Cut_All->Write();
  h1_H_cal_etottracknorm_protons_Cut_All->Write();
  h1_H_cer_npeSum_protons_Cut_All->Write();
  h1_H_RFTime_Dist_protons_Cut_All->Write();
  h1_P_gtr_beta_protons_Cut_All->Write();
  h1_P_gtr_xp_protons_Cut_All->Write();
  h1_P_gtr_yp_protons_Cut_All->Write();
  h1_P_gtr_dp_protons_Cut_All->Write();
  h1_P_gtr_p_protons_Cut_All->Write();
  h1_P_hod_goodscinhit_protons_Cut_All->Write();
  h1_P_hod_goodstarttime_protons_Cut_All->Write();
  h1_P_cal_etotnorm_protons_Cut_All->Write();
  h1_P_cal_etottracknorm_protons_Cut_All->Write();
  h1_P_hgcer_npeSum_protons_Cut_All->Write();
  h1_P_hgcer_xAtCer_protons_Cut_All->Write();
  h1_P_hgcer_yAtCer_protons_Cut_All->Write();
  h1_P_aero_npeSum_protons_Cut_All->Write();
  h1_P_aero_xAtAero_protons_Cut_All->Write();
  h1_P_aero_yAtAero_protons_Cut_All->Write();
  h1_P_kin_MMp_protons_Cut_All->Write();
  h1_P_RFTime_Dist_protons_Cut_All->Write();
  h1_CTime_epCoinTime_ROC1_protons_Cut_All->Write();
  h1_Q2_protons_Cut_All->Write();
  h1_W_protons_Cut_All->Write();
  h1_epsilon_protons_Cut_All->Write();

  h1_CTime_epCoinTime_ROC1_protons_Rndm_Sub->Write();
  h1_P_kin_MMp_protons_Rndm_Sub->Write();

  h2_P_hgcer_npeSum_vs_aero_npeSum_protons_Cut_All->Write();
  h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Cut_All->Write();
  h2_P_hgcer_xAtCer_vs_hgcer_yAtCer_protons_Cut_All->Write();

  d_Cut_Proton_Events_Prompt->cd();
  h1_P_gtr_beta_protons_Cut_Prompt->Write();
  h1_P_RFTime_Dist_protons_Cut_Prompt->Write();
  h1_CTime_epCoinTime_ROC1_protons_Cut_Prompt->Write();
  h1_P_kin_MMp_protons_Cut_Prompt->Write();

//  h2_P_hgcer_npeSum_vs_aero_npeSum_protons_Cut_Prompt->Write();
  h2_CTime_epCoinTime_ROC1_vs_P_kin_MMp_protons_Cut_Prompt->Write();

  d_Cut_Proton_Events_Random->cd();
  h1_P_gtr_beta_protons_Cut_Random->Write();
  h1_P_RFTime_Dist_protons_Cut_Random->Write();
  h1_CTime_epCoinTime_ROC1_protons_Cut_Random->Write();
  h1_P_kin_MMp_protons_Cut_Random->Write();


  OutHisto_file->Close();
}
