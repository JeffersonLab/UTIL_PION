// 20/10/20 - Stephen Kay, University of Regina

// root .c macro to plot output. Reads in a csv file and converts to a tree
#define PlotKinematic_cxx

// Include relevant stuff
#include <TStyle.h>
#include <TGraph.h>
#include <TGraphErrors.h>
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

// Input should be the input root file name (including suffix) and an output file name string (without any suffix)
void PlotKinematic(string Kinematic = "")
{
  TString Hostname = gSystem->HostName();
  TString User = (gSystem->GetUserInfo())->fUser;
  TString Replaypath;
  TString Outpath;
  TString rootFile;
  gStyle->SetPalette(55);

  // Set paths depending on system you're running on
  if(Hostname.Contains("farm")){
    Replaypath = "/group/c-kaonlt/USERS/"+User+"/hallc_replay_lt";
    Outpath = Replaypath+"/UTIL_KAONLT/scripts/CoinTimePeak/OUTPUT";
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
  if(Kinematic == "") {
    cout << "Enter a kinematic setting to plot: ";
    cin >> Kinematic;
  }
  TString TKinematic=Kinematic;
  TString KinCsv = Outpath+"/"+TKinematic+"_Output.csv";
  TString foutpdf = Outpath + "/" + TKinematic + "_CTPeakPlots.pdf";
  TString foutroot = Outpath + "/" + TKinematic + "_CTPeakPlots.root";

  if (gSystem->AccessPathName(KinCsv) == kTRUE){
    cerr << "!!!!! ERROR !!!!! " << endl << KinCsv <<  " not found" << endl <<  "!!!!! ERRROR !!!!!" << endl;
    exit;
  }
  // Slightly labourious but converts csv to tree, then reads all the info from the tree and assigns it to double arrays
  // Can then define our TGraphErrors using the arrays
  TTree *KinData = new TTree("t","tree from"+KinCsv);
  KinData->ReadFile(KinCsv, "RunNumber/D:PionPeak:PionPeakErr:PionWidth:PionWidthErr:KaonPeak:KaonPeakErr:KaonWidth:KaonWidthErr:ProtonPeak:ProtonPeakErr:ProtonWidth:ProtonWidthErr");
  Long64_t nEntries = KinData->GetEntries();
  Double_t RunNum; Double_t RunNumArr[nEntries]; KinData->SetBranchAddress("RunNumber", &RunNum); 
  Double_t PiCTPeak; Double_t PiCTPeakArr[nEntries]; KinData->SetBranchAddress("PionPeak", &PiCTPeak);
  Double_t PiCTPeakErr; Double_t PiCTPeakErrArr[nEntries]; KinData->SetBranchAddress("PionPeakErr", &PiCTPeakErr);
  Double_t PiCTWidth; Double_t PiCTWidthArr[nEntries]; KinData->SetBranchAddress("PionWidth", &PiCTWidth);
  Double_t PiCTWidthErr; Double_t PiCTWidthErrArr[nEntries]; KinData->SetBranchAddress("PionWidthErr", &PiCTWidthErr);
  Double_t KCTPeak; Double_t KCTPeakArr[nEntries]; KinData->SetBranchAddress("KaonPeak", &KCTPeak);
  Double_t KCTPeakErr; Double_t KCTPeakErrArr[nEntries]; KinData->SetBranchAddress("KaonPeakErr", &KCTPeakErr);
  Double_t KCTWidth; Double_t KCTWidthArr[nEntries]; KinData->SetBranchAddress("KaonWidth", &KCTWidth);
  Double_t KCTWidthErr; Double_t KCTWidthErrArr[nEntries]; KinData->SetBranchAddress("KaonWidthErr", &KCTWidthErr);
  Double_t pCTPeak; Double_t pCTPeakArr[nEntries]; KinData->SetBranchAddress("ProtonPeak", &pCTPeak);
  Double_t pCTPeakErr; Double_t pCTPeakErrArr[nEntries]; KinData->SetBranchAddress("ProtonPeakErr", &pCTPeakErr);
  Double_t pCTWidth; Double_t pCTWidthArr[nEntries]; KinData->SetBranchAddress("ProtonWidth", &pCTWidth);
  Double_t pCTWidthErr; Double_t pCTWidthErrArr[nEntries]; KinData->SetBranchAddress("ProtonWidthErr", &pCTWidthErr);

  // We also need x errors annoyingly so just make an array and fill it with 0
  Double_t RunNumErr[nEntries]; // No actual error on this but needed for TGraph errors

  for(Long64_t i = 0; i < nEntries; i++){
    KinData->GetEntry(i);
    RunNumArr[i]=RunNum;
    RunNumErr[i]=0;
    PiCTPeakArr[i]=PiCTPeak;
    PiCTPeakErrArr[i]=PiCTPeakErr;
    PiCTWidthArr[i]=PiCTWidth;
    PiCTWidthErrArr[i]=PiCTWidthErr;
    KCTPeakArr[i]=KCTPeak;
    KCTPeakErrArr[i]=KCTPeakErr;
    KCTWidthArr[i]=KCTWidth;
    KCTWidthErrArr[i]=KCTWidthErr;
    pCTPeakArr[i]=pCTPeak;
    pCTPeakErrArr[i]=pCTPeakErr;
    pCTWidthArr[i]=pCTWidth;
    pCTWidthErrArr[i]=pCTWidthErr;
  }

  auto PiCTPeakPlot = new TGraphErrors(nEntries, RunNumArr, PiCTPeakArr, RunNumErr, PiCTPeakErrArr);
  PiCTPeakPlot->SetTitle(Form("%s #pi CT Peak Position; RunNumber; CT Peak Position/ns", Kinematic.c_str()));
  PiCTPeakPlot->SetMarkerStyle(22); PiCTPeakPlot->SetMarkerSize(1.2); PiCTPeakPlot->SetMarkerColor(2); PiCTPeakPlot->SetLineColor(2);
  auto PiCTWidthPlot = new TGraphErrors(nEntries, RunNumArr, PiCTWidthArr, RunNumErr, PiCTWidthErrArr);
  PiCTWidthPlot->SetTitle(Form("%s #pi CT Peak Width; RunNumber; CT Peak Width/ns", Kinematic.c_str()));
  PiCTWidthPlot->SetMarkerStyle(22); PiCTWidthPlot->SetMarkerSize(1.2); PiCTWidthPlot->SetMarkerColor(2); PiCTWidthPlot->SetLineColor(2);  
  auto KCTPeakPlot = new TGraphErrors(nEntries, RunNumArr, KCTPeakArr, RunNumErr, KCTPeakErrArr);
  KCTPeakPlot->SetTitle(Form("%s K CT Peak Position; RunNumber; CT Peak Position/ns", Kinematic.c_str()));
  KCTPeakPlot->SetMarkerStyle(22); KCTPeakPlot->SetMarkerSize(1.2); KCTPeakPlot->SetMarkerColor(2); KCTPeakPlot->SetLineColor(2);
  auto KCTWidthPlot = new TGraphErrors(nEntries, RunNumArr, KCTWidthArr, RunNumErr, KCTWidthErrArr);
  KCTWidthPlot->SetTitle(Form("%s K CT Peak Width; RunNumber; CT Peak Width/ns", Kinematic.c_str()));
  KCTWidthPlot->SetMarkerStyle(22); KCTWidthPlot->SetMarkerSize(1.2); KCTWidthPlot->SetMarkerColor(2); KCTWidthPlot->SetLineColor(2);
  auto pCTPeakPlot = new TGraphErrors(nEntries, RunNumArr, pCTPeakArr, RunNumErr, pCTPeakErrArr);
  pCTPeakPlot->SetTitle(Form("%s p CT Peak Position; RunNumber; CT Peak Position/ns", Kinematic.c_str()));
  pCTPeakPlot->SetMarkerStyle(22); pCTPeakPlot->SetMarkerSize(1.2); pCTPeakPlot->SetMarkerColor(2); pCTPeakPlot->SetLineColor(2);
  auto pCTWidthPlot = new TGraphErrors(nEntries, RunNumArr, pCTWidthArr, RunNumErr, pCTWidthErrArr);
  pCTWidthPlot->SetTitle(Form("%s p CT Peak Width; RunNumber; CT Peak Width/ns", Kinematic.c_str()));
  pCTWidthPlot->SetMarkerStyle(22); pCTWidthPlot->SetMarkerSize(1.2); pCTWidthPlot->SetMarkerColor(2); pCTWidthPlot->SetLineColor(2);

  TCanvas *c_Summary = new TCanvas("c_Summary", "Summary", 100, 0, 1000, 900);
  c_Summary->Divide(2,3);
  c_Summary->cd(1);
  PiCTPeakPlot->Draw("AEP");
  c_Summary->cd(2);
  PiCTWidthPlot->Draw("AEP");
  c_Summary->cd(3);
  KCTPeakPlot->Draw("AEP");
  c_Summary->cd(4);
  KCTWidthPlot->Draw("AEP");
  c_Summary->cd(5);
  pCTPeakPlot->Draw("AEP");
  c_Summary->cd(6);
  pCTWidthPlot->Draw("AEP");

  c_Summary->Print(foutpdf + '(');

  TCanvas *c_PiPeakPos = new TCanvas("PiPeakPos", "#pi CT Peak Position", 100, 0, 1000, 900);
  PiCTPeakPlot->Draw("AEP");
  c_PiPeakPos->Print(foutpdf);
  TCanvas *c_PiPeakWid = new TCanvas("PiPeakWid", "#pi CT Peak Width", 100, 0, 1000, 900);
  PiCTWidthPlot->Draw("AEP");
  c_PiPeakWid->Print(foutpdf);
  TCanvas *c_KPeakPos = new TCanvas("KPeakPos", "K CT Peak Position", 100, 0, 1000, 900);
  KCTPeakPlot->Draw("AEP");
  c_KPeakPos->Print(foutpdf);
  TCanvas *c_KPeakWid = new TCanvas("KPeakWid", "K CT Peak Width", 100, 0, 1000, 900);
  KCTWidthPlot->Draw("AEP");
  c_KPeakWid->Print(foutpdf);
  TCanvas *c_pPeakPos = new TCanvas("pPeakPos", "p CT Peak Position", 100, 0, 1000, 900);
  pCTPeakPlot->Draw("AEP");
  c_pPeakPos->Print(foutpdf);
  TCanvas *c_pPeakWid = new TCanvas("pPeakWid", "p CT Peak Width", 100, 0, 1000, 900);
  pCTWidthPlot->Draw("AEP");
  c_pPeakWid->Print(foutpdf + ')');

  TFile *OutRoot_file = new TFile(foutroot,"RECREATE");
  PiCTPeakPlot->Write();
  PiCTWidthPlot->Write();
  KCTPeakPlot->Write();
  KCTWidthPlot->Write();
  pCTPeakPlot->Write();
  pCTWidthPlot->Write();
  OutRoot_file->Close();

}
