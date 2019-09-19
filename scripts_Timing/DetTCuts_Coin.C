#define DetTCuts_Coin_cxx

#include "DetTCuts_Coin.h"
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TSpectrum.h>
#include <TList.h>
#include <TPolyMarker.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "TProof.h"

using namespace TMath;

void DetTCuts_Coin::Begin(TTree * /*tree*/)
{
  printf("\n\n");
  TString option = GetOption();
}

void DetTCuts_Coin::SlaveBegin(TTree * /*tree*/)
{
  printf("\n\n");
  TString option = GetOption();

  //HMS Histos
  for (Int_t ipmt = 0; ipmt < 2; ipmt++){
    h1hCerAdcTdcTDiff[ipmt] = new TH1F (Form("hCER%d_timeDiff", ipmt+1), Form("HMS Cer PMT%d AdcTdcTimeDiff", ipmt+1), 200, 0, 200);
    GetOutputList()->Add(h1hCerAdcTdcTDiff[ipmt]);
  }
  for (Int_t i = 0; i < 12; i++){
    h1hdcTdcT[i] = new TH1F( Form("hDC%s_rawTDC", dc_pl_names[i].c_str()), Form("HMS DC Plane %s Raw TDC", dc_pl_names[i].c_str()), 200, -20000, 0);
    GetOutputList()->Add(h1hdcTdcT[i]);
  }
  for (Int_t npl = 0; npl < cal_planes; npl++){ // Loop over all calorimeter planes
    for (Int_t nside = 0; nside < sides; nside++){ //Loop over each side
      for (Int_t ipmt = 0; ipmt < 13; ipmt++){ // Loop over each PMT in a particular plane	
	h1hCalAdcTdcTDiff[npl][nside][ipmt] = new TH1F(Form("hCal%s%d%s_timeDiff", cal_pl_names[npl].c_str(), ipmt+1, nsign[nside].c_str()), Form("HMS Cal %s%d%s AdcTdcTimeDiff", cal_pl_names[npl].c_str(), ipmt+1, nsign[nside].c_str()), 200, -200, 200);
	GetOutputList()->Add(h1hCalAdcTdcTDiff[npl][nside][ipmt]);
      }
    }
  }
  for (Int_t npl = 0; npl < hod_planes; npl++){ // Loop over all hodoscope planes
    for (Int_t nside = 0; nside < sides; nside++){ //Loop over each side
      for (Int_t ipmt = 0; ipmt < hmaxPMT[npl]; ipmt++){ // Loop over each PMT in a particular plane	
	h1hHodoAdcTdcTDiff[npl][nside][ipmt] = new TH1F(Form("hHodo%s%d%s_timeDiff", hod_pl_names[npl].c_str(),ipmt+1,nsign[nside].c_str() ), Form("HMS Hodo %s%d%s AdcTdcTimeDiff", hod_pl_names[npl].c_str(), ipmt+1, nsign[nside].c_str()), 200, -100, 100);
	GetOutputList()->Add(h1hHodoAdcTdcTDiff[npl][nside][ipmt]);
      }
    }
  }
  
  // SHMS Histos
  for (Int_t ipmt = 0; ipmt < 4; ipmt++){
    h1pHGCAdcTdcTDiff[ipmt] = new TH1F (Form("pHGCER%d_timeDiff", ipmt+1), Form("SHMS HGCer PMT%d AdcTdcTimeDiff", ipmt+1), 200, 0, 200);
    GetOutputList()->Add(h1pHGCAdcTdcTDiff[ipmt]);
  }
  for (Int_t nside = 0; nside < sides; nside++){ //Loop over each side
   for (Int_t ipmt = 0; ipmt < 7; ipmt++){ // Loop over PMTs
     h1pAeroAdcTdcTDiff[nside][ipmt] = new TH1F(Form("pAero%d%s_timeDiff", ipmt+1, nsign[nside].c_str()), Form("SHMS Aerogel PMT%d%s AdcTdcTimeDiff", ipmt+1, nsign[nside].c_str()), 200, 0, 400);
     GetOutputList()->Add(h1pAeroAdcTdcTDiff[nside][ipmt]);
   }
  }
  for (Int_t i = 0; i < 12; i++){
    h1pdcTdcT[i] = new TH1F( Form("pDC%s_rawTDC", dc_pl_names[i].c_str()), Form("SHMS DC Plane %s Raw TDC", dc_pl_names[i].c_str()), 200, -20000, 0);
    GetOutputList()->Add(h1pdcTdcT[i]);
  }
  for (Int_t npl = 0; npl < hod_planes; npl++){ // Loop over all hodoscope planes
    for (Int_t nside = 0; nside < sides; nside++){ //Loop over each side
      for (Int_t ipmt = 0; ipmt < pmaxPMT[npl]; ipmt++){ // Loop over each PMT in a particular plane	
	h1pHodoAdcTdcTDiff[npl][nside][ipmt] = new TH1F(Form("pHodo%s%d%s_timeDiff", hod_pl_names[npl].c_str(),ipmt+1,nsign[nside].c_str() ), Form("SHMS Hodo %s%d%s AdcTdcTimeDiff", hod_pl_names[npl].c_str(), ipmt+1,nsign[nside].c_str()), 200, -100, 100);
	GetOutputList()->Add(h1pHodoAdcTdcTDiff[npl][nside][ipmt]);
      }
    }
  }
 for (Int_t nside = 0; nside < sides; nside++){ //Loop over each side
   for (Int_t ipmt = 0; ipmt < 14; ipmt++){ // Loop over PMTs
     h1pPrShAdcTdcTDiff[nside][ipmt] = new TH1F(Form("pPrSh%d%s_timeDiff", ipmt+1, nsign[nside].c_str()), Form("SHMS Pre-Shower PMT%d%s AdcTdcTimeDiff", ipmt+1, nsign[nside].c_str()), 200, -200, 200);
     GetOutputList()->Add(h1pPrShAdcTdcTDiff[nside][ipmt]);
   }
 }
 for(Int_t ipmt = 0; ipmt < 224; ipmt++){
   h1pCalAdcTdcTDiff[ipmt] = new TH1F(Form("pCalPMT%d", ipmt+1), Form("SHMS Calorimeter PMT%d AdcTdcTimeDiff", ipmt+1), 200, -100, 100); 
   GetOutputList()->Add(h1pCalAdcTdcTDiff[ipmt]);
 }

}

Bool_t DetTCuts_Coin::Process(Long64_t entry)
{
  fReader.SetEntry(entry);
  
  // Fill our HMS timing histograms, explicitly select only multiplicity 1 events
  for (Int_t ipmt = 0; ipmt < 2; ipmt++){
    if(H_cer_goodAdcMult[ipmt] == 1) h1hCerAdcTdcTDiff[ipmt]->Fill(H_cer_goodAdcTdcDiffTime[ipmt]);
  }
  // This is a disgustingly bad way of doing this, really need to figure out how to have some ARRAY of readers
  for (Int_t i = 0; i < 12; i++){
    if(i == 0){
      if(H_dc_1u1_nhit[0] == 1) h1hdcTdcT[i]->Fill(H_dc_1u1_rawtdc[0]);
    }
    else if(i == 1){
      if(H_dc_1u2_nhit[0] == 1) h1hdcTdcT[i]->Fill(H_dc_1u2_rawtdc[0]);
    }
    else if(i == 2){
      if(H_dc_1x1_nhit[0] == 1) h1hdcTdcT[i]->Fill(H_dc_1x1_rawtdc[0]);
    }
    else if(i == 3){
      if(H_dc_1x2_nhit[0] == 1) h1hdcTdcT[i]->Fill(H_dc_1x2_rawtdc[0]);
    }
    else if(i == 4){
      if(H_dc_1v1_nhit[0] == 1) h1hdcTdcT[i]->Fill(H_dc_1v1_rawtdc[0]);
    }
    else if(i == 5){
      if(H_dc_1v2_nhit[0] == 1) h1hdcTdcT[i]->Fill(H_dc_1v2_rawtdc[0]);
    }
    else if(i == 6){
      if(H_dc_2u1_nhit[0] == 1) h1hdcTdcT[i]->Fill(H_dc_2u1_rawtdc[0]);
    }
    else if(i == 7){
      if(H_dc_2u2_nhit[0] == 1) h1hdcTdcT[i]->Fill(H_dc_2u2_rawtdc[0]);
    }
    else if(i == 8){
      if(H_dc_2x1_nhit[0] == 1) h1hdcTdcT[i]->Fill(H_dc_2x1_rawtdc[0]);
    }
    else if(i == 9){
      if(H_dc_2x2_nhit[0] == 1) h1hdcTdcT[i]->Fill(H_dc_2x2_rawtdc[0]);
    }
    else if(i == 10){
      if(H_dc_2v1_nhit[0] == 1) h1hdcTdcT[i]->Fill(H_dc_2v1_rawtdc[0]);
    }
    else if(i == 11){
      if(H_dc_2v2_nhit[0] == 1) h1hdcTdcT[i]->Fill(H_dc_2v2_rawtdc[0]);
    }
  }
  
  for (Int_t npl = 0; npl < cal_planes; npl++){ // Loop over all calorimeter planes
    for (Int_t nside = 0; nside < sides; nside++){ //Loop over each side
      for (Int_t ipmt = 0; ipmt < 13; ipmt++){ // Loop over each PMT in a particular plane	
	if (npl == 0 && nside == 0) {
	  if (H_cal_1pr_goodPosAdcMult[ipmt] == 1) h1hCalAdcTdcTDiff[npl][nside][ipmt]->Fill(H_cal_1pr_goodPosAdcTdcDiffTime[ipmt]);
	}
	else if (npl == 1 && nside == 0){
	  if (H_cal_2ta_goodPosAdcMult[ipmt] == 1) h1hCalAdcTdcTDiff[npl][nside][ipmt]->Fill(H_cal_2ta_goodPosAdcTdcDiffTime[ipmt]);
	}
	else if (npl == 2 && nside == 0){
	  if (H_cal_3ta_goodPosAdcMult[ipmt] == 1) h1hCalAdcTdcTDiff[npl][nside][ipmt]->Fill(H_cal_3ta_goodPosAdcTdcDiffTime[ipmt]);
	}  
	else if (npl == 1 && nside == 0){
	  if (H_cal_4ta_goodPosAdcMult[ipmt] == 1) h1hCalAdcTdcTDiff[npl][nside][ipmt]->Fill(H_cal_4ta_goodPosAdcTdcDiffTime[ipmt]);
	}
	else if (npl == 0 && nside == 1){
	  if (H_cal_1pr_goodNegAdcMult[ipmt] == 1) h1hCalAdcTdcTDiff[npl][nside][ipmt]->Fill(H_cal_1pr_goodNegAdcTdcDiffTime[ipmt]);
	}
	else if (npl == 1 && nside == 1){
	  if (H_cal_2ta_goodNegAdcMult[ipmt] == 1) h1hCalAdcTdcTDiff[npl][nside][ipmt]->Fill(H_cal_2ta_goodNegAdcTdcDiffTime[ipmt]);
	}
	else if (npl == 2 && nside == 1){
	  if (H_cal_3ta_goodNegAdcMult[ipmt] == 1) h1hCalAdcTdcTDiff[npl][nside][ipmt]->Fill(H_cal_3ta_goodNegAdcTdcDiffTime[ipmt]);
	}  
	else if (npl == 1 && nside == 1){
	  if (H_cal_4ta_goodNegAdcMult[ipmt] == 1) h1hCalAdcTdcTDiff[npl][nside][ipmt]->Fill(H_cal_4ta_goodNegAdcTdcDiffTime[ipmt]);
	}
      }
    }
  }

  for (Int_t npl = 0; npl < hod_planes; npl++){ // Loop over all hodoscope planes
    for (Int_t nside = 0; nside < sides; nside++){ //Loop over each side
      for (Int_t ipmt = 0; ipmt < hmaxPMT[npl]; ipmt++){ // Loop over each PMT in a particular plane	
	if (npl == 0 && nside == 0){
	  if (H_hod_1x_GoodPosAdcMult[ipmt] == 1) h1hHodoAdcTdcTDiff[npl][nside][ipmt]->Fill(H_hod_1x_GoodPosAdcTdcDiffTime[ipmt]);
	}
	else if (npl == 1 && nside == 0){
	  if (H_hod_1y_GoodPosAdcMult[ipmt] == 1) h1hHodoAdcTdcTDiff[npl][nside][ipmt]->Fill(H_hod_1y_GoodPosAdcTdcDiffTime[ipmt]);
	}
	else if (npl == 2 && nside == 0){
	  if (H_hod_2x_GoodPosAdcMult[ipmt] == 1) h1hHodoAdcTdcTDiff[npl][nside][ipmt]->Fill(H_hod_2x_GoodPosAdcTdcDiffTime[ipmt]);
	}  
	else if (npl == 1 && nside == 0){
	  if (H_hod_2y_GoodPosAdcMult[ipmt] == 1) h1hHodoAdcTdcTDiff[npl][nside][ipmt]->Fill(H_hod_2y_GoodPosAdcTdcDiffTime[ipmt]);
	}
	else if (npl == 0 && nside == 1){
	  if (H_hod_1x_GoodNegAdcMult[ipmt] == 1) h1hHodoAdcTdcTDiff[npl][nside][ipmt]->Fill(H_hod_1x_GoodNegAdcTdcDiffTime[ipmt]);
	}
	else if (npl == 1 && nside == 1){
	  if (H_hod_1y_GoodNegAdcMult[ipmt] == 1) h1hHodoAdcTdcTDiff[npl][nside][ipmt]->Fill(H_hod_1y_GoodNegAdcTdcDiffTime[ipmt]);
	}
	else if (npl == 2 && nside == 1){
	  if (H_hod_2x_GoodNegAdcMult[ipmt] == 1) h1hHodoAdcTdcTDiff[npl][nside][ipmt]->Fill(H_hod_2x_GoodNegAdcTdcDiffTime[ipmt]);
	}  
	else if (npl == 1 && nside == 1){
	  if (H_hod_2y_GoodNegAdcMult[ipmt] == 1) h1hHodoAdcTdcTDiff[npl][nside][ipmt]->Fill(H_hod_2y_GoodNegAdcTdcDiffTime[ipmt]);
	} 	
      }
    }
  }
  
  // Fill our SHMS timing histograms, explicitly select only multiplicity 1 events
  for (Int_t ipmt = 0; ipmt < 4; ipmt++){
    if(P_hgcer_goodAdcMult[ipmt] == 1) h1pHGCAdcTdcTDiff[ipmt]->Fill(P_hgcer_goodAdcTdcDiffTime[ipmt]);
  }
  
  for (Int_t nside = 0; nside < sides; nside++){
    for (Int_t ipmt = 0; ipmt < 7; ipmt++){
      if(nside == 0){
	 if(P_aero_goodPosAdcMult[ipmt] == 1) h1pAeroAdcTdcTDiff[nside][ipmt]->Fill(P_aero_goodPosAdcTdcDiffTime[ipmt]);
      }
      else if(nside == 1){
	 if(P_aero_goodNegAdcMult[ipmt] == 1) h1pAeroAdcTdcTDiff[nside][ipmt]->Fill(P_aero_goodNegAdcTdcDiffTime[ipmt]);  
      }
    }
  }

  for (Int_t i = 0; i < 12; i++){
    if(i == 0){if(P_dc_1u1_nhit[0] == 1) h1pdcTdcT[i]->Fill(P_dc_1u1_rawtdc[0]);}
    else if(i == 1){
      if(P_dc_1u2_nhit[0] == 1) h1pdcTdcT[i]->Fill(P_dc_1u2_rawtdc[0]);
    }
    else if(i == 2){
      if(P_dc_1x1_nhit[0] == 1) h1pdcTdcT[i]->Fill(P_dc_1x1_rawtdc[0]);
    }
    else if(i == 3){
      if(P_dc_1x2_nhit[0] == 1) h1pdcTdcT[i]->Fill(P_dc_1x2_rawtdc[0]);
    }
    else if(i == 4){
      if(P_dc_1v1_nhit[0] == 1) h1pdcTdcT[i]->Fill(P_dc_1v1_rawtdc[0]);
    }
    else if(i == 5){
      if(P_dc_1v2_nhit[0] == 1) h1pdcTdcT[i]->Fill(P_dc_1v2_rawtdc[0]);
    }
    else if(i == 6){
      if(P_dc_2u1_nhit[0] == 1) h1pdcTdcT[i]->Fill(P_dc_2u1_rawtdc[0]);
    }
    else if(i == 7){
      if(P_dc_2u2_nhit[0] == 1) h1pdcTdcT[i]->Fill(P_dc_2u2_rawtdc[0]);
    }
    else if(i == 8){
      if(P_dc_2x1_nhit[0] == 1) h1pdcTdcT[i]->Fill(P_dc_2x1_rawtdc[0]);
    }
    else if(i == 9){
      if(P_dc_2x2_nhit[0] == 1) h1pdcTdcT[i]->Fill(P_dc_2x2_rawtdc[0]);
    }
    else if(i == 10){
      if(P_dc_2v1_nhit[0] == 1) h1pdcTdcT[i]->Fill(P_dc_2v1_rawtdc[0]);
    }
    else if(i == 11){
      if(P_dc_2v2_nhit[0] == 1) h1pdcTdcT[i]->Fill(P_dc_2v2_rawtdc[0]);
    }
  }

  for (Int_t npl = 0; npl < hod_planes; npl++){ // Loop over all hodoscope planes
    for (Int_t nside = 0; nside < sides; nside++){ //Loop over each side
      for (Int_t ipmt = 0; ipmt < pmaxPMT[npl]; ipmt++){ // Loop over each PMT in a particular plane	
	if (npl == 0 && nside == 0){
	  if (P_hod_1x_GoodPosAdcMult[ipmt] == 1) h1pHodoAdcTdcTDiff[npl][nside][ipmt]->Fill(P_hod_1x_GoodPosAdcTdcDiffTime[ipmt]);
	}
	else if (npl == 1 && nside == 0){
	  if (P_hod_1y_GoodPosAdcMult[ipmt] == 1) h1pHodoAdcTdcTDiff[npl][nside][ipmt]->Fill(P_hod_1y_GoodPosAdcTdcDiffTime[ipmt]);
	}
	else if (npl == 2 && nside == 0){
	  if (P_hod_2x_GoodPosAdcMult[ipmt] == 1) h1pHodoAdcTdcTDiff[npl][nside][ipmt]->Fill(P_hod_2x_GoodPosAdcTdcDiffTime[ipmt]);
	}  
	else if (npl == 1 && nside == 0){
	  if (P_hod_2y_GoodPosAdcMult[ipmt] == 1) h1pHodoAdcTdcTDiff[npl][nside][ipmt]->Fill(P_hod_2y_GoodPosAdcTdcDiffTime[ipmt]);
	}
	else if (npl == 0 && nside == 1){
	  if (P_hod_1x_GoodNegAdcMult[ipmt] == 1) h1pHodoAdcTdcTDiff[npl][nside][ipmt]->Fill(P_hod_1x_GoodNegAdcTdcDiffTime[ipmt]);
	}
	else if (npl == 1 && nside == 1){
	  if (P_hod_1y_GoodNegAdcMult[ipmt] == 1) h1pHodoAdcTdcTDiff[npl][nside][ipmt]->Fill(P_hod_1y_GoodNegAdcTdcDiffTime[ipmt]);
	}
	else if (npl == 2 && nside == 1){
	  if (P_hod_2x_GoodNegAdcMult[ipmt] == 1) h1pHodoAdcTdcTDiff[npl][nside][ipmt]->Fill(P_hod_2x_GoodNegAdcTdcDiffTime[ipmt]);
	}  
	else if (npl == 1 && nside == 1){
	  if (P_hod_2y_GoodNegAdcMult[ipmt] == 1) h1pHodoAdcTdcTDiff[npl][nside][ipmt]->Fill(P_hod_2y_GoodNegAdcTdcDiffTime[ipmt]);
	} 	
      }
    }
  }  

  for (Int_t nside = 0; nside < sides; nside++){
    for (Int_t ipmt = 0; ipmt < 14; ipmt++){
      if(nside == 0){
	 if(P_cal_pr_goodPosAdcMult[ipmt] == 1) h1pPrShAdcTdcTDiff[nside][ipmt]->Fill(P_cal_pr_goodPosAdcTdcDiffTime[ipmt]);
      }
      else if(nside == 1){
	 if(P_cal_pr_goodNegAdcMult[ipmt] == 1) h1pPrShAdcTdcTDiff[nside][ipmt]->Fill(P_cal_pr_goodNegAdcTdcDiffTime[ipmt]);  
      }
    }
  }

  for (Int_t ipmt = 0; ipmt < 224; ipmt++){
    if(P_cal_fly_goodAdcMult[ipmt] == 1) h1pCalAdcTdcTDiff[ipmt]->Fill(P_cal_fly_goodAdcTdcDiffTime[ipmt]);
  }

  return kTRUE;

}

void DetTCuts_Coin::SlaveTerminate()
{
}

void DetTCuts_Coin::Terminate()
{
  cout << "Finished processing" << endl;
  printf("\n");
  TString option = GetOption();

  TFile *Histogram_file = new TFile(Form("TimeWindowHistos_Run%i.root",option.Atoi()),"RECREATE");
  TString outputpdf = Form("TimeWindowPlots_Run%i.pdf", option.Atoi()) ; 

  TDirectory *DHMSCER = Histogram_file->mkdir("HMS Cherenkov Timing"); DHMSCER->cd();
  TCanvas *CHMSCER = new TCanvas("CHMSCER", "HMS Cherenkov timing plots", 300,100,1000,900);
  CHMSCER->Divide(2,1);
  for (Int_t ipmt = 0; ipmt < 2; ipmt++){
    TH1F *HMSCER = dynamic_cast<TH1F *>(TProof::GetOutput(Form("hCER%d_timeDiff", ipmt+1), fOutput));
    HMSCER->Write();
    CHMSCER->cd(ipmt+1); HMSCER->Draw();
  }

  TDirectory *DHMSDC = Histogram_file->mkdir("HMS DC Timing"); DHMSDC->cd();
  TCanvas *CHMSDC = new TCanvas("CHMSDC", "HMS DC timing plots", 300,100,1000,900);
  CHMSDC->Divide(4, 3);
  for (Int_t i = 0; i < 12; i++){
    TH1F *HMSDC = dynamic_cast<TH1F *>(TProof::GetOutput(Form("hDC%s_rawTDC", dc_pl_names[i].c_str()), fOutput));
    HMSDC->Write();
    CHMSDC->cd(i+1); HMSDC->Draw();
  }

  TDirectory *DHMSCAL = Histogram_file->mkdir("HMS Calorimeter Timing"); DHMSCAL->cd();  
  TCanvas *CHMSCAL[4][2];
  for (Int_t npl = 0; npl < cal_planes; npl++){ // Loop over all calorimeter planes
    for (Int_t nside = 0; nside < sides; nside++){ //Loop over each side
      CHMSCAL[npl][nside] = new TCanvas(Form("CHMSCAL%s%s", cal_pl_names[npl].c_str(), nsign[nside].c_str()),  Form("HMS Calorimeter %s%s Timing", cal_pl_names[npl].c_str(), nsign[nside].c_str()), 300,100,1000,900);
      CHMSCAL[npl][nside]->Divide(5, 3);
      for (Int_t ipmt = 0; ipmt < 13; ipmt++){ // Loop over each PMT in a particular plane	
	TH1F *HMSCAL = dynamic_cast<TH1F *>(TProof::GetOutput(Form("hCal%s%d%s_timeDiff", cal_pl_names[npl].c_str(), ipmt+1, nsign[nside].c_str()), fOutput));
	HMSCAL->Write();
	CHMSCAL[npl][nside]->cd(ipmt+1); HMSCAL->Draw();
      }
    }
  }

  TDirectory *DHMSHODO = Histogram_file->mkdir("HMS Hodoscope Timing"); DHMSHODO->cd();  
  TCanvas *CHMSHODO[4][2];  
  for (Int_t npl = 0; npl < hod_planes; npl++){ // Loop over all hodoscope planes
    for (Int_t nside = 0; nside < sides; nside++){ //Loop over each side
      CHMSHODO[npl][nside] = new TCanvas(Form("CHMSHODO%s%s", hod_pl_names[npl].c_str(), nsign[nside].c_str()),  Form("HMS Hodoscope %s%s Timing", hod_pl_names[npl].c_str(), nsign[nside].c_str()), 300,100,1000,900);
      CHMSHODO[npl][nside]->Divide(5, 4);
      for (Int_t ipmt = 0; ipmt < hmaxPMT[npl]; ipmt++){ // Loop over each PMT in a particular plane	
	TH1F *HMSHODO = dynamic_cast<TH1F *>(TProof::GetOutput(Form("hHodo%s%d%s_timeDiff", hod_pl_names[npl].c_str(),ipmt+1,nsign[nside].c_str() ), fOutput));
	HMSHODO->Write();
	CHMSHODO[npl][nside]->cd(ipmt+1); HMSHODO->Draw();
      }
    }
  }

  TDirectory *DSHMSHGC = Histogram_file->mkdir("SHMS HGC Timing"); DSHMSHGC->cd();  
  TCanvas *CSHMSHGC = new TCanvas("CSHMSHGC", "SHMS HGC timing plots", 300,100,1000,900);
  CSHMSHGC->Divide(2,2);
  for (Int_t ipmt = 0; ipmt < 4; ipmt++){
    TH1F *SHMSHGC = dynamic_cast<TH1F *>(TProof::GetOutput(Form("pHGCER%d_timeDiff", ipmt+1), fOutput));
    SHMSHGC->Write();
    CSHMSHGC->cd(ipmt+1); SHMSHGC->Draw();
  }

  TDirectory *DSHMSAERO = Histogram_file->mkdir("SHMS Aerogel Cherenkov Timing"); DSHMSAERO->cd();  
  TCanvas *CSHMSAERO[2];
  for (Int_t nside = 0; nside < sides; nside++){ //Loop over each side
      CSHMSAERO[nside] = new TCanvas(Form("CSHMSAERO%s", nsign[nside].c_str()),  Form("SHMS Aerogel  %sPMT Timing", nsign[nside].c_str()), 300,100,1000,900);
      CSHMSAERO[nside]->Divide(2, 4);
    for (Int_t ipmt = 0; ipmt < 7; ipmt++){ // Loop over PMTs
      TH1F *SHMSAERO = dynamic_cast<TH1F *>(TProof::GetOutput(Form("pAero%d%s_timeDiff", ipmt+1, nsign[nside].c_str()), fOutput));
      SHMSAERO->Write();
      CSHMSAERO[nside]->cd(ipmt+1); SHMSAERO->Draw();
    }
  }

  
  TDirectory *DSHMSDC = Histogram_file->mkdir("SHMS DC Timing"); DSHMSDC->cd();  
  TCanvas *CSHMSDC = new TCanvas("CSHMSDC", "SHMS DC timing plots", 300,100,1000,900);
  CSHMSDC->Divide(4, 3);
  for (Int_t i = 0; i < 12; i++){
    TH1F *SHMSDC = dynamic_cast<TH1F *>(TProof::GetOutput(Form("pDC%s_rawTDC", dc_pl_names[i].c_str()), fOutput));
    SHMSDC->Write();
    CSHMSDC->cd(i+1); SHMSDC->Draw();
  }

  TDirectory *DSHMSHODO = Histogram_file->mkdir("SHMS Hodoscope Timing"); DSHMSHODO->cd();  
  TCanvas *CSHMSHODO[4][2];  
  for (Int_t npl = 0; npl < hod_planes; npl++){ // Loop over all hodoscope planes
    for (Int_t nside = 0; nside < sides; nside++){ //Loop over each side
      CSHMSHODO[npl][nside] = new TCanvas(Form("CSHMSHODO%s%s", hod_pl_names[npl].c_str(), nsign[nside].c_str()),  Form("SHMS Hodoscope %s%s Timing", hod_pl_names[npl].c_str(), nsign[nside].c_str()), 300,100,1000,900);
      CSHMSHODO[npl][nside]->Divide(5, 5);
      for (Int_t ipmt = 0; ipmt < pmaxPMT[npl]; ipmt++){ // Loop over each PMT in a particular plane	
	TH1F *SHMSHODO = dynamic_cast<TH1F *>(TProof::GetOutput(Form("pHodo%s%d%s_timeDiff", hod_pl_names[npl].c_str(),ipmt+1,nsign[nside].c_str() ), fOutput));
	SHMSHODO->Write();
	CSHMSHODO[npl][nside]->cd(ipmt+1); SHMSHODO->Draw();
      }
    }
  }

  TDirectory *DSHMSPRSH = Histogram_file->mkdir("SHMS Pre-Shower Timing"); DSHMSPRSH->cd();  
  TCanvas *CSHMSPRSH[2];
  for (Int_t nside = 0; nside < sides; nside++){ //Loop over each side
      CSHMSPRSH[nside] = new TCanvas(Form("CSHMSPRSH%s", nsign[nside].c_str()),  Form("SHMS Pre-Shower  %sPMT Timing", nsign[nside].c_str()), 300,100,1000,900);
      CSHMSPRSH[nside]->Divide(5, 3);    
    for (Int_t ipmt = 0; ipmt < 14; ipmt++){ // Loop over PMTs
      TH1F *SHMSPRSH = dynamic_cast<TH1F *>(TProof::GetOutput(Form("pPrSh%d%s_timeDiff", ipmt+1, nsign[nside].c_str()), fOutput));
      SHMSPRSH->Write();
      CSHMSPRSH[nside]->cd(ipmt+1); SHMSPRSH->Draw();
    }
  }
  
  TDirectory *DSHMSCAL = Histogram_file->mkdir("SHMS Calorimeter Timing"); DSHMSCAL->cd();  
  TCanvas *CSHMSCAL[14]; // 16 histograms per canvas
  for(Int_t row = 0; row < 14; row++){
    CSHMSCAL[row] = new TCanvas(Form("CSHMSCAL%d", row+1),  Form("SHMS Pre-Shower Row %d", row+1), 300,100,1000,900);
     CSHMSCAL[row]->Divide(4, 4);     
     for(Int_t ipmt = 0; ipmt < 16; ipmt++){
      TH1F *SHMSCAL = dynamic_cast<TH1F *>(TProof::GetOutput(Form("pCalPMT%d", (row*16)+ipmt+1), fOutput)); 
      SHMSCAL->Write();
      CSHMSCAL[row]->cd(ipmt+1); SHMSCAL->Draw();
     }
  }

  CHMSCER->Print(outputpdf+"[");
  CHMSCER->Print(outputpdf);
  CHMSDC->Print(outputpdf);
  for (Int_t npl = 0; npl < cal_planes; npl++){ // Loop over all calorimeter planes
    for (Int_t nside = 0; nside < sides; nside++){ //Loop over each side
      CHMSCAL[npl][nside]->Print(outputpdf);
    }
  }
  for (Int_t npl = 0; npl < hod_planes; npl++){ // Loop over all hodoscope planes
    for (Int_t nside = 0; nside < sides; nside++){ //Loop over each side
      CHMSHODO[npl][nside]->Print(outputpdf);
    }
  }
  CSHMSHGC->Print(outputpdf);
  CSHMSAERO[0]->Print(outputpdf);
  CSHMSAERO[1]->Print(outputpdf);
  CSHMSDC->Print(outputpdf);
  for (Int_t npl = 0; npl < hod_planes; npl++){ // Loop over all hodoscope planes
    for (Int_t nside = 0; nside < sides; nside++){ //Loop over each side
      CSHMSHODO[npl][nside]->Print(outputpdf);
    }
  }
  CSHMSPRSH[0]->Print(outputpdf);
  CSHMSPRSH[1]->Print(outputpdf);
  for(Int_t row = 0; row < 14; row++){
    CSHMSCAL[row]->Print(outputpdf);
  }
  CHMSCER->Print(outputpdf+"]");
  
  Histogram_file->Close();

}
