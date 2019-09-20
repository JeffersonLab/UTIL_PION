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
      if(npl == 2 && nside == 1) continue; // Skip 3ta/4ta- since they don't exist!
      if(npl == 3 && nside == 1) continue;
      for (Int_t ipmt = 0; ipmt < 13; ipmt++){ // Loop over each PMT in a particular plane	
	h1hCalAdcTdcTDiff[npl][nside][ipmt] = new TH1F(Form("hCal%s%d%s_timeDiff", cal_pl_names[npl].c_str(), ipmt+1, nsign[nside].c_str()), Form("HMS Cal %s%d%s AdcTdcTimeDiff", cal_pl_names[npl].c_str(), ipmt+1, nsign[nside].c_str()), 200, -100, 100);
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
	else if (npl == 3 && nside == 0){
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
	else if (npl == 3 && nside == 1){
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
    HMSDC_tMin[i] = (HMSDC->GetMean() - (5*HMSDC->GetStdDev()));
    HMSDC_tMax[i] = (HMSDC->GetMean() + (5*HMSDC->GetStdDev()));
    LHMSDC_tMin[i] = new TLine(HMSDC_tMin[i], 0, HMSDC_tMin[i], HMSDC->GetMaximum());
    LHMSDC_tMax[i] = new TLine(HMSDC_tMax[i], 0, HMSDC_tMax[i], HMSDC->GetMaximum());
    LHMSDC_tMin[i]->SetLineColor(kRed); LHMSDC_tMin[i]->SetLineStyle(7); LHMSDC_tMin[i]->SetLineWidth(1);
    LHMSDC_tMax[i]->SetLineColor(kRed); LHMSDC_tMax[i]->SetLineStyle(7); LHMSDC_tMax[i]->SetLineWidth(1);
    CHMSDC->cd(i+1); HMSDC->Draw(); LHMSDC_tMin[i]->Draw("SAME"); LHMSDC_tMax[i]->Draw("SAME");
  }

  TDirectory *DHMSCAL = Histogram_file->mkdir("HMS Calorimeter Timing"); DHMSCAL->cd();  
  TCanvas *CHMSCAL[4][2];
  for (Int_t npl = 0; npl < cal_planes; npl++){ // Loop over all calorimeter planes
    for (Int_t nside = 0; nside < sides; nside++){ //Loop over each side
      if(npl == 2 && nside == 1) continue; // Skip 3ta/4ta- since they don't exist!
      if(npl == 3 && nside == 1) continue;
      CHMSCAL[npl][nside] = new TCanvas(Form("CHMSCAL%s%s", cal_pl_names[npl].c_str(), nsign[nside].c_str()),  Form("HMS Calorimeter %s%s Timing", cal_pl_names[npl].c_str(), nsign[nside].c_str()), 300,100,1000,900);
      CHMSCAL[npl][nside]->Divide(5, 3);
      for (Int_t ipmt = 0; ipmt < 13; ipmt++){ // Loop over each PMT in a particular plane	
	TH1F *HMSCAL = dynamic_cast<TH1F *>(TProof::GetOutput(Form("hCal%s%d%s_timeDiff", cal_pl_names[npl].c_str(), ipmt+1, nsign[nside].c_str()), fOutput));
	HMSCAL->Write();
	HMSCAL_tMin[npl][nside][ipmt] = (HMSCAL->GetMean() - (5*HMSCAL->GetStdDev()));
	HMSCAL_tMax[npl][nside][ipmt] = (HMSCAL->GetMean() + (5*HMSCAL->GetStdDev()));
	LHMSCAL_tMin[npl][nside][ipmt] = new TLine(HMSCAL_tMin[npl][nside][ipmt], 0, HMSCAL_tMin[npl][nside][ipmt], HMSCAL->GetMaximum());
	LHMSCAL_tMax[npl][nside][ipmt] = new TLine(HMSCAL_tMax[npl][nside][ipmt], 0, HMSCAL_tMax[npl][nside][ipmt], HMSCAL->GetMaximum());
	LHMSCAL_tMin[npl][nside][ipmt]->SetLineColor(kRed); LHMSCAL_tMin[npl][nside][ipmt]->SetLineStyle(7); LHMSCAL_tMin[npl][nside][ipmt]->SetLineWidth(1);
	LHMSCAL_tMax[npl][nside][ipmt]->SetLineColor(kRed); LHMSCAL_tMax[npl][nside][ipmt]->SetLineStyle(7); LHMSCAL_tMax[npl][nside][ipmt]->SetLineWidth(1);
	CHMSCAL[npl][nside]->cd(ipmt+1); HMSCAL->Draw(); LHMSCAL_tMin[npl][nside][ipmt]->Draw("SAME"); LHMSCAL_tMax[npl][nside][ipmt]->Draw("SAME");
      }
    }
  }

  TDirectory *DHMSHODO = Histogram_file->mkdir("HMS Hodoscope Timing"); DHMSHODO->cd();  
  TCanvas *CHMSHODO[4][2];  
  for (Int_t npl = 0; npl < hod_planes; npl++){ // Loop over all hodoscope planes
    for (Int_t nside = 0; nside < sides; nside++){ //Loop over each side
      CHMSHODO[npl][nside] = new TCanvas(Form("CHMSHODO%s%s", hod_pl_names[npl].c_str(), nsign[nside].c_str()),  Form("HMS Hodoscope %s%s Timing", hod_pl_names[npl].c_str(), nsign[nside].c_str()), 300,100,1000,900);
      CHMSHODO[npl][nside]->Divide(4, 4);
      for (Int_t ipmt = 0; ipmt < hmaxPMT[npl]; ipmt++){ // Loop over each PMT in a particular plane	
	TH1F *HMSHODO = dynamic_cast<TH1F *>(TProof::GetOutput(Form("hHodo%s%d%s_timeDiff", hod_pl_names[npl].c_str(),ipmt+1,nsign[nside].c_str() ), fOutput));
	HMSHODO->Write();
	HMSHODO_tMin[npl][nside][ipmt] = (HMSHODO->GetMean() - (5*HMSHODO->GetStdDev()));
	HMSHODO_tMax[npl][nside][ipmt] = (HMSHODO->GetMean() + (5*HMSHODO->GetStdDev()));
	LHMSHODO_tMin[npl][nside][ipmt] = new TLine(HMSHODO_tMin[npl][nside][ipmt], 0, HMSHODO_tMin[npl][nside][ipmt], HMSHODO->GetMaximum());
	LHMSHODO_tMax[npl][nside][ipmt] = new TLine(HMSHODO_tMax[npl][nside][ipmt], 0, HMSHODO_tMax[npl][nside][ipmt], HMSHODO->GetMaximum());
	LHMSHODO_tMin[npl][nside][ipmt]->SetLineColor(kRed); LHMSHODO_tMin[npl][nside][ipmt]->SetLineStyle(7); LHMSHODO_tMin[npl][nside][ipmt]->SetLineWidth(1);
	LHMSHODO_tMax[npl][nside][ipmt]->SetLineColor(kRed); LHMSHODO_tMax[npl][nside][ipmt]->SetLineStyle(7); LHMSHODO_tMax[npl][nside][ipmt]->SetLineWidth(1);
	CHMSHODO[npl][nside]->cd(ipmt+1); HMSHODO->Draw(); LHMSHODO_tMin[npl][nside][ipmt]->Draw("SAME"); LHMSHODO_tMax[npl][nside][ipmt]->Draw("SAME");
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
      SHMSAERO_tMin[nside][ipmt] = (SHMSAERO->GetMean() - (5*SHMSAERO->GetStdDev()));
      SHMSAERO_tMax[nside][ipmt] = (SHMSAERO->GetMean() + (5*SHMSAERO->GetStdDev()));
      if(SHMSAERO_tMin[nside][ipmt] < 0) SHMSAERO_tMin[nside][ipmt] = 0;
      LSHMSAERO_tMin[nside][ipmt] = new TLine(SHMSAERO_tMin[nside][ipmt], 0, SHMSAERO_tMin[nside][ipmt], SHMSAERO->GetMaximum());
      LSHMSAERO_tMax[nside][ipmt] = new TLine(SHMSAERO_tMax[nside][ipmt], 0, SHMSAERO_tMax[nside][ipmt], SHMSAERO->GetMaximum());
      LSHMSAERO_tMin[nside][ipmt]->SetLineColor(kRed); LSHMSAERO_tMin[nside][ipmt]->SetLineStyle(7); LSHMSAERO_tMin[nside][ipmt]->SetLineWidth(1);
      LSHMSAERO_tMax[nside][ipmt]->SetLineColor(kRed); LSHMSAERO_tMax[nside][ipmt]->SetLineStyle(7); LSHMSAERO_tMax[nside][ipmt]->SetLineWidth(1);
      CSHMSAERO[nside]->cd(ipmt+1); SHMSAERO->Draw(); LSHMSAERO_tMin[nside][ipmt]->Draw("SAME"); LSHMSAERO_tMax[nside][ipmt]->Draw("SAME");
    }
  }
  
  TDirectory *DSHMSDC = Histogram_file->mkdir("SHMS DC Timing"); DSHMSDC->cd();  
  TCanvas *CSHMSDC = new TCanvas("CSHMSDC", "SHMS DC timing plots", 300,100,1000,900);
  CSHMSDC->Divide(4, 3);
  for (Int_t i = 0; i < 12; i++){
    TH1F *SHMSDC = dynamic_cast<TH1F *>(TProof::GetOutput(Form("pDC%s_rawTDC", dc_pl_names[i].c_str()), fOutput));
    SHMSDC->Write();
    SHMSDC_tMin[i] = (SHMSDC->GetMean() - (5*SHMSDC->GetStdDev()));
    SHMSDC_tMax[i] = (SHMSDC->GetMean() + (5*SHMSDC->GetStdDev()));
    LSHMSDC_tMin[i] = new TLine(SHMSDC_tMin[i], 0, SHMSDC_tMin[i], SHMSDC->GetMaximum());
    LSHMSDC_tMax[i] = new TLine(SHMSDC_tMax[i], 0, SHMSDC_tMax[i], SHMSDC->GetMaximum());
    LSHMSDC_tMin[i]->SetLineColor(kRed); LSHMSDC_tMin[i]->SetLineStyle(7); LSHMSDC_tMin[i]->SetLineWidth(1);
    LSHMSDC_tMax[i]->SetLineColor(kRed); LSHMSDC_tMax[i]->SetLineStyle(7); LSHMSDC_tMax[i]->SetLineWidth(1);
    CSHMSDC->cd(i+1); SHMSDC->Draw(); LSHMSDC_tMin[i]->Draw("SAME"); LSHMSDC_tMax[i]->Draw("SAME");
  }

  TDirectory *DSHMSHODO = Histogram_file->mkdir("SHMS Hodoscope Timing"); DSHMSHODO->cd();  
  TCanvas *CSHMSHODO[4][2];  
  for (Int_t npl = 0; npl < hod_planes; npl++){ // Loop over all hodoscope planes
    for (Int_t nside = 0; nside < sides; nside++){ //Loop over each side
      CSHMSHODO[npl][nside] = new TCanvas(Form("CSHMSHODO%s%s", hod_pl_names[npl].c_str(), nsign[nside].c_str()),  Form("SHMS Hodoscope %s%s Timing", hod_pl_names[npl].c_str(), nsign[nside].c_str()), 300,100,1000,900);
      if (npl != 3) CSHMSHODO[npl][nside]->Divide(5, 3);
      else if (npl == 3) CSHMSHODO[npl][nside]->Divide(7, 3);
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
      SHMSCAL_tMin[row][ipmt] = (SHMSCAL->GetMean() - (5*SHMSCAL->GetStdDev()));
      SHMSCAL_tMax[row][ipmt] = (SHMSCAL->GetMean() + (5*SHMSCAL->GetStdDev()));
      LSHMSCAL_tMin[row][ipmt] = new TLine(SHMSCAL_tMin[row][ipmt], 0, SHMSCAL_tMin[row][ipmt], SHMSCAL->GetMaximum());
      LSHMSCAL_tMax[row][ipmt] = new TLine(SHMSCAL_tMax[row][ipmt], 0, SHMSCAL_tMax[row][ipmt], SHMSCAL->GetMaximum());
      LSHMSCAL_tMin[row][ipmt]->SetLineColor(kRed); LSHMSCAL_tMin[row][ipmt]->SetLineStyle(7); LSHMSCAL_tMin[row][ipmt]->SetLineWidth(1);
      LSHMSCAL_tMax[row][ipmt]->SetLineColor(kRed); LSHMSCAL_tMax[row][ipmt]->SetLineStyle(7); LSHMSCAL_tMax[row][ipmt]->SetLineWidth(1);
      CSHMSCAL[row]->cd(ipmt+1); SHMSCAL->Draw(); SHMSCAL->Draw(); LSHMSCAL_tMin[row][ipmt]->Draw("SAME"); LSHMSCAL_tMax[row][ipmt]->Draw("SAME");
     }
  }

  CHMSCER->Print(outputpdf+"[");
  CHMSCER->Print(outputpdf);
  CHMSDC->Print(outputpdf);
  for (Int_t npl = 0; npl < cal_planes; npl++){ // Loop over all calorimeter planes
    for (Int_t nside = 0; nside < sides; nside++){ //Loop over each side
      if(npl == 2 && nside == 1) continue; // Skip 3ta/4ta- since they don't exist!
      if(npl == 3 && nside == 1) continue;
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
