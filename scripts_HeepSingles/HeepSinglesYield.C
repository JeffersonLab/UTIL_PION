#define HeepSinglesYield_cxx

#include "HeepSinglesYield.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TGaxis.h>

void HeepSinglesYield::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  TString option = GetOption();
}

void HeepSinglesYield::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
  TString option = GetOption();

  EventType   = new TH1F("Event_Type" ,"Counts for each Event Type",20,0,10); 

  TRIG1     = new TH1F("TRIG1","pTRIG1 counts",4000,-1000.0,2000.0);
  TRIG3     = new TH1F("TRIG3","pTRIG3 counts",4000,-1000.0,2000.0);
  TRIG5     = new TH1F("TRIG5","pTRIG5 counts",4000,0.0,1000.0);
  TRIG1_cut = new TH1F("TRIG1_cut","pTRIG1 counts",4000,-1000.0,1000.0);
  TRIG3_cut = new TH1F("TRIG3_cut","pTRIG3 counts",4000,0.0,1000.0);
  TRIG5_cut = new TH1F("TRIG5_cut","pTRIG5 counts",4000,0.0,1000.0);

  HMS_delta            = new TH1F("HMS_delta","HMS #delta;#delta;Counts",100,-50,50);
  HMS_delta_cut        = new TH1F("HMS_delta_cut","HMS #delta Cut;#delta;Counts",100,-50,50);
  HMS_th              = new TH1F("HMS_th","HMS #theta Acceptance;#theta;Counts",100,-0.1,0.1);
  HMS_th_cut          = new TH1F("HMS_th_cut","HMS #theta Acceptance with Cut;#theta;Counts",100,-0.1,0.1);
  HMS_ph              = new TH1F("HMS_ph","HMS #phi Acceptance;#phi;Counts",100,-0.1,0.1);
  HMS_ph_cut          = new TH1F("HMS_ph_cut","HMS #phi Acceptance with Cut;#phi;Counts",100,-0.1,0.1);
  h_etotnorm_before = new TH1F("h_etotnorm_before", "HMS Normalised Calorimeter Energy; Normalised Energy; Counts", 150, 0, 1.5);
  h_etotnorm_cut = new TH1F("h_etotnorm_cut", "HMS Normalised Calorimeter Energy, PID cut; Normalised Energy; Counts", 150, 0, 1.5);
  h_cer_before = new TH1F("h_cer_before", "HMS Cherenkov NPE; NPE; Counts", 50, 0, 25);
  h_cer_cut = new TH1F("h_cer_cut", "HMS Cherenkov NPE; NPE, PID cut; Counts", 50, 0, 25);
  h_cal_cer_before = new TH2F("h_cal_cer_before", "HMS CerNPE(ETotNorm); Normalised Energy; NPE", 150, 0, 1.5, 50, 0, 25);
  h_cal_cer_cut = new TH2F("h_cal_cer_cut", "HMS CerNPE(ETotNorm), PID cut; Normalised Energy; NPE", 150, 0, 1.5, 50, 0, 25);
  HMS_W_Dist = new TH1F("HMS_W_Dist", "HMS Invariant Mass, W, Distribution; W/GeV; Counts", 200, 0, 2);
  HMS_ph_q_Dist = new TH1F("HMS_ph_q_Dist", "HMS #phi_{q} Distribution; #phi_{q}; counts", 200, -0.5, 0.5);

  SHMS_delta           = new TH1F("SHMS_delta","SHMS #delta;#delta;Counts",100,-50,50);
  SHMS_delta_cut       = new TH1F("SHMS_delta_cut","SHMS #delta Cut;#delta;Counts",100,-50,50);
  SHMS_th              = new TH1F("SHMS_th","SHMS #theta Acceptance;#theta;Counts",100,-0.1,0.1);
  SHMS_th_cut          = new TH1F("SHMS_th_cut","SHMS #theta Acceptance with Cut;#theta;Counts",100,-0.1,0.1);
  SHMS_ph              = new TH1F("SHMS_ph","SHMS #phi Acceptance;#phi;Counts",100,-0.1,0.1);
  SHMS_ph_cut          = new TH1F("SHMS_ph_cut","SHMS #phi Acceptance with Cut;#phi;Counts",100,-0.1,0.1);
  p_etotnorm_before = new TH1F("p_etotnorm_before", "SHMS Normalised Calorimeter Energy; Normalised Energy; Counts", 150, 0, 1.5);
  p_etotnorm_cut = new TH1F("p_etotnorm_cut", "SHMS Normalised Calorimeter Energy, PID cut; Normalised Energy; Counts", 150, 0, 1.5);
  
  GetOutputList()->Add(EventType);
  GetOutputList()->Add(TRIG1);
  GetOutputList()->Add(TRIG3);
  GetOutputList()->Add(TRIG5);
  GetOutputList()->Add(TRIG1_cut);
  GetOutputList()->Add(TRIG3_cut);
  GetOutputList()->Add(TRIG5_cut);

  GetOutputList()->Add(HMS_delta);
  GetOutputList()->Add(HMS_delta_cut);
  GetOutputList()->Add(HMS_th);
  GetOutputList()->Add(HMS_th_cut);
  GetOutputList()->Add(HMS_ph);
  GetOutputList()->Add(HMS_ph_cut);
  GetOutputList()->Add(h_etotnorm_before);
  GetOutputList()->Add(h_etotnorm_cut);
  GetOutputList()->Add(h_cer_before);
  GetOutputList()->Add(h_cer_cut);
  GetOutputList()->Add(h_cal_cer_before);
  GetOutputList()->Add(h_cal_cer_cut);
  GetOutputList()->Add(HMS_W_Dist);
  GetOutputList()->Add(HMS_ph_q_Dist);
  GetOutputList()->Add(SHMS_delta);
  GetOutputList()->Add(SHMS_delta_cut);
  GetOutputList()->Add(SHMS_th);
  GetOutputList()->Add(SHMS_th_cut);
  GetOutputList()->Add(SHMS_ph);
  GetOutputList()->Add(SHMS_ph_cut);
  GetOutputList()->Add(p_etotnorm_before);
  GetOutputList()->Add(p_etotnorm_cut);
}

Bool_t HeepSinglesYield::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  fReader.SetEntry(entry);
  
  EventType->Fill(*EvtType); 
  if (*T_coin_pTRIG1_ROC2_tdcTime!=0.0) TRIG1->Fill(*T_coin_pTRIG1_ROC2_tdcTime); // Why only ROC2?
  if (*T_coin_pTRIG3_ROC2_tdcTime!=0.0) TRIG3->Fill(*T_coin_pTRIG3_ROC2_tdcTime);
  if (*T_coin_pTRIG5_ROC2_tdcTime!=0.0) TRIG5->Fill(*T_coin_pTRIG5_ROC2_tdcTime);

  if(*EvtType==1){ // SHMS single
    TRIG1_cut->Fill(*T_coin_pTRIG1_ROC2_tdcTime);
    p_etotnorm_before->Fill(P_cal_etotnorm[0]);
    SHMS_delta->Fill(P_gtr_dp[0]);
    SHMS_th->Fill(P_gtr_th[0]);
    SHMS_ph->Fill(P_gtr_ph[0]);
    
    // Check we actually have a good track
    if(P_hod_goodscinhit[0] != 1) return kTRUE;
    if(P_hod_betanotrack[0] < 0.5 || P_hod_betanotrack[0] > 1.5) return kTRUE;
    if((P_dc_1x1_nhit[0] + P_dc_1u2_nhit[0] + P_dc_1u1_nhit[0] + P_dc_1v1_nhit[0] + P_dc_1x2_nhit[0] + P_dc_1v2_nhit[0]) > 20) return kTRUE;
    if((P_dc_2x1_nhit[0] + P_dc_2u2_nhit[0] + P_dc_2u1_nhit[0] + P_dc_2v1_nhit[0] + P_dc_2x2_nhit[0] + P_dc_2v2_nhit[0]) > 20) return kTRUE; 
    // Cut on delta, theta and phi for the track
    if (P_gtr_dp[0] < -10.0 || P_gtr_dp[0] > 20.0) return kTRUE;
    if (TMath::Abs(P_gtr_th[0]) > 0.080) return kTRUE;
    if (TMath::Abs(P_gtr_ph[0]) > 0.035) return kTRUE;
    
    SHMS_delta_cut->Fill(P_gtr_dp[0]);
    SHMS_th_cut->Fill(P_gtr_th[0]);
    SHMS_ph_cut->Fill(P_gtr_ph[0]);
}

  if(*EvtType==2){ // HMS single
    TRIG3_cut->Fill(*T_coin_pTRIG3_ROC2_tdcTime);
    HMS_delta->Fill(H_gtr_dp[0]);
    HMS_th->Fill(H_gtr_th[0]);
    HMS_ph->Fill(H_gtr_ph[0]);
    h_etotnorm_before->Fill(H_cal_etotnorm[0]);
    h_cer_before->Fill(H_cer_npeSum[0]);
    h_cal_cer_before->Fill(H_cal_etotnorm[0], H_cer_npeSum[0]);

    // Check we have a good track
    if(H_hod_goodscinhit[0] != 1) return kTRUE;
    if(H_hod_betanotrack[0] < 0.8 || H_hod_betanotrack[0] > 1.3) return kTRUE;
    if((H_dc_1x1_nhit[0] + H_dc_1u2_nhit[0] + H_dc_1u1_nhit[0] + H_dc_1v1_nhit[0] + H_dc_1x2_nhit[0] + H_dc_1v2_nhit[0]) > 20) return kTRUE;
    if((H_dc_2x1_nhit[0] + H_dc_2u2_nhit[0] + H_dc_2u1_nhit[0] + H_dc_2v1_nhit[0] + H_dc_2x2_nhit[0] + H_dc_2v2_nhit[0]) > 20) return kTRUE;
    // Cut on delta, theta and phi
    if (TMath::Abs(H_gtr_dp[0]) > 8.0) return kTRUE;
    if (TMath::Abs(H_gtr_th[0]) > 0.080) return kTRUE;
    if (TMath::Abs(H_gtr_ph[0]) > 0.035) return kTRUE;
    // PID Cuts, check for electron
    if (H_cer_npeSum[0] < 1.5) return kTRUE;
    if (H_cal_etotnorm[0] < 0.7) return kTRUE;

    HMS_delta_cut->Fill(H_gtr_dp[0]);
    HMS_th_cut->Fill(H_gtr_th[0]);
    HMS_ph_cut->Fill(H_gtr_ph[0]);
    h_etotnorm_cut->Fill(H_cal_etotnorm[0]);
    h_cer_cut->Fill(H_cer_npeSum[0]);
    h_cal_cer_cut->Fill(H_cal_etotnorm[0], H_cer_npeSum[0]);
    HMS_W_Dist->Fill(HMSW[0]);
    HMS_ph_q_Dist->Fill(HMSph_q[0]);
}

  if(*EvtType==3){ // Single in each
    TRIG5_cut->Fill(*T_coin_pTRIG5_ROC2_tdcTime);
  }  

  return kTRUE;
}

void HeepSinglesYield::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
}

void HeepSinglesYield::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  TString option = GetOption();
  TFile *Histogram_file = new TFile(Form("../HISTOGRAMS/PionLT_HeepSingles_Run%i.root",option.Atoi()),"RECREATE");

  TDirectory *DTiming = Histogram_file->mkdir("Timing and Event Type Summary"); DTiming->cd();
  EventType->Write();
  TRIG1->Write();
  TRIG3->Write();
  TRIG5->Write();
  TRIG1_cut->Write();
  TRIG3_cut->Write();
  TRIG5_cut->Write();

  TDirectory *DHMS = Histogram_file->mkdir("HMS Events"); DHMS->cd();
  HMS_delta->Write();
  HMS_th->Write();
  HMS_ph->Write();
  h_etotnorm_before->Write();  
  h_cer_before->Write();
  h_cal_cer_before->Write();
  HMS_delta_cut->Write();
  HMS_th_cut->Write();
  HMS_ph_cut->Write();
  h_etotnorm_cut->Write();  
  h_cer_cut->Write();
  h_cal_cer_cut->Write();
  HMS_W_Dist->Write();
  HMS_ph_q_Dist->Write();

  TDirectory *DSHMS = Histogram_file->mkdir("SHMS Events"); DSHMS->cd();
  SHMS_delta->Write();
  SHMS_th->Write();
  SHMS_ph->Write();
  SHMS_delta_cut->Write();
  SHMS_th_cut->Write();
  SHMS_ph_cut->Write();
  p_etotnorm_before->Write();

  Histogram_file->Close();

}
