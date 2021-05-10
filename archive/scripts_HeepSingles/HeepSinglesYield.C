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

  TRIG2ROC1     = new TH1F("TRIG2ROC1","ROC1 pTRIG2 counts",4000,-1000.0,2000.0);
  TRIG3ROC1     = new TH1F("TRIG3ROC1","ROC1 pTRIG3 counts",4000,-1000.0,2000.0);
  TRIG5ROC1     = new TH1F("TRIG5ROC1","ROC1 pTRIG5 counts",4000,0.0,1000.0);
  TRIG2ROC1_cut = new TH1F("TRIG2ROC1_cut","ROC1 pTRIG2 counts",4000,-1000.0,1000.0);
  TRIG3ROC1_cut = new TH1F("TRIG3ROC1_cut","ROC1 pTRIG3 counts",4000,0.0,1000.0);
  TRIG5ROC1_cut = new TH1F("TRIG5ROC1_cut","ROC1 pTRIG5 counts",4000,0.0,1000.0);
  TRIG2ROC2     = new TH1F("TRIG2ROC2","ROC2 pTRIG2 counts",4000,-1000.0,2000.0);
  TRIG3ROC2     = new TH1F("TRIG3ROC2","ROC2 pTRIG3 counts",4000,-1000.0,2000.0);
  TRIG5ROC2     = new TH1F("TRIG5ROC2","ROC2 pTRIG5 counts",4000,0.0,1000.0);
  TRIG2ROC2_cut = new TH1F("TRIG2ROC2_cut","ROC2 pTRIG2 counts",4000,-1000.0,1000.0);
  TRIG3ROC2_cut = new TH1F("TRIG3ROC2_cut","ROC2 pTRIG3 counts",4000,0.0,1000.0);
  TRIG5ROC2_cut = new TH1F("TRIG5ROC2_cut","ROC2 pTRIG5 counts",4000,0.0,1000.0);

  HMS_delta            = new TH1F("HMS_delta","HMS #delta;#delta;Counts",100,-20,20);
  HMS_delta_cut        = new TH1F("HMS_delta_cut","HMS #delta Cut;#delta;Counts",100,-20,20);
  HMS_th              = new TH1F("HMS_th","HMS #theta Acceptance;#theta;Counts",100,-0.1,0.1);
  HMS_th_cut          = new TH1F("HMS_th_cut","HMS #theta Acceptance with Cut;#theta;Counts",100,-0.1,0.1);
  HMS_ph              = new TH1F("HMS_ph","HMS #phi Acceptance;#phi;Counts",100,-0.1,0.1);
  HMS_ph_cut          = new TH1F("HMS_ph_cut","HMS #phi Acceptance with Cut;#phi;Counts",100,-0.1,0.1);
  HMS_etotnorm_before = new TH1F("HMS_etotnorm_before", "HMS Normalised Calorimeter Energy; Normalised Energy; Counts", 150, 0, 1.5);
  HMS_etotnorm_cut = new TH1F("HMS_etotnorm_cut", "HMS Normalised Calorimeter Energy, PID cut; Normalised Energy; Counts", 150, 0, 1.5);
  HMS_cer_before = new TH1F("HMS_cer_before", "HMS Cherenkov NPE; NPE; Counts", 50, 0, 25);
  HMS_cer_cut = new TH1F("HMS_cer_cut", "HMS Cherenkov NPE; NPE, PID cut; Counts", 50, 0, 25);
  HMS_cal_cer_before = new TH2F("HMS_cal_cer_before", "HMS CerNPE(ETotNorm); Normalised Energy; NPE", 150, 0, 1.5, 50, 0, 25);
  HMS_cal_cer_cut = new TH2F("HMS_cal_cer_cut", "HMS CerNPE(ETotNorm), PID cut; Normalised Energy; NPE", 150, 0, 1.5, 50, 0, 25);
  HMS_W_Dist = new TH1F("HMS_W_Dist", "HMS W Distribution; W/GeV; Counts", 200, 0.6, 1.6);
  HMS_ph_q_Dist = new TH1F("HMS_ph_q_Dist", "HMS #phi_{q} Distribution; #phi_{q}; Counts", 200, -0.5, 0.5);
  HMS_xpfp_Dist = new TH1F("HMS_xpfp_Dist", "HMS x'_{fp} Distribution; x'_{fp}; Counts", 200, -0.1, 0.1);
  HMS_W_xpfp = new TH2F("HMS_W_xpfp", "HMS x'_{fp}(W); W/GeV; x'_{fP}", 200, 0.6, 1.6, 200, -0.1, 0.1);
  HMS_W_Q2 = new TH2F("HMS_W_Q2", "HMS Q^{2}(W); W/GeV; Q^{2}/(GeV/c)^{2}", 200, 0.2, 1.6, 200, 0, 6);

  SHMS_delta           = new TH1F("SHMS_delta","SHMS #delta;#delta;Counts",100,-20,20);
  SHMS_delta_cut       = new TH1F("SHMS_delta_cut","SHMS #delta Cut;#delta;Counts",100,-20,20);
  SHMS_th              = new TH1F("SHMS_th","SHMS #theta Acceptance;#theta;Counts",100,-0.1,0.1);
  SHMS_th_cut          = new TH1F("SHMS_th_cut","SHMS #theta Acceptance with Cut;#theta;Counts",100,-0.1,0.1);
  SHMS_ph              = new TH1F("SHMS_ph","SHMS #phi Acceptance;#phi;Counts",100,-0.1,0.1);
  SHMS_ph_cut          = new TH1F("SHMS_ph_cut","SHMS #phi Acceptance with Cut;#phi;Counts",100,-0.1,0.1);
  SHMS_etotnorm_before = new TH1F("SHMSS_etotnorm_before", "SHMS Normalised Calorimeter Energy; Normalised Energy; Counts", 150, 0, 1.5);
  SHMS_Cal_HGC_before = new TH2F("SHMS_Cal_HGC_before", "NPE_{HGC}(Normalised Energy); Normalised Energy; HGC NPE", 150, 0, 1.5, 100, 0, 50);
  SHMS_Cal_Aero_before = new TH2F("SHMS_Cal_Aero_before", "NPE_{Aero}(Normalised Energy); Normalised Energy; Aero NPE", 150, 0, 1.5, 100, 0, 50);
  SHMS_Aero_HGC_before = new TH2F("SHMS_Aero_HGC_before", "NPE_{HGC}(NPE_{Aero}); Aero NPE; HGC NPE", 100, 0, 50, 100, 0, 50);
  SHMS_etotnorm_cut = new TH1F("SHMS_etotnorm_cut", "SHMS Normalised Calorimeter Energy, PID cut; Normalised Energy; Counts", 150, 0, 1.5);
  SHMS_Cal_HGC_cut = new TH2F("SHMS_Cal_HGC_cut", "NPE_{HGC}(Normalised Energy), PID cut; Normalised Energy; HGC NPE", 150, 0, 1.5, 100, 0, 50);
  SHMS_Cal_Aero_cut = new TH2F("SHMS_Cal_Aero_cut", "NPE_{Aero}(Normalised Energy), PID cut; Normalised Energy; Aero NPE", 150, 0, 1.5, 100, 0, 50);
  SHMS_Aero_HGC_cut = new TH2F("SHMS_Aero_HGC_cut", "NPE_{HGC}(NPE_{Aero}), PID cut; Aero NPE; HGC NPE", 100, 0, 50, 100, 0, 50);  
  SHMS_W_Dist = new TH1F("SHMS_W_Dist", "SHMS W Distribution; W/GeV; Counts", 200, 0.6, 1.6);
  SHMS_ph_q_Dist = new TH1F("SHMS_ph_q_Dist", "SHMS #phi_{q} Distribution; #phi_{q}; Counts", 200, 2.8, 3.2);
  SHMS_xpfp_Dist = new TH1F("SHMS_xpfp_Dist", "SHMS x'_{fp} Distribution; x'_{fp}; Counts", 200, -0.1, 0.1);
  SHMS_W_xpfp = new TH2F("SHMS_W_xpfp", "SHMS x'_{fp}(W); W/GeV; x'_{fP}", 200, 0.6, 1.6, 200, -0.1, 0.1);
  SHMS_W_Q2 = new TH2F("SHMS_W_Q2", "SHMS Q^{2}(W); W/GeV; Q^{2}/(GeV/c)^{2}", 200, 0.2, 1.6, 200, 0, 6);
  SHMS_W_HGC = new TH2F("SHMS_W_HGC", "SHMS NPE_{HGC}(W); W/GeV; HGC NPE", 200, 0.2, 1.6, 100, 0, 50);

  GetOutputList()->Add(EventType);
  GetOutputList()->Add(TRIG2ROC1);
  GetOutputList()->Add(TRIG3ROC1);
  GetOutputList()->Add(TRIG5ROC1);
  GetOutputList()->Add(TRIG2ROC1_cut);
  GetOutputList()->Add(TRIG3ROC1_cut);
  GetOutputList()->Add(TRIG5ROC1_cut);
  GetOutputList()->Add(TRIG2ROC2);
  GetOutputList()->Add(TRIG3ROC2);
  GetOutputList()->Add(TRIG5ROC2);
  GetOutputList()->Add(TRIG2ROC2_cut);
  GetOutputList()->Add(TRIG3ROC2_cut);
  GetOutputList()->Add(TRIG5ROC2_cut);

  GetOutputList()->Add(HMS_delta);
  GetOutputList()->Add(HMS_delta_cut);
  GetOutputList()->Add(HMS_th);
  GetOutputList()->Add(HMS_th_cut);
  GetOutputList()->Add(HMS_ph);
  GetOutputList()->Add(HMS_ph_cut);
  GetOutputList()->Add(HMS_etotnorm_before);
  GetOutputList()->Add(HMS_etotnorm_cut);
  GetOutputList()->Add(HMS_cer_before);
  GetOutputList()->Add(HMS_cer_cut);
  GetOutputList()->Add(HMS_cal_cer_before);
  GetOutputList()->Add(HMS_cal_cer_cut);
  GetOutputList()->Add(HMS_W_Dist);
  GetOutputList()->Add(HMS_ph_q_Dist);
  GetOutputList()->Add(HMS_xpfp_Dist);
  GetOutputList()->Add(HMS_W_xpfp);
  GetOutputList()->Add(HMS_W_Q2);

  GetOutputList()->Add(SHMS_delta);
  GetOutputList()->Add(SHMS_delta_cut);
  GetOutputList()->Add(SHMS_th);
  GetOutputList()->Add(SHMS_th_cut);
  GetOutputList()->Add(SHMS_ph);
  GetOutputList()->Add(SHMS_ph_cut);
  GetOutputList()->Add(SHMS_etotnorm_before);
  GetOutputList()->Add(SHMS_Cal_HGC_before);
  GetOutputList()->Add(SHMS_Cal_Aero_before);
  GetOutputList()->Add(SHMS_Aero_HGC_before);
  GetOutputList()->Add(SHMS_etotnorm_cut);
  GetOutputList()->Add(SHMS_Cal_HGC_cut);
  GetOutputList()->Add(SHMS_Cal_Aero_cut);
  GetOutputList()->Add(SHMS_Aero_HGC_cut);
  GetOutputList()->Add(SHMS_W_Dist);
  GetOutputList()->Add(SHMS_ph_q_Dist);
  GetOutputList()->Add(SHMS_xpfp_Dist);
  GetOutputList()->Add(SHMS_W_xpfp);
  GetOutputList()->Add(SHMS_W_Q2);
  GetOutputList()->Add(SHMS_W_HGC);
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
  if (*T_coin_pTRIG2_ROC1_tdcTime!=0.0) TRIG2ROC1->Fill(*T_coin_pTRIG2_ROC1_tdcTime);
  if (*T_coin_pTRIG3_ROC1_tdcTime!=0.0) TRIG3ROC1->Fill(*T_coin_pTRIG3_ROC1_tdcTime);
  if (*T_coin_pTRIG5_ROC1_tdcTime!=0.0) TRIG5ROC1->Fill(*T_coin_pTRIG5_ROC1_tdcTime);
  if (*T_coin_pTRIG2_ROC2_tdcTime!=0.0) TRIG2ROC2->Fill(*T_coin_pTRIG2_ROC2_tdcTime); 
  if (*T_coin_pTRIG3_ROC2_tdcTime!=0.0) TRIG3ROC2->Fill(*T_coin_pTRIG3_ROC2_tdcTime);
  if (*T_coin_pTRIG5_ROC2_tdcTime!=0.0) TRIG5ROC2->Fill(*T_coin_pTRIG5_ROC2_tdcTime);

  if(*EvtType==1){ // SHMS single
    SHMS_delta->Fill(P_gtr_dp[0]);
    SHMS_th->Fill(P_gtr_th[0]);
    SHMS_ph->Fill(P_gtr_ph[0]);
    SHMS_etotnorm_before->Fill(P_cal_etotnorm[0]);
    SHMS_Cal_HGC_before->Fill(P_cal_etotnorm[0], P_hgcer_npeSum[0]);
    SHMS_Cal_Aero_before->Fill(P_cal_etotnorm[0], P_aero_npeSum[0]);
    SHMS_Aero_HGC_before->Fill(P_aero_npeSum[0], P_hgcer_npeSum[0]);    
    
    // Cut on delta, theta and phi for the track
    if (TMath::Abs(P_gtr_dp[0]) > 4) return kTRUE;
    if (TMath::Abs(P_gtr_th[0]) > 0.080) return kTRUE;
    if (TMath::Abs(P_gtr_ph[0]) > 0.035) return kTRUE;
    if(P_cal_etotnorm[0] < 0.7) return kTRUE;
    if(P_hgcer_npeSum[0] < 1.5 || P_aero_npeSum[0] < 1.5) return kTRUE;

    SHMS_delta_cut->Fill(P_gtr_dp[0]);
    SHMS_th_cut->Fill(P_gtr_th[0]);
    SHMS_ph_cut->Fill(P_gtr_ph[0]);
    SHMS_etotnorm_cut->Fill(P_cal_etotnorm[0]);
    SHMS_Cal_HGC_cut->Fill(P_cal_etotnorm[0], P_hgcer_npeSum[0]);
    SHMS_Cal_Aero_cut->Fill(P_cal_etotnorm[0], P_aero_npeSum[0]);
    SHMS_Aero_HGC_cut->Fill(P_aero_npeSum[0], P_hgcer_npeSum[0]); 
    SHMS_W_Dist->Fill(SHMSW[0]);
    SHMS_ph_q_Dist->Fill(TMath::Abs(SHMSph_q[0]));
    SHMS_xpfp_Dist->Fill(P_xpfp[0]);
    SHMS_W_xpfp->Fill(SHMSW[0], P_xpfp[0]);
    SHMS_W_Q2->Fill(SHMSW[0], SHMSQ2[0]);
    SHMS_W_HGC->Fill(SHMSW[0], P_hgcer_npeSum[0]);
    TRIG2ROC1_cut->Fill(*T_coin_pTRIG2_ROC1_tdcTime);
    TRIG2ROC2_cut->Fill(*T_coin_pTRIG2_ROC2_tdcTime);
  }
  
  if(*EvtType==2){ // HMS single
    HMS_delta->Fill(H_gtr_dp[0]);
    HMS_th->Fill(H_gtr_th[0]);
    HMS_ph->Fill(H_gtr_ph[0]);
    HMS_etotnorm_before->Fill(H_cal_etotnorm[0]);
    HMS_cer_before->Fill(H_cer_npeSum[0]);
    HMS_cal_cer_before->Fill(H_cal_etotnorm[0], H_cer_npeSum[0]);

    // Cut on delta, theta and phi
    if (TMath::Abs(H_gtr_dp[0]) > 4.0) return kTRUE;
    if (TMath::Abs(H_gtr_th[0]) > 0.080) return kTRUE;
    if (TMath::Abs(H_gtr_ph[0]) > 0.035) return kTRUE;
    // PID Cuts, check for electron
    if (H_cer_npeSum[0] < 1.5) return kTRUE;
    if (H_cal_etotnorm[0] < 0.7) return kTRUE;

    HMS_delta_cut->Fill(H_gtr_dp[0]);
    HMS_th_cut->Fill(H_gtr_th[0]);
    HMS_ph_cut->Fill(H_gtr_ph[0]);
    HMS_etotnorm_cut->Fill(H_cal_etotnorm[0]);
    HMS_cer_cut->Fill(H_cer_npeSum[0]);
    HMS_cal_cer_cut->Fill(H_cal_etotnorm[0], H_cer_npeSum[0]);
    HMS_W_Dist->Fill(HMSW[0]);
    HMS_ph_q_Dist->Fill(HMSph_q[0]);
    HMS_xpfp_Dist->Fill(H_xpfp[0]);
    HMS_W_xpfp->Fill(HMSW[0], H_xpfp[0]);
    HMS_W_Q2->Fill(HMSW[0], HMSQ2[0]); 
    TRIG3ROC1_cut->Fill(*T_coin_pTRIG3_ROC1_tdcTime);
    TRIG3ROC2_cut->Fill(*T_coin_pTRIG3_ROC2_tdcTime);
  }

  if(*EvtType==3){ // Single in each
    Bool_t GoodHMS = kFALSE;
    Bool_t GoodSHMS = kFALSE;
    HMS_delta->Fill(H_gtr_dp[0]);
    HMS_th->Fill(H_gtr_th[0]);
    HMS_ph->Fill(H_gtr_ph[0]);
    HMS_etotnorm_before->Fill(H_cal_etotnorm[0]);
    HMS_cer_before->Fill(H_cer_npeSum[0]);
    HMS_cal_cer_before->Fill(H_cal_etotnorm[0], H_cer_npeSum[0]);
    SHMS_delta->Fill(P_gtr_dp[0]);
    SHMS_th->Fill(P_gtr_th[0]);
    SHMS_ph->Fill(P_gtr_ph[0]);
    SHMS_etotnorm_before->Fill(P_cal_etotnorm[0]);
    SHMS_Cal_HGC_before->Fill(P_cal_etotnorm[0], P_hgcer_npeSum[0]);
    SHMS_Cal_Aero_before->Fill(P_cal_etotnorm[0], P_aero_npeSum[0]);
    SHMS_Aero_HGC_before->Fill(P_aero_npeSum[0], P_hgcer_npeSum[0]);
    
    // Set flag for HMS event passing all cuts
    if(TMath::Abs(H_gtr_dp[0]) < 4.0 && TMath::Abs(H_gtr_th[0]) < 0.080 && TMath::Abs(H_gtr_ph[0]) < 0.035 && H_cer_npeSum[0] > 1.5 && H_cal_etotnorm[0] < 0.7) GoodHMS = kTRUE;
    // Set flag for SHMS events passing all cuts
    if(TMath::Abs(P_gtr_dp[0]) < 4 && TMath::Abs(P_gtr_th[0]) < 0.080 && TMath::Abs(P_gtr_ph[0]) < 0.035 && P_cal_etotnorm[0] > 0.7 && P_aero_npeSum[0] > 1.5 && P_hgcer_npeSum[0] > 1.5) GoodSHMS = kTRUE;
  
    if(GoodHMS == kFALSE && GoodSHMS == kFALSE) return kTRUE;
    else if (GoodHMS == kTRUE && GoodSHMS == kTRUE){
      HMS_delta_cut->Fill(H_gtr_dp[0]);
      HMS_th_cut->Fill(H_gtr_th[0]);
      HMS_ph_cut->Fill(H_gtr_ph[0]);
      HMS_etotnorm_cut->Fill(H_cal_etotnorm[0]);
      HMS_cer_cut->Fill(H_cer_npeSum[0]);
      HMS_cal_cer_cut->Fill(H_cal_etotnorm[0], H_cer_npeSum[0]);
      HMS_W_Dist->Fill(HMSW[0]);
      HMS_ph_q_Dist->Fill(HMSph_q[0]);
      HMS_xpfp_Dist->Fill(H_xpfp[0]);
      HMS_W_xpfp->Fill(HMSW[0], H_xpfp[0]);
      HMS_W_Q2->Fill(HMSW[0], HMSQ2[0]); 
      SHMS_delta_cut->Fill(P_gtr_dp[0]);
      SHMS_th_cut->Fill(P_gtr_th[0]);
      SHMS_ph_cut->Fill(P_gtr_ph[0]);
      SHMS_etotnorm_cut->Fill(P_cal_etotnorm[0]);
      SHMS_Cal_HGC_cut->Fill(P_cal_etotnorm[0], P_hgcer_npeSum[0]);
      SHMS_Cal_Aero_cut->Fill(P_cal_etotnorm[0], P_aero_npeSum[0]);
      SHMS_Aero_HGC_cut->Fill(P_aero_npeSum[0], P_hgcer_npeSum[0]); 
      SHMS_W_Dist->Fill(SHMSW[0]);
      SHMS_ph_q_Dist->Fill(TMath::Abs(SHMSph_q[0]));
      SHMS_xpfp_Dist->Fill(P_xpfp[0]);
      SHMS_W_xpfp->Fill(SHMSW[0], P_xpfp[0]);
      SHMS_W_Q2->Fill(SHMSW[0], SHMSQ2[0]);
      SHMS_W_HGC->Fill(SHMSW[0], P_hgcer_npeSum[0]); 
      TRIG5ROC1_cut->Fill(*T_coin_pTRIG5_ROC1_tdcTime);
      TRIG5ROC2_cut->Fill(*T_coin_pTRIG5_ROC2_tdcTime); 
    }
    else if (GoodHMS == kTRUE && GoodSHMS == kFALSE){
      HMS_delta_cut->Fill(H_gtr_dp[0]);
      HMS_th_cut->Fill(H_gtr_th[0]);
      HMS_ph_cut->Fill(H_gtr_ph[0]);
      HMS_etotnorm_cut->Fill(H_cal_etotnorm[0]);
      HMS_cer_cut->Fill(H_cer_npeSum[0]);
      HMS_cal_cer_cut->Fill(H_cal_etotnorm[0], H_cer_npeSum[0]);
      HMS_W_Dist->Fill(HMSW[0]);
      HMS_ph_q_Dist->Fill(HMSph_q[0]);
      HMS_xpfp_Dist->Fill(H_xpfp[0]);
      HMS_W_xpfp->Fill(HMSW[0], H_xpfp[0]);
      HMS_W_Q2->Fill(HMSW[0], HMSQ2[0]); 
      TRIG5ROC1_cut->Fill(*T_coin_pTRIG5_ROC1_tdcTime);
      TRIG5ROC2_cut->Fill(*T_coin_pTRIG5_ROC2_tdcTime); 
    }
    else if (GoodHMS == kFALSE && GoodSHMS == kTRUE){
      SHMS_delta_cut->Fill(P_gtr_dp[0]);
      SHMS_th_cut->Fill(P_gtr_th[0]);
      SHMS_ph_cut->Fill(P_gtr_ph[0]);
      SHMS_etotnorm_cut->Fill(P_cal_etotnorm[0]);
      SHMS_Cal_HGC_cut->Fill(P_cal_etotnorm[0], P_hgcer_npeSum[0]);
      SHMS_Cal_Aero_cut->Fill(P_cal_etotnorm[0], P_aero_npeSum[0]);
      SHMS_Aero_HGC_cut->Fill(P_aero_npeSum[0], P_hgcer_npeSum[0]); 
      SHMS_W_Dist->Fill(SHMSW[0]);
      SHMS_ph_q_Dist->Fill(TMath::Abs(SHMSph_q[0]));
      SHMS_xpfp_Dist->Fill(P_xpfp[0]);
      SHMS_W_xpfp->Fill(SHMSW[0], P_xpfp[0]);
      SHMS_W_Q2->Fill(SHMSW[0], SHMSQ2[0]);
      SHMS_W_HGC->Fill(SHMSW[0], P_hgcer_npeSum[0]); 
      TRIG5ROC1_cut->Fill(*T_coin_pTRIG5_ROC1_tdcTime);
      TRIG5ROC2_cut->Fill(*T_coin_pTRIG5_ROC2_tdcTime);
    }
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

  //Start of Canvas Painting  
  TCanvas *cHMSCuts = new TCanvas("HMSCuts","Summary of HMS Cuts");
  cHMSCuts->Divide(4,2);
  cHMSCuts->cd(1); HMS_delta->Draw();
  cHMSCuts->cd(2); HMS_delta_cut->Draw();
  cHMSCuts->cd(3); HMS_th->Draw();
  cHMSCuts->cd(4); HMS_th_cut->Draw();
  cHMSCuts->cd(5); HMS_ph->Draw();
  cHMSCuts->cd(6); HMS_ph_cut->Draw();
  cHMSCuts->cd(7); HMS_cal_cer_before->Draw("COLZ");
  cHMSCuts->cd(8); HMS_cal_cer_cut->Draw("COLZ");
  TCanvas *cSHMSCuts = new TCanvas("SHMSCuts","Summary of SHMS Cuts");
  cSHMSCuts->Divide(4,3);
  cSHMSCuts->cd(1); SHMS_delta->Draw();
  cSHMSCuts->cd(2); SHMS_delta_cut->Draw();
  cSHMSCuts->cd(3); SHMS_th->Draw();
  cSHMSCuts->cd(4); SHMS_th_cut->Draw();
  cSHMSCuts->cd(5); SHMS_ph->Draw();
  cSHMSCuts->cd(6); SHMS_ph_cut->Draw();
  cSHMSCuts->cd(7); SHMS_Cal_HGC_before->Draw("COLZ");
  cSHMSCuts->cd(8); SHMS_Cal_HGC_cut->Draw("COLZ");
  cSHMSCuts->cd(9); SHMS_Cal_Aero_before->Draw("COLZ");
  cSHMSCuts->cd(10); SHMS_Cal_Aero_cut->Draw("COLZ");
  cSHMSCuts->cd(11); SHMS_Aero_HGC_before->Draw("COLZ");
  cSHMSCuts->cd(12); SHMS_Aero_HGC_cut->Draw("COLZ");
  //TString foutname = Form("../OUTPUT/HeepSingles_Run%i",option.Atoi());
  // 06/04/21 - SJDK - Horrible hard coded path but just retesting this quickly.
  // PDF output also seems to be broken as well, can't open the produced file
  TString foutname = Form("/group/c-kaonlt/USERS/${USER}/hallc_replay_lt/OUTPUT/Analysis/HeeP/HeepSingles_Run%i",option.Atoi());
  TString outputpdf = foutname + ".pdf";
  TCanvas *cHMS = new TCanvas("cHMS","Summary of HMS Events");
  cHMS->Divide(2,2);
  cHMS->cd(1); HMS_W_Dist->Draw("hist");
  cHMS->Update();
  TLine *Proton_Mass = new TLine(0.938272, 0, 0.938272, gPad->GetUymax()); 
  Proton_Mass->SetLineColor(kRed); Proton_Mass->SetLineWidth(2); Proton_Mass->SetLineStyle(2);
  Proton_Mass->Draw("SAME");
  cHMS->cd(2); HMS_W_xpfp->Draw("COLZ");
  cHMS->Update();
  TLine *Proton_Mass2 = new TLine(0.938272, gPad->GetUymin(), 0.938272, gPad->GetUymax()); // This looks stupid and I hate it, but without it drawing the line screws up. Thanks root.
  Proton_Mass2->SetLineColor(kRed); Proton_Mass2->SetLineWidth(2); Proton_Mass2->SetLineStyle(2);
  Proton_Mass2->Draw("SAME");
  cHMS->cd(3); HMS_ph_q_Dist->Draw();
  cHMS->cd(4); HMS_W_Q2->SetTitleOffset(0.8,"Y"); HMS_W_Q2->Draw("COLZ");
  cHMS->Print(outputpdf + '(');
  TCanvas *cSHMS = new TCanvas("cSHMS","Summary of HMS Events");
  cSHMS->Divide(2,2);
  cSHMS->cd(1); SHMS_W_Dist->Draw("hist");
  cSHMS->Update();
  TLine *Proton_Mass3 = new TLine(0.938272, 0, 0.938272, gPad->GetUymax()); 
  Proton_Mass3->SetLineColor(kRed); Proton_Mass3->SetLineWidth(2); Proton_Mass3->SetLineStyle(2);
  Proton_Mass3->Draw("SAME");
  cSHMS->cd(2); SHMS_W_xpfp->Draw("COLZ");
  cSHMS->Update();
  TLine *Proton_Mass4 = new TLine(0.938272, gPad->GetUymin(), 0.938272,gPad->GetUymax()); // This looks stupid and I hate it, but without it drawing the line screws up. Thanks root.
  Proton_Mass4->SetLineColor(kRed); Proton_Mass4->SetLineWidth(2); Proton_Mass4->SetLineStyle(2);
  Proton_Mass4->Draw("SAME");
  cSHMS->cd(3); SHMS_ph_q_Dist->Draw();
  cSHMS->cd(4); SHMS_W_Q2->SetTitleOffset(0.8,"Y"); SHMS_W_Q2->Draw("COLZ");
  cSHMS->Print(outputpdf + ')');

  //TFile *Histogram_file = new TFile(Form("../HISTOGRAMS/PionLT_HeepSingles_Run%i.root",option.Atoi()),"RECREATE");
  // 06/04/21 - SJDK - Again, horrible hard coded path but just testing
  TFile *Histogram_file = new TFile(Form("/group/c-kaonlt/USERS/${USER}/hallc_replay_lt/HISTOGRAMS/Analysis/HeeP//HeepSingles_Run%i.root",option.Atoi()),"RECREATE");
 
  TDirectory *DTiming = Histogram_file->mkdir("Timing and Event Type Summary"); DTiming->cd();
  EventType->Write();
  TRIG2ROC1->Write();
  TRIG3ROC1->Write();
  TRIG5ROC1->Write();
  TRIG2ROC1_cut->Write();
  TRIG3ROC1_cut->Write();
  TRIG5ROC1_cut->Write();
  TRIG2ROC2->Write();
  TRIG3ROC2->Write();
  TRIG5ROC2->Write();
  TRIG2ROC2_cut->Write();
  TRIG3ROC2_cut->Write();
  TRIG5ROC2_cut->Write();

  TDirectory *DHMS_cuts = Histogram_file->mkdir("HMS Event Selection Cuts"); DHMS_cuts->cd();
  HMS_delta->Write();
  HMS_th->Write();
  HMS_ph->Write();
  HMS_etotnorm_before->Write();  
  HMS_cer_before->Write();
  HMS_cal_cer_before->Write();
  HMS_delta_cut->Write();
  HMS_th_cut->Write();
  HMS_ph_cut->Write();
  HMS_etotnorm_cut->Write();  
  HMS_cer_cut->Write();
  HMS_cal_cer_cut->Write();

  TDirectory *DHMS = Histogram_file->mkdir("HMS Events"); DHMS->cd();
  HMS_W_Dist->Write();
  HMS_ph_q_Dist->Write();
  HMS_xpfp_Dist->Write();
  HMS_W_xpfp->Write();
  HMS_W_Q2->Write();

  TDirectory *DSHMS_cuts = Histogram_file->mkdir("SHMS Event Selection Cuts"); DSHMS_cuts->cd();
  SHMS_delta->Write();
  SHMS_th->Write();
  SHMS_ph->Write();
  SHMS_delta_cut->Write();
  SHMS_th_cut->Write();
  SHMS_ph_cut->Write();
  SHMS_etotnorm_before->Write();
  SHMS_Cal_HGC_before->Write();
  SHMS_Cal_Aero_before->Write();
  SHMS_Aero_HGC_before->Write();
  SHMS_Cal_HGC_cut->Write();
  SHMS_Cal_Aero_cut->Write();
  SHMS_Aero_HGC_cut->Write();

  TDirectory *DSHMS = Histogram_file->mkdir("SHMS Events"); DSHMS->cd();
  SHMS_W_Dist->Write();
  SHMS_ph_q_Dist->Write();
  SHMS_xpfp_Dist->Write();
  SHMS_W_xpfp->Write();
  SHMS_W_Q2->Write();
  SHMS_W_HGC->Write();

  Histogram_file->Close();

}
