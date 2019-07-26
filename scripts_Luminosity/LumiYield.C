#define LumiYield_cxx

#include "LumiYield.h"
#include <fstream>
#include <TStyle.h>


void LumiYield::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  TString option = GetOption();
  Info("Begin", "Options entered: %s", option.Data());
  printf("\n\n");

  fFullShow = kTRUE;
}

void LumiYield::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
  TString PS2_temp = option(0,option.Index("."));
  TString PS3_temp = option(option.Index(".")+1,option.Length());
   
  h_ecut_before = new TH1F("h_ecut_before","HMS CER counts before electron cut",100,0.0,20);
  h_ecut_after  = new TH1F("h_ecut_after" ,"HMS CER counts after electron cut", 100,0.0,20);
  h_ecut_eff    = new TH1F("h_ecut_eff"   ,"HMS CER counts after PID cut",      100,0.0,20);
  p_ecut_before = new TH1F("p_ecut_before","SHMS HGC counts before electron cut",200,0.0,50);
  p_ecut_after  = new TH1F("p_ecut_after" ,"SHMS HGC counts after electron cut", 200,0.0,50);
  p_ecut_eff    = new TH1F("p_ecut_eff"   ,"SHMS HGC counts after PID cut",      200,0.0,50);
  
  h_dp_before   = new TH1F("h_dp_before"  ,"Counts before delta cut",  100,-30,30);
  h_dp_after    = new TH1F("h_dp_after"   ,"Counts after delta cut",   100,-30,30);

  h_th_before   = new TH1F("h_th_before"  ,"Counts before theta cut",  100,-0.1,0.1);
  h_th_after    = new TH1F("h_th_after"   ,"Counts after theta cut",   100,-0.1,0.1);

  h_ph_before   = new TH1F("h_ph_before"  ,"Counts before phi cut",  100,-0.1,0.1);
  h_ph_after    = new TH1F("h_ph_after"   ,"Counts after phi cut",   100,-0.1,0.1);

  h_show_before = new TH1F("h_show_before","Counts in HMS CAL before etotnorm cut",100,0.0,2.0);
  h_show_after  = new TH1F("h_show_after" ,"Counts in HMS CAL after etotnorm cut",100,0.0,2.0);
  p_show_before = new TH1F("p_show_before","Counts in SHMS CAL before etotnorm cut",100,0.0,2.0);
  p_show_after  = new TH1F("p_show_after" ,"Counts in SHMS CAL after etotnorm cut",100,0.0,2.0);
  
  h_track_before   = new TH1F("h_track_before"   ,"Counts before track cut HMS",200,0.0,16.0);
  h_track_after    = new TH1F("h_track_after"    ,"Counts after track cut HMS",200,0.0,16.0);
  h_etrack_before  = new TH1F("h_etrack_before"  ,"Counts before electron track cut HMS",200,0.0,16.0);
  h_etrack_after   = new TH1F("h_etrack_after"   ,"Counts after electron track cut HMS",200,0.0,16.0);
  p_track_before   = new TH1F("p_track_before"   ,"Counts before track cut SHMS",200,0.0,16.0);
  p_track_after    = new TH1F("p_track_after"    ,"Counts after track cut SHMS",200,0.0,16.0);
  p_etrack_before= new TH1F("p_etrack_before","Counts before electron track cut SHMS",200,0.0,16.0);
  p_etrack_after = new TH1F("p_etrack_after" ,"Counts after electron track cut SHMS",200,0.0,16.0);
  EventType   = new TH1F("Event_Type" ,"Counts for each Event Type",20,0,10);
  bcm_before  = new TH1F("bcm_before" ,"Counts before current cut", 100,0,100);
  bcm_after   = new TH1F("bcm_after"  ,"Counts after current cut", 100,0,100);

  EDTM  = new TH1F("EDTM","EDTM counts",4000,-1000.0,2000.0);
  HMS_EDTM  = new TH1F("HMS_EDTM","HMS EDTM counts",4000,-1000.0,1000.0);
  SHMS_EDTM = new TH1F("SHMS_EDTM","SHMS EDTM counts",4000,-1000.0,2000.0);
  TRIG2     = new TH1F("TRIG2","pTRIG2 counts",4000,-1000.0,2000.0);
  TRIG3     = new TH1F("TRIG3","pTRIG3 counts",4000,-1000.0,2000.0);
  TRIG5     = new TH1F("TRIG5","pTRIG5 counts",4000,0.0,1000.0);
  TRIG2_cut = new TH1F("TRIG2_cut","pTRIG1 counts",4000,-1000.0,1000.0);
  TRIG3_cut = new TH1F("TRIG3_cut","pTRIG3 counts",4000,0.0,1000.0);
  TRIG5_cut = new TH1F("TRIG5_cut","pTRIG5 counts",4000,0.0,1000.0);
  
  GetOutputList()->Add(h_ecut_before);
  GetOutputList()->Add(h_ecut_after);
  GetOutputList()->Add(h_ecut_eff);
  GetOutputList()->Add(p_ecut_before);
  GetOutputList()->Add(p_ecut_after);
  GetOutputList()->Add(p_ecut_eff);
  GetOutputList()->Add(h_dp_before);
  GetOutputList()->Add(h_dp_after);
  GetOutputList()->Add(h_th_before);
  GetOutputList()->Add(h_th_after);
  GetOutputList()->Add(h_ph_before);
  GetOutputList()->Add(h_ph_after);
  GetOutputList()->Add(h_show_before);
  GetOutputList()->Add(h_show_after);
  GetOutputList()->Add(p_show_before);
  GetOutputList()->Add(p_show_after);
  GetOutputList()->Add(h_track_before);
  GetOutputList()->Add(h_track_after);
  GetOutputList()->Add(h_etrack_before);
  GetOutputList()->Add(h_etrack_after);
  GetOutputList()->Add(p_track_before);
  GetOutputList()->Add(p_track_after);
  GetOutputList()->Add(p_etrack_before);
  GetOutputList()->Add(p_etrack_after);
 
  GetOutputList()->Add(EventType);
  GetOutputList()->Add(bcm_before);
  GetOutputList()->Add(bcm_after);
  GetOutputList()->Add(EDTM);
  GetOutputList()->Add(HMS_EDTM);
  GetOutputList()->Add(SHMS_EDTM);
  GetOutputList()->Add(TRIG2);
  GetOutputList()->Add(TRIG3);
  GetOutputList()->Add(TRIG5);
  GetOutputList()->Add(TRIG2_cut);
  GetOutputList()->Add(TRIG3_cut);
  GetOutputList()->Add(TRIG5_cut);
}

Bool_t LumiYield::Process(Long64_t entry)
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

  bcm_before->Fill(*H_bcm_bcm4b_AvgCurrent);

  if (*H_bcm_bcm4b_AvgCurrent < 2.5) return kTRUE;

  bcm_after->Fill(*H_bcm_bcm4b_AvgCurrent);

  if (*T_coin_pEDTM_tdcTime!=0.0) EDTM->Fill(*T_coin_pEDTM_tdcTime);
  if (*T_coin_pTRIG2_ROC2_tdcTime!=0.0) TRIG2->Fill(*T_coin_pTRIG2_ROC2_tdcTime);
  if (*T_coin_pTRIG3_ROC2_tdcTime!=0.0) TRIG3->Fill(*T_coin_pTRIG3_ROC2_tdcTime);
  if (*T_coin_pTRIG5_ROC2_tdcTime!=0.0) TRIG5->Fill(*T_coin_pTRIG5_ROC2_tdcTime);
  EventType->Fill(*EvtType);
  
  if ( *EvtType==1 && *T_coin_pEDTM_tdcTime>-5 && *T_coin_pEDTM_tdcTime<5)
    {
      SHMS_EDTM->Fill(*T_coin_pEDTM_tdcTime);
    }

  if (*EvtType==1) // Event was an SHMS Single
    {
      TRIG2_cut->Fill(*T_coin_pTRIG2_ROC2_tdcTime);
    
      //Tracking efficiency calculation, fiducial cut region based off
      //DEF-files/SHMS/PRODUCTION/CUTS/pstackana_reconstruct_cuts.def file
      if (P_hod_goodscinhit[0] == 1                                         &&
	  P_hod_betanotrack[0] > 0.5 && P_hod_betanotrack[0] < 1.4          &&
	  (P_dc_1x1_nhit[0] + P_dc_1u2_nhit[0] + P_dc_1u1_nhit[0] + 
	   P_dc_1v1_nhit[0] + P_dc_1x2_nhit[0] + P_dc_1v2_nhit[0]) < 20     &&	 
	  (P_dc_2x1_nhit[0] + P_dc_2u2_nhit[0] + P_dc_2u1_nhit[0] + 
	   P_dc_2v1_nhit[0] + P_dc_2x2_nhit[0] + P_dc_2v2_nhit[0]) < 20)
	{
	  p_track_before->Fill(P_dc_ntrack[0]);
	  if (P_cal_etotnorm[0] > 0.7 && P_aero_npeSum[0] > 1.5 && P_hgcer_npeSum[0] > 1.5) 
	    {
	      p_etrack_before->Fill(P_dc_ntrack[0]);
	    }
       
	  if (P_dc_ntrack[0] > 0.0) //Requirement that a good track was actually found,, don't require tracking to calculate tracking
	    {
	      p_track_after->Fill(P_dc_ntrack[0]);
	      if (P_cal_etotnorm[0] > 0.7 && P_aero_npeSum[0] > 1.5 && P_hgcer_npeSum[0] > 1.5) 
		{
		  p_etrack_after->Fill(P_dc_ntrack[0]); 
		}
	    }
	}
      
      p_ecut_before->Fill(P_hgcer_npeSum[0]);
      p_show_before->Fill(P_cal_etotnorm[0]);
      if (P_cal_etotnorm[0] < 0.7) return kTRUE;
      if (P_hgcer_npeSum[0] < 1.5) return kTRUE;
      if (P_aero_npeSum[0] < 1.5) return kTRUE;
      if (P_gtr_dp[0] < -10.0 || P_gtr_dp[0] > 20.0) return kTRUE;
      if (TMath::Abs(P_gtr_th[0]) > 0.080) return kTRUE;
      if (TMath::Abs(P_gtr_ph[0]) > 0.035) return kTRUE;
      p_ecut_eff->Fill(P_hgcer_npeSum[0]);
      p_ecut_after->Fill(P_hgcer_npeSum[0]);      
      p_show_after->Fill(P_cal_etotnorm[0]); 
    }
  
  if ( *EvtType==2 &&
      *T_coin_pEDTM_tdcTime>-5 && *T_coin_pEDTM_tdcTime<5)
    {
      HMS_EDTM->Fill(*T_coin_pEDTM_tdcTime);
    }

  if (*EvtType==2) // Event was an HMS Single
    {
      TRIG3_cut->Fill(*T_coin_pTRIG3_ROC2_tdcTime);
      //Tracking efficiency calculation, fiducial cut region based off
      //DEF-files/HMS/PRODUCTION/CUTS/hstackana_reconstruct_cuts.def file
      if (H_hod_goodscinhit[0] == 1                                         &&
	  H_hod_betanotrack[0] > 0.8 && H_hod_betanotrack[0] < 1.3          &&
	  (H_dc_1x1_nhit[0] + H_dc_1u2_nhit[0] + H_dc_1u1_nhit[0] + 
	   H_dc_1v1_nhit[0] + H_dc_1x2_nhit[0] + H_dc_1v2_nhit[0]) < 20     &&	 
	  (H_dc_2x1_nhit[0] + H_dc_2u2_nhit[0] + H_dc_2u1_nhit[0] + 
	   H_dc_2v1_nhit[0] + H_dc_2x2_nhit[0] + H_dc_2v2_nhit[0]) < 20)     
	{
	  h_track_before->Fill(H_dc_ntrack[0]);
	  if (H_cer_npeSum[0] > 1.5 && H_cal_etotnorm[0] > 0.6) {
	    h_etrack_before->Fill(H_dc_ntrack[0]);
	  }
       
	  if (H_dc_ntrack[0] > 0.0) //Requirement that a good track was actually found, don't require tracking to calculate tracking
	    {
	      h_track_after->Fill(H_dc_ntrack[0]);
	      if (H_cer_npeSum[0] > 1.5 && H_cal_etotnorm[0] > 0.6) {
		h_etrack_after->Fill(H_dc_ntrack[0]);
	      }
	    }
	}
      
      h_ecut_before->Fill(H_cer_npeSum[0]);
      h_dp_before->Fill(H_gtr_dp[0]);
      h_th_before->Fill(H_tr_tg_th[0]);
      h_ph_before->Fill(H_tr_tg_ph[0]);
      h_show_before->Fill(H_cal_etotnorm[0]);

      if (H_cal_etotnorm[0] < 0.6) return kTRUE; //.9
      if (H_cer_npeSum[0] < 1.5) return kTRUE;
      if (TMath::Abs(H_gtr_dp[0]) > 8.0) return kTRUE;
      if (TMath::Abs(H_tr_tg_th[0]) > 0.080) return kTRUE;
      if (TMath::Abs(H_tr_tg_ph[0]) > 0.035) return kTRUE;

      h_ecut_after->Fill(H_cer_npeSum[0]);
      h_dp_after->Fill(H_gtr_dp[0]);
      h_th_after->Fill(H_tr_tg_th[0]);
      h_ph_after->Fill(H_tr_tg_ph[0]);
      h_show_after->Fill(H_cal_etotnorm[0]);
      h_ecut_eff->Fill(H_cer_npeSum[0]);
    }

  return kTRUE;
}

void LumiYield::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}

void LumiYield::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  printf("\n\n");
  
  TString option = GetOption();
  TString PS2_temp = option(0,option.Index("."));
  TString PS3_temp = option(option.Index(".")+1,option.Length());
  PS2 = PS2_temp.Atof();
  PS3 = PS3_temp.Atof();

  Info("Terminate","Selection rules have been applied, plotting results");
  cout << Form("Using prescale factors: PS2 %.0f, PS3 %.0f\n",PS2,PS3);
  if(PS2 == 0) PS2 = 1;
  if(PS3 == 0) PS3 = 1;
  cout << Form("Total number of events: %.0f\n",EventType->GetEntries());
  cout << Form("Number of EDTM  Events: %.0f\n",(PS2*SHMS_EDTM->Integral() + PS3*HMS_EDTM->Integral()));
  cout << Form("Number of TRIG2 Events: %.0f\n",(PS2*TRIG2_cut->Integral()));
  cout << Form("Number of TRIG3 Events: %.0f\n",(PS3*TRIG3_cut->Integral()));
  //cout << Form("Number of TRIG5 Events: %.0f\n\n",TRIG5->Integral());

  cout << Form("Number of HMS good events: %.0f",(PS3*h_ecut_eff->GetEntries())) << "  +/- " << sqrt(PS3*h_ecut_eff->GetEntries()) << endl;
  cout << Form("Calculated tracking efficiency: %f +/- %f\n",h_track_after->GetEntries()/h_track_before->GetEntries(),(h_track_after->GetEntries()/h_track_before->GetEntries())*sqrt((1/h_track_after->GetEntries()) + (1/h_track_before->GetEntries())));
  cout << Form("Calculated electron tracking efficiency: %f +/- %f\n",h_etrack_after->GetEntries()/h_etrack_before->GetEntries(),(h_etrack_after->GetEntries()/h_etrack_before->GetEntries())*sqrt((1/h_etrack_after->GetEntries()) + (1/h_etrack_before->GetEntries())));
  cout << Form("Calculated HMS Cherenkov efficiency: %f +/- %f\n\n",h_ecut_eff->GetEntries()/h_ecut_after->GetEntries(),(h_ecut_eff->GetEntries()/h_ecut_after->GetEntries())*sqrt((1/h_ecut_eff->GetEntries()) + (1/h_ecut_after->GetEntries())));

  cout << Form("Number of SHMS good events: %.0f",(PS2*p_ecut_eff->GetEntries())) << "  +/- " << sqrt(PS2*p_ecut_eff->GetEntries()) << endl;
  cout << Form("Calculated tracking efficiency: %f +/- %f\n",p_track_after->GetEntries()/p_track_before->GetEntries(),(p_track_after->GetEntries()/p_track_before->GetEntries())*sqrt((1/p_track_after->GetEntries()) + (1/p_track_before->GetEntries())));
  cout << Form("Calculated electron tracking efficiency: %f +/- %f\n",p_etrack_after->GetEntries()/p_etrack_before->GetEntries(),(p_etrack_after->GetEntries()/p_etrack_before->GetEntries())*sqrt((1/p_etrack_after->GetEntries()) + (1/p_etrack_before->GetEntries())));
   
  TCanvas *c_ID_cut;
  c_ID_cut = new TCanvas("c_ID_cut","Particle ID Information");
  c_ID_cut->Divide(3,2);
  c_ID_cut->cd(1);
  h_ecut_before->Draw();
  c_ID_cut->cd(2);
  h_ecut_after->Draw();
  c_ID_cut->cd(3);
  h_ecut_eff->Draw();
  c_ID_cut->cd(4);
  gPad->SetLogy();
  p_ecut_before->Draw();
  c_ID_cut->cd(5);
  gPad->SetLogy();
  p_ecut_after->Draw();
  c_ID_cut->cd(6);
  gPad->SetLogy();
  p_ecut_eff->Draw();

  TCanvas *c_etot_cut;
  c_etot_cut = new TCanvas("c_etot_cut","Calorimeter Energy Information");
  c_etot_cut->Divide(2,2);
  c_etot_cut->cd(1);
  h_show_before->Draw();
  c_etot_cut->cd(2);
  h_show_after->Draw();
  c_etot_cut->cd(3);
  p_show_before->Draw();
  c_etot_cut->cd(4);
  p_show_after->Draw();
  
  TCanvas *c_track_cut;
  c_track_cut = new TCanvas("c_track_cut","Good Track Information");
  c_track_cut->Divide(2,2);
  c_track_cut->cd(1);
  gPad->SetLogy();
  h_track_before->Draw();
  c_track_cut->cd(2);
  gPad->SetLogy();
  h_track_after->Draw();
  c_track_cut->cd(3);
  gPad->SetLogy();
  p_etrack_before->Draw();
  c_track_cut->cd(4);
  gPad->SetLogy();
  p_etrack_after->Draw();

  TCanvas *c_bcm_cut;
  c_bcm_cut = new TCanvas("c_bcm_cut","Beam Current Monitor Information");
  c_bcm_cut->Divide(2,1);
  c_bcm_cut->cd(1);
  bcm_before->Draw();
  c_bcm_cut->cd(2);
  bcm_after->Draw();

  TCanvas *c_EventType = new TCanvas("Event Type","Event Type Information");
  c_EventType->Divide(2,2);
  c_EventType->cd(1);
  EventType->Draw();
  c_EventType->cd(2);
  EDTM->Draw();
  c_EventType->cd(3);
  TRIG2->Draw();
  c_EventType->cd(4);
  TRIG3->Draw();
 
  ofstream myfile1;
  myfile1.open ("LumiYield_Data.dat", fstream::app);
  myfile1 << Form("%.0f %.0f %.0f %.0f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f",
		  //HMS Evts
		  (h_ecut_eff->GetEntries()),sqrt(h_ecut_eff->GetEntries()),
		  //SHMS Evts
  		  (p_ecut_eff->GetEntries()),sqrt(p_ecut_eff->GetEntries()),
		  //HMS Track
		  h_track_after->GetEntries()/h_track_before->GetEntries(),(h_track_after->GetEntries()/h_track_before->GetEntries())*sqrt((1/h_track_after->GetEntries()) + (1/h_track_before->GetEntries())),
		  //e Track
		  h_etrack_after->GetEntries()/h_etrack_before->GetEntries(),(h_etrack_after->GetEntries()/h_etrack_before->GetEntries())*sqrt((1/h_etrack_after->GetEntries()) + (1/h_etrack_before->GetEntries())),
		  //SHMS Track
		  p_track_after->GetEntries()/p_track_before->GetEntries(),(p_track_after->GetEntries()/p_track_before->GetEntries())*sqrt((1/p_track_after->GetEntries()) + (1/p_track_before->GetEntries())),
		  //h Track
		  p_etrack_after->GetEntries()/p_etrack_before->GetEntries(),(p_etrack_after->GetEntries()/p_etrack_before->GetEntries())*sqrt((1/p_etrack_after->GetEntries()) + (1/p_etrack_before->GetEntries())),
		  //Accept EDTM
		  (SHMS_EDTM->Integral() + HMS_EDTM->Integral()),
		  //PS2
		  PS2,
		  //PS3
		  PS3);;
  myfile1.close();
}
