#ifndef HeepSinglesYield_h
#define HeepSingesYield_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TArc.h>

// Headers needed by this particular selector

class HeepSinglesYield : public TSelector {
 public :
  TTreeReader     fReader;  //!the tree reader
  TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

  //Declare Histograms

  TH1F           *EventType;

  TH1F           *TRIG1ROC1; 
  TH1F           *TRIG3ROC1; 
  TH1F           *TRIG5ROC1;
  TH1F           *TRIG1ROC1_cut; 
  TH1F           *TRIG3ROC1_cut; 
  TH1F           *TRIG5ROC1_cut; 
  TH1F           *TRIG1ROC2; 
  TH1F           *TRIG3ROC2; 
  TH1F           *TRIG5ROC2;
  TH1F           *TRIG1ROC2_cut; 
  TH1F           *TRIG3ROC2_cut; 
  TH1F           *TRIG5ROC2_cut; 

  TH1F           *HMS_delta;
  TH1F           *HMS_delta_cut;
  TH1F           *HMS_th;
  TH1F           *HMS_th_cut;
  TH1F           *HMS_ph;
  TH1F           *HMS_ph_cut;
  TH1F           *HMS_etotnorm_before;
  TH1F           *HMS_etotnorm_cut;
  TH1F           *HMS_cer_before;
  TH1F           *HMS_cer_cut;
  TH2F           *HMS_cal_cer_before;
  TH2F           *HMS_cal_cer_cut;

  TH1F           *HMS_W_Dist;
  TH1F           *HMS_ph_q_Dist;
  TH1F           *HMS_xpfp_Dist;
  TH2F           *HMS_W_xpfp;
  TH2F           *HMS_W_Q2;

  TH1F           *SHMS_delta;
  TH1F           *SHMS_delta_cut;
  TH1F           *SHMS_th;
  TH1F           *SHMS_th_cut;
  TH1F           *SHMS_ph;
  TH1F           *SHMS_ph_cut;
  TH1F           *SHMS_etotnorm_before;
  TH2F           *SHMS_Cal_HGC_before;
  TH2F           *SHMS_Cal_Aero_before;
  TH2F           *SHMS_Aero_HGC_before;
  TH1F           *SHMS_etotnorm_cut;
  TH2F           *SHMS_Cal_HGC_cut;
  TH2F           *SHMS_Cal_Aero_cut;
  TH2F           *SHMS_Aero_HGC_cut;

  TH1F           *SHMS_W_Dist;
  TH1F           *SHMS_ph_q_Dist;
  TH1F           *SHMS_xpfp_Dist;
  TH2F           *SHMS_W_xpfp;
  TH2F           *SHMS_W_Q2;
  TH2F           *SHMS_W_HGC;
  TH1F           *SHMS_delta_test;

  // Readers to access the data (delete the ones you do not need).
  TTreeReaderArray<Double_t> P_gtr_beta         = {fReader, "P.gtr.beta"};
  TTreeReaderArray<Double_t> P_gtr_th           = {fReader, "P.gtr.th"};
  TTreeReaderArray<Double_t> P_gtr_ph           = {fReader, "P.gtr.ph"};
  TTreeReaderArray<Double_t> H_gtr_beta         = {fReader, "H.gtr.beta"};
  TTreeReaderArray<Double_t> H_gtr_th           = {fReader, "H.gtr.th"};
  TTreeReaderArray<Double_t> H_gtr_ph           = {fReader, "H.gtr.ph"};
  TTreeReaderArray<Double_t> H_cal_etotnorm     = {fReader, "H.cal.etotnorm"}; 
  TTreeReaderArray<Double_t> H_cer_npeSum       = {fReader, "H.cer.npeSum"};
  TTreeReaderArray<Double_t> P_cal_etotnorm     = {fReader, "P.cal.etotnorm"};
  TTreeReaderArray<Double_t> P_aero_npeSum      = {fReader, "P.aero.npeSum"};
  TTreeReaderArray<Double_t> P_hgcer_npeSum     = {fReader, "P.hgcer.npeSum"};
  TTreeReaderArray<Double_t> H_gtr_dp           = {fReader, "H.gtr.dp"};
  TTreeReaderArray<Double_t> P_gtr_dp           = {fReader, "P.gtr.dp"};
  TTreeReaderArray<Double_t> H_gtr_p            = {fReader, "H.gtr.p"};
  TTreeReaderArray<Double_t> P_gtr_p            = {fReader, "P.gtr.p"}; 
  TTreeReaderArray<Double_t> H_xpfp             = {fReader, "H.dc.xp_fp"};
  TTreeReaderArray<Double_t> P_xpfp             = {fReader, "P.dc.xp_fp"};

  // Kinematics
  TTreeReaderArray<Double_t> HMSQ2                 = {fReader, "H.kin.primary.Q2"};
  TTreeReaderArray<Double_t> HMSW                  = {fReader, "H.kin.primary.W"};
  TTreeReaderArray<Double_t> HMSepsilon            = {fReader, "H.kin.primary.epsilon"};
  TTreeReaderArray<Double_t> HMSph_q               = {fReader, "H.kin.primary.ph_q"};  
  TTreeReaderArray<Double_t> SHMSQ2                = {fReader, "P.kin.primary.Q2"};
  TTreeReaderArray<Double_t> SHMSW                 = {fReader, "P.kin.primary.W"};
  TTreeReaderArray<Double_t> SHMSepsilon           = {fReader, "P.kin.primary.epsilon"};
  TTreeReaderArray<Double_t> SHMSph_q              = {fReader, "P.kin.primary.ph_q"};  
  
  TTreeReaderValue<Int_t>    EvtType            = {fReader, "fEvtHdr.fEvtType"};
 
  // HMS Tracking
  TTreeReaderArray<Double_t> H_tr_beta               = {fReader, "H.tr.beta"};
  TTreeReaderArray<Double_t> H_tr_chi2               = {fReader, "H.tr.chi2"};
  TTreeReaderArray<Double_t> H_tr_ndof               = {fReader, "H.tr.ndof"};
  TTreeReaderArray<Double_t> H_hod_goodscinhit       = {fReader, "H.hod.goodscinhit"};
  TTreeReaderArray<Double_t> H_hod_betanotrack       = {fReader, "H.hod.betanotrack"};
  TTreeReaderArray<Double_t> H_dc_ntrack             = {fReader, "H.dc.ntrack"};
  TTreeReaderArray<Double_t> H_dc_1x1_nhit           = {fReader, "H.dc.1x1.nhit"};
  TTreeReaderArray<Double_t> H_dc_1u2_nhit           = {fReader, "H.dc.1u2.nhit"};
  TTreeReaderArray<Double_t> H_dc_1u1_nhit           = {fReader, "H.dc.1u1.nhit"};
  TTreeReaderArray<Double_t> H_dc_1v1_nhit           = {fReader, "H.dc.1v1.nhit"};
  TTreeReaderArray<Double_t> H_dc_1x2_nhit           = {fReader, "H.dc.1x2.nhit"};
  TTreeReaderArray<Double_t> H_dc_1v2_nhit           = {fReader, "H.dc.1v2.nhit"};
  TTreeReaderArray<Double_t> H_dc_2x1_nhit           = {fReader, "H.dc.2x1.nhit"};
  TTreeReaderArray<Double_t> H_dc_2u2_nhit           = {fReader, "H.dc.2u2.nhit"};
  TTreeReaderArray<Double_t> H_dc_2u1_nhit           = {fReader, "H.dc.2u1.nhit"};
  TTreeReaderArray<Double_t> H_dc_2v1_nhit           = {fReader, "H.dc.2v1.nhit"};
  TTreeReaderArray<Double_t> H_dc_2x2_nhit           = {fReader, "H.dc.2x2.nhit"};
  TTreeReaderArray<Double_t> H_dc_2v2_nhit           = {fReader, "H.dc.2v2.nhit"};

  // SHMS Tracking
  TTreeReaderArray<Double_t> P_tr_beta               = {fReader, "P.tr.beta"};
  TTreeReaderArray<Double_t> P_tr_chi2               = {fReader, "P.tr.chi2"};
  TTreeReaderArray<Double_t> P_tr_ndof               = {fReader, "P.tr.ndof"};
  TTreeReaderArray<Double_t> P_hod_goodscinhit       = {fReader, "P.hod.goodscinhit"};
  TTreeReaderArray<Double_t> P_hod_betanotrack       = {fReader, "P.hod.betanotrack"};
  TTreeReaderArray<Double_t> P_dc_ntrack             = {fReader, "P.dc.ntrack"};
  TTreeReaderArray<Double_t> P_dc_1x1_nhit           = {fReader, "P.dc.1x1.nhit"};
  TTreeReaderArray<Double_t> P_dc_1u2_nhit           = {fReader, "P.dc.1u2.nhit"};
  TTreeReaderArray<Double_t> P_dc_1u1_nhit           = {fReader, "P.dc.1u1.nhit"};
  TTreeReaderArray<Double_t> P_dc_1v1_nhit           = {fReader, "P.dc.1v1.nhit"};
  TTreeReaderArray<Double_t> P_dc_1x2_nhit           = {fReader, "P.dc.1x2.nhit"};
  TTreeReaderArray<Double_t> P_dc_1v2_nhit           = {fReader, "P.dc.1v2.nhit"};
  TTreeReaderArray<Double_t> P_dc_2x1_nhit           = {fReader, "P.dc.2x1.nhit"};
  TTreeReaderArray<Double_t> P_dc_2u2_nhit           = {fReader, "P.dc.2u2.nhit"};
  TTreeReaderArray<Double_t> P_dc_2u1_nhit           = {fReader, "P.dc.2u1.nhit"};
  TTreeReaderArray<Double_t> P_dc_2v1_nhit           = {fReader, "P.dc.2v1.nhit"};
  TTreeReaderArray<Double_t> P_dc_2x2_nhit           = {fReader, "P.dc.2x2.nhit"};
  TTreeReaderArray<Double_t> P_dc_2v2_nhit           = {fReader, "P.dc.2v2.nhit"};

  // Timing
  TTreeReaderValue<Double_t> T_coin_pTRIG1_ROC1_tdcTime   = {fReader, "T.coin.pTRIG1_ROC1_tdcTime"}; 
  TTreeReaderValue<Double_t> T_coin_pTRIG3_ROC1_tdcTime   = {fReader, "T.coin.pTRIG3_ROC1_tdcTime"}; 
  TTreeReaderValue<Double_t> T_coin_pTRIG5_ROC1_tdcTime   = {fReader, "T.coin.pTRIG5_ROC1_tdcTime"}; 
  TTreeReaderValue<Double_t> T_coin_pTRIG1_ROC2_tdcTime   = {fReader, "T.coin.pTRIG1_ROC2_tdcTime"}; 
  TTreeReaderValue<Double_t> T_coin_pTRIG3_ROC2_tdcTime   = {fReader, "T.coin.pTRIG3_ROC2_tdcTime"}; 
  TTreeReaderValue<Double_t> T_coin_pTRIG5_ROC2_tdcTime   = {fReader, "T.coin.pTRIG5_ROC2_tdcTime"}; 

  TTreeReaderValue<Double_t> pEDTM              = {fReader, "T.coin.pEDTM_tdcTime"};

  HeepSinglesYield(TTree * /*tree*/ =0) {EventType=0, TRIG1ROC1=0, TRIG3ROC1=0, TRIG5ROC1=0, TRIG1ROC1_cut=0, TRIG3ROC1_cut=0, TRIG5ROC1_cut=0, TRIG1ROC2=0, TRIG3ROC2=0, TRIG5ROC2=0, TRIG1ROC2_cut=0, TRIG3ROC2_cut=0, TRIG5ROC2_cut=0, HMS_delta=0, HMS_delta_cut=0, HMS_th=0, HMS_th_cut=0, HMS_ph=0, HMS_ph_cut=0, HMS_etotnorm_before=0, HMS_etotnorm_cut=0, HMS_cer_before=0, HMS_cer_cut=0, HMS_cal_cer_before =0, HMS_cal_cer_cut=0, HMS_W_Dist=0, HMS_ph_q_Dist=0, HMS_xpfp_Dist=0, HMS_W_xpfp=0, HMS_W_Q2=0, SHMS_delta=0, SHMS_delta_cut=0, SHMS_th=0, SHMS_th_cut=0, SHMS_ph=0, SHMS_ph_cut=0, SHMS_etotnorm_before=0, SHMS_Cal_HGC_before=0, SHMS_Cal_Aero_before=0, SHMS_Aero_HGC_before=0, SHMS_etotnorm_cut=0, SHMS_Cal_HGC_cut=0, SHMS_Cal_Aero_cut=0, SHMS_Aero_HGC_cut=0, SHMS_W_Dist=0, SHMS_ph_q_Dist=0, SHMS_xpfp_Dist=0, SHMS_W_xpfp=0, SHMS_W_Q2=0, SHMS_W_Dist=0, SHMS_ph_q_Dist=0, SHMS_xpfp_Dist=0, SHMS_W_xpfp=0, SHMS_W_Q2=0, SHMS_W_HGC=0;}
  virtual ~HeepSinglesYield() { }
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();

  ClassDef(HeepSinglesYield,0);
};

#endif

#ifdef HeepSinglesYield_cxx
void HeepSinglesYield::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the reader is initialized.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  fReader.SetTree(tree);
}

Bool_t HeepSinglesYield::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  return kTRUE;
}

#endif // #ifdef HeepSinglesYield_cxx
