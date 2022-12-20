#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TBox.h>
#include <TPolyLine.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
void TLT_scaler() {
    
  //ofstream myfile("EDTM_3_Study.csv");
  int EDTMlo[6] = {3100,2755,3250,2590,3650,3810};
  int EDTMhi[6] = {3400,3060,3550,2900,3950,4110};
  int x = 3;
  int y = 3;
  int PS[17] = {1,2,3,5,9,17,33,65,129,257,513,1025,2049,4097,8193,16385,32769};
  TString inputroot;
  int runs[20] = {35,36,37,38,39,40,41,46,47,48,52,53,54,55,56,57,58,59,60,61};
  int psv[20] = {0,1,2,3,4,5,6,11,0,1,5,6,7,8,9,10,11,12,13,14};
  for(int i = 0; i < 20; i++){
    int nrun = 14100 + runs[i];
    int pres = psv[i];
    int psf = PS[pres];
  inputroot=Form("../../../ROOTfiles/Analysis/Lumi/Pion_replay_luminosity_%d_%d.root",nrun,-1);

  cout << " INfile = " << inputroot << endl;
  TFile *fsimc = new TFile(inputroot);
  //  TTree *tsimc = (TTree*)fsimc->Get("T");
  // Define branches
  Double_t  tcurrent;
  // tsimc->SetBranchAddress("H.bcm.bcm4a.AvgCurrent",&tcurrent);
  Double_t  tEDTM;
  // tsimc->SetBranchAddress("T.coin.pEDTM_tdcTimeRaw",&tEDTM);
  double currLim = 53;
  if(nrun>13906) currLim = 8;
  int EDTMacc24 = 0;
  int EDTMacc2 = 0;
  int EDTMacc4 = 0;
  int EDTMacc13 = 0;
  int EDTMacc0 = 0;
  // loop over entries
  // Long64_t nentries = tsimc->GetEntries();
  //for (int i = 0; i < nentries; i++) {
    //if(i%10000==0) cout << endl << i << endl;
    //tsimc->GetEntry(i);
    //if(tcurrent>=currLim){
  //if(tEDTM>EDTMlo[3])
  //	if(tEDTM<EDTMhi[1])
  //	  EDTMacc24++;
  //  if(tEDTM>EDTMlo[1])
  //	if(tEDTM<EDTMhi[1])
  //	  EDTMacc2++;
  //  if(tEDTM>EDTMlo[3])
  //	if(tEDTM<EDTMhi[3])
  //	  EDTMacc4++;
  //  if(tEDTM>EDTMlo[0])
  //	if(tEDTM<EDTMhi[2])
  //	  EDTMacc13++;
  //  if(tEDTM>0)
  //	EDTMacc0++;
      //}
  // }
    

  // SCALER TREES
  TTree *Psimc = (TTree*)fsimc->Get("TSP");
  // Define branches
  Double_t  Pcurrent;
  Psimc->SetBranchAddress("P.BCM4A.scalerCurrent",&Pcurrent);
  Double_t  PEDTM;
  Psimc->SetBranchAddress("P.EDTM.scaler",&PEDTM);

  double EDTMSum = 0;
  double EDTMCur = 0;
  double EDTMPre = 0;
  // loop over entries
  Long64_t Pnentries = Psimc->GetEntries();
  for (int i = 0; i < Pnentries; i++) {
    //if(i%10000==0) cout << i;
    Psimc->GetEntry(i);
    EDTMCur = PEDTM;
    //if(Pcurrent>currLim)
      EDTMSum += EDTMCur-EDTMPre;
    EDTMPre = PEDTM;
  }

  TTree *Hsimc = (TTree*)fsimc->Get("TSH");
  // Define branches
  Double_t  Hcurrent;
  Hsimc->SetBranchAddress("H.BCM4A.scalerCurrent",&Hcurrent);
  Double_t  HEDTM;
  Hsimc->SetBranchAddress("H.EDTM.scaler",&HEDTM);
  Double_t hTime;
  Hsimc->SetBranchAddress("H.1MHz.scalerTime",&hTime);
  Double_t HpTRIG3;
  Hsimc->SetBranchAddress("H.pTRIG3.scaler",&HpTRIG3);
  Double_t HL1;
  Hsimc->SetBranchAddress("H.hL1ACCP.scaler",&HL1);
  double HEDTMSum = 0;
  double HEDTMCur = 0;
  double HEDTMPre = 0;
  double HtimeSum = 0;
  double HtimeCur = 0;
  double HtimePre = 0;
  double HTSum = 0;
  double HTCur = 0;
  double HTPre = 0;
  double HL1Sum = 0;
  double HL1Cur = 0;
  double HL1Pre = 0;
  // loop over entries
  Long64_t Hnentries = Hsimc->GetEntries();
  for (int i = 0; i < Hnentries; i++) {
    //if(i%10000==0) cout << i;
    Hsimc->GetEntry(i);
    HEDTMCur = HEDTM;
    HtimeCur = hTime;
    HTCur = HpTRIG3;
    HL1Cur = HL1;
    //if(Hcurrent>currLim){
      HEDTMSum += HEDTMCur-HEDTMPre;
      HtimeSum += HtimeCur-HtimePre;
      HTSum += HTCur-HTPre;
      HL1Sum += HL1Cur-HL1Pre;
      // }
    HtimePre = hTime;
    HEDTMPre = HEDTM;
    HTPre = HpTRIG3;
    HL1Pre = HL1;
  }


  //OUTPUT
  if(HEDTMSum>EDTMSum) 
    EDTMSum = HEDTMSum;
  int EDTMps = EDTMSum/psf; // EDTMSum/129+EDTMSum/17-EDTMSum/(129*17);
  int HTps = HTSum/psf;
  cout << "\nRun " << nrun
       << "\nTime: " << HtimeSum
       << "\nTrig Count: " << HTSum
       << "\nTrig Count after Prescale: " << HTps
       << "\nL1 Count: " << HL1Sum
       << "\nTrig Ratiot: " << HTps/HL1Sum
       << "\nEDTM Sent: " << EDTMSum
       << "\nEDTM Sent after Prescale: " << EDTMps << endl;
    // << "\nEDTM Accepted: " << EDTMacc24
    //   << "\n TLT: " << EDTMacc24*1.0/EDTMps*1.0 << endl;

  //myfile << nrun << "," << pres << "," << psf << "," << setprecision(9) << HtimeSum << "," << HL1Sum << "," << HTSum << "," << EDTMSum << endl;
  }
  //myfile.close();
}
