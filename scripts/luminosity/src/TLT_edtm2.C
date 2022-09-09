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
#include <TChain.h>
#include <TFile.h>
void TLT_edtm2() {
    
  ofstream myfile("EDTM_2_Study.csv");
  // tdc Time Raw spectra ranges
  int EDTMlo[6] = {3100,2755,3250,2600,3650,3810};
  int EDTMhi[6] = {3400,3055,3550,2900,3950,4110};
  int pTRIG1lo[6] = {3100,2755,3250,2600,3650,3810};
  int pTRIG1hi[6] = {3400,3055,3550,2900,3950,4110};
  // Prescale factors
  int PS[18] = {0,1,2,3,5,9,17,33,65,129,257,513,1025,2049,4097,8193,16385,32769};

  TString inputroot;
  
  //Runs for EDTM Study #2 (13000 added later) along with SHMS(p) and HMS (h) prescale values
  int runs[13] = {899,900,901,902,903,905,912,913,914,915,916,918,919};
  int psh[13] = {5,-1,6,5,5,5,3,-1,1,1,2,2,2};
  int psp[13] = {-1,8,8,8,8,8,6,5,-1,-1,5,5,5};

  for(int i = 0; i < 13; i++){
    // Set correct run number and prescale factor
    int nrun = 13000 + runs[i];
    int psfh = PS[psh[i]+1];
    int psfp = PS[psp[i]+1];
    
    //Adjusting spectrometer trigers to 3/4 instead of EL_REAL for 2 runs
    int pha = 0;
    if(nrun==13912||nrun==13901) pha = -1;

    inputroot=Form("../../../ROOTfiles/Analysis/Lumi/Pion_replay_luminosity_%d_%d.root",nrun,-1);

    cout << " INfile = " << inputroot << endl;
    TFile *fsimc = new TFile(inputroot);





    //Load T Tree-------------------------------------------------------------------------------------

    TTree *tsimc = (TTree*)fsimc->Get("T");
    tsimc->SetMakeClass(1);
    // Define branches
    Double_t  tcurrent;
    tsimc->SetBranchAddress("H.bcm.bcm4a.AvgCurrent",&tcurrent);
    Double_t  tEDTM;
    tsimc->SetBranchAddress("T.coin.pEDTM_tdcTimeRaw",&tEDTM);
    //TBranch *EvtB = tsimc->GetBranch("Event_Branch");
    TBranch *EvtH;
    Int_t  t;
    tsimc->SetBranchAddress("fEvtHdr.fEvtType",&t,&EvtH);

    // Set Current Limit for runs
    double currLim = 50;
    if(nrun>13906) currLim = 7;

    // EDTM accepted for hms, shms, and all non-zero
    int EDTMacch = 0;
    int EDTMaccp = 0;
    int EDTMacca = 0;


    //cout << "\n" << EDTMlo[3+pha] << "\t" << EDTMhi[3+pha] << "\t" << currLim << endl;
    // loop over entries
    Long64_t nentries = tsimc->GetEntries();
    for (int i = 0; i < nentries; i++) {
      //break;
      //if(i%10000==0) cout << endl << i;
      
      tsimc->GetEntry(i);
      // Set Current Cut
      if(tcurrent>=currLim){
	// Cut at EDTM events for SHMS
	if(tEDTM>EDTMlo[1+pha]&&tEDTM<EDTMhi[1+pha]&&(t==1||t==3)){
	  EDTMaccp++;
	}
	// Cut at EDTM events for HMS
	if(tEDTM>EDTMlo[3+pha]&&tEDTM<EDTMhi[3+pha]&&t==2){
	  EDTMacch++;
	}
	// Cut at EDTM events, non-specific
	if(tEDTM>0){
	  EDTMacca++;
	}
      }
    }
    

    // SCALER TREES-------------------------------------------------------------------------------------
   

    Double_t  Hcurrent;
    Double_t  HEDTM;
    Double_t hTime;
    Double_t HpTRIG1;
    Double_t HpTRIG2;
    Double_t HpTRIG3;
    Double_t HpTRIG4;
    Double_t HL1;

    Double_t HEDTMSum = 0, HEDTMCur = 0, HEDTMPre = 0,
             HtimeSum = 0, HtimeCur = 0, HtimePre = 0,
             HT1Sum = 0, HT1Cur = 0, HT1Pre = 0,
             HT2Sum = 0, HT2Cur = 0, HT2Pre = 0,
             HT3Sum = 0, HT3Cur = 0, HT3Pre = 0,
             HT4Sum = 0, HT4Cur = 0, HT4Pre = 0,
             HL1Sum = 0, HL1Cur = 0, HL1Pre = 0;

    if(psfh>psfp){
      TTree *Hsimc = (TTree*)fsimc->Get("TSH");
      // Define branches
      Hsimc->SetBranchAddress("H.BCM4A.scalerCurrent",&Hcurrent);
      Hsimc->SetBranchAddress("H.EDTM.scaler",&HEDTM);
      Hsimc->SetBranchAddress("H.1MHz.scalerTime",&hTime);
      Hsimc->SetBranchAddress("H.pTRIG1.scaler",&HpTRIG1);
      Hsimc->SetBranchAddress("H.pTRIG2.scaler",&HpTRIG2);
      Hsimc->SetBranchAddress("H.pTRIG3.scaler",&HpTRIG3);
      Hsimc->SetBranchAddress("H.pTRIG4.scaler",&HpTRIG4);
      Hsimc->SetBranchAddress("H.hL1ACCP.scaler",&HL1);      
    // loop over entries
    Long64_t Hnentries = Hsimc->GetEntries();
    for (int i = 0; i < Hnentries; i++) {
      Hsimc->GetEntry(i);
      HEDTMCur = HEDTM;
      HtimeCur = hTime;
      HT1Cur = HpTRIG1;
      HT2Cur = HpTRIG2;
      HT3Cur = HpTRIG3;
      HT4Cur = HpTRIG4;
      HL1Cur = HL1;
      if(Hcurrent>currLim){
	HEDTMSum += HEDTMCur-HEDTMPre;
	HtimeSum += HtimeCur-HtimePre;
	HT1Sum += HT1Cur-HT1Pre;
	HT2Sum += HT2Cur-HT2Pre;
	HT3Sum += HT3Cur-HT3Pre;
	HT4Sum += HT4Cur-HT4Pre;
	HL1Sum += HL1Cur-HL1Pre;
      }
      HtimePre = hTime;
      HEDTMPre = HEDTM;
      HT1Pre = HpTRIG1;
      HT2Pre = HpTRIG2;
      HT3Pre = HpTRIG3;
      HT4Pre = HpTRIG4;
      HL1Pre = HL1;
    }
    }
    else{
      TTree *Hsimc = (TTree*)fsimc->Get("TSP");
      // Define branches
      Hsimc->SetBranchAddress("P.BCM4A.scalerCurrent",&Hcurrent);
      Hsimc->SetBranchAddress("P.EDTM.scaler",&HEDTM);
      Hsimc->SetBranchAddress("P.1MHz.scalerTime",&hTime);
      Hsimc->SetBranchAddress("P.pTRIG1.scaler",&HpTRIG1);
      Hsimc->SetBranchAddress("P.pTRIG2.scaler",&HpTRIG2);
      Hsimc->SetBranchAddress("P.pTRIG3.scaler",&HpTRIG3);
      Hsimc->SetBranchAddress("P.pTRIG4.scaler",&HpTRIG4);
      Hsimc->SetBranchAddress("P.pL1ACCP.scaler",&HL1);      
    // loop over entries
    Long64_t Hnentries = Hsimc->GetEntries();
    for (int i = 0; i < Hnentries; i++) {
      Hsimc->GetEntry(i);
      HEDTMCur = HEDTM;
      HtimeCur = hTime;
      HT1Cur = HpTRIG1;
      HT2Cur = HpTRIG2;
      HT3Cur = HpTRIG3;
      HT4Cur = HpTRIG4;
      HL1Cur = HL1;
      if(Hcurrent>currLim){
	HEDTMSum += HEDTMCur-HEDTMPre;
	HtimeSum += HtimeCur-HtimePre;
	HT1Sum += HT1Cur-HT1Pre;
	HT2Sum += HT2Cur-HT2Pre;
	HT3Sum += HT3Cur-HT3Pre;
	HT4Sum += HT4Cur-HT4Pre;
	HL1Sum += HL1Cur-HL1Pre;
      }
      HtimePre = hTime;
      HEDTMPre = HEDTM;
      HT1Pre = HpTRIG1;
      HT2Pre = HpTRIG2;
      HT3Pre = HpTRIG3;
      HT4Pre = HpTRIG4;
      HL1Pre = HL1;
    }
    }



    //OUTPUT----------------------------------------------------------------------------------------------------------
    
    Double_t HT, PT, HLTF, PLTF, HCLTF, PCLTF;

    // Select EL_REAL (pTRIG2/4) or 3/4 (pTRIG 1/3) by pha value
    if(pha==0){ // pTRIG 2 and 4
      //Set scaler counts
      HT = HT4Sum;
      PT = HT2Sum;
      //Set Lifetime Factors to account for prescale
      // For pTRIG 2 and 4, 2 is the earlier trigger (SHMS)
      PLTF = psfp*1.0;
      //Alternate conditions given for single-arm runs
      if(psfp>0){
	HLTF = psfh*1.0/(1-1/PLTF);
	  }	
      else{
	HLTF = psfh*1.0;
      }	
      if(psfh==0) 
	HCLTF=0;
      else HCLTF = 1.0/psfh;
      if(psfp==0) 
	PCLTF=0;
      else PCLTF = 1.0/psfp;
    }
    else{// pTRIG 1 and 3
      //Set scaler counts
      HT = HT3Sum;
      PT = HT1Sum;

      //Set Lifetime Factors to account for prescale
      // For pTRIG 1 and 3, 3 is the earlier trigger (HMS)
      HLTF = psfh*1.0;
      //Alternate conditions given for single-arm runs
      if(psfh>0){
	PLTF = psfp*1.0/(1-1/HLTF);
	  }	
      else{
	PLTF = psfp*1.0;
      }	  
      if(psfh==0) 
	HCLTF=0;
      else HCLTF = 1.0/psfh;
      if(psfp==0) 
	PCLTF=0;
      else PCLTF = 1.0/psfp;    
    }
    cout << "\nRun " << nrun << setprecision(9)
	 << "\tTime: " << HtimeSum
	 << "\tSHMS Prescale: " << psfp
	 << "\tHMS Prescale: " << psfh

	 << "\nSHMS Trig Count: " << PT
	 << "\tHMS Trig Count: " << HT
	 << "\tL1 Acc Count: " << HL1Sum

	 << "\nSHMS Trig Count PS: " << PT*PCLTF
	 << "\tHMS Trig Count PS: " << HT*HCLTF
	 << "\tCPULT: " << HL1Sum*1.0/(PT*PCLTF+HT*HCLTF)

	 << "\nEDTM Sent: " << HEDTMSum
         << "\tEDTM Accepted SHMS: " << EDTMaccp
         << "\tEDTM Accepted HMS: " << EDTMacch
         << "\t TLT SHMS: " << PLTF*EDTMaccp/HEDTMSum
         << "\t TLT HMS: " << HLTF*EDTMacch/HEDTMSum << endl;

    myfile << nrun << "," << HtimeSum << "," << psfp << "," << psfh << "," << setprecision(9) << PT << "," << HT << "," << HL1Sum << "," << HEDTMSum << "," << EDTMaccp << "," << EDTMacch << endl;
  }
  myfile.close();
}
