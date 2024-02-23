/*
	plotSHMSCerEfficiency.C
	Author: Nathan Heinrich
	
	A short script for compairing the efficency of SHMS Cerenkov detectors over delta.
	Includes plots of SHMS HGC, Aero, and NGC.

*/

#include <time.h>
#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include <TH2.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TObjArray.h>
#include <TF1.h>

//Files and containers and histograms
TFile *input1;
TTree *tree1;

//Cut Names
string cutNames[4] = {"Electron", "Pion", "Kaon", "Proton"};

//histograms
TH1D *th1_cal, *th1_calCut, *th1_hgcer, *th1_hgcerCut, *th1_aero, *th1_aeroCut, *th1_ngcer, *th1_ngcerCut;
TH1D *th1_cal_eff, *th1_hgcer_eff, *th1_aero_eff, *th1_ngcer_eff;

TH2D *th2_aeroXhgcer, *th2_ngcerXcal, *th2_ngcerXaero, *th2_ngcerXhgcer;

TH2D *th2_fpXhgcer, *th2_fpXngcer, *th2_fpXaero;
TH2D *th2_fpXhgcer_cut, *th2_fpXngcer_cut, *th2_fpXaero_cut;
TH3D *th3_fpXhgcer_eff, *th3_fpXngcer_eff, *th3_fpXaero_eff;
TH2D *th2_fpXhgcer_eff2D, *th2_fpXngcer_eff2D, *th2_fpXaero_eff2D;

//variables for cutting trees and plotting
Double_t delta, calEtot, hgcerNpeSum, aeroNpeSum, ngcerNpeSum;
Double_t dc_ntrack, InsideDipoleExit, goodStartTime;
Double_t xfp, yfp;

//cuts
const Double_t calEtotLowPi = 0.1; //normaized energy
const Double_t calEtotLowe = 0.7; //normaized energy
const Double_t hgcerNpeSumLow = 1.5; //unit NPE
const Double_t aeroNpeSumLow = 1.5; //unit NPE
const Double_t ngcerNpeSumLow = 1.5; //unit NPE

Bool_t calCut, hgcerCut, aeroCut, ngcerCut, fpcut;

const Int_t INILENGTH = 64;
Int_t NumEvents = -1;


void makePlots ( TString rootFile, Int_t runNum, int NumEvents, int cutType ) 
{

    // make empty histograms	
	th1_hgcer = new TH1D("hgcerShould", "hgcerShould", 120, -30.0, 30.0);
	th1_hgcerCut = new TH1D("hgcerDid", "hgcerDid", 120, -30.0, 30.0);
	th1_aero = new TH1D("aeroShould", "aeroShould", 120, -30.0, 30.0);
	th1_aeroCut = new TH1D("aeroDid", "aeroDid", 120, -30.0, 30.0);
	th1_ngcer = new TH1D("aeroShould", "aeroShould", 120, -30.0, 30.0);
	th1_ngcerCut = new TH1D("aeroDid", "aeroDid", 120, -30.0, 30.0);
	
	
	th2_aeroXhgcer = new TH2D("aeroNpeSumVhgcerNpeSum","aeroNpeSumVhgcerNpeSum", 100, 0.0, 35, 100, 0.0, 35);
	th2_aeroXhgcer->GetXaxix()->SetNameTitle("Aerogel NPE");
	th2_aeroXhgcer->GetYaxix()->SetNameTitle("HGC NPE");
	th2_aeroXhgcer->SetStats(0);
	
	th2_ngcerXcal = new TH2D("ngcerNpeSumVP.cal.etottracknorm","ngcerNpeSumVP.cal.etottracknorm", 100, 0.0, 35, 100, 0.0, 1.6);
	th2_ngcerXcal->GetXaxix()->SetNameTitle("NGC NPE");
	th2_ngcerXcal->GetYaxix()->SetNameTitle("Normalized Calorimeter Energy");
	th2_ngcerXcal->SetStats(0);
	
	th2_ngcerXaero = new TH2D("ngcerNpeSumVaeroNpeSum","ngcerNpeSumVaeroNpeSum", 100, 0.0, 35, 100, 0.0, 35);
	th2_ngcerXaero->GetXaxix()->SetNameTitle("NGC NPE");
	th2_ngcerXaero->GetYaxix()->SetNameTitle("Aerogel NPE");
	th2_ngcerXaero->SetStats(0);
	
	th2_ngcerXhgcer = new TH2D("ngcerNpeSumVhgcerNpeSum","ngcerNpeSumVhgcerNpeSum", 100, 0.0, 35, 100, 0.0, 35);
	th2_ngcerXhgcer->GetXaxix()->SetNameTitle("NGC NPE");
	th2_ngcerXhgcer->GetYaxix()->SetNameTitle("HGC NPE");
	th2_ngcerXhgcer->SetStats(0);
	
	th2_fpXhgcer = new TH2D("fpVhgcereff_should", "fpVhgcereff_should", 80, -40.0, 40.0, 80, -40.0, 40.0);
    th2_fpXngcer = new TH2D("fpVngcereff_should", "fpVngcereff_should", 80, -40.0, 40.0, 80, -40.0, 40.0);
    th2_fpXaero = new TH2D( "fpVaeroeff_should", "fpVaeroeff_should",   80, -40.0, 40.0, 80, -40.0, 40.0);
    th2_fpXhgcer_cut = new TH2D("fpVhgcereff_did", "fpVhgcereff_did",   80, -40.0, 40.0, 80, -40.0, 40.0);
    th2_fpXngcer_cut = new TH2D("fpVngcereff_did", "fpVngcereff_did",   80, -40.0, 40.0, 80, -40.0, 40.0);
    th2_fpXaero_cut = new TH2D( "fpVaeroeff_did", "fpVaeroeff_did",     80, -40.0, 40.0, 80, -40.0, 40.0);

	
	input1 = new TFile(rootFile, "READ");
	cout << "\n";
	
	if(!input1)
	{
		cout << "File not open properly!\nTried to open:\n     " << rootFile;
		return;
	}
	
	tree1 = dynamic_cast <TTree*> (input1->Get("T")); //get T tree from root files
	
	tree1->SetBranchAddress("P.cal.etottracknorm", &calEtot);
	tree1->SetBranchAddress("P.hgcer.npeSum", &hgcerNpeSum);
	tree1->SetBranchAddress("P.ngcer.npeSum", &ngcerNpeSum);
	tree1->SetBranchAddress("P.aero.npeSum", &aeroNpeSum);
	tree1->SetBranchAddress("P.gtr.dp", &delta);
	tree1->SetBranchAddress("P.dc.x_fp", &xfp);
	tree1->SetBranchAddress("P.dc.y_fp", &yfp);
	tree1->SetBranchAddress("P.dc.ntrack", &dc_ntrack); 
	tree1->SetBranchAddress("P.dc.InsideDipoleExit", &InsideDipoleExit); 
	tree1->SetBranchAddress("P.hod.goodstarttime", &goodStartTime);
	
	Int_t nEntries;
	if (NumEvents == -1)
	{
	    nEntries = tree1->GetEntries();
	}else{
	    nEntries = NumEvents;
	}
	cout << "****************************\n" << nEntries << " Entries to be processed in part 1\n";
	for(Int_t iEntry = 0; iEntry < nEntries; iEntry++)
	{
		tree1->GetEntry(iEntry);
		if (iEntry % 10000 == 0) cout << iEntry << endl;
		
        if((dc_ntrack > 0) & (InsideDipoleExit == 1) & (goodStartTime == 1) ) //basic Cuts
        {
            fpcut = (delta < 25.0) && (delta > -15.5);
            // Fill 2D PID plots
            if(fpcut){
                th2_aeroXhgcer->Fill(aeroNpeSum, hgcerNpeSum);
                th2_ngcerXcal->Fill(ngcerNpeSum, calEtot);
                th2_ngcerXaero->Fill(ngcerNpeSum, aeroNpeSum);
                th2_ngcerXhgcer->Fill(ngcerNpeSum, hgcerNpeSum);
            }
        
            if(cutType == 0) // electron
            {
                calCut = (calEtotLowe < calEtot);
                hgcerCut = hgcerNpeSumLow < hgcerNpeSum;
                aeroCut = aeroNpeSumLow < aeroNpeSum;
                ngcerCut = ngcerNpeSumLow < ngcerNpeSum;
            }else if (cutType == 1) // Pion
            {
                calCut = (calEtotLowPi < calEtot);
                hgcerCut = hgcerNpeSumLow < hgcerNpeSum;
                aeroCut = aeroNpeSumLow < aeroNpeSum;
                ngcerCut = ngcerNpeSumLow > ngcerNpeSum;
            }else if (cutType == 2) // Kaon
            {
                calCut = (calEtotLowPi < calEtot);
                hgcerCut = hgcerNpeSumLow > hgcerNpeSum;
                aeroCut = aeroNpeSumLow < aeroNpeSum;
                ngcerCut = ngcerNpeSumLow > ngcerNpeSum;
            }else if (cutType == 3) // Proton
            {
                calCut = (calEtotLowPi < calEtot);
                hgcerCut = hgcerNpeSumLow > hgcerNpeSum;
                aeroCut = aeroNpeSumLow > aeroNpeSum;
                ngcerCut = ngcerNpeSumLow > ngcerNpeSum;
            }
            if (calCut & aeroCut & ngcerCut & fpcut) // should for hgcer
            {
                th1_hgcer->Fill(delta);
                th2_fpXhgcer->Fill(xfp,yfp);
                if (hgcerCut)
                {
                    th1_hgcerCut->Fill(delta);
                    th2_fpXhgcer_cut->Fill(xfp,yfp); // change to fpAthgcer
                }
            }
            if (calCut & aeroCut & hgcerCut & fpcut) // should for ngcer
            {
                th1_ngcer->Fill(delta);
                th2_fpXngcer->Fill(xfp,yfp);
                if (ngcerCut)
                {
                    th1_ngcerCut->Fill(delta);
                    th2_fpXngcer_cut->Fill(xfp,yfp);
                }
            }
            if(calCut & hgcerCut & ngcerCut & fpcut) // should for aero
            {
                th1_aero->Fill(delta);
                th2_fpXaero->Fill(xfp,yfp);
                if(aeroCut)
                {
                    th1_aeroCut->Fill(delta);
                    th2_fpXaero_cut->Fill(xfp,yfp);
                }
            }   
        }	
	}
	
	// do division of plots to get efficiency plots
	Bool_t junk; // for holding return of TH1->Divide() 
	th1_hgcer_eff = new TH1D("hgcer_eff", "hgcer_eff", 120, -30.0, 30.0);
	th1_hgcer_eff->GetXaxix()->SetNameTitle("#delta");
	th1_hgcer_eff->GetYaxix()->SetNameTitle("Efficiency");
	
	th1_aero_eff = new TH1D("aero_eff", "aero_eff", 120, -30.0, 30.0);
	th1_aero_eff->GetXaxix()->SetNameTitle("#delta");
	th1_aero_eff->GetYaxix()->SetNameTitle("Efficiency");
	th1_ngcer_eff = new TH1D("ngcer_eff", "ngcer_eff", 120, -30.0, 30.0);
	th1_ngcer_eff->GetXaxix()->SetNameTitle("#delta");
	th1_ngcer_eff->GetYaxix()->SetNameTitle("Efficiency");
	junk = th1_hgcer_eff->Divide(th1_hgcerCut,th1_hgcer);
    junk = th1_ngcer_eff->Divide(th1_ngcerCut,th1_ngcer);
    junk = th1_aero_eff->Divide(th1_aeroCut,th1_aero);
    
    th2_fpXhgcer_eff2D = new TH2D("fpVhgcereff_2Deff", "fpVhgcereff_2Deff", 80, -40.0, 40.0, 80, -40.0, 40.0);
    junk = th2_fpXhgcer_eff2D->Divide(th2_fpXhgcer_cut, th2_fpXhgcer);
    th2_fpXhgcer_eff2D->GetXaxix()->SetNameTitle("Focal Plane X");
    th2_fpXhgcer_eff2D->GetYaxix()->SetNameTitle("Focal Plane Y");
    th2_fpXhgcer_eff2D->SetStats(0);
    
    th2_fpXngcer_eff2D = new TH2D("fpVngcereff_2Deff", "fpVngcereff_2Deff", 80, -40.0, 40.0, 80, -40.0, 40.0);
    junk = th2_fpXngcer_eff2D->Divide(th2_fpXngcer_cut, th2_fpXngcer);
    th2_fpXngcer_eff2D->GetXaxix()->SetNameTitle("Focal Plane X");
    th2_fpXngcer_eff2D->GetYaxix()->SetNameTitle("Focal Plane Y");
    th2_fpXngcer_eff2D->SetStats(0);
    
    th2_fpXaero_eff2D = new TH2D("fpVaeroeff_2Deff", "fpVaeroeff_2Deff", 80, -40.0, 40.0, 80, -40.0, 40.0);
    junk = th2_fpXaero_eff2D->Divide(th2_fpXaero_cut, th2_fpXaero);
    th2_fpXaero_eff2D->GetXaxix()->SetNameTitle("Focal Plane X");
    th2_fpXaero_eff2D->GetYaxix()->SetNameTitle("Focal Plane Y");
    th2_fpXaero_eff2D->SetStats(0);
    
    cout << "Finished making plots, saving to pdf.\n";
    //make and print canvas output to pdf
    TCanvas *c1 = new TCanvas (Form("SHMS_%s_PID_Plots_%d", cutNames[cutType].c_str(), runNum), Form("SHMS_%s_PID_Plots_%d", cutNames[cutType].c_str(), runNum), 2400, 2400);
    c1->Divide(2,2);
    
    c1->cd(1);
    gPad->SetLogz();
    th2_aeroXhgcer->Draw("colz");

    c1->cd(2);
    gPad->SetLogz();
    th2_ngcerXcal->Draw("colz");

    c1->cd(3);
    gPad->SetLogz();
    th2_ngcerXhgcer->Draw("colz");

    c1->cd(4);
    gPad->SetLogz();
    th2_ngcerXaero->Draw("colz");
    
    TCanvas *c2 = new TCanvas (Form("SHMS_%s_EffvDelta_Plots_%d", cutNames[cutType].c_str(), runNum), Form("SHMS_%s_EffvDelta_Plots_%d", cutNames[cutType].c_str(), runNum), 2400, 2400);
    c2->Divide(1,3);
    c2->cd(1);
    th1_hgcer_eff->Draw("E0");
    
    c2->cd(2);
    th1_ngcer_eff->Draw("E0E1");
    
    c2->cd(3);
    th1_aero_eff->Draw("E0E2");
    
    TCanvas *c3 = new TCanvas (Form("SHMS_%s_2DEff_Plots_%d", cutNames[cutType].c_str(), runNum), Form("SHMS_%s_2DEff_Plots_%d", cutNames[cutType].c_str(), runNum), 2400, 2400);
    c3->Divide(1,3);
    c3->cd(1);
    th2_fpXhgcer_eff2D->Draw("colz");
    
    c3->cd(2);
    th2_fpXngcer_eff2D->Draw("colz");
    
    c3->cd(3);
    th2_fpXaero_eff2D->Draw("colz");
    
    c1->Print(Form("SHMS_%s_PIDeffPlots_%d.pdf(", cutNames[cutType].c_str(), runNum));
    c2->Print(Form("SHMS_%s_PIDeffPlots_%d.pdf", cutNames[cutType].c_str(), runNum));
    c3->Print(Form("SHMS_%s_PIDeffPlots_%d.pdf)", cutNames[cutType].c_str(), runNum));

    return;
}

void plotSHMSCerEfficiency (TString pathToRootFile, Int_t runNum, int NumEventsInput, int cutType)
{
    gROOT->SetBatch(1);
  	cout << "\n\n";

	cout << "Running Run: '"<<runNum<<"' for " << NumEventsInput << " Events\n";
	
	makePlots(pathToRootFile, runNum, NumEventsInput, cutType);
	
	return;
    
}



