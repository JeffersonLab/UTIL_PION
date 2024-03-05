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
Double_t xfp, yfp, xpfp, ypfp;

//cuts
const Double_t calEtotLowPi = 0.1; //normaized energy
const Double_t calEtotLowe = 0.7; //normaized energy
const Double_t hgcerNpeSumLow = 1.5; //unit NPE
const Double_t aeroNpeSumLow = 1.5; //unit NPE
const Double_t ngcerNpeSumLow = 1.5; //unit NPE

//z pos of cerenkovs - in cm
const Double_t HGCER_z = 156.27;
const Double_t NGCER_z = -89.1;
const Double_t AERO_z = 231.0;

Bool_t calCut, hgcerCut, aeroCut, ngcerCut, fpcut;

const Int_t INILENGTH = 64;
Int_t NumEvents = -1;


void makePlots ( TString rootFile, Int_t runNum, int NumEvents, int cutType ) 
{
    gStyle->SetOptTitle(0);
    gStyle->SetLabelSize(0.04,"X");
    gStyle->SetLabelSize(0.04,"Y");
    gStyle->SetLabelSize(0.04,"Z");
    // make empty histograms	
	th1_hgcer = new TH1D("hgcerShould", "hgcerShould", 120, -10.0, 20.0);
	th1_hgcerCut = new TH1D("hgcerDid", "hgcerDid", 120, -10.0, 20.0);
	th1_aero = new TH1D("aeroShould", "aeroShould", 120, -10.0, 20.0);
	th1_aeroCut = new TH1D("aeroDid", "aeroDid", 120, -10.0, 20.0);
	th1_ngcer = new TH1D("aeroShould", "aeroShould", 120, -10.0, 20.0);
	th1_ngcerCut = new TH1D("aeroDid", "aeroDid", 120, -10.0, 20.0);
	
	
	th2_aeroXhgcer = new TH2D("aeroNpeSumVhgcerNpeSum","aeroNpeSumVhgcerNpeSum", 100, 0.0, 35, 100, 0.0, 35);
	th2_aeroXhgcer->GetXaxis()->SetNameTitle("Aerogel NPE","Aerogel NPE");
	th2_aeroXhgcer->GetYaxis()->SetNameTitle("HGC NPE","HGC NPE");
	th2_aeroXhgcer->GetXaxis()->SetLabelSize(0.04);
	th2_aeroXhgcer->GetYaxis()->SetLabelSize(0.04);
	th2_aeroXhgcer->SetStats(0);
	
	th2_ngcerXcal = new TH2D("ngcerNpeSumVP.cal.etottracknorm","ngcerNpeSumVP.cal.etottracknorm", 100, 0.0, 35, 100, 0.0, 1.6);
	th2_ngcerXcal->GetXaxis()->SetNameTitle("NGC NPE","NGC NPE");
	th2_ngcerXcal->GetYaxis()->SetNameTitle("Normalized Calorimeter Energy","Normalized Calorimeter Energy");
	th2_ngcerXcal->GetXaxis()->SetLabelSize(0.04);
	th2_ngcerXcal->GetYaxis()->SetLabelSize(0.04);
	th2_ngcerXcal->SetStats(0);
	
	th2_ngcerXaero = new TH2D("ngcerNpeSumVaeroNpeSum","ngcerNpeSumVaeroNpeSum", 100, 0.0, 35, 100, 0.0, 35);
	th2_ngcerXaero->GetXaxis()->SetNameTitle("NGC NPE","NGC NPE");
	th2_ngcerXaero->GetYaxis()->SetNameTitle("Aerogel NPE","Aerogel NPE");
	th2_ngcerXaero->GetXaxis()->SetLabelSize(0.04);
	th2_ngcerXaero->GetYaxis()->SetLabelSize(0.04);
	th2_ngcerXaero->SetStats(0);
	
	th2_ngcerXhgcer = new TH2D("ngcerNpeSumVhgcerNpeSum","ngcerNpeSumVhgcerNpeSum", 100, 0.0, 35, 100, 0.0, 35);
	th2_ngcerXhgcer->GetXaxis()->SetNameTitle("NGC NPE","NGC NPE");
	th2_ngcerXhgcer->GetYaxis()->SetNameTitle("HGC NPE","HGC NPE");
	th2_ngcerXhgcer->GetXaxis()->SetLabelSize(0.04);
	th2_ngcerXhgcer->GetYaxis()->SetLabelSize(0.04);
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
	tree1->SetBranchAddress("P.dc.xp_fp", &xpfp);
	tree1->SetBranchAddress("P.dc.yp_fp", &ypfp);
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
            
            // Fill 2D PID plots
            if(fpcut){
                th2_aeroXhgcer->Fill(aeroNpeSum, hgcerNpeSum);
                th2_ngcerXaero->Fill(ngcerNpeSum, aeroNpeSum);
                th2_ngcerXhgcer->Fill(ngcerNpeSum, hgcerNpeSum);
            }
            
            if (fpcut & hgcerCut & aeroCut) th2_ngcerXcal->Fill(ngcerNpeSum, calEtot);
            
            if (calCut & aeroCut & ngcerCut & fpcut) // should for hgcer
            {
                th1_hgcer->Fill(delta);
                th2_fpXhgcer->Fill((xfp+xpfp*HGCER_z),(yfp+ypfp*HGCER_z));
                if (hgcerCut)
                {
                    th1_hgcerCut->Fill(delta);
                    th2_fpXhgcer_cut->Fill((xfp+xpfp*HGCER_z),(yfp+ypfp*HGCER_z)); // change to fpAthgcer
                }
            }
            if (calCut & aeroCut & hgcerCut & fpcut) // should for ngcer
            {
                th1_ngcer->Fill(delta);
                th2_fpXngcer->Fill((xfp+xpfp*NGCER_z),(yfp+ypfp*NGCER_z));
                if (ngcerCut)
                {
                    th1_ngcerCut->Fill(delta);
                    th2_fpXngcer_cut->Fill((xfp+xpfp*NGCER_z),(yfp+ypfp*NGCER_z));
                }
            }
            if(calCut & hgcerCut & ngcerCut & fpcut) // should for aero
            {
                th1_aero->Fill(delta);
                th2_fpXaero->Fill((xfp+xpfp*AERO_z),(yfp+ypfp*AERO_z));
                if(aeroCut)
                {
                    th1_aeroCut->Fill(delta);
                    th2_fpXaero_cut->Fill((xfp+xpfp*AERO_z),(yfp+ypfp*AERO_z));
                }
            }   
        }	
	}
	
	// do division of plots to get efficiency plots
	Bool_t junk; // for holding return of TH1->Divide() 
	th1_hgcer_eff = new TH1D("hgcer_eff", "hgcer_eff", 120, -10.0, 20.0);
	junk = th1_hgcer_eff->Divide(th1_hgcerCut,th1_hgcer);
	th1_hgcer_eff->GetXaxis()->SetNameTitle("#delta (%)","#delta (%)");
	th1_hgcer_eff->GetYaxis()->SetNameTitle("HGCER Efficiency","HGCER Efficiency");
	th1_hgcer_eff->GetXaxis()->SetLabelSize(0.04);
	th1_hgcer_eff->GetYaxis()->SetLabelSize(0.04);
	th1_hgcer_eff->SetStats(0);
	
	th1_aero_eff = new TH1D("aero_eff", "aero_eff", 120, -10.0, 20.0);
	junk = th1_aero_eff->Divide(th1_aeroCut,th1_aero);
	th1_aero_eff->GetXaxis()->SetNameTitle("#delta (%)","#delta (%)");
	th1_aero_eff->GetYaxis()->SetNameTitle("Aerogel Efficiency","Aerogel Efficiency");
	th1_aero_eff->GetXaxis()->SetLabelSize(0.04);
	th1_aero_eff->GetYaxis()->SetLabelSize(0.04);
	th1_aero_eff->SetStats(0);
	
	th1_ngcer_eff = new TH1D("ngcer_eff", "ngcer_eff", 120, -10.0, 20.0);
	junk = th1_ngcer_eff->Divide(th1_ngcerCut,th1_ngcer);
	th1_ngcer_eff->GetXaxis()->SetNameTitle("#delta (%)","#delta (%)");
	th1_ngcer_eff->GetYaxis()->SetNameTitle("NGCER Efficiency","NGCER Efficiency");
	th1_ngcer_eff->GetXaxis()->SetLabelSize(0.04);
	th1_ngcer_eff->GetYaxis()->SetLabelSize(0.04);
	th1_ngcer_eff->SetStats(0);
	
    
    //calculate binamial errors manually
    Int_t BINS = th1_hgcer_eff->GetNbinsX();
    Double_t should, did, err;
    for(int i = 0; i < BINS; i++)
    {
        should = th1_hgcer->GetBinContent(i);
        did = th1_hgcerCut->GetBinContent(i);
        if (should != 0)
        {
            err = TMath::Sqrt(((did*should) - (did*did))/(should*should*should));
        }else{
            err = 0;
        }
        
        //cout << "Error in bin: " << th1_hgcer_eff->GetBinError(i) << ", Calculated Error: " << err << '\n';
        
        th1_hgcer_eff->SetBinError(i , err);
    }
    
    BINS = th1_ngcer_eff->GetNbinsX();
    for(int i = 0; i < BINS; i++)
    {
        should = th1_ngcer->GetBinContent(i);
        did = th1_ngcerCut->GetBinContent(i);
        if (should != 0)
        {
            err = TMath::Sqrt(((did*should) - (did*did))/(should*should*should));
        }else{
            err = 0;
        }
        
        //cout << "Error in bin: " << th1_ngcer_eff->GetBinError(i) << ", Calculated Error: " << err << '\n';
        
        th1_ngcer_eff->SetBinError(i , err);
    }
    
    BINS = th1_aero_eff->GetNbinsX();
    for(int i = 0; i < BINS; i++)
    {
        should = th1_aero->GetBinContent(i);
        did = th1_aeroCut->GetBinContent(i);
        if (should != 0)
        {
            err = TMath::Sqrt(((did*should) - (did*did))/(should*should*should));
        }else{
            err = 0;
        }
        
        //cout << "Error in bin: " << th1_aero_eff->GetBinError(i) << ", Calculated Error: " << err << '\n';
        
        th1_aero_eff->SetBinError(i , err);
    }

    //auto set axis ranges
    cout << "th1_hgcer_eff->GetBinContent(th1_hgcer_eff->GetMaximumBin())" << th1_hgcer_eff->GetBinContent(th1_hgcer_eff->GetMaximumBin()) << '\n';
    th1_hgcer_eff->SetMinimum(0.8*(th1_hgcer_eff->GetBinContent(th1_hgcer_eff->GetMinimumBin())));
	th1_hgcer_eff->SetMaximum(1.2*(th1_hgcer_eff->GetBinContent(th1_hgcer_eff->GetMaximumBin())));
	th1_ngcer_eff->SetMinimum(0.8*(th1_ngcer_eff->GetBinContent(th1_ngcer_eff->GetMinimumBin())));
	th1_ngcer_eff->SetMaximum(1.2*(th1_ngcer_eff->GetBinContent(th1_ngcer_eff->GetMaximumBin())));
	th1_aero_eff->SetMinimum(0.8*(th1_aero_eff->GetBinContent(th1_aero_eff->GetMinimumBin())));
	th1_aero_eff->SetMaximum(1.2*(th1_aero_eff->GetBinContent(th1_aero_eff->GetMaximumBin())));    
    
    th2_fpXhgcer_eff2D = new TH2D("fpVhgcereff_2Deff", "fpVhgcereff_2Deff", 80, -40.0, 40.0, 80, -40.0, 40.0);
    junk = th2_fpXhgcer_eff2D->Divide(th2_fpXhgcer_cut, th2_fpXhgcer);
    th2_fpXhgcer_eff2D->GetXaxis()->SetNameTitle("X at HGCER (cm)","X at HGCER (cm)");
    th2_fpXhgcer_eff2D->GetYaxis()->SetNameTitle("Y at HGCER (cm)","Y at HGCER (cm)");
    th2_fpXhgcer_eff2D->GetYaxis()->SetLabelSize(0.04);
    th2_fpXhgcer_eff2D->GetXaxis()->SetLabelSize(0.04);
    th2_fpXhgcer_eff2D->SetStats(0);
    
    th2_fpXngcer_eff2D = new TH2D("fpVngcereff_2Deff", "fpVngcereff_2Deff", 80, -40.0, 40.0, 80, -40.0, 40.0);
    junk = th2_fpXngcer_eff2D->Divide(th2_fpXngcer_cut, th2_fpXngcer);
    th2_fpXngcer_eff2D->GetXaxis()->SetNameTitle("X at NGCER (cm)","X at NGCER (cm)");
    th2_fpXngcer_eff2D->GetYaxis()->SetNameTitle("Y at NGCER (cm)","Y at NGCER (cm)");
    th2_fpXngcer_eff2D->GetYaxis()->SetLabelSize(0.04);
    th2_fpXngcer_eff2D->GetXaxis()->SetLabelSize(0.04);
    th2_fpXngcer_eff2D->SetStats(0);
    
    th2_fpXaero_eff2D = new TH2D("fpVaeroeff_2Deff", "fpVaeroeff_2Deff", 80, -40.0, 40.0, 80, -40.0, 40.0);
    junk = th2_fpXaero_eff2D->Divide(th2_fpXaero_cut, th2_fpXaero);
    th2_fpXaero_eff2D->GetXaxis()->SetNameTitle("X at Aerogel (cm)","X at Aerogel (cm)");
    th2_fpXaero_eff2D->GetYaxis()->SetNameTitle("Y at Aerogel (cm)","Y at Aerogel (cm)");
    th2_fpXaero_eff2D->GetYaxis()->SetLabelSize(0.04);
    th2_fpXaero_eff2D->GetXaxis()->SetLabelSize(0.04);
    th2_fpXaero_eff2D->SetStats(0);
    
    cout << "Finished making plots, saving to pdf.\n";
    //make and print canvas output to pdf
    TCanvas *c1_1 = new TCanvas (Form("SHMS_%s_PID_Plots_%d_1", cutNames[cutType].c_str(), runNum), Form("SHMS_%s_PID_Plots_%d", cutNames[cutType].c_str(), runNum), 600, 600);
    c1_1->SetMargin(0.15,0.15,0.15,0.15);
    gPad->SetLogz();
    th2_aeroXhgcer->Draw("colz");

    TCanvas *c1_2 = new TCanvas (Form("SHMS_%s_PID_Plots_%d_2", cutNames[cutType].c_str(), runNum), Form("SHMS_%s_PID_Plots_%d", cutNames[cutType].c_str(), runNum), 600, 600);
    c1_2->SetMargin(0.15,0.15,0.15,0.15);
    gPad->SetLogz();
    th2_ngcerXcal->Draw("colz");

    TCanvas *c1_3 = new TCanvas (Form("SHMS_%s_PID_Plots_%d_3", cutNames[cutType].c_str(), runNum), Form("SHMS_%s_PID_Plots_%d", cutNames[cutType].c_str(), runNum), 600, 600);
    c1_3->SetMargin(0.15,0.15,0.15,0.15);
    gPad->SetLogz();
    th2_ngcerXhgcer->Draw("colz");

    TCanvas *c1_4 = new TCanvas (Form("SHMS_%s_PID_Plots_%d_4", cutNames[cutType].c_str(), runNum), Form("SHMS_%s_PID_Plots_%d", cutNames[cutType].c_str(), runNum), 600, 600);
    c1_4->SetMargin(0.15,0.15,0.15,0.15);
    gPad->SetLogz();
    th2_ngcerXaero->Draw("colz");
    
    TCanvas *c2_1 = new TCanvas (Form("SHMS_%s_Effv#delta_Plots_%d_1", cutNames[cutType].c_str(), runNum), Form("SHMS_%s_Effv#delta_Plots_%d", cutNames[cutType].c_str(), runNum), 600, 600);
    c2_1->SetMargin(0.15,0.15,0.15,0.15);
    th1_hgcer_eff->Draw("E0");
    
    TCanvas *c2_2 = new TCanvas (Form("SHMS_%s_Effv#delta_Plots_%d_2", cutNames[cutType].c_str(), runNum), Form("SHMS_%s_Effv#delta_Plots_%d", cutNames[cutType].c_str(), runNum), 600, 600);
    c2_2->SetMargin(0.15,0.15,0.15,0.15);
    th1_ngcer_eff->Draw("E0");
    
    TCanvas *c2_3 = new TCanvas (Form("SHMS_%s_Effv#delta_Plots_%d_3", cutNames[cutType].c_str(), runNum), Form("SHMS_%s_Effv#delta_Plots_%d", cutNames[cutType].c_str(), runNum), 600, 600);
    c2_3->SetMargin(0.15,0.15,0.15,0.15);
    th1_aero_eff->Draw("E0");
    
    TCanvas *c3_1 = new TCanvas (Form("SHMS_%s_2DEff_Plots_%d_1", cutNames[cutType].c_str(), runNum), Form("SHMS_%s_2DEff_Plots_%d", cutNames[cutType].c_str(), runNum), 600, 600);
    c2_1->SetMargin(0.15,0.15,0.15,0.15);
    th2_fpXhgcer_eff2D->Draw("colz");
    
    TCanvas *c3_2 = new TCanvas (Form("SHMS_%s_2DEff_Plots_%d_2", cutNames[cutType].c_str(), runNum), Form("SHMS_%s_2DEff_Plots_%d", cutNames[cutType].c_str(), runNum), 600, 600);
    c3_2->SetMargin(0.15,0.15,0.15,0.15);
    th2_fpXngcer_eff2D->Draw("colz");
    
    TCanvas *c3_3 = new TCanvas (Form("SHMS_%s_2DEff_Plots_%d_3", cutNames[cutType].c_str(), runNum), Form("SHMS_%s_2DEff_Plots_%d", cutNames[cutType].c_str(), runNum), 600, 600);
    c3_3->SetMargin(0.15,0.15,0.15,0.15);
    th2_fpXaero_eff2D->Draw("colz");
    
    c1_1->Print(Form("SHMS_%s_PIDeffPlots_%d.pdf(", cutNames[cutType].c_str(), runNum));
    c1_2->Print(Form("SHMS_%s_PIDeffPlots_%d.pdf", cutNames[cutType].c_str(), runNum));
    c1_3->Print(Form("SHMS_%s_PIDeffPlots_%d.pdf", cutNames[cutType].c_str(), runNum));
    c1_4->Print(Form("SHMS_%s_PIDeffPlots_%d.pdf", cutNames[cutType].c_str(), runNum));
    c2_1->Print(Form("SHMS_%s_PIDeffPlots_%d.pdf", cutNames[cutType].c_str(), runNum));
    c2_2->Print(Form("SHMS_%s_PIDeffPlots_%d.pdf", cutNames[cutType].c_str(), runNum));
    c2_3->Print(Form("SHMS_%s_PIDeffPlots_%d.pdf", cutNames[cutType].c_str(), runNum));
    c3_1->Print(Form("SHMS_%s_PIDeffPlots_%d.pdf", cutNames[cutType].c_str(), runNum));
    c3_2->Print(Form("SHMS_%s_PIDeffPlots_%d.pdf", cutNames[cutType].c_str(), runNum));
    c3_3->Print(Form("SHMS_%s_PIDeffPlots_%d.pdf)", cutNames[cutType].c_str(), runNum));

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




