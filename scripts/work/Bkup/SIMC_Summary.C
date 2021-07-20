/*
    SIMC_Summary.C 
    Author: Nathan Heinrich
    
    Code produces plots of variables from simc, as well as calculate the expected rate for running at beam current of 70uA, and total charge of the simulation is 1mC.
    All Histograms will be properly weighted (unless they're cut summaries). 
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
#include <fstream>

const Bool_t Debug_Flag = false;

//values we get from the tree
Float_t Weight, Em, Pm; 
Float_t hsdelta, hsxptar, hsyptar, ssdelta, ssxptar, ssyptar; // acceptance variables

// The range for integrating the missing mass peak
const Double_t MMIntegralLow  = -0.01;
const Double_t MMIntegralHigh = 0.01;

//acceptance cuts
const Float_t hsdeltaCutHigh =  8.0;
const Float_t hsdeltaCutLow  = -8.0;

const Float_t hsxptarCutHigh =  0.080;
const Float_t hsxptarCutLow  = -0.080;

const Float_t hsyptarCutHigh =  0.035;
const Float_t hsyptarCutLow  = -0.035;

const Float_t ssdeltaCutHigh =  15;
const Float_t ssdeltaCutLow  = -15;

const Float_t ssxptarCutHigh =  0.040;
const Float_t ssxptarCutLow  = -0.040;

const Float_t ssyptarCutHigh =  0.0240;
const Float_t ssyptarCutLow  = -0.0240;

const Double_t BeamCurrent = 70.0; //   uA
const Double_t TotCharge = 1.0; //      mC

// Converts the line that I get from the 
Double_t GetNormFac(TString normfacString)
{

    if(Debug_Flag) cout << normfacString << endl;
    
    Int_t StartIndex = 0;
    
    //find the first number in the string
    /*for(Int_t i = 0; i < normfacString.Length() && StartIndex == 0; i++)
    {
        if (normfacString[i] > 47 && normfacString < 58) // if ascii is between 0 (ascii value is 48) and 9 (ascii value is 57) 
        {
            StartIndex = i;
        }
    }*/
    StartIndex = normfacString.First('0'); // fortran sci notation should always start with 0
    
    Double_t normfac = 0;
    
    if(Debug_Flag) cout << normfacString[StartIndex] << endl;
    // I'm just going to hard code this because the output is always in the form: 0.123456E+01
    normfac += (normfacString[StartIndex]-'0');         // should be 0 always
    normfac += (normfacString[StartIndex+2]-'0')/10.0;  // get tenths
    normfac += (normfacString[StartIndex+3]-'0')/100.0; // get hundreths
    normfac += (normfacString[StartIndex+4]-'0')/1000.0;// you get this now right?
    normfac += (normfacString[StartIndex+5]-'0')/10000.0;
    normfac += (normfacString[StartIndex+6]-'0')/100000.0;
    normfac += (normfacString[StartIndex+7]-'0')/1000000.0;
    
    // now we need to get the exponent
    Int_t exponent = 0;
    exponent += (normfacString[StartIndex+10]-'0')*10; // tens place
    exponent += (normfacString[StartIndex+11]-'0');   // ones place
    // record if need to flip sign
    if (normfacString[StartIndex+9] == '-')
    {
        exponent = -1*exponent;
    }
    
    return normfac*TMath::Power(10.0, exponent);
}

//do not include the file extension of input. 
//spectrometer input is: 'C' for coin, 'S' for SHMS, 'H' for HMS, 'N' for no cuts
void SIMC_Summary (TString Filename, TString SPEC_FLAG) 
{

    //root file, for getting the data
    TString ROOTFileName = "worksim/"+Filename+".root";
    //Sim input file could be used for ngen and total charge 
    //TString InfileName = "infiles/"+Filename+".inp"; //shouldn't need this
    // Out file, need to get normfac
    TString OutFileName = "outfiles/"+Filename+".hist";
    
    //start by tring to open ROOTfile
    TFile *RootFile = new TFile(ROOTFileName, "READ");
    if (!RootFile->IsOpen())
    {
        cout << "Root File failed to open!!!\nPlease Check that '"+ROOTFileName+"' exists\n---Ending---\n";
        return;
    }
    
    TTree *DataTree = dynamic_cast <TTree*> (RootFile->Get("h10")); //get h10 tree from root file
    
    //Assign all of the desired tree values to variables 
    DataTree->SetBranchAddress("Weight", &Weight);
    DataTree->SetBranchAddress("Em", &Em); // missing energy
    DataTree->SetBranchAddress("Pm", &Pm); // missing momentum
    
    //acceptance variables
    DataTree->SetBranchAddress("hsdelta", &hsdelta);
    DataTree->SetBranchAddress("hsxptar", &hsxptar);
    DataTree->SetBranchAddress("hsyptar", &hsyptar);
    DataTree->SetBranchAddress("ssdelta", &ssdelta);
    DataTree->SetBranchAddress("ssxptar", &ssxptar);
    DataTree->SetBranchAddress("ssyptar", &ssyptar);
    
    TCanvas *MMCanvas = new TCanvas(Filename, "Missing Mass "+Filename, 1200, 800);
    TH1D *MMHist = new TH1D ("MMHist", "Missing Mass", 100, -1, 1);
    TH1D *MMHistIntegral = new TH1D ("MMHistIntegral", "Integrated Missing Mass", 100, -1, 1);
    
    //acceptance cuts
    Bool_t hsdeltaCut, hsxptarCut, hsyptarCut, ssdeltaCut, ssxptarCut, ssyptarCut;
    Bool_t HMSAcceptanceCut, SHMSAcceptanceCut;
    
    // Loop over all events in the Root file
    Int_t nEntries = DataTree->GetEntries();
    for(Int_t iEntry = 0; iEntry < nEntries; iEntry++) // event loop
    {
        DataTree->GetEntry(iEntry);
        
        //acceptance cuts
        //Bool_t AcceptanceCuts = true; //for testing
        hsdeltaCut = (hsdelta < hsdeltaCutHigh && hsdelta > hsdeltaCutLow); 
        hsxptarCut = (hsxptar < hsxptarCutHigh && hsxptar > hsxptarCutLow); 
        hsyptarCut = (hsyptar < hsyptarCutHigh && hsyptar > hsyptarCutLow);
        
        ssdeltaCut = (ssdelta < ssdeltaCutHigh && ssdelta > ssdeltaCutLow);
        ssxptarCut = (ssxptar < ssxptarCutHigh && ssxptar > ssxptarCutLow);
        ssyptarCut = (ssyptar < ssyptarCutHigh && ssyptar > ssyptarCutLow);
        
        //Statement for picking which 
        if ( SPEC_FLAG == "C"){
            HMSAcceptanceCut = hsdeltaCut && hsxptarCut && hsyptarCut;
            SHMSAcceptanceCut = ssdeltaCut && ssxptarCut && ssyptarCut;
        } else if ( SPEC_FLAG == "H") {
            HMSAcceptanceCut = hsdeltaCut && hsxptarCut && hsyptarCut;
            SHMSAcceptanceCut = true;
        } else if ( SPEC_FLAG == "S") {
            HMSAcceptanceCut = true;
            SHMSAcceptanceCut = ssdeltaCut && ssxptarCut && ssyptarCut;
        } else {
            HMSAcceptanceCut = true;
            SHMSAcceptanceCut = true;
        }
        
        
        Bool_t AcceptanceCuts = HMSAcceptanceCut && SHMSAcceptanceCut;
        
        Double_t MM = Em*Em-Pm*Pm; // calculate missing mass
        
        //cut on acceptance variables
        if(AcceptanceCuts)
        {
            
            MMHist->Fill(MM, Weight); // fill weigthed histogram for missing mass
            
            // create a histogram that will be filled in on plot to represent the integral
            if(MM < MMIntegralHigh && MM > MMIntegralLow)
            {
                MMHistIntegral->Fill(MM, Weight);
            }
        }
    } // end Event loop
    
    Double_t counts = MMHistIntegral->Integral(); // get the total counts under the peak  
    
    ifstream OutFile;
    OutFile.open(OutFileName);
    
    if (!OutFile.good())
    {
        cout << ".hist file not found!!!\nPlease Check that '"+OutFileName+"' exists\n---Ending---\n";
        return;
    }
    
    TString normfacString;
    Double_t normfac;
    Bool_t normfacFound = false;
    
    //find the line that contains "normfac"
    while (!OutFile.eof() && !normfacFound)
    {
        char temp[256];
        OutFile.getline(temp, 256, '\n');
        
        normfacString = temp;
        normfacFound = normfacString.Contains("normfac");
    }
    
    if(!normfacFound)
    {
        cout << "File read failed!!!\nDid not find 'normfac' in: "+OutFileName+"\n---Ending---\n";
        return;
    }
    
    
    //get the value out of the string.
    normfac = GetNormFac(normfacString);
    
    // output all the relevent quantaties
    cout << "counts: " << counts << "\n";
    cout << "ngen: " << nEntries << "\n"; //ngen should always be the total number of events in the ROOTfile
    cout << "normfac: " << normfac << endl;
    
    Double_t BeamTime = TotCharge/(BeamCurrent/1000); // factor of 1000 to convert from uA to mA
    cout << "time: " << BeamTime << " s\n";
    
    Double_t Rate = counts*normfac/(nEntries*BeamTime); //should be the count rate in Hz
    cout << "Count rate: " << Rate << " Hz\n";
    
    //draw histograms
    MMHist->Draw("BAR");
    MMHistIntegral->SetFillStyle(3144);
    MMHistIntegral->SetFillColor(kBlue);
    MMHistIntegral->Draw("SAME BAR");
    
    // make legend
    TLegend *legend = new TLegend(0.1, 0.75, 0.30, 0.9);
    legend->AddEntry("MMHistIntegral", "Intgrated Peak", "f");
    legend->AddEntry((TObject*)0, Form("Value of integral: %f", counts), "");
    
    legend->Draw();
    //MMCanvas->WriteObject();
    
    //close files
    //OutFile.close();
    //RootFile->Close();
    return;
}

