/****************************************
    plotCombinedLumi.C
    Author: Nathan Heinrich

    Quick root script to plot the lumi combined over diff. settings.

*****************************************/

TH1D *resHist;

//sets graph points color and style to something
void PlotColors (TGraphErrors *g, int i)
{
    g->SetMarkerStyle(20+(i%4));
    if(i > 8){
        g->SetMarkerColor(i+21);
    }else{
        g->SetMarkerColor(i+1);
    }
    //g->GetXaxis()->SetTitle("Current (uA)");
    g->GetYaxis()->SetTitle("EDTM TLT");
}

void LineColors (TF1 *TFi, int i)
{
    //TFi->SetLineStyle(1+(i%10));
    if(i > 8){
        TFi->SetLineColor(i+21);
    }else{
        TFi->SetLineColor(i+1);
    }
}

//does fit and resets intercept to 1 at x=0 in linear fit
// Yes this causes a memory leak, I don't care
TGraphErrors* ReBase(TGraphErrors *ge1, TF1 *&TFi, bool rate)
{
    Double_t *reBase1, *reX, *reEX;
    TFi = new TF1("log1", "[0]+[1]*log(x)",0.001,5);
    reBase1 = ge1->GetY();
    int N = ge1->GetN();
    reX = ge1->GetX();
    reEX = ge1->GetEX();
    Double_t *reBase2 = new Double_t[N];
    Double_t *reX2 = new Double_t[N];
    Double_t *reEX2 = new Double_t[N];
    for(int i = 0; i < N; i++)
    {
        //reBase2[i]= reBase1[i] - (reIntercept - 1);
        if(rate){ // conver to kHz
            reX2[i] = reX[i]/1000;
            reEX2[i] = reEX[i]/1000;
        }else{
            reX2[i] = reX[i];
            reEX2[i] = reEX[i];
        }
    }
    TGraphErrors *ge3 = new TGraphErrors(N, reX2, reBase1, reEX2, ge1->GetEY());
    if(rate){
        ge3->GetXaxis()->SetTitle("Rate (kHz)");
    }else{
        ge3->GetXaxis()->SetTitle("Current (uA)");
    }
    PlotColors(ge3,0);
    ge3->SetTitle(ge1->GetTitle());
    ge3->SetName(ge1->GetName());
    ge3->Fit(TFi);
    ge3->Draw("AP");
    TLegend* l1 = new TLegend(0.5, 0.85, 0.9, 0.9);
    Double_t reIntercept = TFi->GetParameter(0);
    for(int i = 0; i < N; i++)
    {
        reBase2[i]= reBase1[i];// - (TFi->Eval(reX2[i]));
        resHist->Fill(reBase2[i]);
    }

    l1->AddEntry(TFi, Form("y = (%f #pm %f)log(x) + (%f #pm %f)",TFi->GetParameter(1), TFi->GetParError(1),TFi->GetParameter(0), TFi->GetParError(0)));
    l1->Draw();
    
    
    TGraphErrors *ge2 = new TGraphErrors(N, reX2, reBase2, reEX2, ge1->GetEY());
    return ge2;
}

void PlotCombinedTLT (TString dataType)
{
    gROOT->SetBatch(kTRUE);
    int NFILES;
    TString *filenames;
    TString OutFileName;
    dataType.ToLower();
    //cout << dataType;
    bool rate = false;
    
    // This is my nifty way of handling the amount of inputfiles, so that I need only write stuff once,
    // danamically allocated arrays, a marvel to behold 
    cout << "running TLT\n";
    NFILES = 7;
    filenames = new TString[NFILES];
    filenames[0] = "../OUTPUTS/ExclusiveLumi/CBlock_exc1.csv";
    filenames[1] = "../OUTPUTS/ExclusiveLumi/CBlock_exc2.csv";
    filenames[2] = "../OUTPUTS/ExclusiveLumi/CBlock_exc3.csv";
    filenames[3] = "../OUTPUTS/ExclusiveLumi/CBlock_exc4.csv";
    filenames[4] = "../OUTPUTS/ExclusiveLumi/CBlock_exc5.csv";
    filenames[5] = "../OUTPUTS/sidis/CBlock_sidis1.csv";
    filenames[6] = "../OUTPUTS/sidis/CBlock_sidis2.csv";
    OutFileName = "TLT";
    rate = true;
    

    //yes yes, needlessly complicated. At least I need only rewrite the above everytime the files change!
    TCanvas **C = new TCanvas*[NFILES]; 
    TGraphErrors **Gi = new TGraphErrors*[NFILES];
    TGraphErrors **Gf = new TGraphErrors*[NFILES];
    TF1 **TFi = new TF1*[NFILES];
    TMultiGraph *mg = new TMultiGraph(Form("%s Residuals Combined", OutFileName.Data()),Form("%s Data Combined", OutFileName.Data()));
    
    resHist = new TH1D("residuals","residuals",50,-0.03,0.02);

    for(int i = 0; i < NFILES; i++)
    {
        C[i] = new TCanvas(Form("c%d", i),Form("c%d", i),10, 10, 1000, 800);
        Gi[i] = new TGraphErrors(filenames[i], "%lg, %lg, %lg");
        Gi[i]->SetName(filenames[i]);
        Gi[i]->SetTitle(filenames[i]);
        Gf[i] = ReBase(Gi[i], TFi[i], rate);
        int lastSlash = filenames[i].Last('/');
        int firstDot = filenames[i].First('.');
        Gf[i]->SetName(filenames[i]);
        Gf[i]->SetTitle(filenames[i]);
        if (!rate)
            Gf[i]->GetXaxis()->SetTitle("Current (uA)");
        else
            Gf[i]->GetXaxis()->SetTitle("Rate (kHz)");
        //Gf[i]->SetName((filenames[i](lastSlash,firstDot-lastSlash)).Data());
        //Gf[i]->SetTitle((filenames[i](lastSlash,firstDot-lastSlash)).Data());
        PlotColors(Gf[i], i);
        LineColors(TFi[i],i);
        //Gi[i]->Draw();
        if(i == 0) C[i]->Print(Form("../OUTPUTS/CombinedPlot%s.pdf(", OutFileName.Data()));
        else C[i]->Print(Form("../OUTPUTS/CombinedPlot%s.pdf", OutFileName.Data()));
        mg->Add(Gf[i]);
    }
    
    TCanvas *cf = new TCanvas("cf","cf",10, 10, 1000, 800);
    if (!rate)
        mg->GetXaxis()->SetTitle("Current (uA)");
    else
        mg->GetXaxis()->SetTitle("Rate (kHz)");
        
    mg->GetYaxis()->SetTitle("Coin Blocking Correction");
    
    mg->Draw("AP");
    for(int i = 0; i < NFILES; i++)
    {
        //TFi[i]->Draw("SAME C");
    }

    //TLegend* l2 = new TLegend(0.5, 0.85, 0.9, 0.9);
    cf->Print(Form("../OUTPUTS/CombinedPlot%s.pdf", OutFileName.Data()));
    
    TLegend* l1 = new TLegend(0.1, 0.1, 0.4, 0.22);
    //l1->AddEntry(lin2, Form("y = (%f #pm %f)log(x) + (%f #pm %f)",lin2->GetParameter(1), lin2->GetParError(1),lin2->GetParameter(0), lin2->GetParError(0)));
    for(int i = 0; i < NFILES; i++)
    {
        l1->AddEntry(Gf[i], Gf[i]->GetTitle());
    }
    l1->Draw();
    
    
    cf->Print(Form("../OUTPUTS/CombinedPlot%s.pdf", OutFileName.Data()));

    TCanvas *cf2 = new TCanvas("cf2","cf2",10, 10, 1000, 800);
    resHist->Draw();
    cf2->Print(Form("../OUTPUTS/CombinedPlot%s.pdf)", OutFileName.Data()));
    return;
}








