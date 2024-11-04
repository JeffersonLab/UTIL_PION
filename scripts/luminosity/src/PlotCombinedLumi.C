/****************************************
    plotCombinedLumi.C
    Author: Nathan Heinrich

    Quick root script to plot the lumi combined over diff. settings.

*****************************************/

//sets graph points color and style to something
void PlotColors (TGraphErrors *g, int i)
{
    g->SetMarkerStyle(20+(i%4));
    g->SetMarkerColor(i+1);
    g->GetXaxis()->SetTitle("Current (uA)");
    g->GetYaxis()->SetTitle("Norm. Track Yield");
}

//does fit and resets intercept to 1 at x=0 in linear fit
// Yes this causes a memory leak, I don't care
TGraphErrors* ReBase(TGraphErrors *ge1)
{
    Double_t *reBase1;
    TF1 *lin1 = new TF1("lin1", "[0]+[1]*x",0,80);

    PlotColors(ge1,0);
    ge1->Fit(lin1);
    ge1->Draw("AP");
    
    TLegend* l1 = new TLegend(0.5, 0.7, 0.9, 0.9);
    Double_t reIntercept = lin1->GetParameter(0);

    l1->AddEntry(lin1, "y = mx + b");
    l1->AddEntry((TObject*)0, Form("m = %f #pm %f", lin1->GetParameter(1), lin1->GetParError(1)));
    l1->AddEntry((TObject*)0, Form("b = %f #pm %f", lin1->GetParameter(0), lin1->GetParError(0)));
    l1->Draw();
    
    reBase1 = ge1->GetY();
    int N = ge1->GetN();
    Double_t *reBase2 = new Double_t[N];
    for(int i = 0; i < N; i++)
    {
        reBase2[i]= reBase1[i] - (reIntercept - 1);
    }
    TGraphErrors *ge2 = new TGraphErrors(N, ge1->GetX(), reBase2, ge1->GetEX(), ge1->GetEY());
    return ge2;
}

void PlotCombinedLumi (TString dataType)
{
    int NFILES;
    TString *filenames;
    
    TString fileName1 = "../OUTPUTS/ExclusiveLumi/yield_data_exc1_clean.csv";
    TString fileName2 = "../OUTPUTS/ExclusiveLumi/yield_data_exc2_clean.csv";
    TString fileName3 = "../OUTPUTS/ExclusiveLumi/yield_data_exc3_clean.csv";
    TString fileName4 = "../OUTPUTS/ExclusiveLumi/yield_data_exc4_clean.csv";
    TString fileName5 = "../OUTPUTS/ExclusiveLumi/yield_data_exc5_clean.csv";
    TString fileName6 = "../OUTPUTS/sidis/yield_data_sidis1_clean.csv";
    TString fileName7 = "../OUTPUTS/sidis/yield_data_sidis2_clean.csv";
    NFILES = 7;
    
    filenames = new TString[NFILES];
    filenames[0] = fileName1;
    filenames[1] = fileName2;
    filenames[2] = fileName3;
    filenames[3] = fileName4;
    filenames[4] = fileName5;
    filenames[5] = fileName6;
    filenames[6] = fileName7;

    TString OutFileName = "Coin";
    TCanvas **C = new TCanvas*[NFILES]; 
    TGraphErrors **Gi = new TGraphErrors*[NFILES];
    TGraphErrors **Gf = new TGraphErrors*[NFILES];
    TMultiGraph *mg = new TMultiGraph(Form("%s Data Combined", OutFileName.Data()),Form("%s Data Combined", OutFileName.Data()));
    
    for(int i = 0; i < NFILES; i++)
    {
        C[i] = new TCanvas(Form("c%d", i),Form("c%d", i),10, 10, 1000, 800);
        Gi[i] = new TGraphErrors(filenames[i], "%lg, %lg, %lg");
        Gi[i]->SetName(filenames[i]);
        Gi[i]->SetTitle(filenames[i]);
        Gf[i] = ReBase(Gi[i]);
        int lastSlash = filenames[i].Last('/');
        int firstDot = filenames[i].First('.');
        Gf[i]->SetName(filenames[i]);
        Gf[i]->SetTitle(filenames[i]);
        //Gf[i]->SetName((filenames[i](lastSlash,firstDot-lastSlash)).Data());
        //Gf[i]->SetTitle((filenames[i](lastSlash,firstDot-lastSlash)).Data());
        PlotColors(Gf[i], i);
        
        if(i == 0) C[i]->Print(Form("../OUTPUTS/CombinedPlot%s.pdf(", OutFileName.Data()));
        else C[i]->Print(Form("../OUTPUTS/CombinedPlot%s.pdf", OutFileName.Data()));
        mg->Add(Gf[i]);
    }
    
    TCanvas *cf = new TCanvas("cf","cf",10, 10, 1000, 800);
    mg->GetXaxis()->SetTitle("Current (uA)");
    mg->GetYaxis()->SetTitle("Renormalized Yield");
    
    TF1 *lin2 = new TF1("lin2", "[0]+[1]*x",0,80);
    mg->Fit(lin2);
    mg->Draw("AP");
    
    TLegend* l1 = new TLegend(0.5, 0.7, 0.9, 0.9);
    l1->AddEntry(lin2, Form("y = (%f #pm %f)x + (%f #pm %f)",lin2->GetParameter(1), lin2->GetParError(1),lin2->GetParameter(0), lin2->GetParError(0)));
    for(int i = 0; i < NFILES; i++)
    {
        l1->AddEntry(Gf[i], Gf[i]->GetTitle());
    }
    l1->Draw();
    
    cf->Print(Form("../OUTPUTS/CombinedPlot%s.pdf)", OutFileName.Data()));
}








