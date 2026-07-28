/****************************************
    plotCombinedLumi.C
    Author: Nathan Heinrich

    Quick root script to plot the lumi combined over diff. settings.

*****************************************/

//sets graph points color and style to something
void PlotColors (TGraphErrors *g, int i)
{
    g->SetMarkerStyle(20+(i%4));
    if(i > 8){
        g->SetMarkerColor(i+21);
    }else{
        g->SetMarkerColor(i+1);
    }
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
    
    TLegend* l1 = new TLegend(0.5, 0.85, 0.9, 0.9);
    Double_t reIntercept = lin1->GetParameter(0);

    l1->AddEntry(lin1, Form("y = (%f #pm %f)x + (%f #pm %f)",lin1->GetParameter(1), lin1->GetParError(1),lin1->GetParameter(0), lin1->GetParError(0)));
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

void PlotCombinedLumi_Scalar (TString dataType)
{
    gROOT->SetBatch(kTRUE);
    int NFILES;
    TString *filenames;
    TString OutFileName;
    dataType.ToLower();
    //cout << dataType;
    TString coin = "coin";
    TString coin_rate = "coin_rate";
    TString carbon = "carbon";
    TString singles = "singles";
    TString singles_rate = "singles_rate";
    TString shms = "shms";
    TString hms = "hms";
    TString shms_rate = "shms_rate";
    TString hms_rate = "hms_rate";
    TString lh2_all = "lh2_all";
    bool rate = false;
    
    // This is my nifty way of handling the amount of inputfiles, so that I need only write stuff once,
    // danamically allocated arrays, a marvel to behold 
    if (!dataType.CompareTo(coin) ) {
        cout << "running Coin\n";
        NFILES = 7;
        filenames = new TString[NFILES];
        filenames[0] = "../OUTPUTS/ExclusiveLumi/yield_data_exc1_ScYvCur.csv";
        filenames[1] = "../OUTPUTS/ExclusiveLumi/yield_data_exc2_ScYvCur.csv";
        filenames[2] = "../OUTPUTS/ExclusiveLumi/yield_data_exc3_ScYvCur.csv";
        filenames[3] = "../OUTPUTS/ExclusiveLumi/yield_data_exc4_ScYvCur.csv";
        filenames[4] = "../OUTPUTS/ExclusiveLumi/yield_data_exc5_ScYvCur.csv";
        filenames[5] = "../OUTPUTS/sidis/yield_data_sidis1_ScvCur.csv";
        filenames[6] = "../OUTPUTS/sidis/yield_data_sidis2_ScvCur.csv";
        OutFileName = "Coin";
    } else if (!dataType.CompareTo(coin_rate) ) {
        cout << "running Coin vrs Rate\n";
        rate = true;
        NFILES = 7;
        filenames = new TString[NFILES];
        filenames[0] = "../OUTPUTS/ExclusiveLumi/yield_data_exc1_ScYvRate.csv";
        filenames[1] = "../OUTPUTS/ExclusiveLumi/yield_data_exc2_ScYvRate.csv";
        filenames[2] = "../OUTPUTS/ExclusiveLumi/yield_data_exc3_ScYvRate.csv";
        filenames[3] = "../OUTPUTS/ExclusiveLumi/yield_data_exc4_ScYvRate.csv";
        filenames[4] = "../OUTPUTS/ExclusiveLumi/yield_data_exc5_ScYvRate.csv";
        filenames[5] = "../OUTPUTS/sidis/yield_data_sidis1_ScvRate.csv";
        filenames[6] = "../OUTPUTS/sidis/yield_data_sidis2_ScvRate.csv";
        OutFileName = "CoinRate";
    }else if (!dataType.CompareTo(carbon)){
        cout << "running Carbon\n";
        NFILES = 6;
        filenames = new TString[NFILES];
        filenames[0] = "../OUTPUTS/Lumi6-4/SHMS/yield_data_Carbon_clean.csv";
        filenames[1] = "../OUTPUTS/Lumi6-4/HMS/yield_data_Carbon_clean.csv";
        filenames[2] = "../OUTPUTS/Lumi9-2/pt1/yield_data_Carbon_HMS.csv";
        filenames[3] = "../OUTPUTS/Lumi9-2/pt1/yield_data_Carbon_SHMS.csv";
        filenames[4] = "../OUTPUTS/Lumi9-2/pt2/yield_data_Carbon_HMS.csv";
        filenames[5] = "../OUTPUTS/Lumi9-2/pt2/yield_data_Carbon_SHMS.csv";
        OutFileName = "Carbon";
    }else if (!dataType.CompareTo(singles)){
        cout << "running Singles\n";
        NFILES = 6;
        filenames = new TString[NFILES];
        filenames[0] = "../OUTPUTS/Lumi6-4/SHMS/yield_data_LH2_Sc.csv";
        filenames[1] = "../OUTPUTS/Lumi6-4/HMS/yield_data_LH2_Sc.csv";
        filenames[2] = "../OUTPUTS/Lumi9-2/pt1/yield_data_LH2_HMS_Sc.csv";
        filenames[3] = "../OUTPUTS/Lumi9-2/pt1/yield_data_LH2_SHMS_Sc.csv";
        filenames[4] = "../OUTPUTS/Lumi9-2/pt2/yield_data_LH2_HMS_Sc.csv";
        filenames[5] = "../OUTPUTS/Lumi9-2/pt2/yield_data_LH2_SHMS_Sc.csv";
        OutFileName = "LH2_Singles";
    }else if (!dataType.CompareTo(singles_rate)){
        cout << "running Singles\n";
        NFILES = 6;
        rate = true;
        filenames = new TString[NFILES];
        filenames[0] = "../OUTPUTS/Lumi6-4/SHMS/yield_data_LH2_ScVr.csv";
        filenames[1] = "../OUTPUTS/Lumi6-4/HMS/yield_data_LH2_ScVr.csv";
        filenames[2] = "../OUTPUTS/Lumi9-2/pt1/yield_data_LH2_HMS_ScVr.csv";
        filenames[3] = "../OUTPUTS/Lumi9-2/pt1/yield_data_LH2_SHMS_ScVr.csv";
        filenames[4] = "../OUTPUTS/Lumi9-2/pt2/yield_data_LH2_HMS_ScVr.csv";
        filenames[5] = "../OUTPUTS/Lumi9-2/pt2/yield_data_LH2_SHMS_ScVr.csv";
        OutFileName = "LH2_Singles_Rate";
    }else if (!dataType.CompareTo(shms)){
        cout << "running Singles\n";
        NFILES = 3;
        filenames = new TString[NFILES];
        filenames[0] = "../OUTPUTS/Lumi6-4/SHMS/yield_data_LH2_Sc.csv";
        filenames[1] = "../OUTPUTS/Lumi9-2/pt1/yield_data_LH2_SHMS_Sc.csv";
        filenames[2] = "../OUTPUTS/Lumi9-2/pt2/yield_data_LH2_SHMS_Sc.csv";
        OutFileName = "SHMS_LH2";
    }else if (!dataType.CompareTo(hms)){
        cout << "running Singles\n";
        NFILES = 3;
        filenames = new TString[NFILES];
        filenames[0] = "../OUTPUTS/Lumi6-4/HMS/yield_data_LH2_Sc.csv";
        filenames[1] = "../OUTPUTS/Lumi9-2/pt1/yield_data_LH2_HMS_Sc.csv";
        filenames[2] = "../OUTPUTS/Lumi9-2/pt2/yield_data_LH2_HMS_Sc.csv";
        OutFileName = "HMS_LH2";
    }else if (!dataType.CompareTo(shms_rate)){
        cout << "running Singles\n";
        NFILES = 3;
        rate = true;
        filenames = new TString[NFILES];
        filenames[0] = "../OUTPUTS/Lumi6-4/SHMS/yield_data_LH2_ScVr.csv";
        filenames[1] = "../OUTPUTS/Lumi9-2/pt1/yield_data_LH2_SHMS_ScVr.csv";
        filenames[2] = "../OUTPUTS/Lumi9-2/pt2/yield_data_LH2_SHMS_ScVr.csv";
        OutFileName = "SHMS_rate_LH2";
    }else if (!dataType.CompareTo(hms_rate)){
        cout << "running Singles\n";
        NFILES = 3;
        rate = true;
        filenames = new TString[NFILES];
        filenames[0] = "../OUTPUTS/Lumi6-4/HMS/yield_data_LH2_ScVr.csv";
        filenames[1] = "../OUTPUTS/Lumi9-2/pt1/yield_data_LH2_HMS_ScVr.csv";
        filenames[2] = "../OUTPUTS/Lumi9-2/pt2/yield_data_LH2_HMS_ScVr.csv";
        OutFileName = "HMS_rate_LH2";
    }else if (!dataType.CompareTo(lh2_all)){
        cout << "running LH2_ALL\n";
        NFILES = 13;
        filenames = new TString[NFILES];
        filenames[0] = "../OUTPUTS/Lumi6-4/SHMS/yield_data_LH2_ScvCur.csv";
        filenames[1] = "../OUTPUTS/Lumi6-4/HMS/yield_data_LH2_ScvCur.csv";
        filenames[2] = "../OUTPUTS/Lumi9-2/pt1/yield_data_LH2_HMS_Sc.csv";
        filenames[3] = "../OUTPUTS/Lumi9-2/pt1/yield_data_LH2_SHMS_Sc.csv";
        filenames[4] = "../OUTPUTS/Lumi9-2/pt2/yield_data_LH2_HMS_Sc.csv";
        filenames[5] = "../OUTPUTS/Lumi9-2/pt2/yield_data_LH2_SHMS_Sc.csv";
        filenames[6] = "../OUTPUTS/ExclusiveLumi/yield_data_exc1_ScYvCur.csv";
        filenames[7] = "../OUTPUTS/ExclusiveLumi/yield_data_exc2_ScYvCur.csv";
        filenames[8] = "../OUTPUTS/ExclusiveLumi/yield_data_exc3_ScYvCur.csv";
        filenames[9] = "../OUTPUTS/ExclusiveLumi/yield_data_exc4_ScYvCur.csv";
        filenames[10] = "../OUTPUTS/ExclusiveLumi/yield_data_exc5_ScYvCur.csv";
        filenames[11] = "../OUTPUTS/sidis/yield_data_sidis1_ScvCur.csv";
        filenames[12] = "../OUTPUTS/sidis/yield_data_sidis2_ScvCur.csv";
        OutFileName = "LH2_ALL";
    }else{
        cout << "choose one of: 'coin', 'carbon', 'singles', or 'lh2_all'. \nShutting down!\n";
        return;
    }

    //yes yes, needlessly complicated. At least I need only rewrite the above everytime the files change!
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
        
        if(i == 0) C[i]->Print(Form("../OUTPUTS/CombinedPlot%s_Scalar.pdf(", OutFileName.Data()));
        else C[i]->Print(Form("../OUTPUTS/CombinedPlot%s_Scalar.pdf", OutFileName.Data()));
        mg->Add(Gf[i]);
    }
    
    TCanvas *cf = new TCanvas("cf","cf",10, 10, 1000, 800);
    if (!rate)
        mg->GetXaxis()->SetTitle("Current (uA)");
    else
        mg->GetXaxis()->SetTitle("3/4 Rate (MHz)");
        
    mg->GetYaxis()->SetTitle("Renormalized Scalar Yield");
    
    TF1 *lin2 = new TF1("lin2", "[0]+[1]*x",0,80);
    mg->Fit(lin2);
    mg->Draw("AP");
    
    cf->Print(Form("../OUTPUTS/CombinedPlot%s_Scalar.pdf", OutFileName.Data()));
    TLegend* l1 = new TLegend(0.1, 0.1, 0.4, 0.22);
    l1->AddEntry(lin2, Form("y = (%f #pm %f)x + (%f #pm %f)",lin2->GetParameter(1), lin2->GetParError(1),lin2->GetParameter(0), lin2->GetParError(0)));
    for(int i = 0; i < NFILES; i++)
    {
        l1->AddEntry(Gf[i], Gf[i]->GetTitle());
    }
    l1->Draw();
    
    cf->Print(Form("../OUTPUTS/CombinedPlot%s_Scalar.pdf)", OutFileName.Data()));
}








