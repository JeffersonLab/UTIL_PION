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
    //g->GetXaxis()->SetTitle("Current (uA)");
    g->GetYaxis()->SetTitle("Norm. Track Yield");
}

//does fit and resets intercept to 1 at x=0 in linear fit
// Yes this causes a memory leak, I don't care
TGraphErrors* ReBase(TGraphErrors *ge1, TGraphErrors *TLTi, bool rate)
{
    Double_t *reBase1, *reBaseErr1, *reX, *reEX, *TLT;
    TF1 *lin1 = new TF1("lin1", "[0]+[1]*x",0,80);

    reBase1 = ge1->GetY();
    reBaseErr1 = ge1->GetEY();
    int N = ge1->GetN();
    reX = ge1->GetX();
    reEX = ge1->GetEX();
    TLT = TLTi->GetY();
    lin1->FixParameter(1,0);
    Double_t *reBase2 = new Double_t[N];
    Double_t *reBaseErr2 = new Double_t[N];
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
    ge3->Fit(lin1);
    ge3->Draw("AP");
    TLegend* l1 = new TLegend(0.5, 0.85, 0.9, 0.9);
    Double_t reIntercept = lin1->GetParameter(0);
    Double_t errweightsum = 0;
    Double_t errsum = 0;
    for(int i = 0; i < N; i++)
    {
        errweightsum += reBase1[i]*reBaseErr1[i];
        errsum += reBaseErr1[i];
    }

    double additionalError = 0.0095;
    for(int i = 0; i < N; i++)
    {
        reBase2[i]= reBase1[i] - (errweightsum/errsum - 1);
        //reBase2[i]= reBase1[i] - (reIntercept - 1);
        //reBaseErr2[i] = reBaseErr1[i];// + 0.095;
        reBaseErr2[i] = sqrt(reBaseErr1[i]*reBaseErr1[i] + (additionalError*additionalError)/(TLT[i]*TLT[i]));
        //reBaseErr2[i] = reBaseErr1[i];
    }

   

    l1->AddEntry(lin1, Form("y = (%f #pm %f)x + (%f #pm %f)",lin1->GetParameter(1), lin1->GetParError(1),lin1->GetParameter(0), lin1->GetParError(0)));
    l1->Draw();
    
    
    TGraphErrors *ge2 = new TGraphErrors(N, reX2, reBase2, reEX2, reBaseErr2);
    return ge2;
}

void PlotCombinedLumi (TString dataType)
{
    gROOT->SetBatch(kTRUE);
    int NFILES;
    TString *filenames;
    TString *filenamesTLT;
    TString OutFileName;
    dataType.ToLower();
    //cout << dataType;
    TString coin = "coin";
    TString coin_rate = "coin_rate";
    TString coin_daq = "coin_daq";
    TString coin_Detrate = "coin_detrate";
    TString coinSHMSrate = "coinshmsrate";
    TString coinHMSrate = "coinhmsrate";
    TString carbon = "carbon";
    TString singles = "singles";
    TString shms = "shms";
    TString hms = "hms";
    TString lh2_all = "lh2_all";
    bool rate = false;
    
    filenamesTLT = new TString[NFILES];
    filenamesTLT[0] = "../OUTPUTS/ExclusiveLumi/TLT_exc1.csv";
    filenamesTLT[1] = "../OUTPUTS/ExclusiveLumi/TLT_exc2.csv";
    filenamesTLT[2] = "../OUTPUTS/ExclusiveLumi/TLT_exc3.csv";
    filenamesTLT[3] = "../OUTPUTS/ExclusiveLumi/TLT_exc4.csv";
    filenamesTLT[4] = "../OUTPUTS/ExclusiveLumi/TLT_exc5.csv";
    filenamesTLT[5] = "../OUTPUTS/sidis/TLT_sidis1.csv";
    filenamesTLT[6] = "../OUTPUTS/sidis/TLT_sidis2.csv";

    // This is my nifty way of handling the amount of inputfiles, so that I need only write stuff once,
    // danamically allocated arrays, a marvel to behold 
    if (!dataType.CompareTo(coin) ) {
        cout << "running Coin\n";
        NFILES = 7;
        filenames = new TString[NFILES];
        filenames[0] = "../OUTPUTS/ExclusiveLumi/yield_data_exc1_clean.csv";
        filenames[1] = "../OUTPUTS/ExclusiveLumi/yield_data_exc2_clean.csv";
        filenames[2] = "../OUTPUTS/ExclusiveLumi/yield_data_exc3_clean.csv";
        filenames[3] = "../OUTPUTS/ExclusiveLumi/yield_data_exc4_clean.csv";
        filenames[4] = "../OUTPUTS/ExclusiveLumi/yield_data_exc5_clean.csv";
        filenames[5] = "../OUTPUTS/sidis/yield_data_sidis1_clean.csv";
        filenames[6] = "../OUTPUTS/sidis/yield_data_sidis2_clean.csv";
        OutFileName = "Coin";
    } else if (!dataType.CompareTo(coin_rate) ) {
        cout << "running Coin vrs Rate\n";
        rate = true;
        NFILES = 7;
        filenames = new TString[NFILES];
        filenames[0] = "../OUTPUTS/ExclusiveLumi/yield_data_exc1_TrYvRate.csv";
        filenames[1] = "../OUTPUTS/ExclusiveLumi/yield_data_exc2_TrYvRate.csv";
        filenames[2] = "../OUTPUTS/ExclusiveLumi/yield_data_exc3_TrYvRate.csv";
        filenames[3] = "../OUTPUTS/ExclusiveLumi/yield_data_exc4_TrYvRate.csv";
        filenames[4] = "../OUTPUTS/ExclusiveLumi/yield_data_exc5_TrYvRate.csv";
        filenames[5] = "../OUTPUTS/sidis/yield_data_sidis1_TrvRate.csv";
        filenames[6] = "../OUTPUTS/sidis/yield_data_sidis2_TrvRate.csv";
        OutFileName = "CoinRate";
    } else if (!dataType.CompareTo(coin_daq) ) {
        cout << "running Coin vrs DAQ Rate\n";
        rate = true;
        NFILES = 7;
        filenames = new TString[NFILES];
        filenames[0] = "../OUTPUTS/ExclusiveLumi/yield_data_exc1_Daq.csv";
        filenames[1] = "../OUTPUTS/ExclusiveLumi/yield_data_exc2_Daq.csv";
        filenames[2] = "../OUTPUTS/ExclusiveLumi/yield_data_exc3_Daq.csv";
        filenames[3] = "../OUTPUTS/ExclusiveLumi/yield_data_exc5_Daq.csv";
        filenames[4] = "../OUTPUTS/sidis/yield_data_sidis1_Daq.csv";
        filenames[5] = "../OUTPUTS/sidis/yield_data_sidis2_Daq.csv";
        filenames[6] = "../OUTPUTS/ExclusiveLumi/yield_data_exc4_Daq.csv";
        OutFileName = "CoinDAQRate";
    } else if (!dataType.CompareTo(coin_Detrate) ) {
        cout << "running Coin vrs Rate\n";
        rate = true;
        NFILES = 7;
        filenames = new TString[NFILES];
        filenames[0] = "../OUTPUTS/ExclusiveLumi/yield_data_exc1_TrYvDRate.csv";
        filenames[1] = "../OUTPUTS/ExclusiveLumi/yield_data_exc2_TrYvDRate.csv";
        filenames[2] = "../OUTPUTS/ExclusiveLumi/yield_data_exc3_TrYvDRate.csv";
        filenames[3] = "../OUTPUTS/ExclusiveLumi/yield_data_exc5_TrYvDRate.csv";
        filenames[4] = "../OUTPUTS/sidis/yield_data_sidis1_TrvDRate.csv";
        filenames[5] = "../OUTPUTS/sidis/yield_data_sidis2_TrvDRate.csv";
        filenames[6] = "../OUTPUTS/ExclusiveLumi/yield_data_exc4_TrYvDRate.csv";
        OutFileName = "CoinDetRate";
    } else if (!dataType.CompareTo(coinSHMSrate) ) {
        cout << "running Coin vrs Rate\n";
        rate = true;
        NFILES = 7;
        filenames = new TString[NFILES];
        filenames[0] = "../OUTPUTS/ExclusiveLumi/yield_data_exc1_TrYvPRate.csv";
        filenames[1] = "../OUTPUTS/ExclusiveLumi/yield_data_exc2_TrYvPRate.csv";
        filenames[2] = "../OUTPUTS/ExclusiveLumi/yield_data_exc3_TrYvPRate.csv";
        filenames[3] = "../OUTPUTS/ExclusiveLumi/yield_data_exc5_TrYvPRate.csv";
        filenames[4] = "../OUTPUTS/sidis/yield_data_sidis1_TrvPRate.csv";
        filenames[5] = "../OUTPUTS/sidis/yield_data_sidis2_TrvPRate.csv";
        filenames[6] = "../OUTPUTS/ExclusiveLumi/yield_data_exc4_TrYvPRate.csv";
        OutFileName = "CoinSHMSRate";
    } else if (!dataType.CompareTo(coinHMSrate) ) {
        cout << "running Coin vrs Rate\n";
        rate = true;
        NFILES = 7;
        filenames = new TString[NFILES];
        filenames[0] = "../OUTPUTS/ExclusiveLumi/yield_data_exc1_TrYvHRate.csv";
        filenames[1] = "../OUTPUTS/ExclusiveLumi/yield_data_exc2_TrYvHRate.csv";
        filenames[2] = "../OUTPUTS/ExclusiveLumi/yield_data_exc3_TrYvHRate.csv";
        filenames[3] = "../OUTPUTS/ExclusiveLumi/yield_data_exc5_TrYvHRate.csv";
        filenames[4] = "../OUTPUTS/sidis/yield_data_sidis1_TrvHRate.csv";
        filenames[5] = "../OUTPUTS/sidis/yield_data_sidis2_TrvHRate.csv";
        filenames[6] = "../OUTPUTS/ExclusiveLumi/yield_data_exc4_TrYvHRate.csv";
        OutFileName = "CoinHMSRate";
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
        filenames[0] = "../OUTPUTS/Lumi6-4/SHMS/yield_data_LH2_clean.csv";
        filenames[1] = "../OUTPUTS/Lumi6-4/HMS/yield_data_LH2_clean.csv";
        filenames[2] = "../OUTPUTS/Lumi9-2/pt1/yield_data_LH2_HMS.csv";
        filenames[3] = "../OUTPUTS/Lumi9-2/pt1/yield_data_LH2_SHMS.csv";
        filenames[4] = "../OUTPUTS/Lumi9-2/pt2/yield_data_LH2_HMS.csv";
        filenames[5] = "../OUTPUTS/Lumi9-2/pt2/yield_data_LH2_SHMS.csv";
        OutFileName = "LH2_Singles";
    }else if (!dataType.CompareTo(shms)){
        cout << "running Singles\n";
        NFILES = 3;
        filenames = new TString[NFILES];
        filenames[0] = "../OUTPUTS/Lumi6-4/SHMS/yield_data_LH2_clean.csv";
        filenames[1] = "../OUTPUTS/Lumi9-2/pt1/yield_data_LH2_SHMS.csv";
        filenames[2] = "../OUTPUTS/Lumi9-2/pt2/yield_data_LH2_SHMS.csv";
        OutFileName = "SHMS_LH2";
    }else if (!dataType.CompareTo(hms)){
        cout << "running Singles\n";
        NFILES = 3;
        filenames = new TString[NFILES];
        filenames[0] = "../OUTPUTS/Lumi6-4/HMS/yield_data_LH2_clean.csv";
        filenames[1] = "../OUTPUTS/Lumi9-2/pt1/yield_data_LH2_HMS.csv";
        filenames[2] = "../OUTPUTS/Lumi9-2/pt2/yield_data_LH2_HMS.csv";
        OutFileName = "HMS_LH2";
    }else if (!dataType.CompareTo(lh2_all)){
        cout << "running LH2_ALL\n";
        NFILES = 13;
        filenames = new TString[NFILES];
        filenames[0] = "../OUTPUTS/Lumi6-4/SHMS/yield_data_LH2_clean.csv";
        filenames[1] = "../OUTPUTS/Lumi6-4/HMS/yield_data_LH2_clean.csv";
        filenames[2] = "../OUTPUTS/Lumi9-2/pt1/yield_data_LH2_HMS.csv";
        filenames[3] = "../OUTPUTS/Lumi9-2/pt1/yield_data_LH2_SHMS.csv";
        filenames[4] = "../OUTPUTS/Lumi9-2/pt2/yield_data_LH2_HMS.csv";
        filenames[5] = "../OUTPUTS/Lumi9-2/pt2/yield_data_LH2_SHMS.csv";
        filenames[6] = "../OUTPUTS/ExclusiveLumi/yield_data_exc1_clean.csv";
        filenames[7] = "../OUTPUTS/ExclusiveLumi/yield_data_exc2_clean.csv";
        filenames[8] = "../OUTPUTS/ExclusiveLumi/yield_data_exc3_clean.csv";
        filenames[9] = "../OUTPUTS/ExclusiveLumi/yield_data_exc4_clean.csv";
        filenames[10] = "../OUTPUTS/ExclusiveLumi/yield_data_exc5_clean.csv";
        filenames[11] = "../OUTPUTS/sidis/yield_data_sidis1_clean.csv";
        filenames[12] = "../OUTPUTS/sidis/yield_data_sidis2_clean.csv";
        OutFileName = "LH2_ALL";
    }else{
        cout << "choose one of: 'coin', 'carbon', 'singles', or 'lh2_all'. \nShutting down!\n";
        return;
    }

    //yes yes, needlessly complicated. At least I need only rewrite the above everytime the files change!
    TCanvas **C = new TCanvas*[NFILES]; 
    TGraphErrors **Gi = new TGraphErrors*[NFILES];
    TGraphErrors **TLTi = new TGraphErrors*[NFILES];
    TGraphErrors **Gf = new TGraphErrors*[NFILES];
    TMultiGraph *mg = new TMultiGraph(Form("%s Data Combined", OutFileName.Data()),Form("%s Data Combined", OutFileName.Data()));
    
    for(int i = 0; i < NFILES; i++)
    {
        C[i] = new TCanvas(Form("c%d", i),Form("c%d", i),10, 10, 1000, 800);
        Gi[i] = new TGraphErrors(filenames[i], "%lg, %lg, %lg");
        TLTi[i] = new TGraphErrors(filenamesTLT[i], "%lg, %lg, %lg");
        Gi[i]->SetName(filenames[i]);
        Gi[i]->SetTitle(filenames[i]);
        Gf[i] = ReBase(Gi[i], TLTi[i], rate);
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
        
    mg->GetYaxis()->SetTitle("Renormalized Yield");
    
    TF1 *lin2 = new TF1("lin2", "[0]+[1]*x",0,80);
    //TF1 *lin2 = new TF1("lin2", "[0]",0,80);
    //lin2->FixParameter(1,0);
    mg->Fit(lin2);
    mg->Draw("AP");
    /*
    TF1 *lin3 = new TF1("lin3", "[0]+[1]*x",0,80);
    lin3->SetParameter(1,lin2->GetParameter(1));
    lin3->SetParameter(0,lin2->GetParameter(0)+lin2->GetParError(0));
    lin3->SetLineStyle(2);
    lin3->SetLineWidth(10);
    //lin3->Draw("sames AP");
    
    TF1 *lin4 = new TF1("lin4", "[0]+[1]*x",0,80);
    lin4->SetParameter(1,lin2->GetParameter(1));
    lin4->SetParameter(0,lin2->GetParameter(0)-lin2->GetParError(0));
    lin4->SetLineStyle(4);
    lin4->SetLineWidth(6);
    //lin4->Draw("sames AP");
    */
    cf->Print(Form("../OUTPUTS/CombinedPlot%s.pdf", OutFileName.Data()));
    TLegend* l1 = new TLegend(0.1, 0.1, 0.4, 0.22);
    l1->AddEntry(lin2, Form("y = (%f #pm %f)x + (%f #pm %f)",lin2->GetParameter(1), lin2->GetParError(1),lin2->GetParameter(0), lin2->GetParError(0)));
    //l1->AddEntry(lin2, Form("y = (%f #pm %f)",lin2->GetParameter(0), lin2->GetParError(0)));
    l1->AddEntry("", Form("Reduced #Chi^2 = %f", lin2->GetChisquare()/lin2->GetNDF()));
    for(int i = 0; i < NFILES; i++)
    {
        l1->AddEntry(Gf[i], Gf[i]->GetTitle());
    }
    l1->Draw();
    
    cf->Print(Form("../OUTPUTS/CombinedPlot%s.pdf)", OutFileName.Data()));
}








