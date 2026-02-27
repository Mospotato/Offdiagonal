#include <iostream>
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "../../PlotHeaders/include.h"
#include "../../PlotHeaders/CompactPads.h"
void drawRatio(Int_t Energy = 7)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(111);
    gStyle->SetLineWidth(2);
    TLegend *leg = new TLegend(0.7, 0.70, 0.9, 0.90);
    SetLegend(leg);
    leg->SetTextSize(0.035);
    TCanvas *cas = new TCanvas("cas", "", 600, 500);
    TFile *fRatio = new TFile(Form("File/%d/glauber.root", Energy));
    TH1D *hUrQMD = (TH1D *)fRatio->Get("RefMult");
    leg->AddEntry(hUrQMD, "UrQMD", "l");
    hUrQMD->SetLineWidth(2);
    hUrQMD->SetTitle("");
    hUrQMD->GetXaxis()->CenterTitle(1);
    hUrQMD->GetXaxis()->SetTitle("RefMultPi");
    hUrQMD->GetXaxis()->SetTitleOffset(1.1);
    TH1D *hSimulation = (TH1D *)fRatio->Get("RefMultSim");
    leg->AddEntry(hSimulation, "Glauber Fit", "l");
    hSimulation->SetLineColor(kRed);
    hSimulation->SetLineWidth(2);
    cas->cd();
    gPad->SetLogy(1);
    gPad->SetTicks(1, 1);
    gPad->SetMargin(0.06, 0.01, 0.1, 0.01);
    hUrQMD->Draw("hist");
    hUrQMD->GetXaxis()->CenterTitle(1);
    hUrQMD->GetXaxis()->SetTitleFont(62);
    hUrQMD->GetXaxis()->SetRange(0, 3200);
    hSimulation->Draw("histsame");
    cas->cd();
    leg->Draw();
    cas->Print(Form("File/%d/Glauber.pdf", Energy));
}