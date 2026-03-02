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
void drawRatio(Int_t Energy = 218)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(111);
    gStyle->SetLineWidth(2);
    TLegend *leg = new TLegend(0.4, 0.70, 0.7, 0.90);
    SetLegend(leg);
    leg->SetTextSize(0.035);
    TCanvas *cas = new TCanvas("cas", "", 800, 700);
    gPad->SetMargin(0.12, 0.04, 0.15, 0.1);
    TFile *fRatio = new TFile(Form("File/%d/glauber.root", Energy));
    TH1D *hUrQMD = (TH1D *)fRatio->Get("RefMult");
    leg->AddEntry(hUrQMD, "UrQMD, Ru+Ru, #sqrt{s_{NN}} = 200 GeV", "l");
    hUrQMD->SetLineWidth(2);
    hUrQMD->SetTitle("");
    TAxis *xAxis{hUrQMD->GetXaxis()}, *yAxis{hUrQMD->GetYaxis()};
    xAxis->CenterTitle(1);
    xAxis->SetTitle("EPD Multiplicity #Sigma_{i=1}^{744}#xi_{i}");
    xAxis->SetTitleOffset(1.1);
    xAxis->SetTitleSize(0.05);
    xAxis->SetLabelSize(0.04);
    yAxis->SetTitleOffset(1.2);
    yAxis->SetTitleSize(0.05);
    yAxis->SetLabelSize(0.04);
    yAxis->SetTitle("Number of Events");
    TH1D *hSimulation = (TH1D *)fRatio->Get("RefMultSim");
    leg->AddEntry(hSimulation, "Glauber Fit", "l");
    hSimulation->SetLineColor(kRed);
    hSimulation->SetLineWidth(2);
    cas->cd();
    gPad->SetLogy(1);
    gPad->SetTicks(1, 1);
    hUrQMD->Draw("hist");
    hSimulation->Draw("histsame");
    cas->cd();
    leg->Draw();
    cas->Print(Form("File/%d/Glauber.pdf", Energy));
}