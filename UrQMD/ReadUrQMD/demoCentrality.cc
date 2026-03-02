#include "../../PlotHeaders/include.h"
void demoCentrality(int Energy = 218)
{
    auto fInput = TFile::Open(Form("File/%d.root", Energy));
    if (!fInput || fInput->IsZombie())
    {
        std::cerr << "Error: Could not open file File" << Energy << ".root" << std::endl;
        return;
    }
    auto hEPDNpart = (TH2D *)fInput->Get("Mult2Npart");
    if (!hEPDNpart)    {
        std::cerr << "Error: Could not find histogram Mult2Npart in file File" << Energy << ".root" << std::endl;
        return;
    }
    auto hEPDImpact = (TH2D *)fInput->Get("Mult2Impact");
    if (!hEPDImpact)    {
        std::cerr << "Error: Could not find histogram Mult2Impact in file File" << Energy << ".root" << std::endl;
        return;
    }
    auto drawCanvas = [&](TH2D *hist, const char *xTitle, const char *yTitle, const char *title, double ymax = 0)
    {
        gStyle->SetOptStat(0);
        TCanvas *canvas = new TCanvas(Form("cas_%s", hist->GetName()), title, 800, 700);
        gPad->SetMargin(0.12, 0.04, 0.15, 0.1);
        gPad->SetLogz();
        TAxis *xAxis{hist->GetXaxis()}, *yAxis{hist->GetYaxis()};
        xAxis->SetTitleSize(0.05);
        xAxis->CenterTitle();
        xAxis->SetLabelSize(0.04);
        yAxis->SetTitleSize(0.05);
        yAxis->SetLabelSize(0.04);
        yAxis->SetRangeUser(0, ymax > 0 ? ymax : yAxis->GetXmax());
        hist->SetXTitle(xTitle);
        hist->SetYTitle(yTitle);
        hist->SetTitle(title);
        hist->Draw("COL");
        gSystem->mkdir("pdf", true);
        canvas->SaveAs(Form("pdf/%d_%s.pdf", Energy, hist->GetName()));
    };
    drawCanvas(hEPDNpart, "EPD Multiplicity #Sigma_{i=1}^{744}#xi_{i}", "Average Number of Participants #LTN_{part}#GT", Form("UrQMD, %s, #sqrt{s_{NN}} = 200 GeV", Energy != 218 ? "Zr+Zr" : "Ru+Ru"), 200);
    drawCanvas(hEPDImpact, "EPD Multiplicity #Sigma_{i=1}^{744}#xi_{i}", "Impact Parameter b [fm]", Form("UrQMD, %s, #sqrt{s_{NN}} = 200 GeV", Energy != 218 ? "Zr+Zr" : "Ru+Ru"), 14);
}