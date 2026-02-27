#include <iostream>
#include <sstream>
#include <fstream>
#include "TF1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TFile.h"
#include "TSystem.h"
void GetCentrality(Int_t Energy, Bool_t is32Cent = false)
{
    gSystem->Exec(Form("cp Configuration.txt %s/%d/Configuration.txt", is32Cent ? "32Cent" : "File", Energy));
    TFile *file = new TFile(Form("File/%d/glauber.root", Energy));
    TH1D *sim = (TH1D *)file->Get("RefMultSim");
    TH1D *data = (TH1D *)file->Get("RefMult");
    double integral = sim->Integral();
    Int_t centflag = 0;
    std::vector<Int_t> CentVec;
    Double_t sum{0}, content{0}, Threshold{0};
    std::vector<Double_t> CountsDivision;
    CountsDivision = is32Cent ? std::vector<Double_t>{0., 0.025, 0.05, 0.075, 0.10, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8} : std::vector<Double_t>{0., 0.05, 0.10, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
    Int_t nCentrality = CountsDivision.size() - 1;
    for (Int_t ibin = data->GetNbinsX(); ibin > 0; ibin--)
    {
        TH1D *&hist = CentVec.size() > 1 ? sim : data;
        content = hist->GetBinContent(ibin);
        if (CentVec.empty() && content > 1.)
        {
            CentVec.push_back(ibin);
            centflag++;
            Threshold = integral * CountsDivision[centflag];
        }
        sum += content;
        if (sum > Threshold)
        {
            CentVec.push_back(ibin);
            if (++centflag > nCentrality)
            {
                break;
            }
            Threshold = integral * CountsDivision[centflag];
        }
    }
    gSystem->mkdir(Form("%s/%d", is32Cent ? "32Cent" : "File", Energy), kTRUE);
    std::ofstream CentDef(Form("%s/%d/Definition.txt",  is32Cent ? "32Cent" : "File", Energy));
    for (auto &cent : CentVec)
    {
        std::cout << cent;
        CentDef << cent;
        if (cent != CentVec.back())
        {
            std::cout << ',';
            CentDef << ',';
        }
    }
    std::cout << '\n';
    CentDef << '\n';
    CentDef.close();
}
