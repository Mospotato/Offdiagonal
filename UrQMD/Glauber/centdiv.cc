#include <fstream>
#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
void centdiv(Int_t sNN)
{
    TFile *inf = new TFile(Form("File/%d/glauber.root", sNN), "READ");
    TH1D *hRefMult = (TH1D *)inf->Get("RefMult");
    std::vector<Double_t> CountsDivision;
    std::vector<int> CentVec;
    Int_t nCent = 9;
    std::vector<float> CentDivision = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
    Double_t entries = hRefMult->GetEntries();
    for (auto &division : CentDivision)
    {
        CountsDivision.emplace_back(division * entries);
    }
    Int_t binx = hRefMult->GetNbinsX();
    Double_t sum = 0;
    Double_t content;
    Int_t centFlag = 0;
    Int_t MaxBin = hRefMult->GetNbinsX();
    for (Int_t ibin = MaxBin; ibin > 0; ibin--)
    {
        content = hRefMult->GetBinContent(ibin);
        if (!content)
        {
            continue;
        }
        if (CentVec.empty() && content != 0.)
        {
            CentVec.push_back(ibin);
        }
        sum += content;
        if (sum > CountsDivision[centFlag])
        {
            CentVec.push_back(ibin);
            std::cout << round(CentDivision[centFlag] * 100) << "% :" << ibin << std::endl;
            if (++centFlag >= nCent)
            {
                break;
            }
        }
    }
    std::ofstream CentDef(Form("File/%d/RawDividion.txt", sNN));
    for(auto &cent: CentVec){
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
    return;
}