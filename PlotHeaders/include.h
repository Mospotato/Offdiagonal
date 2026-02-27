#ifndef PLOTHEADERS_INCLUDE_H
#define PLOTHEADERS_INCLUDE_H
#include <map>
#include <iostream>
#include "TF1.h"
#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TPad.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TError.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TMultiGraph.h"
#include "TLegendEntry.h"
#include "TPaletteAxis.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

#include "GraphTools.h"
#include "CompactPads.h"
#include "CloneWithEdge.h"
std::vector<std::string> CentDef{"0-5%", "5-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "0-80%"};
std::map<Int_t, Double_t> EnergyFloatMap{
    {3, 3.0},
    {7, 7.7},
    {9, 9.2},
    {11, 11.5},
    {14, 14.6},
    {17, 17.3},
    {19, 19.6},
    {27, 27.0},
    {39, 39.0},
    {54, 54.4},
    {62, 62.4},
    {200, 200.0}
};
std::map<Int_t, std::string> EnergySyntax{
    {3, "3"},
    {7, "7.7"},
    {9, "9.2"},
    {11, "11.5"},
    {14, "14.6"},
    {17, "17.3"},
    {19, "19.6"},
    {27, "27"},
    {39, "39"},
    {54, "54.4"},
    {62, "62.4"},
    {200, "200"}};
std::map<Int_t, Double_t> muB{
    {3, 760},
    {7, 420},
    {9, 355},
    {11, 315},
    {14, 264},
    {17, 230},
    {19, 206},
    {27, 156},
    {39, 112},
    {54, 83},
    {62, 73}};
std::map<std::string, std::string> SpeciesSyntax{
    {"Proton", "p"},
    {"Pbar", "#bar{p}"},
    {"Kplus", "K^{+}"},
    {"Kminus", "K^{-}"},
    {"Lambda", "#Lambda"},
    {"LambdaBar", "#bar{#Lambda}"},
    {"Xi", "#Xi"},
    {"XiBar", "#bar{#Xi}"},
};
void SetLegend(TLegend *leg)
{
    leg->SetLineColor(0);
    leg->SetFillColor(10);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    return;
}
TGraphErrors *ScaleXaxis(TGraph *gBenchmarck, TGraph *gTarget)
{
    TGraphErrors *gScaled = (TGraphErrors *)gTarget->Clone();
    Double_t x;
    for (Int_t i = 0; i < gScaled->GetN(); i++)
    {
        x = gBenchmarck->GetPointX(i);
        gScaled->SetPointX(i, x);
    }
    return gScaled;
}
std::string GetSuffix(const Int_t &ipt, const Int_t &npt)
{
    if (npt == 1)
        return "";
    if (ipt == 0)
        return "(";
    if (ipt == npt - 1)
        return ")";
    return "";
}
double modifyNumber(double input, int sign)
{
    double absInput = std::fabs(input);

    int scale = 0;
    while (absInput < 1.0)
    {
        absInput *= 10;
        scale--;
    }
    int MaxDigit = 4;
    std::vector<int> digits(MaxDigit, 0);
    int position = 0;
    while (absInput > 0)
    {
        int digit = static_cast<int>(absInput);
        if (digit != 0)
        {
            if (digits[position] == 0)
            {
                digits[position] = digit; // Assign the first non-zero digit
                position++;
                if (position == MaxDigit)
                {
                    break;
                }
            }
        }
        absInput = (absInput - digit) * 10; // Move to the next digit
    }
    digits[1] += sign;
    if(digits[1]==5) digits[1] += sign;
    double modifiedNumber = 0.;
    for (int i = 0; i < 3; i++)
    {
        modifiedNumber += digits[i] * std::pow(10, scale - i);
    }
    return input < 0 ? -modifiedNumber : modifiedNumber;
}

TH1 *Cumulate(TH1 *hRaw, Bool_t IsGreater)
{
    TH1D *hCumulative = (TH1D *)hRaw->Clone(Form("%s_Cumulative", hRaw->GetName()));
    hCumulative->Reset();
    Double_t Total = hRaw->Integral();
    Double_t Integral = 0;
    for (Int_t i = 1; i <= hRaw->GetNbinsX(); i++)
    {
        Integral += hRaw->GetBinContent(i);
        hCumulative->SetBinContent(i, IsGreater ? 1. - (Integral / Total) : Integral / Total);
    }
    return hCumulative;
}
#endif