#include "TFile.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TClass.h"
#include "TString.h"
#include "TROOT.h"

#include "../../PlotHeaders/include.h"

#include <array>
#include <string>
#include <vector>

namespace
{
    void CollectObjectPaths(TDirectory *dir, const std::string &prefix, std::vector<std::string> &out)
    {
        if (!dir)
            return;

        TIter nextKey(dir->GetListOfKeys());
        TKey *key = nullptr;
        while ((key = static_cast<TKey *>(nextKey())))
        {
            const std::string keyName = key->GetName();
            const std::string path = prefix.empty() ? keyName : (prefix + "/" + keyName);

            TClass *cls = gROOT->GetClass(key->GetClassName());
            if (cls && cls->InheritsFrom(TDirectory::Class()))
            {
                TDirectory *subdir = dynamic_cast<TDirectory *>(dir->Get(keyName.c_str()));
                CollectObjectPaths(subdir, path, out);
            }
            else
            {
                out.push_back(path);
            }
        }
    }

    std::string SanitizeName(std::string s)
    {
        for (char &c : s)
        {
            if (c == '/' || c == ':' || c == ' ')
                c = '_';
        }
        return s;
    }
    void DrawPair(TGraph *objRu, TGraph *objZr, const std::string &path)
    {
        if (!objRu || !objZr)
            return;
        const std::string canvasName = "c_" + SanitizeName(path);
        auto *c = new TCanvas(canvasName.c_str(), path.c_str(), 800, 700);
        c->cd();
        gPad->SetTicks(1, 1);
        gPad->SetMargin(0.14, 0.02, 0.13, 0.08); // left, right, bottom, top
        // Intentionally keep draw options empty so you can control styles externally.
        auto SetGeneralStyle = [](TGraph *obj, Color_t color, Marker_t marker) {
            if (auto *graph = dynamic_cast<TGraphErrors *>(obj))
            {
                graph->SetMarkerColor(color);
                graph->SetMarkerStyle(marker);
                graph->SetMarkerSize(2);
            }
        };
        SetGeneralStyle(objZr, kBlue, 21);
        auto objZrShifted = ShiftXaxis(dynamic_cast<TGraphErrors*>(objZr), 5);
        SetGeneralStyle(objRu, kRed, 20);
        SetGeneralStyle(objZrShifted, kBlue, 21);
        std::vector<TGraph*> gVec{objRu, objZrShifted};
        auto yRange = FindRange(gVec);
        auto yDiff = yRange.second - yRange.first;
        auto frame = gPad->DrawFrame(0, yRange.first - 0.1 * yDiff, 190, yRange.second + 0.1 * yDiff);
        frame->Draw("0");
        frame->SetTitle(Form("UrQMD, #sqrt{s_{NN}}=200 GeV, 0.4 < p_{T} < 1.6 GeV/c, |y| < 0.5;%s;%s", objRu->GetXaxis()->GetTitle(), objRu->GetYaxis()->GetTitle()));
        objRu->Draw("P");
        // DrawEdge(objRu, "P");
        // auto frame = objRu->GetHistogram();
        if (frame)
        {
            frame->SetTitleSize(0.05);
            TAxis *xAxis{frame->GetXaxis()}, *yAxis{frame->GetYaxis()};
            xAxis->SetTitleSize(0.05);
            xAxis->SetLabelSize(0.04);
            xAxis->SetTitleOffset(1.2);
            yAxis->SetTitleSize(0.05);
            yAxis->SetLabelSize(0.04);
            yAxis->SetTitleOffset(1.3);
        }
        objZrShifted->Draw("p");
        // DrawEdge(objZrShifted, "P");
        TLegend *leg = new TLegend(0.2, 0.3, 0.4, 0.4);
        SetLegend(leg);
        leg->AddEntry(objRu, "Ru+Ru", "p");
        leg->AddEntry(objZrShifted, "Zr+Zr", "p");
        leg->Draw();
        GetEdgeFromLegend(leg, true);
        gSystem->mkdir("Compare", true);
        c->Print(Form("Compare/%s.pdf", path.c_str()));
    }
}

void plot()
{
    gStyle->SetOptStat(0);
    std::array<TFile *, 2> fInput;
    fInput[0] = TFile::Open("Collected/218.root", "READ");
    if (!fInput[0] || fInput[0]->IsZombie())
    {
        Error("plot", "Failed to open input file: Collected/218.root");
        return;
    }

    fInput[1] = TFile::Open("Collected/219.root", "READ");
    if (!fInput[1] || fInput[1]->IsZombie())
    {
        Error("plot", "Failed to open input file: Collected/219.root");
        return;
    }

    std::vector<std::string> paths;
    CollectObjectPaths(fInput[0], "", paths);

    if (paths.empty())
    {
        Warning("plot", "No objects found in Collected/218.root");
        return;
    }

    for (const auto &path : paths)
    {
        auto objRu = (TGraph*)fInput[0]->Get(path.c_str());
        auto objZr = (TGraph*)fInput[1]->Get(path.c_str());

        if (!objRu || !objZr)
        {
            Warning("plot", "Missing object '%s' in one of the files", path.c_str());
            continue;
        }
        DrawPair(objRu, objZr, path);
    }
}