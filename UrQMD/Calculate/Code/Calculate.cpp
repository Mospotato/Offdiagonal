#define URQMD_CALCULATE_CXX
#include <iostream>
#include "TStyle.h"
#include "TSystem.h"
#include "Calculate.h"
#include "../../Dependency/Configuration.h"

void Calculate::Process()
{
    TProfile3D *pMoment = (TProfile3D *)fInput->Get("Moment");
    if (!pMoment)
    {
        Warning("Systematic::Calculate", "Abort due to Missing Profile for Moment!");
        return;
    }
    auto &Evaluator = KochRatio::getInstance();
    Evaluator.SetProfiles(pMoment);
    Evaluator.Process();
}

void Calculate::generateGraphs(const std::string &object, TList *&list, std::vector<TGraphErrors *> &graph)
{
    auto &Evaluator = KochRatio::getInstance();
    list = new TList();
    size_t nResult = Evaluator.ResultName[object].size();
    graph.resize(nResult);
    for (size_t iR = 0; iR < nResult; iR++)
    {
        graph[iR] = new TGraphErrors();
        graph[iR]->SetName(Evaluator.ResultName[object][iR].c_str());
        list->Add(graph[iR]);
    }
}

Bool_t Calculate::Init(Int_t Energy)
{
    gStyle->SetMarkerStyle(20);
    auto &Config = Configuration::getInstance();
    Config.SetCentrality(Energy);
    if (gSystem->AccessPathName(Form("../ReadUrQMD/File/%d.root", Energy)))
    {
        std::cerr << "Abort due to Missing Input File!\n";
        return kTRUE;
    }
    fInput = TFile::Open(Form("../ReadUrQMD/File/%d.root", Energy), "READ");
    if (!fInput || fInput->IsZombie())
    {
        Warning("Calculate::Init", "Abort due to Missing Input File!");
        return kTRUE;
    }
    if (ReadRefMult())
    {
        Warning("Calculate::Init", "Abort due to RefMult Distribution!");
        return kTRUE;
    }
    gSystem->mkdir("Result", kTRUE);
    fOutput = new TFile(Form("Result/%d.root", Energy), "RECREATE");
    if (!fOutput || fOutput->IsZombie())
    {
        Warning("Calculate::Init", "Abort due to Missing Output File!");
        return kTRUE;
    }
    for (auto &entry : ProxySet::getInstance().Entries)
    {
        auto &object = entry.Name;
        for (Int_t icent = 0; icent < Config.nCent; icent++)
        {
            for (Int_t IsRapidity = 0; IsRapidity < 2; IsRapidity++)
            {
                generateGraphs(object, lAcceptance[object][icent][IsRapidity], gAcceptance[object][icent][IsRapidity]);
            }
        }
        for (Int_t iAcc = 1; iAcc <= Config.nAcceptance; iAcc++)
        {
            generateGraphs(object, List[iAcc][object], gDefault[iAcc][object]);
            GenerateHistogram(iAcc, object);
        }
    }
    Info("Calculate::Init", "Initialization Completed!");
    return kFALSE;
}

void Calculate::GenerateHistogram(Int_t iAcc, const std::string &object)
{
    auto &Config = Configuration::getInstance();
    ChangeDirectory(Form("Numerator/%d/%s", iAcc - 1, object.c_str()));
    for (Int_t icent = 0; icent < Config.nCent; icent++)
    {
        hNumerator[iAcc][object][icent] = new TH1D(Form("%s", Config.CentDef[icent].c_str()), "", 100, 0, 100);
    }
    ChangeDirectory(Form("Denominator/%d/%s", iAcc - 1, object.c_str()));
    for (Int_t icent = 0; icent < Config.nCent; icent++)
    {
        hDenominator[iAcc][object][icent] = new TH1D(Form("%s", Config.CentDef[icent].c_str()), "", 100, 0, 100);
    }
    return;
}

void Calculate::LogResult(Int_t iAcc, ProxyEntry *entry)
{
    auto &Config = Configuration::getInstance();
    bool IsRapidity = iAcc <= Config.RapiditySize + 1;
    auto &Evaluator = KochRatio::getInstance();
    auto &object = entry->Name;
    size_t nResult = Evaluator.ResultName[object].size();
    auto &gdefault = gDefault[iAcc][object];
    for (Int_t icent = 0; icent < Config.nCent; icent++)
    {
        auto &gacceptance = gAcceptance[object][icent];
        auto valueMap = Evaluator.GetResult(false, iAcc, entry, icent, 2);
        auto errorMap = Evaluator.GetResult(true, iAcc, entry, icent, 2);
        for (size_t ig = 0; ig < nResult; ig++)
        {
            gdefault[ig]->SetPoint(icent, Config.Npart[icent], valueMap[ig]);
            Double_t error = sqrt(errorMap[ig]);
            gdefault[ig]->SetPointError(icent, 0, error);
            if (iAcc > 1)
            {
                int ip = IsRapidity ? iAcc - 2 : iAcc - 2 - Config.RapiditySize;
                double x = IsRapidity ? Config.SystematicRapidityVec[ip] : Config.SystematicPtVec[ip];
                gacceptance[IsRapidity][ig]->SetPoint(ip, x, valueMap[ig]);
                gacceptance[IsRapidity][ig]->SetPointError(ip, 0, error);
            }
            else
            {
                int ip = IsRapidity ? Config.RapiditySize : Config.PtSize;
                gacceptance[0][ig]->SetPoint(ip, Config.ptMax, valueMap[ig]);
                gacceptance[0][ig]->SetPointError(ip, 0, error);
                gacceptance[1][ig]->SetPoint(ip, Config.yMax, valueMap[ig]);
                gacceptance[1][ig]->SetPointError(ip, 0, error);
            }
        }
    }
    FillPercentage(iAcc, object, true, 2);
    FillPercentage(iAcc, object, false, 2);
}

void Calculate::Terminate()
{
    auto &Config = Configuration::getInstance();
    for (auto &entry : ProxySet::getInstance().Entries)
    {
        auto &object = entry.Name;
        for (Int_t iAcc = 1; iAcc <= Config.nAcceptance; iAcc++)
        {
            LogResult(iAcc, &entry);
            ChangeDirectory(Form("%s/%d", object.c_str(), iAcc - 1));
            List[iAcc][object]->Write();
        }
        WriteGraphAcceptance(object);
    }
    delete fInput;
    delete fOutput;
    return;
}

void Calculate::FillPercentage(Int_t iAcc, const std::string &object, const Bool_t &IsNumerator, const Int_t &Order)
{
    auto &Evaluator = KochRatio::getInstance();
    auto &document = Document::getInstance();
    std::pair<int, int> pair = std::pair<int, int>(IsNumerator, Order - IsNumerator);
    Component key(pair);
    auto &component = Evaluator.SeriesHolder[object][key].Components;
    for (Int_t icent = 0; icent < Configuration::getInstance().nCent; icent++)
    {
        auto &brick = Evaluator.CumulantBrick[iAcc][object][icent];
        TH1D *hPercentage = IsNumerator ? hNumerator[iAcc][object][icent] : hDenominator[iAcc][object][icent];
        Double_t &Normalization = Evaluator.Cumulant[iAcc][object][icent][key];
        Int_t Sign = IsNumerator ? 1 : -1;
        for (auto &pair : component)
        {
            std::vector<int> &Expand = ComponentHelper::getInstance()->Expand(pair.first);
            Double_t &Coefficient = pair.second;
            std::string object;
            if (Coefficient * Sign > 0)
            {
                object += "-";
            }
            if (Expand[0] == Expand[1])
            {
                object += "#color[0]{|}#LT" + document.SyntaxMap[Expand[0]] + "^{2}#GT_{c}#color[0]{|}";
            }
            else
            {
                object += "#color[0]{|}#LT" + document.SyntaxMap[Expand[0]] + document.SyntaxMap[Expand[1]] + "#GT_{c}#color[0]{|}";
            }
            Double_t Value = brick[pair.first];
            hPercentage->Fill(object.c_str(), Coefficient * Value * 100 / Normalization);
        }
        hPercentage->LabelsDeflate();
    }
}

void Calculate::WriteGraphAcceptance(const std::string &object)
{
    auto &Config = Configuration::getInstance();
    auto &lacceptance = lAcceptance[object];
    for (Int_t icent = 0; icent < Config.nCent; icent++)
    {
        for (Int_t IsRapidity = 0; IsRapidity < 2; IsRapidity++)
        {
            std::string AccDir = IsRapidity ? "Rapidity" : "Momentum";
            std::string directory = Form("%s/%s/%s", AccDir.c_str(), object.c_str(), Config.CentDef[icent].c_str());
            ChangeDirectory(directory);
            lacceptance[icent][IsRapidity]->Write();
        }
    }
    return;
}

void Calculate::ChangeDirectory(const std::string &Directory)
{
    if (!fOutput->GetDirectory(Directory.c_str()))
        fOutput->mkdir(Directory.c_str());
    fOutput->cd(Directory.c_str());
    return;
}

Bool_t Calculate::ReadRefMult()
{
    auto &Config = Configuration::getInstance();
    auto &Evaluator = KochRatio::getInstance();
    Evaluator.MinRefMult = Config.CentVec.back();
    Evaluator.MaxRefMult = Config.CentVec.front();
    Evaluator.nEvent.resize(Config.nCent, 0);
    TH1D *hRefMult = (TH1D *)fInput->Get("RefMult");
    if (hRefMult == nullptr)
    {
        Warning("Systematic::ReadRefMult", "Abort due to Missing RefMult Distribution!");
        return kTRUE;
    }
    for (Int_t im = hRefMult->GetNbinsX(); im > 0; im--)
    {
        Int_t Centrality = Config.GetCentrality(im);
        if (Centrality < 0)
        {
            continue;
        }
        Double_t Content = hRefMult->GetBinContent(im);
        Evaluator.sEvent[im].first = Content;
        Evaluator.sEvent[im].second = Centrality;
        Evaluator.nEvent[Centrality] += Content;
    }
    Config.Npart.reserve(Config.nCent);

    TH2F *hMult2Npart = (TH2F *)fInput->Get("Mult2Npart");
    if (!hMult2Npart)
    {
        std::cerr << "Error: hMult2Npart not found in file." << std::endl;
        return kTRUE;
    }
    for (Int_t iCent = 0; iCent < Config.nCent; iCent++)
    {
        Double_t Npart = hMult2Npart->ProjectionY("hTemp", Config.CentVec[iCent + 1], Config.CentVec[iCent])->GetMean();
        Config.Npart.push_back(Npart);
    }
    Info("Systematic::ReadRefMult", "Finish Reading Refmult Distribution");
    return kFALSE;
}