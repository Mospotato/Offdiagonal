#include "TFile.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TError.h"
#include "../../PlotHeaders/include.h"
void collect(int Energy = 218)
{
    TFile *fInput = new TFile(Form("Result/%d.root", Energy), "READ");
    if (!fInput || fInput->IsZombie()) {
        Error("collect", "Failed to open input file.");
        return;
    }
    std::map<std::string, std::array<std::string, 2>> ObjectMap{
        {"BS_Lambda/0/B1S1", {"BS", "#LTBS#GT_{c}"}},
        {"BS_Lambda/0/S2", {"S2", "#LTS^{2}#GT_{c}"}},
        {"BS_Lambda/0/B2", {"B2", "#LTB^{2}#GT_{c}"}},
        {"BS_Lambda/0/R11", {"CBS", "C_{BS}=#LTBS#GT_{c}/#LTS^{2}#GT_{c}"}},
        {"BQ_Lambda/0/B1Q1", {"BQ", "#LTBQ#GT_{c}"}},
        {"BQ_Lambda/0/Q2", {"Q2", "#LTQ^{2}#GT_{c}"}},
        {"BQ_Lambda/0/R11", {"CBQ", "#LTBQ#GT_{c}/#LTQ^{2}#GT_{c}"}},
        {"QS_Lambda/0/Q1S1", {"QS", "#LTQS#GT_{c}"}},
    };
    std::map<std::string, TGraphErrors*> graphsMap;
    for (const auto& [objectPath, labels] : ObjectMap) {
        TGraphErrors* graph = dynamic_cast<TGraphErrors*>(fInput->Get(objectPath.c_str()));
        if (!graph) {
            Error("collect", "Failed to retrieve object: %s", objectPath.c_str());
            continue;
        }
        graph->SetTitle(Form("%s vs N_{part};Average Number of Participanting Nucleons #LTN_{part}#GT;%s", labels[0].c_str(), labels[1].c_str()));
        graphsMap[labels[0]] = graph;
    }
    auto gRBQQS = GetRatio(graphsMap["BQ"], graphsMap["QS"]);
    if (gRBQQS) {
        gRBQQS->SetTitle("CBQ/CQS vs N_{part};Average Number of Participanting Nucleons #LTN_{part}#GT;#LTBQ#GT_{c}/#LTQS#GT_{c}");
    }
    auto gRB2Q2 = GetRatio(graphsMap["B2"], graphsMap["Q2"]);
    if (gRB2Q2) {
        gRB2Q2->SetTitle("B2/Q2 vs N_{part};Average Number of Participanting Nucleons #LTN_{part}#GT;#LTB^{2}#GT_{c}/#LTQ^{2}#GT_{c}");
    }
    gSystem->mkdir("Collected", true);
    TFile *fOutput = new TFile(Form("Collected/%d.root", Energy), "RECREATE");
    if (!fOutput || fOutput->IsZombie()) {
        Error("collect", "Failed to create output file.");
        return;
    }
    for (const auto& [label, graph] : graphsMap) {
        if (graph) {
            graph->Write(label.c_str());
        }
    }
    if (gRBQQS) {
        gRBQQS->Write("BQQS");
    }
    if (gRB2Q2) {
        gRB2Q2->Write("B2Q2");
    }
    fOutput->Close();
    fInput->Close();
}