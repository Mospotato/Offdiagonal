#ifndef URQMD_CALCULATE_HH
#define URQMD_CALCULATE_HH
// STL
#include <map>
#include <string>
#include <vector>
#include <unordered_map>
// ROOT Library
#include "TH1D.h"
#include "TH2.h"
#include "TList.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "../../Dependency/KochRatio.h"
struct Calculate
{
    TFile *fInput;
    TFile *fOutput;
    std::map<Int_t, std::map<std::string, std::map<Int_t, TH1D *>>> hNumerator;
    std::map<Int_t, std::map<std::string, std::map<Int_t, TH1D *>>> hDenominator;
    std::map<Int_t, std::unordered_map<std::string, TList*>> List;
    std::map<Int_t, std::unordered_map<std::string, std::vector<TGraphErrors *>>> gDefault;
    std::unordered_map<std::string, std::unordered_map<Int_t, std::array<TList*, 2>>> lAcceptance;
    std::unordered_map<std::string, std::unordered_map<Int_t, std::array<std::vector<TGraphErrors *>, 2>>> gAcceptance;
    Calculate(){};
    ~Calculate() {};
    Bool_t Init(Int_t energy);
    Bool_t ReadRefMult();
    void Process();
    void ChangeDirectory(const std::string &);
    void generateGraphs(const std::string &, TList *&list, std::vector<TGraphErrors *> &graph);
    void Terminate();

private:
    void LogResult(Int_t iAcc, ProxyEntry *);
    void FillPercentage(Int_t iAcc, const std::string &, const Bool_t &, const Int_t &Order);
    void GenerateHistogram(Int_t iAcc, const std::string &);
    void WriteGraphAcceptance(const std::string &);
};
#endif