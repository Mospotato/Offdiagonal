#ifndef KOCHRATIO_HH
#define KOCHRATIO_HH
// C++ STL
#include <unordered_set>
#include <unordered_map>
#include <array>
#include <fstream>
// ROOT Library
#include "TH1D.h"
#include "TFile.h"
#include "TError.h"
#include "TProfile3D.h"
#include "TGraphErrors.h"
// User Defined Library
#include "Dictionary.h"
#include "Configuration.h"
struct ProxyEntry
{
    ProxyEntry() {}
    ~ProxyEntry() {}
    ProxyEntry(const std::string &name,
               const std::unordered_set<int> &baryonSet,
               const std::unordered_set<int> &strangeSet,
               const std::unordered_set<int> &chargeSet,
               const std::pair<int, int> &conservedPair,
               Double_t factor) : Name(name), BaryonSet(baryonSet), StrangeSet(strangeSet), ChargeSet(chargeSet), ConservedPair(conservedPair), Factor(factor) {}
    ProxyEntry(const ProxyEntry &entry)
    {
        this->Name = entry.Name;
        this->BaryonSet = entry.BaryonSet;
        this->StrangeSet = entry.StrangeSet;
        this->ChargeSet = entry.ChargeSet;
        this->ConservedPair = entry.ConservedPair;
        this->Factor = entry.Factor;
    }
    ProxyEntry &operator=(const ProxyEntry &entry)
    {
        this->Name = entry.Name;
        this->BaryonSet = entry.BaryonSet;
        this->StrangeSet = entry.StrangeSet;
        this->ChargeSet = entry.ChargeSet;
        this->ConservedPair = entry.ConservedPair;
        this->Factor = entry.Factor;
        return *this;
    }
    // Move Constructor
    ProxyEntry(ProxyEntry &&entry) noexcept
    {
        this->Name = std::move(entry.Name);
        this->BaryonSet = std::move(entry.BaryonSet);
        this->StrangeSet = std::move(entry.StrangeSet);
        this->ChargeSet = std::move(entry.ChargeSet);
        this->ConservedPair = entry.ConservedPair;
        this->Factor = entry.Factor;
    }
    // Move Assignment
    ProxyEntry &operator=(ProxyEntry &&entry) noexcept
    {
        this->Name = std::move(entry.Name);
        this->BaryonSet = std::move(entry.BaryonSet);
        this->StrangeSet = std::move(entry.StrangeSet);
        this->ChargeSet = std::move(entry.ChargeSet);
        this->ConservedPair = entry.ConservedPair;
        this->Factor = entry.Factor;
        return *this;
    }
    std::string Name;
    std::unordered_set<int> BaryonSet;
    std::unordered_set<int> StrangeSet;
    std::unordered_set<int> ChargeSet;
    std::pair<int, int> ConservedPair = {1, 2};
    Double_t Factor;
};
struct ProxySet
{
    static ProxySet &getInstance()
    {
        static ProxySet instance;
        return instance;
    }
    std::vector<ProxyEntry> Entries;
    ProxySet();
};
struct KochRatio
{
    using Pool = std::unordered_map<Component, Double_t, ComponentHash>;
    using Element = std::unordered_map<std::string, std::map<Int_t, Pool>>;
    KochRatio();
    ~KochRatio() {};
    static KochRatio &getInstance()
    {
        static KochRatio instance;
        return instance;
    }
    void SetProfiles(TProfile3D *moment)
    {
        Ratio.clear();
        Error.clear();
        Cumulant.clear();
        this->pMoment = moment;
        return;
    }
    void Process();
    void LoadMoments(Int_t iAcc, const Int_t &RefMult);
    Double_t GetCumulant(const Pool &, const Component &);
    Double_t GetError(const Pool &, const Component &);
    Double_t GetRatioError(Int_t Order, const Pool &, const Pool &, Int_t, Int_t);
    Double_t getError(const Pool &, std::unordered_map<Component, double, ComponentHash>&);
    Double_t getDerivative(const Pool &, const Component &, const Component &);
    void Calculate(Int_t iAcc, const std::string &name);
    void Task(Int_t iAcc, Int_t RefMult);
    void SetStatusOff(int type);
    void SetStatusOff(std::vector<int> &type)
    {
        for (auto tmp : type)
        {
            SetStatusOff(tmp);
        }
    }
    std::vector<Double_t> GetResult(Bool_t IsError, Int_t iAcc, ProxyEntry *, const Int_t &Centrality, const Int_t &Order);
    //=====================================================
    TProfile3D *pMoment;
    //=====================================================
    // Initialize depending on the dataset
    Int_t nProxy;
    std::unordered_map<std::string, std::unordered_map<Component, Series, ComponentHash>> SeriesHolder;
    std::unordered_map<std::string, std::pair<int, int>> ConservedPairMap;

    std::unordered_map<std::string, std::set<Component>> Result;
    std::unordered_map<std::string, std::set<Component>> Collection;

    std::unordered_set<int> OffSet;

    std::unordered_map<std::string, std::vector<std::string>> ResultName;

    Int_t MinRefMult;
    Int_t MaxRefMult;
    Double_t Weight;
    Double_t ErrorWeight;
    Int_t Centrality;
    std::vector<Double_t> nEvent;
    std::unordered_map<Int_t, std::pair<Double_t, Int_t>> sEvent;
    std::map<Int_t, Pool> MomentPool;
    std::map<Int_t, Pool> CumulantPool;
    //=====================================================
    // Value need to refresh when TProfile is changed.
    std::map<Int_t, Element> Ratio;
    std::map<Int_t, Element> RatioError;
    std::map<Int_t, Element> Error;
    std::map<Int_t, Element> Cumulant;
    std::map<Int_t, Element> CumulantBrick;
    std::map<Int_t, Element> ErrorBrick;
    //=======================================================
};
#endif