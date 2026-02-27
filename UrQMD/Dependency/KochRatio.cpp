#define KOCHRATIO_CXX
#include <iostream>
#include <sstream>
#include <algorithm>
#include "KochRatio.h"
#include "ThreadPool.h"
Component BuildKey(int leftType, Int_t leftPower, int rightType, Int_t rightPower)
{
    Component key;
    if (leftPower)
    {
        key[leftType] = leftPower;
    }
    if (rightPower)
    {
        key[rightType] = rightPower;
    }
    return key;
}

ProxySet::ProxySet()
{
    Entries.reserve(8);
    Entries.emplace_back(ProxyEntry("Lambda", {2212, 3122}, {3122, 321}, -3., OffdiagonalType::kBS));
    Entries.emplace_back(ProxyEntry("All", {1}, {2}, -3., OffdiagonalType::kBS));
    Entries.emplace_back(ProxyEntry("Lambda", {2212, 3122}, {2212, 211, 321}, 1., OffdiagonalType::kBQ));
    Entries.emplace_back(ProxyEntry("All", {1}, {3}, 1., OffdiagonalType::kBQ));
    Entries.emplace_back(ProxyEntry("Lambda", {2212, 211, 321}, {3122, 321}, 1., OffdiagonalType::kQS));
    Entries.emplace_back(ProxyEntry("All", {2}, {3}, 1., OffdiagonalType::kQS));
}
void KochRatio::LoadMoments(Int_t iAcc, const Int_t &RefMult)
{
    Dictionary *Dict = Dictionary::getInstance();
    auto &momentPool = MomentPool[iAcc];
    Int_t iBin = 1;
    for (auto &component : Dict->Constituent)
    {
        momentPool[component] = pMoment->GetBinContent(RefMult, iBin, iAcc);
        iBin++;
    }
    auto &cumulantPool = CumulantPool[iAcc];
    for (auto &component : Dict->Constituent)
    {
        cumulantPool[component] = GetCumulant(momentPool, component);
    }
}

void KochRatio::Calculate(Int_t iAcc, const std::string &object)
{
    Pool &cumulantPool = CumulantPool[iAcc];
    Pool &momentPool = MomentPool[iAcc];
    Pool &brickCumulant = CumulantBrick[iAcc][object][Centrality];
    for (const auto &component : Collection[object])
    {
        brickCumulant[component] += Weight * cumulantPool[component];
    }
    Pool &error = ErrorBrick[iAcc][object][Centrality];
    for (const auto &component : Result[object])
    {
        error[component] += ErrorWeight * GetError(momentPool, component);
    }
    auto *helper = ComponentHelper::getInstance();
    auto getComposition = [&](const Series &series, const Pool &pool) -> Double_t
    {
        Double_t value = 0;
        for (const auto &pair : series.Components)
        {
            value += pair.second * pool.at(pair.first);
        }
        return value;
    };
    Pool moment, cumulant;
    auto &cumulantComposition = Cumulant[iAcc][object][Centrality];
    for (auto &Obs : SeriesHolder[object])
    {
        moment[Obs.first] = getComposition(Obs.second, momentPool);
        cumulant[Obs.first] = getComposition(Obs.second, cumulantPool);
        cumulantComposition[Obs.first] += Weight * cumulant[Obs.first];
    }
    auto &errorComposition = Error[iAcc][object][Centrality];
    for (auto &Obs : SeriesHolder[object])
    {
        Int_t Order = helper->GetOrder(Obs.first);
        if (Order > 2)
        {
            continue;
        }
        errorComposition[Obs.first] += ErrorWeight * GetError(moment, Obs.first);
    }
    auto &ratio = Ratio[iAcc][object][Centrality];
    auto &ratioError = RatioError[iAcc][object][Centrality];
    auto Pair = GetConservedTypes(GetOffdiagonalType(object));
    for (Int_t Order = 0; Order < 2; Order++)
    {
        Component numerator = BuildKey(Pair.first, 1, Pair.second, Order);
        Component denominator = BuildKey(Pair.first, 0, Pair.second, Order + 1);
        ratio[numerator] += Weight * cumulant[numerator] / cumulant[denominator];
        ratioError[numerator] += ErrorWeight * GetRatioError(Order + 1, moment, cumulant, GetOffdiagonalType(object));
    }
    return;
}

Double_t KochRatio::GetCumulant(const Pool &moment, const Component &component)
{
    auto *helper = ComponentHelper::getInstance();
    Int_t Order = helper->GetOrder(component);
    if (Order == 1)
        return moment.at(component);
    std::vector<int> &Expand = helper->Expand(component);
    auto &partition = Partition::instance();
    Double_t Cumulant = 0.;
    std::vector<std::vector<std::vector<int>>> &sets = partition.getPartitions(Order);
    for (const auto &set : sets)
    {
        Double_t subset = 1;
        for (const auto &power : set)
        {
            std::vector<int> expansion(power.size(), 0);
            std::transform(power.begin(), power.end(), expansion.begin(), [&](int idx)
                           { return Expand[idx]; });
            Component component(expansion);
            subset *= moment.at(component);
        }
        Cumulant += partition.getCoefficient(set.size()) * subset;
    }
    return Cumulant;
}

Double_t KochRatio::GetError(const Pool &moment, const Component &component)
{
    auto *helper = ComponentHelper::getInstance();
    std::unordered_map<Component, double, ComponentHash> componentSet;
    auto &sets = helper->Subset(component);
    for (const auto &differential : sets)
    {
        Double_t value = getDerivative(moment, component, differential);
        componentSet.insert({differential, value});
    }
    Double_t error = getError(moment, componentSet);
    return error;
}

Double_t KochRatio::GetRatioError(Int_t Order, const Pool &moment, const Pool &cumulant, OffdiagonalType type)
{
    auto *helper = ComponentHelper::getInstance();
    auto Pair = GetConservedTypes(type);
    Component numerator = BuildKey(Pair.first, 1, Pair.second, Order - 1);
    Component denominator = BuildKey(Pair.first, 0, Pair.second, Order);
    Component MaxIndex = BuildKey(Pair.first, 1, Pair.second, Order);
    std::unordered_map<Component, double, ComponentHash> componentSet;
    auto &sets = helper->Subset(MaxIndex);
    for (const auto &differential : sets)
    {
        if (helper->GetOrder(differential) > Order)
        {
            continue;
        }
        Double_t numeValue = getDerivative(moment, numerator, differential);
        numeValue /= cumulant.at(denominator);
        Double_t denoValue = getDerivative(moment, denominator, differential);
        denoValue *= cumulant.at(numerator) / std::pow(cumulant.at(denominator), 2);
        componentSet.insert({differential, numeValue - denoValue});
    }
    Double_t error = getError(moment, componentSet);
    return error;
}

Double_t KochRatio::getDerivative(const Pool &moment, const Component &component, const Component &differential)
{
    auto *derive = DeriveHelper::getInstance();
    for (auto &pair : differential.Internal)
    {
        if (component.Internal.find(pair.first) == component.Internal.end())
        {
            return 0;
        }
        else if (component.Internal.at(pair.first) < pair.second)
        {
            return 0;
        }
    }
    Product product = derive->Derive(component, differential);
    Double_t value = product.Internal.second;
    for (const auto &pair : product.Internal.first.Components)
    {
        value *= std::pow(moment.at(pair.first), pair.second);
    }
    return value;
}

Double_t KochRatio::getError(const Pool &moment, std::unordered_map<Component, double, ComponentHash> &componentSet)
{
    Double_t error = 0;
    for (const auto &iPair : componentSet)
    {
        for (const auto &jPair : componentSet)
        {
            error += iPair.second * jPair.second * (moment.at(iPair.first * jPair.first) - moment.at(iPair.first) * moment.at(jPair.first));
        }
    }
    return error;
}

void KochRatio::Task(Int_t iAcc, Int_t RefMult)
{
    LoadMoments(iAcc, RefMult);
    for (const auto &entry : ProxySet::getInstance().Entries)
    {
        Calculate(iAcc, entry.Name);
    }
}

void KochRatio::Process()
{

    auto &Config = Configuration::getInstance();
    for (Int_t mult = MinRefMult; mult < MaxRefMult; mult++)
    {
        if (sEvent[mult].first < 10)
        {
            continue;
        }
        Centrality = sEvent[mult].second;
        Weight = sEvent[mult].first / nEvent[Centrality];
        ErrorWeight = Weight / nEvent[Centrality];
#if 0
        for (int iAcc = 1; iAcc <= Config.nAcceptance; iAcc++)
        {
            pool.enqueue([this, iAcc, mult]
                                   { Task(iAcc, mult); });
        }
        pool.waitForCompletion();
#else
        for (int iAcc = 1; iAcc <= Config.nAcceptance; iAcc++)
        {
            Task(iAcc, mult);
        }
#endif
    }
    Info("KochRatio::Process", "Finish Processing");
}

void KochRatio::SetStatusOff(int type)
{
    Info("KochRatio::SetStatusOff", "Set %s Off!", Document::getInstance().StringMap[type].c_str());
    OffSet.insert(type);
    return;
}

KochRatio::KochRatio()
{
    auto *helper = ComponentHelper::getInstance();
    auto *Dict = Dictionary::getInstance();
    auto &PDG = PDGData::getInstance();
    for (const auto &entry : ProxySet::getInstance().Entries)
    {
        auto &object = entry.Name;
        OffdiagTypeMap[object] = entry.Type;
        auto Pair = GetConservedTypes(entry.Type);
        std::vector<Particle> LeftArray, RightArray;
        if (object.find("All") == std::string::npos)
        {
            LeftArray.reserve(2 * entry.LeftSet.size());
            RightArray.reserve(2 * entry.RightSet.size());
            for (auto &pdg : entry.LeftSet)
            {
                Int_t value = GetConservedValue(PDG.PDGMap.at(pdg), Pair.first);
                LeftArray.push_back(std::make_pair(pdg, value));
                LeftArray.push_back(std::make_pair(-pdg, -value));
            }
            for (auto &pdg : entry.RightSet)
            {
                Int_t value = GetConservedValue(PDG.PDGMap.at(pdg), Pair.second);
                RightArray.push_back(std::make_pair(pdg, value));
                RightArray.push_back(std::make_pair(-pdg, -value));
            }
        }
        else
        {
            LeftArray.push_back(std::make_pair(Pair.first, 1));
            RightArray.push_back(std::make_pair(Pair.second, 1));
        }
        for (const auto &power : Dict->Powers)
        {
            Series Bi(power[0], LeftArray);
            Bi *= Series(power[1], RightArray);
            for (auto &component : Bi.Components)
            {
                if (component.first.size() == 0)
                {
                    continue;
                }
                Collection[object].insert(component.first);
                if (helper->GetOrder(component.first) > 2)
                {
                    continue;
                }
                Result[object].insert(component.first);
            }
            Component key = BuildKey(Pair.first, power[0], Pair.second, power[1]);
            SeriesHolder[object][key] = std::move(Bi);
        }
        auto &document = Document::getInstance();
        ResultName[object] = {
            "R11", Form("%s1%s1", document.SyntaxMap[Pair.first].c_str(), document.SyntaxMap[Pair.second].c_str()), Form("%s2", document.SyntaxMap[Pair.second].c_str()), Form("%s2", document.SyntaxMap[Pair.first].c_str())};
        for (const auto &component : Result[object])
        {
            ResultName[object].push_back(component.GetString());
        }
        Info("KochRatio::AddObject", "%s Added. %lu Components Generated.", object.c_str(), Collection[object].size());
    }
}

std::vector<Double_t> KochRatio::GetResult(Bool_t IsError, Int_t iAcc, ProxyEntry *entry, const Int_t &Centrality, const Int_t &Order)
{
    std::vector<Double_t> value;
    auto &object = entry->Name;
    auto Pair = GetConservedTypes(GetOffdiagonalType(object));
    Component LeftSqure = BuildKey(Pair.first, Order, Pair.second, 0);
    Component numerator = BuildKey(Pair.first, 1, Pair.second, Order - 1);
    Component denominator = BuildKey(Pair.first, 0, Pair.second, Order);
    value.reserve(ResultName[object].size());
    if (IsError)
    {
        value = {fabs(entry->Factor) * RatioError[iAcc][object][Centrality][numerator], Error[iAcc][object][Centrality][numerator], Error[iAcc][object][Centrality][denominator], Error[iAcc][object][Centrality][LeftSqure]};
        auto &error = ErrorBrick[iAcc][object][Centrality];
        for (const auto &component : Result[object])
        {
            value.push_back(error[component]);
        }
    }
    else
    {
        value = {entry->Factor * Ratio[iAcc][object][Centrality][numerator], Cumulant[iAcc][object][Centrality][numerator], Cumulant[iAcc][object][Centrality][denominator], Cumulant[iAcc][object][Centrality][LeftSqure]};
        auto &brick = CumulantBrick[iAcc][object][Centrality];
        for (const auto &component : Result[object])
        {
            value.push_back(brick[component]);
        }
    }
    return value;
}