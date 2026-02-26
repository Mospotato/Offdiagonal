#define KOCHRATIO_CXX
#include <iostream>
#include <sstream>
#include "KochRatio.h"
#include "ThreadPool.h"
ProxySet::ProxySet()
{
    Entries.reserve(12);
    Entries.emplace_back(ProxyEntry("pK", {2212}, {321}, {}, {1, 2}, 1.));
    Entries.emplace_back(ProxyEntry("KLambda", {3122}, {321}, {}, {1, 2}, 1.));
    Entries.emplace_back(ProxyEntry("KXi", {3312}, {321}, {}, {1, 2}, 1.));
    Entries.emplace_back(ProxyEntry("pLambda", {2212}, {3122}, {}, {1, 2}, 1.));
    Entries.emplace_back(ProxyEntry("pXi", {2212}, {3312}, {}, {1, 2}, 1.));

    Entries.emplace_back(ProxyEntry("Lambda", {2212, 3122}, {3122, 321}, {}, {1, 2}, -3.));
    Entries.emplace_back(ProxyEntry("Xi", {2212, 3122, 3312}, {3122, 3312, 321}, {}, {1, 2}, -3.));
    Entries.emplace_back(ProxyEntry("AllBS", {1}, {2}, {}, {1, 2}, -3.));

    Entries.emplace_back(ProxyEntry("pPi", {2212}, {}, {211}, {1, 3}, 1.));
    Entries.emplace_back(ProxyEntry("LambdaPi", {3122}, {}, {211}, {1, 3}, 1.));
    Entries.emplace_back(ProxyEntry("KPi", {}, {321}, {211}, {2, 3}, 1.));
    Entries.emplace_back(ProxyEntry("AllBQ", {1}, {}, {3}, {1, 3}, 1.));
    Entries.emplace_back(ProxyEntry("AllSQ", {}, {2}, {3}, {2, 3}, 1.));
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
    const auto &conservedPair = ConservedPairMap[object];
    Int_t primaryType = conservedPair.first;
    Int_t secondaryType = conservedPair.second;

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
    for (Int_t Order = 0; Order < 2; Order++)
    {
        Component numerator, denominator, ratioKey;
        numerator[primaryType] = 1;
        if (Order)
        {
            numerator[secondaryType] = Order;
        }
        denominator[secondaryType] = Order + 1;
        ratioKey[primaryType] = 1;
        if (Order)
        {
            ratioKey[secondaryType] = Order;
        }
        ratio[ratioKey] += Weight * -3. * cumulant[numerator] / cumulant[denominator];
        ratioError[ratioKey] += ErrorWeight * GetRatioError(Order + 1, moment, cumulant, primaryType, secondaryType);
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

Double_t KochRatio::GetRatioError(Int_t Order, const Pool &moment, const Pool &cumulant, Int_t primaryType, Int_t secondaryType)
{
    auto *helper = ComponentHelper::getInstance();
    Component numerator, denominator, MaxIndex;
    numerator[primaryType] = 1;
    denominator[secondaryType] = Order;
    MaxIndex[primaryType] = 1;
    MaxIndex[secondaryType] = Order;
    if (Order > 1)
    {
        numerator[secondaryType] = Order - 1;
    }
    std::unordered_map<Component, double, ComponentHash> componentSet;
    auto &sets = helper->Subset(MaxIndex);
    for (const auto &differential : sets)
    {
        if (helper->GetOrder(differential) > Order)
        {
            continue;
        }
        Double_t numeValue = getDerivative(moment, numerator, differential);
        numeValue *= -3. / cumulant.at(denominator);
        Double_t denoValue = getDerivative(moment, denominator, differential);
        denoValue *= 3. * cumulant.at(numerator) / std::pow(cumulant.at(denominator), 2);
        componentSet.insert({differential, numeValue + denoValue});
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
    auto &document = Document::getInstance();
    for (const auto &entry : ProxySet::getInstance().Entries)
    {
        auto &object = entry.Name;
        ConservedPairMap[object] = entry.ConservedPair;

        std::vector<Particle> BaryonArray, StrangeArray, ChargeArray;
        if (object.find("All") != 0)
        {
            BaryonArray.reserve(2 * entry.BaryonSet.size());
            StrangeArray.reserve(2 * entry.StrangeSet.size());
            ChargeArray.reserve(2 * entry.ChargeSet.size());
            for (auto &pdg : entry.BaryonSet)
            {
                Int_t Baryon = PDG.PDGMap.at(pdg)->Baryon;
                BaryonArray.push_back(std::make_pair(pdg, Baryon));
                BaryonArray.push_back(std::make_pair(-pdg, -Baryon));
            }
            for (auto &pdg : entry.StrangeSet)
            {
                Int_t Strangeness = PDG.PDGMap.at(pdg)->Strangeness;
                StrangeArray.push_back(std::make_pair(pdg, Strangeness));
                StrangeArray.push_back(std::make_pair(-pdg, -Strangeness));
            }
            for (auto &pdg : entry.ChargeSet)
            {
                Int_t Charge = PDG.PDGMap.at(pdg)->Charge;
                ChargeArray.push_back(std::make_pair(pdg, Charge));
                ChargeArray.push_back(std::make_pair(-pdg, -Charge));
            }
        }
        else
        {
            BaryonArray.push_back(std::make_pair(1, 1));  // Hold All Baryons
            StrangeArray.push_back(std::make_pair(2, 1)); // Hold All Strange Particles
            ChargeArray.push_back(std::make_pair(3, 1));  // Hold All Charge Particles
        }

        std::unordered_map<int, std::vector<Particle> *> chargeBase = {
            {1, &BaryonArray},
            {2, &StrangeArray},
            {3, &ChargeArray}};

        Int_t primaryType = entry.ConservedPair.first;
        Int_t secondaryType = entry.ConservedPair.second;

        for (const auto &power : Dict->Powers)
        {
            Series Bi(power[0], *chargeBase[primaryType]);
            Bi *= Series(power[1], *chargeBase[secondaryType]);
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
            Component key;
            if (power[0])
            {
                key[primaryType] = power[0];
            }
            if (power[1])
            {
                key[secondaryType] = power[1];
            }
            SeriesHolder[object][key] = std::move(Bi);
        }

        std::string firstName = document.SyntaxMap[primaryType] + "1" + document.SyntaxMap[secondaryType] + "1";
        std::string secondName = document.SyntaxMap[secondaryType] + "2";
        ResultName[object] = {"R11", firstName, secondName};
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
    Int_t primaryType = entry->ConservedPair.first;
    Int_t secondaryType = entry->ConservedPair.second;
    Component numerator, denominator;
    numerator[primaryType] = 1;
    if (Order > 1)
    {
        numerator[secondaryType] = Order - 1;
    }
    denominator[secondaryType] = Order;

    value.reserve(ResultName[object].size());
    if (IsError)
    {
        value = {RatioError[iAcc][object][Centrality][numerator], fabs(entry->Factor) * Error[iAcc][object][Centrality][numerator], Error[iAcc][object][Centrality][denominator]};
        auto &error = ErrorBrick[iAcc][object][Centrality];
        for (const auto &component : Result[object])
        {
            value.push_back(error[component]);
        }
    }
    else
    {
        value = {Ratio[iAcc][object][Centrality][numerator], entry->Factor * Cumulant[iAcc][object][Centrality][numerator], Cumulant[iAcc][object][Centrality][denominator]};
        auto &brick = CumulantBrick[iAcc][object][Centrality];
        for (const auto &component : Result[object])
        {
            value.push_back(brick[component]);
        }
    }
    return value;
}