// STL
#include <set>
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_set>
// ROOT Master Header
#include "TROOT.h"
// User Defined Library
#include "KochRatio.h"
#include "Helper.h"
#include "Series.h"
#include "PDGData.h"
R__LOAD_LIBRARY(libDependency.so)
void bsdictionary()
{
    auto *combination = Combination::getInstance();
    std::vector<std::vector<int>> Powers;
    for (int i = 1; i <= 4; i++)
    {
        PowerVec vec = combination->getCombinations(i, 2);
        Powers.insert(Powers.end(), vec.begin(), vec.end());
    }
    std::ofstream LogFile("log.txt");
    auto &BSMap = BSDocument::getInstance();
    Int_t nSpecies = static_cast<int>(BSMap.IndexMap.size());
    auto &PDG = PDGData::getInstance();
    std::vector<Particle> BaryonArrary{{2212, 1}, {-2212, -1}};
    BaryonArrary.reserve(2 * nSpecies);
    std::vector<Particle> StrangeArrary{{321, 1}, {-321, -1}};
    StrangeArrary.reserve(2 * nSpecies);
    for (const auto &pdg : BSMap.SpeciesSet)
    {
        if (PDG.PDGMap[pdg]->Baryon)
        {
            BaryonArrary.emplace_back(std::make_pair(pdg, PDG.PDGMap[pdg]->Baryon));
            BaryonArrary.emplace_back(std::make_pair(-pdg, -PDG.PDGMap[pdg]->Baryon));
        }
        if (PDG.PDGMap[pdg]->Strangeness)
        {
            StrangeArrary.emplace_back(std::make_pair(pdg, PDG.PDGMap[pdg]->Strangeness));
            StrangeArrary.emplace_back(std::make_pair(-pdg, -PDG.PDGMap[pdg]->Strangeness));
        }
    }
    std::set<Component> Individual;
    for (const auto &power : Powers)
    {
        Series B(power[0], BaryonArrary);
        Series S(power[1], StrangeArrary);
        Series Bi(B, S);
        LogFile << Form("(%d, %d) Number of Components: %zu\n", power[0], power[1], Bi.Components.size());
        for (auto &component : Bi.Components)
        {
            Individual.insert(component.first);
            LogFile << "Coeff: " << component.second << " Component: " << component.first.GetSyntax() << '\n';
        }
    }
    int iComponent = 1;
    for (const auto &Component : Individual)
    {
        LogFile << iComponent++ << " " << Component.GetSyntax() << '\n';
    }
    LogFile << Individual.size() << " Components Generated!\n";
    BaryonArrary.clear();
    StrangeArrary.clear();
    BaryonArrary.emplace_back(std::make_pair(100, 1));
    StrangeArrary.emplace_back(std::make_pair(200, 1));
    std::set<Component> Conserved;
    for (const auto &power : Powers)
    {
        Series B(power[0], BaryonArrary);
        Series S(power[1], StrangeArrary);
        Series Bi(B, S);
        LogFile << Form("(%d, %d) Number of Components: %zu\n", power[0], power[1], Bi.Components.size());
        for (const auto &pair : Bi.Components)
        {
            Conserved.insert(pair.first);
            LogFile << "Coeff: " << pair.second << " Component: " << pair.first.GetSyntax() << '\n';
        }
    }
    iComponent = 1;
    for (const auto &component : Conserved)
    {
        LogFile << iComponent++ << " " << component.GetSyntax() << '\n';
    }
    LogFile << Conserved.size() << " Components Generated!\n";
    LogFile.close();
}

void demoPartition()
{
    int n = 2;
    auto &helper = Partition::instance();
    auto partitions = helper.getPartitions(n);
    for (const auto &partition : partitions)
    {
        std::cout << "Coeff: " << helper.getCoefficient(static_cast<int>(partition.size())) << " ";
        for (const auto &part : partition)
        {
            std::cout << "{ ";
            for (const auto &elem : part)
            {
                std::cout << elem << " ";
            }
            std::cout << "} ";
        }
        std::cout << "\n";
    }
}

void demostrate()
{
    // bsdictionary();
    KochRatio calculate;
    KochRatio::Pool moment;
    moment[{321}] = 10;
    moment[{321, 321}] = 100;
    moment[{321, 321, 321}] = 900;
    Component component({321, 321, 321});
    Double_t cumulant = calculate.GetCumulant(moment, component);
    Info("KochRatio::GetCumulant", "%s: Moment: %f Cumulant: %f", component.GetSyntax().c_str(), moment.at(component), cumulant);
    return;
}