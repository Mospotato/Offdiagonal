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
void dictionary()
{
    auto *combination = Combination::getInstance();
    std::vector<std::vector<int>> Powers;
    for (int i = 1; i <= 4; i++)
    {
        PowerVec vec = combination->getCombinations(i, 2);
        Powers.insert(Powers.end(), vec.begin(), vec.end());
    }
    std::ofstream LogFile("log.txt");
    auto &BSMap = Document::getInstance();
    auto &PDG = PDGData::getInstance();
    auto GenerateArray = [&](std::vector<int> PDGVec, int type)
    {
        std::vector<Particle> Array;
        Array.reserve(2 * PDGVec.size());
        for (const auto &pdg : PDGVec)
        {
            int value = GetConservedValue(PDG.PDGMap.at(pdg), type);
            Array.push_back(std::make_pair(pdg, value));
            Array.push_back(std::make_pair(-pdg, -value));
        }
        return Array;
    };
    //====================================================================
    // Generate the collection of all possible components
    auto BaryonArrary = GenerateArray(std::vector<int>{2212, 3122}, 1);
    auto StrangeArrary = GenerateArray(std::vector<int>{3122, 321}, 2);
    auto ChargeArrary = GenerateArray(std::vector<int>{2212, 211, 321}, 3);
    std::vector<std::vector<Particle> *> BaseArrays{&BaryonArrary, &StrangeArrary, &ChargeArrary};
    std::vector<std::array<int, 2>> CombinationVec{{0, 1}, {0, 2}, {2, 1}};
    std::set<Component> Constituent;
    auto GetConstituents = [&]()
    {
        for (const auto &power : Powers)
        {
            for (int ic = 0; ic < CombinationVec.size(); ic++)
            {
                auto pair = CombinationVec[ic];
                Series Left(power[0], *BaseArrays[pair[0]]);
                Series Right(power[1], *BaseArrays[pair[1]]);
                Series Product(Left, Right);
                LogFile << Form("%s_(%d, %d) Number of Components: %zu\n", GetConservedString(static_cast<OffdiagonalType>(ic)).c_str(), power[0], power[1], Product.Components.size());
                for (const auto &Component : Product.Components)
                {
                    Constituent.insert(Component.first);
                    LogFile << "Coeff: " << Component.second << " Component: " << Component.first.GetSyntax() << '\n';
                }
            }
        }
    };
    GetConstituents();
    int iComponent = 1;
    for (const auto &Component : Constituent)
    {
        LogFile << iComponent++ << " " << Component.GetSyntax() << '\n';
    }
    LogFile << Constituent.size() << " Components Generated!\n";
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
    dictionary();
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