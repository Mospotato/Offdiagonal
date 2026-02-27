#define DEP_DICTIONARY_CXX
#include <iostream>
#include <array>
#include "Helper.h"
#include "Dictionary.h"
Dictionary::Dictionary()
{
    auto *combination = Combination::getInstance();
    for (int i = 1; i <= 4; i++)
    {
        PowerVec vec = combination->getCombinations(i, 2);
        Powers.insert(Powers.end(), vec.begin(), vec.end());
    }
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
    auto GetConstituents = [&]()
    {
        for (const auto &power : Powers)
        {
            for (auto &pair : CombinationVec)
            {
                Series Left(power[0], *BaseArrays[pair[0]]);
                Series Right(power[1], *BaseArrays[pair[1]]);
                Series Product(Left, Right);
                for (const auto &Component : Product.Components)
                {
                    Constituent.insert(Component.first);
                }
            }
        }
    };
    GetConstituents();
    BaryonArrary.clear();
    StrangeArrary.clear();
    ChargeArrary.clear();
    BaryonArrary = {{1, 1}};
    StrangeArrary = {{2, 1}};
    ChargeArrary = {{3, 1}};
    GetConstituents();
    nConstituent = static_cast<int>(Constituent.size());
    printf("Info in <Dictionary::Instance>: %d Conserved Components will be Loaded!\n", nConstituent);
}