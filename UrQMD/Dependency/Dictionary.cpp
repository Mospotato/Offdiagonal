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
    auto &document = Document::getInstance();
    nSpecies = static_cast<int>(document.IndexMap.size());
    auto &PDG = PDGData::getInstance();
    std::vector<Particle> BaryonArrary{{2212, 1}, {-2212, -1}};
    BaryonArrary.reserve(2 * nSpecies);
    std::vector<Particle> StrangeArrary{{321, 1}, {-321, -1}};
    StrangeArrary.reserve(2 * nSpecies);
    std::vector<Particle> ChargeArrary{{211, 1}, {-211, -1}};
    ChargeArrary.reserve(2 * nSpecies);
    for (const auto &pdg : document.SpeciesSet)
    {
        if (PDG.PDGMap[pdg]->Baryon)
        {
            BaryonArrary.push_back(std::make_pair(pdg, PDG.PDGMap[pdg]->Baryon));
            BaryonArrary.push_back(std::make_pair(-pdg, -PDG.PDGMap[pdg]->Baryon));
        }
        if (PDG.PDGMap[pdg]->Strangeness)
        {
            StrangeArrary.push_back(std::make_pair(pdg, PDG.PDGMap[pdg]->Strangeness));
            StrangeArrary.push_back(std::make_pair(-pdg, -PDG.PDGMap[pdg]->Strangeness));
        }
        if (PDG.PDGMap[pdg]->Charge)
        {
            ChargeArrary.push_back(std::make_pair(pdg, PDG.PDGMap[pdg]->Charge));
            ChargeArrary.push_back(std::make_pair(-pdg, -PDG.PDGMap[pdg]->Charge));
        }
    }
    //====================================================================
    // Generate the collection of all possible components
    std::vector<std::vector<Particle> *> BaseArrays{&BaryonArrary, &StrangeArrary, &ChargeArrary};
    std::vector<std::array<int, 2>> CombinationVec{{0, 1}, {0, 2}, {1, 2}};
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