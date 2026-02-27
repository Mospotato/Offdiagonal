#define PDGDATA_CXX_
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "PDGData.h"
PDGData::PDGData()
{
    std::ifstream fin("list-all.dat");
    if (!fin.is_open())
    {
        std::cerr << "Error: Particle File not found!\n";
        exit(EXIT_FAILURE);
    }
    std::string token;
    std::getline(fin, token);
    int nParticle = 0;
    EntryVec.reserve(337);
    while (std::getline(fin, token))
    {
        std::istringstream iss(token);
        int pdgid, stable, stat, str, bary, charge, charm;
        double mass, degeneracy;
        std::string name;
        if (iss >> pdgid >> name >> stable >> mass >> degeneracy >> stat >> bary >> charge >> str >> charm)
        {
            EntryVec.emplace_back(ParticleEntry(bary, str, charge, mass));
            PDGMap[pdgid] = &EntryVec[nParticle];
            nParticle++;
            // cout << pdgid << " " << name << " " << mass << " " << bary << " " << str << "\n";
        }
    }
    fin.close();
}
int GetConservedValue(const ParticleEntry *entry, int type)
{
    switch (type)
    {
    case 1:
        return entry->Baryon;
    case 2:
        return entry->Strangeness;
    case 3:
        return entry->Charge;
    default:
        return 0;
    }
}
std::pair<int, int> GetConservedTypes(OffdiagonalType type)
{
    switch (type)
    {
    case OffdiagonalType::kBQ:
        return {1, 3};
    case OffdiagonalType::kQS:
        return {3, 2};
    case OffdiagonalType::kBS:
        return {1, 2};
    default:
        throw std::invalid_argument("Invalid OffdiagonalType");
    }
}
std::string GetConservedString(OffdiagonalType type)
{
    switch (type)
    {
    case OffdiagonalType::kBS:
        return "BS";
    case OffdiagonalType::kQS:
        return "QS";
    case OffdiagonalType::kBQ:
        return "BQ";
    default:
        throw std::invalid_argument("Invalid OffdiagonalType");
    }
}