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
