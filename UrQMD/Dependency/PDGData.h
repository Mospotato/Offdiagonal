#ifndef PDGDATA_HH_
#define PDGDATA_HH_
#include <map>
#include <vector>
#include <string>
struct ParticleEntry
{
    int Baryon;
    int Strangeness;
    int Charge;
    double Mass;
    ParticleEntry() : Baryon(0), Strangeness(0), Charge(0), Mass(0.0) {};
    // Constructor
    ParticleEntry(int b, int s, int q, double m) : Baryon(b), Strangeness(s), Charge(q), Mass(m) {};
    // Copy Constructor
    ParticleEntry(const ParticleEntry &rhs) : Baryon(rhs.Baryon), Strangeness(rhs.Strangeness), Charge(rhs.Charge), Mass(rhs.Mass) {};
    // Copy assignment operator
    ParticleEntry &operator=(const ParticleEntry &rhs)
    {
        Baryon = rhs.Baryon;
        Strangeness = rhs.Strangeness;
        Charge = rhs.Charge;
        Mass = rhs.Mass;
        return *this;
    }
    // Move assignment operator
    ParticleEntry &operator=(ParticleEntry &&rhs)
    {
        Baryon = rhs.Baryon;
        Strangeness = rhs.Strangeness;
        Charge = rhs.Charge;
        Mass = rhs.Mass;
        return *this;
    }
    // Move constructor
    ParticleEntry(ParticleEntry &&rhs)
    {
        *this = std::move(rhs);
    }
    ~ParticleEntry() {};
};

struct PDGData
{
    static PDGData &getInstance()
    {
        static PDGData instance;
        return instance;
    }
    PDGData();
    ~PDGData() {}
    std::vector<ParticleEntry> EntryVec;
    std::map<int, ParticleEntry *> PDGMap;
};
#endif
