#ifndef EPD_CENTRALITY_ESTIMATOR_HH
#define EPD_CENTRALITY_ESTIMATOR_HH
#include "TF1.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include <array>
struct EPDSimulator
{
    static EPDSimulator& getInstance()
    {
        static EPDSimulator instance;
        return instance;
    }
    void Reset()
    {
        MIP.fill(0);
    }
    void FillEPD(const TLorentzVector &vector);
    int GetTileIndex(const TLorentzVector &vector) const;
    double GetnMIP() const;
private:
    EPDSimulator(){};
    ~EPDSimulator(){};
    std::array<double, 744> MIP; // Fixed size for EPD tiles
};


#endif