#define EPD_CENTRALITY_ESTIMATOR_CXX
#include "EPDSimulator.h"
#include "TRandom3.h"
#include <cmath>
#include <numeric>
int EPDSimulator::GetTileIndex(const TLorentzVector &vector) const
{
    static std::array<double, 17> EtaBoundary{5.09, 4.42, 4.03, 3.74, 3.47, 3.26, 3.08, 2.94, 2.81, 2.69, 2.59, 2.50, 2.41, 2.34, 2.27, 2.20, 2.14};
    double abeta = fabs(vector.Eta());
    if (abeta > EtaBoundary[0] || abeta < EtaBoundary[16])
    {
        return -1; // Out of EPD acceptance
    }
    int RingIndex{0};
    for (int ir = 0; ir < 16; ++ir)
    {
        if (abeta > EtaBoundary[ir + 1] && abeta < EtaBoundary[ir])
        {
            RingIndex = ir;
            break;
        }
    }
    double PhiWidth = TMath::Pi() / (RingIndex ? 24 : 12);
    double phi = vector.Phi();
    if (phi < 0)
    {
        phi += 2 * TMath::Pi();
    }
    int SegmentIndex = std::floor(phi / PhiWidth);
    int TileIndex = 0;
    if (RingIndex)
    {
        TileIndex = SegmentIndex + RingIndex * 48 - 24;
    }
    else
    {
        TileIndex = SegmentIndex;
    }
    if (TileIndex < 0 || TileIndex >= 744)
    {
        throw std::out_of_range("Calculated tile index is out of range: " + std::to_string(TileIndex));
    }
    return TileIndex;
}

void EPDSimulator::FillEPD(const TLorentzVector &vector)
{
    int tileIndex = GetTileIndex(vector);
    if (tileIndex == -1)
    {
        return; // Particle is out of EPD acceptance
    }
    MIP[tileIndex] += gRandom->Landau(1, 0.2); // Simulate MIP response with a Landau distribution
    int MaxMIP = 7; // Maximum MIP value to improve correlation with RefMult
    if (MIP[tileIndex] > MaxMIP)
    {
        MIP[tileIndex] = MaxMIP;
    }
    return;
}

double EPDSimulator::GetnMIP() const
{
    return std::accumulate(MIP.begin(), MIP.end(), 0.0);
}