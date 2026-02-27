#define URQMD_CONFIGURATION_CXX
#include <map>
#include <cmath>
#include "Configuration.h"
#include "TError.h"
Configuration::Configuration() : yMax(0.5), ptMin(0.4), ptMax(1.6)
{
    BeamRapidityMap = {{3, 1.05}, {3.2, 1.13}};
    SystematicPtVec = {1.2, 1.3, 1.4, 1.5};
    PtSize = static_cast<int>(SystematicPtVec.size());
    SystematicRapidityVec = {0.1, 0.2, 0.3, 0.4};
    RapiditySize = static_cast<int>(SystematicRapidityVec.size());
    // nAcceptance = PtSize + RapiditySize + 1;
    nAcceptance = 1;
}

bool Configuration::IsAccepted(const double &pt, const double &y)
{
    return pt > ptMin && pt < ptMax && fabs(y) < yMax;
}

int Configuration::GetCentrality(const int &Refmult)
{
    if (Refmult >= CentVec.front() || Refmult < CentVec.back())
    {
        return -1;
    }
    for (int iCent = 0; iCent < nCent; iCent++)
    {
        if (Refmult < CentVec[iCent] && Refmult >= CentVec[iCent + 1])
        {
            return iCent;
        }
    }
    return -1;
}

void Configuration::SetCentrality(const int &Energy)
{
    static std::map<int, std::vector<int>> CentVecMap{
        {3, {118, 64, 53, 37, 26, 18, 12, 7, 4, 2}},
        {7, {409, 283, 233, 158, 105, 67, 41, 24, 13, 6}},
        {11, {587, 382, 313, 212, 140, 89, 54, 31, 16, 8}},
        {14, {640, 442, 362, 244, 161, 102, 61, 35, 18, 9}},
        {19, {802, 526, 431, 289, 190, 121, 72, 40, 21, 10}},
        {27, {903, 597, 488, 327, 215, 135, 80, 45, 23, 11}},
        {39, {1000, 682, 557, 372, 243, 152, 90, 50, 25, 12}},
        {62, {1194, 784, 641, 427, 277, 173, 102, 56, 28, 13}},
        {200, {1413, 945, 770, 508, 328, 202, 119, 65, 32, 15}},
        {218, {3100, 2466, 2167, 1639, 1195, 836, 558, 353, 210, 117}}, // RuRu 200 GeV
        {219, {3100, 2466, 2167, 1639, 1195, 836, 558, 353, 210, 117}}, // RuRu 200 GeV
    };
    CentVec = CentVecMap[Energy];
    nCent = CentVec.size() - 1;
    CentDef = nCent == 9 ? std::vector<std::string>{"0-5%", "5-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "0-80%"} : std::vector<std::string>{"0-2.5%", "2.5-5%", "5-7.5%", "7.5-10%", "10-12.5%", "12.5-15%", "15-17.5%", "17.5-20%", "20-22.5%", "22.5-25%", "25-27.5%", "27.5-30%", "30-32.5%", "32.5-35%", "35-37.5%", "37.5-40%", "40-42.5%", "42.5-45%", "45-47.5%", "47.5-50%", "50-52.5%", "52.5-55%", "55-57.5%", "60-62.5%", "62.5-65%", "65-67.5%", "67.5-70%", "70-72.5%", "72.5-75%", "75-77.5%", "77.5-80%"};
    Info("Configuration", "Centrality Definition for Energy %d GeV Loaded!", Energy);
}
std::vector<int> Configuration::GetRapidityVec(const double &y)
{
    std::vector<int> RapidityVec;
    RapidityVec.reserve(RapiditySize + 1);
    RapidityVec.push_back(0);
    for (int iy = 0; iy < RapiditySize; iy++)
    {
        if (std::fabs(y) < SystematicRapidityVec[iy])
        {
            for (int iSubY = iy + 1; iSubY <= RapiditySize; iSubY++)
            {
                RapidityVec.push_back(iSubY);
            }
            break;
        }
    }
    return RapidityVec;
}

std::vector<int> Configuration::GetAcceptanceVec(const double &pt, const double &y)
{
    const double absY = std::fabs(y);
    std::vector<int> AcceptanceVec;
    AcceptanceVec.reserve(nAcceptance);
    AcceptanceVec.push_back(0);
    return AcceptanceVec;
    for (int iy = 0; iy < RapiditySize; iy++)
    {
        if (absY < SystematicRapidityVec[iy])
        {
            for (int iSubY = iy + 1; iSubY <= RapiditySize; iSubY++)
            {
                AcceptanceVec.push_back(iSubY);
            }
            break;
        }
    }
    for (int iPt = 0; iPt < PtSize; iPt++)
    {
        if (pt < SystematicPtVec[iPt])
        {
            for (int iSubPt = iPt + RapiditySize + 1; iSubPt < nAcceptance; iSubPt++)
            {
                AcceptanceVec.push_back(iSubPt);
            }
            break;
        }
    }
    return AcceptanceVec;
}