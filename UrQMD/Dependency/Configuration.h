#ifndef URQMD_CONFIGURATION_HH
#define URQMD_CONFIGURATION_HH
#include <map>
#include <vector>
#include <string>
struct Configuration
{
    double yMax;
    double ptMin;
    double ptMax;

    std::vector<double> Npart;
    std::vector<int> CentVec;
    std::vector<std::string> CentDef;
    int nCent;
    int PtSize;
    std::vector<float> SystematicPtVec;
    int RapiditySize;
    std::vector<float> SystematicRapidityVec;
    int nAcceptance;
    std::map<int, double> BeamRapidityMap;
    
    static Configuration &getInstance()
    {
        static Configuration instance;
        return instance;
    }
    Configuration();
    ~Configuration() {};
    int GetCentrality(const int &);
    void SetCentrality(const int & Energy);
    std::vector<int> GetRapidityVec(const double &y);
    std::vector<int> GetAcceptanceVec(const double &pt, const double &y);
    bool IsAccepted(const double &pt, const double &y);
};

#endif