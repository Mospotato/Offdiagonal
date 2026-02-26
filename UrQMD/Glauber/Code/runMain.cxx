#include <string>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TSystem.h"
#include "StNbdFitMaker.h"
#include "StNegativeBinomial.h"
void MakeConfiguration(const std::string filename)
{
    std::ofstream file(filename);
    file << "nppBin 10\n";
    file << "nppMin 0.0\n";
    file << "nppMax 10.0\n";
    file << "kBin 10\n";
    file << "kMin 0.0\n";
    file << "kMax 10.0\n";
    file << "xBin 10\n";
    file << "xMin 0.0\n";
    file << "xMax 1.0\n";
    file << "eBin 10\n";
    file << "eMin 0.5\n";
    file << "eMax 1.0\n";
    file.close();
    printf("Configuration file created\n");
    return;
}
std::map<std::string, double> readConfig(const std::string& filename) {
    std::map<std::string, double> config;
    std::ifstream file(filename);
    if (!file.is_open())
    {
        MakeConfiguration(filename);
        file.open(filename);
    }
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key;
        double value;
        if (iss >> key >> value) {
            config[key] = value;
        }
    }
    return config;
}
int main(int argc, char **argv)
{
    if (argc < 2)
    {
        printf("Usage: %s [Energy]\n", argv[0]);
        return 1;
    }
    Int_t Energy = std::stoi(argv[1]);
    Bool_t IsConstantEfficiency = argc >= 3 ? std::stoi(argv[2]) : true;
    gSystem->mkdir(Form("File/%d", Energy), kTRUE);
    StNbdFitMaker maker;
    maker.SetOutputName(Form("File/%d/glauber.root", Energy));
    maker.SetMinimumMultiplicityCut(120);
    //=================================================
    TFile *fData = TFile::Open(Form("../ReadUrQMD/File/%d.root", Energy));
    TFile *fGlauber = TFile::Open(Form("NcollNpart/%d.root", Energy));
    if (!fData || !fGlauber)
    {
        printf("Cannot open input files\n");
        delete fData;
        delete fGlauber;
        return 1;
    }
    auto config = readConfig("Configuration.txt");
    maker.ReadData(fData, fGlauber);
    maker.Scan(1000000,
                static_cast<int>(config["nppBin"]), config["nppMin"], config["nppMax"],
                static_cast<int>(config["kBin"]), config["kMin"], config["kMax"],
                static_cast<int>(config["xBin"]), config["xMin"], config["xMax"],
                static_cast<int>(config["eBin"]), config["eMin"], config["eMax"], 
                1.0, IsConstantEfficiency);
    delete fData;
    delete fGlauber;
    return 0;
}