#include <iostream>
#include "getopt.h"
#include "Calculate.h"
int main(int argc, char *argv[])
{
    Int_t opt{0}, Energy{0};
    std::string USAGE{"Usage: ./run -e Energy"};
    while ((opt = getopt(argc, argv, "e:")) != -1)
    {
        switch (opt)
        {
        case 'e':
            Energy = std::stoi(optarg);
            break;
        default:
            std::cerr << "Usage: " << argv[0] << " -e [Energy]\n";
            return 1;
        }
    }
    if (!Energy)
    {
        std::cerr << "Abort! " << USAGE << std::endl;
        return 1;
    }
    Calculate Calc;
    if (Calc.Init(Energy))
    {
        return 1;
    }
    Calc.Process();
    Calc.Terminate();
    return 0;
}