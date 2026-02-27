#ifndef DEP_BSDICTIONARY_HH
#define DEP_BSDICTIONARY_HH
#include <set>
#include <unordered_set>
#include "Series.h"
#include "PDGData.h"
struct Dictionary 
{
    static Dictionary *getInstance()
    {
        static Dictionary instance;
        return &instance;
    }
    Dictionary();
    ~Dictionary() {};
    int nConstituent;
    std::vector<std::vector<int>> Powers;
    std::set<Component> Constituent;
};
#endif