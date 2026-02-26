#define URQMD_COMPONENT_CXX
#include <algorithm>
#include "Helper.h"
#include "Component.h"
#include "ThreadPool.h"
Document::Document() : SpeciesSet{3122, 3312}
{
    int index = 0;
    for (const auto &element : SpeciesSet)
    {
        IndexMap[element] = index++;
    }
    StringMap = {
        {2212, "Proton"},
        {-2212, "Pbar"},
        {321, "Kplus"},
        {-321, "Kminus"},
        {221, "Pi+"},
        {-221, "Pi-"},
        {3122, "Lambda"},
        {-3122, "LambdaBar"},
        {3212, "Sigma0"},
        {-3212, "Sigma0Bar"},
        {3312, "Xi"},
        {-3312, "XiBar"},
        {3334, "Omega-"},
        {-3334, "AOmega+"},
        {3322, "Xi0"},
        {-3322, "AXi0"},
        {3222, "Sigma+"},
        {-3222, "ASigma+"},
        {3112, "Sigma-"},
        {-3112, "ASigma-"},
        {2112, "Neutron"},
        {-2112, "ANeutron"},
        {311, "K0"},
        {-311, "AK0"},
        {1, "Baryon"},
        {2, "Strangeness"},
        {3, "Charge"}};
    SyntaxMap= {
            {2212, "p"},
            {-2212, "#bar{p}"},
            {321, "K^{+}"},
            {-321, "K^{-}"},
            {221, "#pi^{+}"},
            {-221, "#pi^{-}"},
            {3122, "#Lambda"},
            {-3122, "#bar{#Lambda}"},
            {3212, "#Sigma^{0}"},
            {-3212, "#bar{#Sigma^{0}}"},
            {3312, "#Xi"},
            {-3312, "#bar{#Xi}"},
            {3334, "#Omega^{-}"},
            {-3334, "#bar{#Omega}^{+}"},
            {3322, "#Xi^{0}"},
            {-3322, "#bar{#Xi}^{0}"},
            {3222, "#Sigma^{+}"},
            {-3222, "#bar{#Sigma}^{+}"},
            {3112, "#Sigma^{-}"},
            {-3112, "#bar{#Sigma}^{-}"},
            {2112, "n"},
            {-2112, "#bar{n}"},
            {311, "K^{0}"},
            {-311, "#bar{K}^{0}"},
            {1, "B"},
            {2, "S"},
            {3, "Q"}};
}
Component::Component(const std::pair<int, int> &pair)
{
    if (pair.first)
    {
        Internal.insert({1, pair.first});
    }
    if (pair.second)
    {
        Internal.insert({2, pair.second});
    }
}
Component::Component(const std::vector<int> &species)
{
    for (const auto &element : species)
    {
        if (!element)
        {
            Internal = {};
            continue;
        }
        Internal[element]++;
    }
}
Component::Component(const std::initializer_list<int> &list)
{
    for (const auto &element : list)
    {
        if (!element)
        {
            Internal = {};
            continue;
        }
        Internal[element]++;
    }
}

Component &Component::operator*=(const Component &rhs)
{
    for (const auto &element : rhs.Internal)
    {
        if (Internal.find(element.first) == Internal.end())
        {
            if (element.second)
                Internal[element.first] = element.second;
        }
        else
        {
            Internal[element.first] += element.second;
            if (Internal[element.first] == 0)
            {
                auto it = Internal.find(element.first);
                if (it != Internal.end())
                {
                    Internal.erase(it);
                }
            }
        }
    }
    return *this;
}

Component &Component::operator/=(const Component &rhs)
{
    Component temp;
    for (const auto &element : rhs.Internal)
    {
        if (Internal.find(element.first) == Internal.end())
        {
            Internal[element.first] = -1. * element.second;
        }
        else
        {
            Internal[element.first] -= element.second;
            if (Internal[element.first] == 0)
            {
                auto it = Internal.find(element.first);
                if (it != Internal.end())
                {
                    Internal.erase(it);
                }
            }
        }
    }
    return *this;
}

std::string Component::GetString() const
{
    std::string syntax = "";
    auto &BSMap = Document::getInstance();
    for (const auto &pair : Internal)
    {
        if (pair.second)
        {
            syntax += BSMap.StringMap[pair.first];
            syntax += std::to_string(pair.second);
        }
    }
    return syntax;
}

std::string Component::GetSyntax() const
{
    std::string syntax = "";
    auto &BSMap = Document::getInstance();
    for (const auto &pair : Internal)
    {
        if (!pair.second)
            continue;
        syntax += pair.second > 1 ? (BSMap.SyntaxMap[pair.first] + "^{" + std::to_string(pair.second) + "}") : BSMap.SyntaxMap[pair.first];
    }
    return syntax;
}

int ComponentHelper::GetOrder(const Component &component)
{
    if (OrderMap.find(component) != OrderMap.end())
    {
        return OrderMap[component];
    }
    int order = 0;
    for (const auto &element : component.Internal)
    {
        order += element.second;
    }
    OrderMap[component] = order;
    return OrderMap[component];
}

std::vector<int> &ComponentHelper::Expand(const Component &component)
{
    if (ExpandMap.find(component) != ExpandMap.end())
    {
        return ExpandMap[component];
    }
    ExpandMap[component].resize(6, 0);
    size_t Index = 0;
    for (const auto &element : component.Internal)
    {
        for (int ip = 0; ip < element.second; ip++)
        {
            ExpandMap[component][Index++] = element.first;
        }
    }
    return ExpandMap[component];
}

std::unordered_set<Component, ComponentHash> &ComponentHelper::Subset(const Component &component)
{
    if (SubsetMap.find(component) != SubsetMap.end())
    {
        return SubsetMap[component];
    }
    auto &Expansion = Expand(component);
    auto &partition = Partition::instance();
    auto &sets = partition.getPartitions(GetOrder(component));
    auto &subset = SubsetMap[component];
    for (const auto &set : sets)
    {
        for (const auto &power : set)
        {
            std::vector<int> expansion(power.size());
            std::transform(power.begin(), power.end(), expansion.begin(), [&](int idx)
                           { return Expansion[idx]; });
            subset.insert(Component(expansion));
        }
    }
    return SubsetMap[component];
}