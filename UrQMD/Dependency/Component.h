#ifndef URQMD_COMPONENT_HH
#define URQMD_COMPONENT_HH
#include <map>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
struct Document
{
    std::vector<int> SpeciesSet;
    static Document &getInstance()
    {
        static Document instance;
        return instance;
    }
    Document();
    ~Document() {};
    std::map<int, std::string> StringMap;
    std::map<int, std::string> SyntaxMap;
    std::map<int, std::string> IndexMap;
};

struct Component
{
    std::map<int, int> Internal;
    Component() : Internal() {}
    Component(const Component &rhs) : Internal(rhs.Internal) {}
    Component(Component &&rhs) : Internal(std::move(rhs.Internal)) {}
    Component(const std::vector<int> &species);
    Component(const std::pair<int, int> &);
    Component(const std::initializer_list<int> &list);
    size_t GetHash() const
    {
        std::size_t seed = Internal.size();
        for (const auto &pair : Internal)
        {
            seed ^= std::hash<int>{}(pair.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= std::hash<int>{}(pair.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

    std::size_t size() const
    {
        return Internal.size();
    }
    std::string GetString() const;
    std::string GetSyntax() const;

    Component &operator=(const Component &rhs)
    {
        if (this != &rhs)
        {
            Internal = rhs.Internal;
        }
        return *this;
    }
    Component &operator=(Component &&rhs)
    {
        Internal = std::move(rhs.Internal);
        return *this;
    }

    bool operator<(const Component &rhs) const
    {
        return this->GetHash() < rhs.GetHash();
    }
    bool operator==(const Component &rhs) const
    {
        return this->GetHash() == rhs.GetHash();
    }
    bool operator!=(const Component &rhs) const
    {
        return this->GetHash() != rhs.GetHash();
    }
    int operator[](int type) const
    {
        return Internal.find(type) != Internal.end() ? Internal.at(type) : 0;
    }
    int &operator[](int type)
    {
        return Internal[type];
    }
    Component &operator*=(const Component &rhs);
    Component operator*(const Component &rhs) const
    {
        Component result = *this;
        result *= rhs;
        return result;
    }

    Component &operator/=(const Component &rhs);
    Component operator/(const Component &rhs) const
    {
        Component result = *this;
        result /= rhs;
        return result;
    }
};
struct ComponentHash
{
    size_t operator()(const Component &component) const
    {
        return component.GetHash();
    }
};

struct ComponentPairHash
{
    size_t operator()(const std::pair<Component, Component> &pair) const
    {
        size_t hash1 = pair.first.GetHash();
        size_t hash2 = pair.second.GetHash();
        return hash1 ^ (hash2 + 0x9e3779b9 + (hash1 << 6) + (hash1 >> 2));
    }
};

struct ComponentHelper
{
    static ComponentHelper *getInstance()
    {
        static ComponentHelper instance;
        return &instance;
    }
    ComponentHelper() {};
    ~ComponentHelper() {};
    int GetOrder(const Component &component);
    std::vector<int> &Expand(const Component &component);
    std::unordered_set<Component, ComponentHash> &Subset(const Component &component);

    std::unordered_map<Component, int, ComponentHash> OrderMap;
    std::unordered_map<Component, std::vector<int>, ComponentHash> ExpandMap;
    std::unordered_map<Component, std::unordered_set<Component, ComponentHash>, ComponentHash> SubsetMap;
};
#endif