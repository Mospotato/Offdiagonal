#ifndef SERIES_HH
#define SERIES_HH
#include <vector>
#include <cmath>
#include <unordered_map>
#include "Helper.h"
#include "Component.h"
using Particle = std::pair<int, int>;
struct Series
{
    Series() {};
    Series(const Component &component);
    Series(int n, const std::vector<Particle> &);
    Series(Series &, Series &);
    Series(int, const std::vector<Particle> &, int, const std::vector<Particle> &);
    Series(const Series &rhs) : Components(rhs.Components) {};
    // Move constructor
    Series(Series &&rhs) { *this = std::move(rhs); };
    // Copy assignment operator
    Series &operator=(const Series &rhs)
    {
        if (this != &rhs)
        {
            Components = rhs.Components;
        }
        return *this;
    }
    // Move assignment operator
    Series &operator=(Series &&rhs)
    {
        if (this != &rhs)
        {
            Components = std::move(rhs.Components);
        }
        return *this;
    }
    Series &operator+=(const Series &rhs);
    Series operator+(const Series &rhs) const
    {
        Series result = *this;
        result += rhs;
        return result;
    }
    Series &operator*=(const Series &);
    Series operator*(const Series &rhs) const
    {
        Series result = *this;
        result *= rhs;
        return result;
    }
    using Constituent = std::unordered_map<Component, double, ComponentHash>;
    Constituent Components;
};

struct Product
{
    std::pair<Series, double> Internal;
    Product() {};
    ~Product() {};
    Product(const Component &, const Component &, double);
    Product(const Product &rhs) : Internal(rhs.Internal) {};
    Product(Product &&rhs) { *this = std::move(rhs); };
    Product &operator=(const Product &rhs)
    {
        if (this != &rhs)
        {
            Internal = rhs.Internal;
        }
        return *this;
    }
    Product &operator=(Product &&rhs)
    {
        if (this != &rhs)
        {
            Internal = std::move(rhs.Internal);
        }
        return *this;
    }
};

struct DeriveHelper
{
    std::unordered_map<std::pair<Component, Component>, double, ComponentPairHash> CoefficientMap;
    std::unordered_map<std::pair<Component, Component>, Product, ComponentPairHash> ProductMap;
    static DeriveHelper *getInstance()
    {
        static DeriveHelper instance;
        return &instance;
    }
    DeriveHelper() {};
    ~DeriveHelper() {};
    Product &Derive(const Component &component, const Component &);
    double &Coefficient(const std::pair<Component, Component>&);
};
#endif