#define SERIES_CXX
#include <numeric>
#include <algorithm>
#include "Series.h"
Series::Series(int Order, const std::vector<Particle> &particleArray)
{
    int nElements = static_cast<int>(particleArray.size());
    auto *combination = Combination::getInstance();
    auto &Powers = combination->getCombinations(Order, nElements);
    for (size_t ip = 0; ip < Powers.size(); ip++)
    {
        std::vector<int> &power = Powers[ip];
        Component component;
        double sign = 1;
        for (int i = 0; i < nElements; ++i)
        {
            if (power[i])
            {
                sign *= std::pow(particleArray[i].second, power[i]);
                component[particleArray[i].first] += power[i];
            }
        }
        Components[component] = sign * combination->Multinomial(Order, power);
    }
}

Series::Series(Series &Left, Series &Right)
{
    if (Left.Components.empty() && Right.Components.empty())
    {
        *this = Series();
    }
    if (Left.Components.empty())
    {
        *this = Right;
    }
    if (Right.Components.empty())
    {
        *this = Left;
    }
    for (const auto &iComponent : Left.Components)
    {
        for (const auto &jComponent : Right.Components)
        {
            Component component = iComponent.first * jComponent.first;
            Components[component] += iComponent.second * jComponent.second;
        }
    }
}

Series::Series(const Component &component)
{
    static ComponentHelper *helper = ComponentHelper::getInstance();
    int Order = helper->GetOrder(component);
    if (Order)
    {
        Components.insert({component, 1});
    }
}

Series::Series(int iOrder, const std::vector<Particle> &bArray, int jOrder, const std::vector<Particle> &sArray)
{
    Series Left(iOrder, bArray);
    Series Right(jOrder, sArray);
    *this = Left * Right;
}

Series &Series::operator*=(const Series &rhs)
{
    if (rhs.Components.empty())
    {
        return *this;
    }
    if (this->Components.empty())
    {
        *this = rhs;
        return *this;
    }
    Series lhs;
    lhs.Components = std::move(this->Components);
    for (const auto &icomponent : lhs.Components)
    {
        for (const auto &jcomponent : rhs.Components)
        {
            Component newComponent = icomponent.first * jcomponent.first;
            Components[newComponent] += icomponent.second * jcomponent.second;
        }
    }
    return *this;
}

Series &Series::operator+=(const Series &rhs)
{
    if (rhs.Components.empty())
    {
        return *this;
    }
    for (const auto &jComponent : rhs.Components)
    {
        if (Components.find(jComponent.first) == Components.end())
        {
            Components.insert(jComponent);
        }
        else
        {
            Components[jComponent.first] += jComponent.second;
        }
    }
    return *this;
}

Product::Product(const Component &variable, const Component &subscript, double value)
{
    if (variable.Internal.empty())
    {
        *this = Product();
    }
    Series product;
    Component divided = variable / subscript;
    for (auto &pair : divided.Internal)
    {
        if (pair.second)
        {
            product.Components[{pair.first}] = pair.second;
        }
    }
    Internal = std::make_pair(product, value);
}

Product &DeriveHelper::Derive(const Component &variable, const Component &subscript)
{
    std::pair<Component, Component> key = std::make_pair(variable, subscript);
    if (ProductMap.find(key) != ProductMap.end())
    {
        return ProductMap[key];
    }
    ProductMap[key] = Product(variable, subscript, Coefficient(key));
    return ProductMap[key];
}

double &DeriveHelper::Coefficient(const std::pair<Component, Component> &key)
{
    if (CoefficientMap.find(key) != CoefficientMap.end())
    {
        return CoefficientMap[key];
    }
    auto *binomial = Binomial::getInstance();
    auto *helper = ComponentHelper::getInstance();
    Component divided = key.first / key.second;
    CoefficientMap[key] = std::pow(-1, helper->GetOrder(divided));
    for (const auto &element : key.first.Internal)
    {
        CoefficientMap[key] *= binomial->getBinomial(element.second, key.second[element.first]);
    }
    return CoefficientMap[key];
}