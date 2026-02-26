#define STATISTICS_HELPER_CXX
#include <math.h>
#include "Helper.h"
double Factorial::getFactorial(int n)
{
    if (Memory.find(n) != Memory.end())
    {
        return Memory[n];
    }
    double result = 1;
    for (int i = 1; i <= n; ++i)
    {
        result *= i;
    }
    Memory[n] = result;
    return result;
}

double Binomial::getBinomial(int n, int k)
{
    std::pair<int, int> key = std::make_pair(n, k);
    if (Memory.find(key) != Memory.end())
    {
        return Memory[key];
    }
    auto *factorial = Factorial::getInstance();
    double result = factorial->getFactorial(n) / (factorial->getFactorial(k) * factorial->getFactorial(n - k));
    Memory[key] = result;
    return result;
}

double Combination::Multinomial(int n, const std::vector<int> &powers)
{
    double result = 1;
    for (int power : powers)
    {
        result *= Factorial::getInstance()->getFactorial(power);
    }
    result = Factorial::getInstance()->getFactorial(n) / result;
    return result;
}

void Combination::generateCombinations(int n, int k, std::vector<int> &current, std::vector<std::vector<int>> &result)
{
    if (k == 0)
    {
        if (n == 0)
        {
            result.push_back(current);
        }
        return;
    }

    for (int i = 0; i <= n; ++i)
    {
        current.push_back(i);
        generateCombinations(n - i, k - 1, current, result);
        current.pop_back();
    }
}

std::vector<std::vector<int>> &Combination::getCombinations(int Order, int nElements)
{
    std::vector<int> current;
    std::pair<int, int> key = std::make_pair(Order, nElements);
    if (Memory.find(key) != Memory.end())
    {
        return Memory[key];
    }
    std::vector<std::vector<int>> combinations;
    generateCombinations(Order, nElements, current, combinations);
    Memory[key] = combinations;
    return Memory[key];
}

std::vector<std::vector<int>> Partition::concat(const std::vector<std::vector<int>> &first, const std::vector<int> &second)
{
    std::vector<std::vector<int>> result = first;
    result.push_back(second);
    return result;
}

std::vector<std::vector<std::vector<int>>> Partition::getAllPartitions(
    const std::vector<std::vector<int>> &fixedParts, std::vector<int> &suffixElements)
{
    std::vector<std::vector<std::vector<int>>> result;
    result.push_back(concat(fixedParts, suffixElements));
    auto suffixPartitions = getTuplePartitions(suffixElements);
    for (auto &suffixPartition : suffixPartitions)
    {
        auto subPartitions = getAllPartitions(concat(fixedParts, suffixPartition.first), suffixPartition.second);
        result.insert(result.end(), subPartitions.begin(), subPartitions.end());
    }
    return result;
}

std::vector<std::pair<std::vector<int>, std::vector<int>>> Partition::getTuplePartitions(std::vector<int> &elements)
{
    std::vector<std::pair<std::vector<int>, std::vector<int>>> result;
    if (elements.size() < 2)
        return result;
    size_t n = elements.size();
    for (size_t pattern = 1; pattern < static_cast<size_t>(1 << (n - 1)); ++pattern)
    {
        std::vector<int> set1 = {elements[0]};
        std::vector<int> set2;
        for (size_t index = 1; index < n; ++index)
        {
            if ((pattern >> (index - 1)) & 1)
            {
                set2.push_back(elements[index]);
            }
            else
            {
                set1.push_back(elements[index]);
            }
        }
        result.push_back({set1, set2});
    }
    return result;
}

std::vector<std::vector<std::vector<int>>> &Partition::getPartitions(size_t n)
{
    if (Memory.find(n) != Memory.end())
    {
        return Memory[n];
    }
    std::vector<int> elements;
    for (int i = 0; i < static_cast<int>(n); i++)
    {
        elements.push_back(i);
    }
    Memory[n] = getAllPartitions(elements);
    return Memory[n];
}

double &Partition::getCoefficient(size_t nv)
{
    if (CoefficientMap.find(nv) != CoefficientMap.end())
    {
        return CoefficientMap[nv];
    }
    auto *factorial = Factorial::getInstance();
    CoefficientMap[nv] = pow(-1, nv - 1) * factorial->getFactorial(nv - 1);
    return CoefficientMap[nv];
}