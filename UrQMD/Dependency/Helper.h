#ifndef STATISTICS_HELPER_HH
#define STATISTICS_HELPER_HH
#include <map>
#include <vector>
#include <unordered_map>
struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &pair) const
    {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};
using PowerVec = std::vector<std::vector<int>>;
using MemoryMap = std::unordered_map<std::pair<int, int>, std::vector<std::vector<int>>, pair_hash>;

struct Combination
{
    static Combination *getInstance()
    {
        static Combination instance;
        instance.Memory[std::make_pair(0, 0)] = {{}};
        return &instance;
    }
    MemoryMap Memory;
    PowerVec &getCombinations(int Order, int nElements);
    double Multinomial(int n, const std::vector<int> &powers);
    void generateCombinations(int n, int k, std::vector<int> &current, std::vector<std::vector<int>> &result);
};

struct Factorial
{
    static Factorial *getInstance()
    {
        static Factorial instance;
        return &instance;
    }
    std::unordered_map<int, double> Memory;
    double getFactorial(int n);
};

struct Binomial
{
    static Binomial *getInstance()
    {
        static Binomial instance;
        return &instance;
    }
    std::unordered_map<std::pair<int, int>, double, pair_hash> Memory;
    double getBinomial(int n, int k);
};

struct Partition
{
private:
    Partition() = default;
    std::map<size_t, std::vector<std::vector<std::vector<int>>>> Memory;
    std::map<size_t, double> CoefficientMap;
    std::vector<std::vector<int>> concat(const std::vector<std::vector<int>> &first, const std::vector<int> &second);
    std::vector<std::vector<std::vector<int>>> getAllPartitions(std::vector<int> &elements)
    {
        return getAllPartitions({}, elements);
    }
    std::vector<std::vector<std::vector<int>>> getAllPartitions(
        const std::vector<std::vector<int>> &fixedParts, std::vector<int> &suffixElements);

    std::vector<std::pair<std::vector<int>, std::vector<int>>> getTuplePartitions(std::vector<int> &elements);
public:
    static Partition &instance()
    {
        static Partition instance;
        return instance;
    }
    std::vector<std::vector<std::vector<int>>> &getPartitions(size_t n);
    double &getCoefficient(size_t nv);
    Partition(const Partition &) = delete;
    Partition &operator=(const Partition &) = delete;
};
#endif