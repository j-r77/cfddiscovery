#include "setutil.h"
#include <algorithm>
#include <numeric>

SubsetIterator::SubsetIterator(int size, int subSize)
:fSeed((1 << subSize) - 1), fMaxSeed(1 << size), fNrSubs(binomialCoeff(size, subSize)) {

}

std::bitset<64> SubsetIterator::next() {
    int x = fSeed;
    int u = x & (-x);
    int v = x + u;
    fSeed = v + (((v ^ x) / u) >> 2);
    return std::bitset<64>(x);
}

long int SubsetIterator::nrSubs() {
    return fNrSubs;
}

std::vector<int> range(int min, int max, int step) {
    std::vector<int> res;
    res.reserve((max - min)/step);
    for (int i = min; i != max; i += step) {
        res.push_back(i);
    }
    return res;
}

std::vector<int> iota(int max) {
    std::vector<int> iotas(max);
    std::iota(iotas.begin(), iotas.end(), 0);
    return iotas;
}

int binomialCoeff(int n, int k)
{
    int res = 1;
 
    // Since C(n, k) = C(n, n-k)
    if ( k > n - k )
        k = n - k;
 
    // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
    for (int i = 0; i < k; ++i)
    {
        res *= (n - i);
        res /= (i + 1);
    }
 
    return res;
}

std::vector<std::vector<int>> cartesianProduct(const std::vector<int>& sizes) {
    std::vector<int> counters(sizes.size());
    std::vector<std::vector<int>> cartprod;
    cartprod.reserve(product(sizes));

    while (counters[0] < sizes[0]) {
        std::vector<int> row(sizes.size());
        for (unsigned i = 0; i < sizes.size(); i++) {
            row[i] = counters[i];
        }
        int incIx = counters.size() - 1;
        counters[incIx]++;
        while (incIx > 0 && counters[incIx] == sizes[incIx]) {
            counters[incIx] = 0;
            counters[--incIx]++;
        }
        cartprod.push_back(row);
    }
    return cartprod;
}

void subsetsLengthK(const int size, const int k, std::vector<std::bitset<32>>& subs) {
    const int maxB = pow(2, size);
    int x = pow(2, k) - 1;
    while (x < maxB) {
        subs.push_back(std::bitset<32>(x));
        int u = x & (-x);
        int v = x + u;
        x = v + (((v ^ x) / u) >> 2);
    }
}

std::vector<std::bitset<32>> allSubsets(const int size) {
    std::vector<std::bitset<32>> subs;
    subs.reserve(pow(2, size)-2);
    for (int i = 1; i < size; i++) {
        subsetsLengthK(size, i, subs);
    }
    return subs;
}

std::vector<std::bitset<32>> allSubsetsIncl(const int size) {
    std::vector<std::bitset<32>> subs;
    subs.reserve(pow(2, size)-2);
    for (int i = 1; i <= size; i++) {
        subsetsLengthK(size, i, subs);
    }
    return subs;
}
