#ifndef UTIL_SETUTIL_H_
#define UTIL_SETUTIL_H_

#include <cmath>
#include <ctime>
#include <bitset>
#include <vector>
#include <algorithm>
#include <iostream>
#include <map>
#include <random>
#include <functional>

struct SubsetIterator {
    SubsetIterator(int, int);
    std::bitset<64> next();
    long int nrSubs();

    long int fSeed;
    const long int fMaxSeed;
    const long int fNrSubs;
};

int binomialCoeff(int n, int k);
std::vector<std::vector<int>> cartesianProduct(const std::vector<int>&);
std::vector<std::bitset<32>> allSubsets(const int);
std::vector<std::bitset<32>> allSubsetsIncl(const int);
void subsetsLengthK(const int, const int, std::vector<std::bitset<32>>&);
std::vector<int> range(int, int, int=1);
std::vector<int> iota(int);

template <typename T>
typename T::value_type product(const T& items) {
   typename T::value_type prod{};
   for (const typename T::value_type& i : items) {
       prod *= i;
   }
   return prod;
}

template <typename T>
typename T::value_type implode(const T& collection) {
    typename T::value_type result;
    for (const auto& c : collection) {
        result = join(result, c);
    }
    return result;
}

template<typename T>
bool isSubsetOf(const T& sub, const T& super) {
    return std::includes(super.begin(), super.end(), sub.begin(), sub.end());
}

template<typename T>
bool isStrictSubsetOf(const T& sub, const T& super) {
    return sub.size() < super.size() && std::includes(super.begin(), super.end(), sub.begin(), sub.end());
}

template <typename T>
bool containsStrictSubsetOf(const T& collection, const typename T::value_type& item) {
    for (const typename T::value_type& s : collection) {
        if (s.size() < item.size() && isSubsetOf(s, item)) return true;
    }
    return false;
}

template <typename T>
bool containsSubsetOf(const T& collection, const typename T::value_type& item) {
    for (const typename T::value_type& s : collection) {
        if (isSubsetOf(s, item)) return true;
    }
    return false;
}

template <typename T>
bool containsSupersetOf(const T& collection, const typename T::value_type& item) {
    for (const typename T::value_type& s : collection) {
        if (isSubsetOf(item, s)) return true;
    }
    return false;
}

template<typename T>
T subset(const T& items, const std::bitset<64> mask) {
    T sub;
    sub.reserve(mask.count());
    for (unsigned mi = 0; mi < items.size(); mi++) {
        if (mask[mi]) {
            sub.push_back(items[mi]);
        }
    }
    return sub;
}

template<typename T>
T subset(const T& items, const typename T::value_type leaveOut) {
    if (items.size() == 1) return T();
    T sub;
    sub.reserve(items.size() - 1);
    for (unsigned mi = 0; mi < items.size(); mi++) {
        if (items[mi] != leaveOut) {
            sub.push_back(items[mi]);
        }
    }
    return sub;
}

template<typename T>
void insertSorted(const T& in, const typename T::value_type item, T& out) {
    auto it = in.begin();
    while (it != in.end() && *it < item) { it++; }
    if (it != in.begin()) {
        out.insert(out.begin(), in.begin(), it);
    }
    int offset = (int)(it - in.begin());
    out.insert(out.begin() + offset, item);
    if (in.size() && it != in.end()) {
        out.insert(out.begin() + offset + 1, it, in.end());
    }
}

template<typename T>
T* intersection(const T* lhs, const T* rhs) {
    T* isect = new T(lhs->size());
    auto it = std::set_intersection(lhs->begin(), lhs->end(), rhs->begin(), rhs->end(), isect->begin());
    isect->resize((int)(it - isect->begin()));
    return isect;
}

template<typename T>
T intersection(const T& lhs, const T& rhs) {
    T isect(std::min(lhs.size(), rhs.size()));
    auto it = std::set_intersection(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), isect.begin());
    isect.resize((int)(it - isect.begin()));
    return isect;
}

template<typename T>
T* setdiff(const T* lhs, const T* rhs) {
    T* res = new T();
    std::set_difference(lhs->begin(), lhs->end(), rhs->begin(), rhs->end(), std::inserter(*res, res->begin()));
    return res;
}

template<typename T>
T setdiff(const T& lhs, const T& rhs) {
    T res;
    std::set_difference(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), std::inserter(res, res.begin()));
    return res;
}

template<typename T>
T* join(const T* lhs, const T* rhs) {
    T* uni = new T(lhs->size() + rhs->size());
    auto it = std::set_union(lhs->begin(), lhs->end(), rhs->begin(), rhs->end(), uni->begin());
    uni->resize((int)(it - uni->begin()));
    return uni;
}

template<typename T>
T join(const T& lhs, const T& rhs) {
    T uni(lhs.size() + rhs.size());
    auto it = std::set_union(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), uni.begin());
    uni.resize((int)(it - uni.begin()));
    return uni;
}

template<typename T>
T join(const T& lhs, const typename T::value_type& rhsItem) {
    T res;
    res.reserve(lhs.size() + 1);
    insertSorted(lhs, rhsItem, res);
    return res;
}

template<typename T>
T* join(const T* lhs, const typename T::value_type& rhsItem) {
    T* res = new T();
    res->reserve(lhs->size() + 1);
    insertSorted(*lhs, rhsItem, *res);
    return res;
}

template<typename T>
T remove(const T& lhs, const typename T::value_type& rhsItem) {
    T res;
    res.reserve(lhs.size() - 1);
    const auto& lower = std::lower_bound(lhs.begin(), lhs.end(), rhsItem);
    res.insert(res.begin(), lhs.begin(), lower);
    res.insert(res.begin() + (lower - lhs.begin()), lower + 1, lhs.end());
    return res;
}

template <typename T>
bool contains(const T& collection, const typename T::value_type& item) {
    return std::find(collection.begin(), collection.end(), item) != collection.end();
}

template <typename T>
bool containsKey(const T& collection, const typename T::key_type& item) {
    return collection.find(item) != collection.end();
}

template <typename T>
typename T::value_type randElem(const T& collection) {
    std::mt19937 gen(time(0));
    std::uniform_int_distribution<int> dis(0, collection.size()-1);
    return collection[dis(gen)];
}

template <typename T>
T randSubset(const T& collection, int size) {
    std::vector<int> indices = iota(collection.size());
    std::mt19937 gen(time(0));
    std::shuffle(indices.begin(), indices.end(), gen);
    T res(size);
    for (int i = 0; i < size; i++) {
        res[i] = collection[indices[i]];
    }
    return res;
}

template <typename T>
std::vector<T> randPartitions(const T& collection, int nr) {
    std::vector<int> indices = iota(collection.size());
    std::mt19937 gen(time(0));
    std::shuffle(indices.begin(), indices.end(), gen);
    std::vector<T> res(nr);
    int size = collection.size() / nr;
    int j = 0;
    int block = 0;
    for (int i = 0; i < collection.size(); i++) {
        res[j].push_back(collection[indices[i]]);
        block++;
        if (block == size) {
            std::sort(res[j].begin(), res[j].end());
            block = 0;
            j++;
        }
    }
    return res;
}

template <typename T>
bool has(const T& collection, std::function<bool(typename T::value_type)> f) {
    for (const auto& i : collection) {
        if (f(i)) {
            return true;
        }
    }
    return false;
}

template <typename T>
T where(const T& collection, std::function<bool(typename T::value_type)> f) {
    T res;
    for (const auto& i : collection) {
        if (f(i)) {
            res.push_back(i);
        }
    }
    return res;
}

template<typename T>
void shuffle(T& collection) {
    std::mt19937 gen(time(0));
    std::shuffle(collection.begin(), collection.end(), gen);
    return;
}

template<typename T>
void shuffle(T& collection, unsigned seed) {
    std::mt19937 gen(seed);
    std::shuffle(collection.begin(), collection.end(), gen);
    return;
}

template <typename T, typename S>
T projection(const T& collection, const S& indices) {
    T res;
    for (int i : indices) {
        res.push_back(collection[i]);
    }
    return res;
}

template<typename T>
T sorted(const T& collection) {
    T copy = collection;
    std::sort(copy.begin(), copy.end());
    return copy;
}

template<typename T>
double precisionAt(const T& found, const T& expected, int at) {
    int prec = 0;
    for (int i = 0; i < at; i++) {
        int expPos = 0;
        for (; expPos <= at; expPos++) {
            if (expPos < at && found[i] == expected[expPos]) break;
        }
        prec += std::abs(i - expPos);
    }
    return prec;
}

template<typename T>
std::vector<T> split(const T& collection, const typename T::value_type& spl) {
    std::vector<T> res(1);
    for (int i : collection) {
        if (i == spl) {
            res.push_back(T());
        }
        else {
            res.back().push_back(i);
        }
    }
    res.pop_back();
    return res;
}

template <typename T>
T getMaxElem(const std::vector<std::pair<T, int> >& collection) {
    int max = -1;
    int maxI = -1;
    for (int i = 0; i < collection.size(); i++) {
        const auto& elem = collection[i];
        if (elem.second > max) {
            max = elem.second;
            maxI = i;
        }
    }
    return collection[maxI].first;
}

#endif //UTIL_SETUTIL_H_
