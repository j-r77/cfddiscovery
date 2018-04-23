#ifndef CFDD_GENERATORSTORE_H
#define CFDD_GENERATORSTORE_H

#include <map>
#include "../util/prefixtree.h"
#include "../util/setutil.h"

template<typename T>
class GeneratorStore {
public:
    bool addMinGen(const Itemset& newset, int supp, const T& hash) {
        const auto& genIt = fGenMap.find(hash);
        if (genIt != fGenMap.end()) {
            if (genIt->second.hasStrictSubset(newset, supp)) {
                return false;
            }
            genIt->second.insert(newset, supp);
        }
        else {
            fGenMap[hash].insert(newset, supp);
        }
        return true;
    }

    bool addMinGen(const Itemset& newset, int supp, const T& hash, std::vector<Itemset>& subs) {
        const auto& genIt = fGenMap.find(hash);
        if (genIt != fGenMap.end()) {
            if (genIt->second.hasSubset(newset, supp)) {
                const auto& ptSubs = genIt->second.getSets();
                subs.insert(subs.begin(), ptSubs.begin(), ptSubs.end());
                return false;
            }
            genIt->second.insert(newset, supp);
        }
        else {
            fGenMap[hash].insert(newset, supp);
        }
        return true;
    }

    bool isMinGen(const Itemset& newset, int supp, const T& hash) {
        const auto& genIt = fGenMap.find(hash);
        if (genIt != fGenMap.end()) {
            if (genIt->second.hasStrictSubset(newset, supp)) {
                return false;
            }
        }
        return true;
    }

    std::vector<Itemset> getMinGens(const Itemset& newset, int supp, const T& hash) {
        const auto& genIt = fGenMap.find(hash);
        if (genIt != fGenMap.end()) {
            return genIt->second.getSubsets(newset, supp);
        }
    }
private:
    std::map<T, PrefixTree<Itemset,int>> fGenMap;
};

#endif //CFDD_GENERATORSTORE_H
