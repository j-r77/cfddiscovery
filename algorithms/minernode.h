#ifndef ALGORITHMS_MINERNODE_H_
#define ALGORITHMS_MINERNODE_H_

#include "../data/types.h"
#include <numeric>

int support(const SimpleTidList& tids);
int support(const PartitionTidList& tids);
int hash(const SimpleTidList& tids);
int hash(const PartitionTidList& tids);

template<typename T>
struct MinerNode {
    MinerNode()
        : fItem(-1), fSupp(-1), fHash(-1) {

    }

    MinerNode(const Item& item, int supp)
        : fItem(item), fSupp(supp), fHash(0) {
            //tidmap.reserve(supp);
    }

    MinerNode(const Item& item)
        : fItem(item), fSupp(0), fHash(0) {

    }

    MinerNode(const Item& item, const T& tids)
        : fItem(item), fTids(tids), fSupp(0), fHash(0) {
        hashTids();
        fSupp = supp();
    }

    MinerNode(const Item& item, const T& tids, int supp)
        :fItem(item), fTids(tids), fSupp(supp), fHash(0), fPrefix() {
        hashTids();
    }

    MinerNode(const Item& item, const T& tids, int supp, const Itemset& prefix)
        :fItem(item), fTids(tids), fSupp(supp), fHash(0), fPrefix(prefix) {
        hashTids();
    }

    MinerNode(const Item& item, const T& tids, const Itemset& prefix)
        :fItem(item), fTids(tids), fSupp(0), fHash(0), fPrefix(prefix) {
        hashTids();
        fSupp = supp();
    }

    MinerNode(const Item& item, const T& tids, int supp, int hash)
        :fItem(item), fTids(tids), fSupp(supp), fHash(hash) {

    }

    bool operator< (const MinerNode<T>& rhs) const {
        //return supp() < rhs.supp() || (supp() == rhs.supp() && tidmap.size() > rhs.tidmap.size());
        //return fCands.size() > rhs.fCands.size();
        return lessthan(fTids, rhs.fTids);
    }

    bool operator== (const MinerNode<T>& rhs) const {
        return fItem == rhs.fItem && fPrefix == rhs.fPrefix && equals(fTids, rhs.fTids);
    }

    int supp() const {
        if (fSupp) return fSupp;
        return support(fTids);
    }

    int resupp() {
        fSupp = support(fTids);
        return fSupp;
    }

    void hashTids() {
        fHash = hash(fTids);
    }

    Item fItem;
    T fTids;
    int fSupp;
    int fHash;
    Itemset fPrefix;
    Itemset fCands;
};

typedef std::vector<MinerNode<PartitionTidList> >::const_iterator MNIter;

#endif //ALGORITHMS_MINERNODE_H_
