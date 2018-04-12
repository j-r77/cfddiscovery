#include "types.h"

const int PartitionTidList::SEP = -1;

bool PartitionTidList::operator==(const PartitionTidList& b) const {
    return fNrSets == b.fNrSets && fTids == b.fTids;
}

bool PartitionTidList::operator!=(const PartitionTidList& b) const {
    return fNrSets != b.fNrSets || fTids != b.fTids;
}

bool PartitionTidList::operator<(const PartitionTidList& b) const {
    return lessthan(*this, b);
}

Itemset itemset(int i) {
    Itemset res(1);
    res[0] = i;
    return res;
}

PartitionTidList convert(const SimpleTidList& tids) {
    return { tids, 1 };
}

SimpleTidList convert(const PartitionTidList& tids) {
    auto res = tids.fTids;
    if (tids.fNrSets > 1) {
        std::sort(res.begin(), res.end());
        res.erase(res.begin(), res.begin()+tids.fNrSets-1);
    }
    return res;
}

bool lessthan(const PartitionTidList& lhs, const PartitionTidList& rhs) {
    return lhs.fNrSets < rhs.fNrSets || (lhs.fNrSets == rhs.fNrSets && lhs.fTids < rhs.fTids);
}

bool lessthan(const SimpleTidList& lhs, const SimpleTidList& rhs) {
    return lhs.size() > rhs.size();
}

bool equals(const PartitionTidList& lhs, const PartitionTidList& rhs) {
    return lhs.fNrSets == rhs.fNrSets && lhs.fTids == rhs.fTids;
}