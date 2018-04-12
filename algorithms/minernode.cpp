#include "minernode.h"

int support(const SimpleTidList& tids) {
    return tids.size();
}

int support(const PartitionTidList& tids) {
    if (tids.fTids.size())
        return tids.fTids.size() + 1 - tids.fNrSets;
    return 0;
}

int hash(const SimpleTidList& tids) {
    return std::accumulate(tids.begin(), tids.end(), 1);
}

int hash(const PartitionTidList& tids) {
    return std::accumulate(tids.fTids.begin(), tids.fTids.end(), 1) + (tids.fNrSets - 1);
}
