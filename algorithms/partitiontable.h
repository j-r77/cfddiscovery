#ifndef ALGORITHMS_PARTITIONTABLE_H_
#define ALGORITHMS_PARTITIONTABLE_H_

#include <vector>
#include <utility>
#include <unordered_map>
#include <map>
#include <set>
#include <string>
#include "../data/types.h"

class PartitionTable {
public:
    static int fDbSize;

    static PartitionTidList intersection(const PartitionTidList &lhs, const PartitionTidList &rhs);
    static std::vector<PartitionTidList> intersection(const PartitionTidList &lhs, const std::vector<PartitionTidList*> rhses);
    static std::vector<PartitionTidList> intersection(const PartitionTidList &lhs, const std::vector<const PartitionTidList*> rhses);
    static std::vector<std::pair<int, std::vector<int> > > partitionMap(const PartitionTidList& x, const PartitionTidList& xa);
    static std::vector<int> partitionMap2(const PartitionTidList& x, const PartitionTidList& xa);
    static int partitionError(const std::vector<std::vector<int> >& x, const std::vector<std::vector<int> >& xa);
    static int partitionError(const PartitionTidList&, const PartitionTidList&);
    static int constify(const PartitionTidList&, const PartitionTidList&);
    static SimpleTidList violations(const PartitionTidList&, const PartitionTidList&);
    static SimpleTidList violations(const PartitionTidList&, const PartitionTidList&, int& e);
    static bool violatedInCleaned(const PartitionTidList& x, const PartitionTidList& xa, int);
};

#endif //ALGORITHMS_PARTITIONTABLE_H_