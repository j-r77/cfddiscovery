#ifndef DATA_TYPES_H_
#define DATA_TYPES_H_

#include <vector>
#include <utility>
#include <unordered_map>
#include <map>
#include <set>
#include <string>
#include "../util/setutil.h"

#ifndef uint
typedef unsigned int uint;
#endif
typedef std::vector<int> Transaction;
typedef int Item;
typedef std::vector<Item> SimpleTidList;
typedef std::vector<Item> Itemset;
typedef std::vector<int> Diffset;

struct PartitionTidList {
    PartitionTidList()
        :fNrSets(0){

    }
    PartitionTidList(const SimpleTidList& tids, int nrSets)
            :fTids(tids), fNrSets(nrSets) {

    }
    SimpleTidList fTids;
    int fNrSets;
    static const int SEP;// = -1;
    bool operator==(const PartitionTidList&) const;
    bool operator!=(const PartitionTidList&) const;
    bool operator<(const PartitionTidList&) const;
};

Itemset itemset(int);
PartitionTidList convert(const SimpleTidList&);
SimpleTidList convert(const PartitionTidList&);
bool lessthan(const PartitionTidList&, const PartitionTidList&);
bool equals(const PartitionTidList&, const PartitionTidList&);
bool lessthan(const SimpleTidList&, const SimpleTidList&);

struct pairhash {
public:
	template <typename T1, typename T2> 
	size_t operator()(const std::pair<T1,T2>& p) const
	{
		return std::hash<T1>()(p.first) ^ std::hash<T2>()(p.second);
	}
};

namespace std {
template <typename T1, typename T2>  struct hash<pair<T1,T2> >
{
    size_t operator()(const pair<T1,T2>& p) const
    {
        return hash<T1>()(p.first) ^ hash<T2>()(p.second);
    }
};
}

namespace std {
template <typename T> struct hash<vector<T> >
{
    size_t operator()(const vector<T>& xs) const
    {
        size_t res = 0;
        for (const T& x : xs) {
            res ^= hash<T>()(x) + 0x9e3779b9 + (res << 6) + (res >> 2);
        }
        return res;
    }
};
}
#endif //DATA_TYPES_H_