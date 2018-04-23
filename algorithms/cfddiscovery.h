#ifndef ALGORITHMS_CTANE_H_
#define ALGORITHMS_CTANE_H_

#include <map>
#include <vector>
#include "minernode.h"
#include "../util/prefixtree.h"
#include "partitiontable.h"
#include "../data/database.h"
#include "../util/setutil.h"
#include "../data/cfd.h"
#include "generatorstore.h"

enum SUBSTRATEGY {
    DFS, BFS
};

class CFDDiscovery {
public:
    CFDDiscovery(Database&);
    int nrCFDs() const;
    CFDList getCFDs() const;
    void ctane(int minsup, int maxsize, double conf = 1);
    void integratedDFS(int minsup, int maxsize, double conf = 1);
    void itemsetsFirstBFS(int minsup, int maxsize, SUBSTRATEGY= SUBSTRATEGY::BFS, double conf = 1);
    void itemsetsFirstDFS(int minsup, int maxsize, SUBSTRATEGY= SUBSTRATEGY::BFS, double conf = 1);
    void fdsFirstBFS(int minsup, int maxsize, SUBSTRATEGY= SUBSTRATEGY::BFS, double conf = 1);
    void fdsFirstDFS(int minsup, int maxsize, SUBSTRATEGY= SUBSTRATEGY::BFS, double conf = 1);

    void integratedDFS(const Itemset &, std::vector<MinerNode<PartitionTidList> > &);
    void itemsetsFirstDFS(const Itemset &, std::vector<MinerNode<SimpleTidList> > &, SUBSTRATEGY= SUBSTRATEGY::BFS);
    void mineFDs(const SimpleTidList&, const Itemset&, const Itemset&);
    void mineFDsDFS(const SimpleTidList &, const Itemset &, const Itemset &);
    void mineFDsDFS(std::vector<MinerNode<PartitionTidList> > &, const Itemset &, const Itemset &);
    void fdsFirstDFS(const Itemset &, std::vector<MinerNode<PartitionTidList> > &, SUBSTRATEGY= SUBSTRATEGY::BFS);
    void minePatternsBFS(const Itemset &lhs, int rhs, const PartitionTidList &allTids);
    void minePatternsDFS(const Itemset &lhs, int rhs, const PartitionTidList &allTids);
    void minePatternsDFS(const Itemset &, std::vector<MinerNode<SimpleTidList> > &, const Itemset &, int,
                         std::vector<std::vector<std::pair<int, int> > > &,
                         std::vector<std::pair<Itemset, std::vector<int> > > &, std::vector<int> &);

    std::vector<MinerNode<PartitionTidList> > getPartitionSingletons();
    std::vector<MinerNode<PartitionTidList> > getPartitionSingletons(const SimpleTidList&, const Itemset&);
    std::vector<MinerNode<PartitionTidList> > getAllSingletons(int);
    std::vector<MinerNode<SimpleTidList> > getSingletons(int);

protected:
    bool precedes(const Itemset& a, const Itemset& b);
    void pruneCands(std::vector<MinerNode<PartitionTidList> >& items, const Itemset& sub, int out);
    void pruneCands(std::vector<MinerNode<SimpleTidList> >& items, const Itemset& sub, int out);
    bool isConstRule(const PartitionTidList& items, int rhsA);
    bool isConstRulePartition(const SimpleTidList &items, const std::vector<std::vector<std::pair<int, int> > > &rhses);
    int getPartitionSupport(const SimpleTidList& pids, const std::vector<int>& partitions);
    int getPartitionError(const SimpleTidList& pids, const std::vector<std::pair<Itemset, std::vector<int> > >& partitions);

private:
    int fMinSup;
    int fMaxSize;
    double fMinConf;
    Database& fDb;
    CFDList fCFDs;

    GeneratorStore<int> fGens;
    std::map<Itemset,PartitionTidList> fStore;
    std::unordered_map<Itemset,SimpleTidList> fTidStore;
    PrefixTree<Itemset, Itemset> fCandStore;
    Itemset fAllAttrs;
    std::vector<MinerNode<SimpleTidList> > fItemsLayer;
    std::map<std::pair<int,int>, std::vector<Itemset> > fFreeMap;
    std::set<Itemset> fFreeItemsets;
    std::unordered_map<int, std::vector<Itemset> > fRules;
    Itemset fAttOrder;
};

#endif // ALGORITHMS_CTANE_H_