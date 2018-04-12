#ifndef ALGORITHMS_CTANE_H_
#define ALGORITHMS_CTANE_H_

#include "baseminer.h"

enum SUBSTRATEGY {
    DFS, BFS
};

class CTane : public BaseMiner {
public:
    CTane(Database&);
    int nrCFDs() const;
    CFDList getCFDs() const;
    void mine(int, int, double=1);
    void mineFree(int, double=1);
    void mineFreeDepth(int, int, double=1);
    void mineFreeDepth(const Itemset&, std::vector<MinerNode<PartitionTidList> >&);
    void mineItemsets(const Itemset&, std::vector<MinerNode<SimpleTidList> >, const std::vector<MinerNode<PartitionTidList> >&);
    Itemset getConstantCFDs(const Itemset&, const SimpleTidList&, const std::vector<MinerNode<PartitionTidList> >&);
    void getCFDs(const Itemset&, const std::vector<MinerNode<PartitionTidList> >&);
    void getVariableCFDs(const Itemset&, std::vector<MinerNode<PartitionTidList> >&, std::vector<std::pair<Itemset, int> >&);
    void pruneCands(std::vector<MinerNode<PartitionTidList> >& items, const Itemset& sub, int out);
    void pruneCands(std::vector<MinerNode<SimpleTidList> >& items, const Itemset& sub, int out);
        std::vector<std::pair<Itemset,SimpleTidList> > getItemsetLayer();
    bool isConstRule(const PartitionTidList& items, int rhsA);
    bool isConstRulePart(const SimpleTidList& items, const Itemset& rhses);
    bool isConstRulePart2(const SimpleTidList& items, const std::vector<std::vector<std::pair<int,int> > >& rhses);
    bool precedes(const Itemset& a, const Itemset& b);
    void mineItemsetsFirst(int, int, SUBSTRATEGY=SUBSTRATEGY::BFS, double=1);
    void mineItemsetsFirstDepth(int, int, SUBSTRATEGY=SUBSTRATEGY::BFS, double=1);
    void mineItemsetsFirstDepth(const Itemset&, std::vector<MinerNode<SimpleTidList> >&, SUBSTRATEGY=SUBSTRATEGY::BFS);
    void mineFDsDepth(const SimpleTidList&, const Itemset&, const Itemset&);
    void mineFDsDepth(std::vector<MinerNode<PartitionTidList> >&, const Itemset&, const Itemset&);
    void mineFDsFirst(int, int, SUBSTRATEGY=SUBSTRATEGY::BFS, double=1);
    void mineFDsFirstDepth(int, int, SUBSTRATEGY=SUBSTRATEGY::BFS, double=1);
    void mineFDsFirstDepth(const Itemset&, std::vector<MinerNode<PartitionTidList> >&, SUBSTRATEGY=SUBSTRATEGY::BFS);
    void mineFDs(const SimpleTidList&, const Itemset&, const Itemset&);
    void mineTPs(const Itemset& lhs, int rhs, const PartitionTidList& lhsTids, const PartitionTidList& allTids);
    void mineTPsDepth(const Itemset& lhs, int rhs, const PartitionTidList& lhsTids, const PartitionTidList& allTids);
    void mineTPsDepth(const Itemset&, std::vector<MinerNode<SimpleTidList> >&, const Itemset&, int,
                      std::vector<std::vector<std::pair<int,int> > >&, std::vector<std::pair<Itemset, std::vector<int> > >&, std::vector<int>&);
    int getPartitionSupport(const SimpleTidList& pids, const std::vector<int>& partitions);
    int getPartitionError(const SimpleTidList& pids, const std::vector<std::pair<Itemset, std::vector<int> > >& partitions);

private:
    GeneratorStore<int> fGens;
    std::map<Itemset,PartitionTidList> fStore;
    std::unordered_map<Itemset,SimpleTidList> fTidStore;
    PrefixTree<Itemset, Itemset> fCandStore;
    Itemset fAllAttrs;
    std::vector<MinerNode<SimpleTidList> > fItemsLayer;
    double fMinConf;
    std::map<std::pair<int,int>, std::vector<Itemset> > fFreeMap;
    std::set<Itemset> fFreeItemsets;
    std::unordered_map<int, std::vector<Itemset> > fRules;
    Itemset fAttOrder;
public:
    int fVisited;
    int fValidated;
};

#endif // ALGORITHMS_CTANE_H_