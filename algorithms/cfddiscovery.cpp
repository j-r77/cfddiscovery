#include "cfddiscovery.h"
#include "../util/output.h"
#include "partitiontable.h"
#include <iterator>

CFDDiscovery::CFDDiscovery(Database& db)
    :fDb(db) {
}

int CFDDiscovery::nrCFDs() const {
    return fCFDs.size();
}

CFDList CFDDiscovery::getCFDs() const {
    return fCFDs;
}

bool CFDDiscovery::precedes(const Itemset& a, const Itemset& b) {
    if (a.size() > b.size()) return false;
    if (a == b) return false;
    Itemset pattern(fDb.nrAttrs());
    for (int v : b) {
        if (v > 0) pattern[fDb.getAttrIndex(v)] = v;
        else pattern[-1-v] = v;
    }
    for (int i : a) {
        if (i > 0) {
            int bv = pattern[fDb.getAttrIndex(i)];
            if (bv != i) {
                return false;
            }
        }
        else {
            if (pattern[-1-i] == 0) return false;
        }
    }
    return true;
}

void CFDDiscovery::pruneCands(std::vector<MinerNode<PartitionTidList> >& items, const Itemset& sub, int out) {
    Itemset cset = join(sub, out);
    for (int i = 0; i < items.size(); i++) {
        MinerNode<PartitionTidList> &inode = items[i];
        const Itemset iset = join(inode.fPrefix, inode.fItem);
        bool applies = false;
        auto nAttrs = fDb.getAttrVectorItems(iset);
        if (out < 0) {
            applies = std::binary_search(iset.begin(), iset.end(), out);
        }
        else {
            applies = std::binary_search(iset.begin(), iset.end(), out) || std::binary_search(iset.begin(), iset.end(), -1-fDb.getAttrIndex(out));
        }
        if (applies && precedes(sub, iset)) {
            if (out < 0)
                inode.fCands = where(inode.fCands, [this,out](int i){return fDb.getAttrIndex(i) != fDb.getAttrIndex(out);});
            else
                inode.fCands = where(inode.fCands, [this,out](int i){return i < 0 || fDb.getAttrIndex(i) != fDb.getAttrIndex(out);});
            if (inode.fCands.size() && sub.size()) inode.fCands = intersection(inode.fCands, sub);
        }
    }
}

void CFDDiscovery::pruneCands(std::vector<MinerNode<SimpleTidList> >& items, const Itemset& sub, int out) {
    for (int i = 0; i < items.size(); i++) {
        MinerNode<SimpleTidList> &inode = items[i];
        const Itemset iset = join(inode.fPrefix, inode.fItem);
        if (isSubsetOf(sub, iset)) {
            inode.fCands = where(inode.fCands, [this,out](int i){return i < 0 || fDb.getAttrIndex(i) != fDb.getAttrIndex(out);});
        }
    }
}

void CFDDiscovery::ctane(int minsup, int maxSize, double minconf) {
    fMaxSize = maxSize;
    fMinSup = minsup;
    fMinConf = minconf;
    fAllAttrs = range(-fDb.nrAttrs(), 0);
    PartitionTable::fDbSize = fDb.size();
    std::vector<MinerNode<PartitionTidList> > items;
    items = getPartitionSingletons();
    auto constants = getSingletons(fMinSup);
    for (const auto& cp : constants) {
        items.emplace_back(cp.fItem, convert(cp.fTids));
        fAllAttrs.push_back(cp.fItem);
    }
    for (auto& a : items) {
        Itemset at;
        for (int cat : fAllAttrs) {
            if (a.fItem < 0 && cat >= 0) continue;
            if (a.fItem != cat && fDb.getAttrIndex(a.fItem) == fDb.getAttrIndex(cat)) continue;
            at.push_back(cat);
        }
        a.fCands = at;
    }
    fCandStore = PrefixTree<Itemset, Itemset>();
    fGens.addMinGen(Itemset(), support(convert(iota(fDb.size()))), 0);
    fStore[Itemset()] = convert(iota(fDb.size()));
    fCandStore.insert(Itemset(), fAllAttrs);

    while (!items.empty()) {
        std::map<Itemset,PartitionTidList> newStore;
        PrefixTree<Itemset, Itemset> newCandStore;
        std::vector<std::vector<std::pair<Itemset,int> > > next(items.size());

        for (int i = 0; i < items.size(); i++) {
            MinerNode<PartitionTidList>& inode = items[i];
            const Itemset iset = join(inode.fPrefix, inode.fItem);
            auto insect = intersection(iset, inode.fCands);
            if (iset.size() > 1) {
                for (int out : insect) {
                    Itemset sub = subset(iset, out);
                    // Discard Variable -> Constant and Constant -> Variable CFDs
                    if (out < 0) {
                        if (sub.size() && !has(sub, [](int si) -> bool { return si < 0; })) continue;
                        if (inode.fTids.fNrSets == 1) continue;
                    }
                    else {
                        if (!sub.size() || has(sub, [](int si) -> bool { return si < 0; })) continue;
                    }
                    auto storedSub = fStore.find(sub);
                    if (storedSub == fStore.end()) {
                        continue;
                    }
                    if (support(storedSub->second) < fMinSup) continue;
                    double e = (out < 0) ? PartitionTable::partitionError(storedSub->second, inode.fTids)
                                         : setdiff(storedSub->second.fTids, inode.fTids.fTids).size();
                    double conf = 1 - (e / support(storedSub->second));
                    if (conf >= fMinConf) {
                        fCFDs.emplace_back(sub, out);
                    }
                    if (conf >= 1) {
                        inode.fCands = intersection(inode.fCands, sub);
                        pruneCands(items, sub, out);
                    }
                }
            }
            if (inode.fCands.empty()) continue;
            if (iset.size() == fMaxSize) continue;
            newStore[iset] = inode.fTids;

            auto nodeAttrs = fDb.getAttrVector(iset);
            for (int j = i+1; j < items.size(); j++) {
                if (j == i) continue;
                const auto& jnode = items[j];
                if (jnode.fPrefix != inode.fPrefix) continue;
                if (std::binary_search(nodeAttrs.begin(), nodeAttrs.end(), fDb.getAttrIndex(jnode.fItem))) continue;
                next[i].emplace_back(iset, j);
            }
        }
        for (int i = 0; i < items.size(); i++) {
            MinerNode<PartitionTidList> &inode = items[i];
            const Itemset iset = join(inode.fPrefix, inode.fItem);
            newCandStore.insert(iset, inode.fCands);
        }

        std::vector<MinerNode<PartitionTidList> > suffix;
        for (int i = 0; i < items.size(); i++) {
            std::vector<PartitionTidList*> expands;
            std::vector<MinerNode<PartitionTidList> > tmpSuffix;
            for (auto& newsetTup : next[i]) {
                int j = newsetTup.second;
                Itemset newset = join(newsetTup.first, items[j].fItem);
                auto c = items[i].fCands;
                for (int zz : newset) {
                    auto zsub = subset(newset, zz);
                    auto storedSub = newStore.find(zsub);
                    if (storedSub == newStore.end()) {
                        c.clear();
                        break;
                    }
                    auto subCandsPtr = newCandStore.find(zsub);
                    if (subCandsPtr) {
                        const Itemset& subCands = *subCandsPtr;
                        if (subCands.size()) c = intersection(c, subCands);
                        else c.clear();
                    }
                    else {
                        c.clear();
                    }
                    if (c.empty()) {
                        break;
                    }
                }
                if (c.size()) {
                    expands.push_back(&items[j].fTids);
                    tmpSuffix.emplace_back(items[j].fItem);
                    tmpSuffix.back().fCands = c;
                    tmpSuffix.back().fPrefix = newsetTup.first;
                }
            }
            const auto exps = PartitionTable::intersection(items[i].fTids, expands);
            for (int e = 0; e < exps.size(); e++) {
                Itemset newset = join(tmpSuffix[e].fPrefix, tmpSuffix[e].fItem);
                if (support(exps[e]) >= (fMinSup*fMinConf)) {
                    suffix.emplace_back(tmpSuffix[e].fItem, exps[e]);
                    suffix.back().fCands = tmpSuffix[e].fCands;
                    suffix.back().fPrefix = tmpSuffix[e].fPrefix;
                }
            }
        }
        fStore.swap(newStore);
        fCandStore = newCandStore;
        items.swap(suffix);
    }
}

bool CFDDiscovery::isConstRule(const PartitionTidList& items, int rhsA) {
    int rhsVal;
    bool first = true;
    for (int pi = 0; pi <= items.fTids.size(); pi++) {
        if (pi == items.fTids.size() || items.fTids[pi] == PartitionTidList::SEP) {
            auto tup = fDb.getRow(items.fTids[pi-1]);
            if (first) {
                first = false;
                rhsVal = tup[rhsA];
            }
            else {
                if (tup[rhsA] != rhsVal) return false;
            }
        }
    }
    return true;
}

void CFDDiscovery::integratedDFS(const Itemset &prefix, std::vector<MinerNode<PartitionTidList> > &items) {
    for (int ix = items.size()-1; ix >= 0; ix--) {
        MinerNode<PartitionTidList>& inode = items[ix];
        const Itemset iset = join(prefix, inode.fItem);
        std::vector<std::vector<std::pair<Itemset,int> > > next(items.size());
        std::vector<const PartitionTidList*> expands;
        std::vector<MinerNode<PartitionTidList> > tmpSuffix;
        std::vector<MinerNode<PartitionTidList> > suffix;
        auto insect = intersection(iset, inode.fCands);
        for (int out : insect) {
            Itemset sub = subset(iset, out);
            if (out < 0) {
                if (sub.size() && !has(sub, [](int si) -> bool { return si < 0; })) continue;
                if (inode.fTids.fNrSets == 1 || isConstRule(inode.fTids, -1-out)) continue;
            }
            else {
                if (!sub.size() || has(sub, [](int si) -> bool { return si < 0; })) continue;
            }
            auto storedSub = fStore.find(sub);
            if (storedSub == fStore.end()) {
                continue;
            }
            double e = out < 0 ? PartitionTable::partitionError(storedSub->second, inode.fTids)
                               : setdiff(storedSub->second.fTids, inode.fTids.fTids).size();
            double conf = 1 - (e / support(storedSub->second));
            Itemset cSub = sub;
            bool lhsGen = true;
            auto rsPtr = fRules.find(out);
            if (rsPtr != fRules.end()) {
                for (const auto& subRule : rsPtr->second) {
                    if (out < 0 && !has(subRule, [](int si) -> bool { return si < 0; })) continue;
                    if (precedes(subRule, cSub)) {
                        lhsGen = false;
                        break;
                    }
                }
            }
            if (!lhsGen) continue;


            if (conf >= fMinConf) {
                fCFDs.emplace_back(cSub, out);
            }
            if (conf >= 1) {
                fRules[out].push_back(cSub);
                if (out > 0) fRules[-1 - fDb.getAttrIndex(out)].push_back(cSub);
                inode.fCands = intersection(inode.fCands, cSub);
                pruneCands(items, sub, out);
            }
        }

        if (iset.size() == fMaxSize) continue;
        fStore[iset] = inode.fTids;
        fCandStore.insert(iset, inode.fCands);

        auto nodeAttrs = fDb.getAttrVector(iset);
        for (int jx = items.size()-1; jx > ix; jx--) {
            const auto& jnode = items[jx];
            if (std::binary_search(nodeAttrs.begin(), nodeAttrs.end(), fDb.getAttrIndex(jnode.fItem))) continue;
            Itemset newset = join(iset, jnode.fItem);
            auto c = intersection(inode.fCands, jnode.fCands);
            for (int zz : newset) {
                auto zsub = subset(newset, zz);
                auto storedSub = fStore.find(zsub);
                if (storedSub == fStore.end()) {
                    c.clear();
                    break;
                }
                auto subCandsPtr = fCandStore.find(zsub);
                if (subCandsPtr) {
                    const Itemset& subCands = *subCandsPtr;
                    if (subCands.size()) c = intersection(c, subCands);
                    else c.clear();
                }
                else c.clear();
                if (c.empty()) break;
            }
            if (c.size()) {
                tmpSuffix.emplace_back(jnode.fItem);
                tmpSuffix.back().fCands = c;
                expands.push_back(&jnode.fTids);
            }
        }

        const auto exps = PartitionTable::intersection(inode.fTids, expands);
        for (int e = 0; e < exps.size(); e++) {
            if (support(exps[e]) >= fMinSup) {
                bool gen = true;
                auto c = tmpSuffix[e].fCands;
                auto newset = join(iset, tmpSuffix[e].fItem);
                auto sp = std::make_pair(support(exps[e]), exps[e].fNrSets);
                auto fm2 = fFreeMap.find(sp);
                if (fm2 != fFreeMap.end()) {
                    auto freeCands = fm2->second;
                    for (const auto &subCand : freeCands) {
                        if (isSubsetOf(subCand, newset)) {
                            gen = false;
                            break;
                        }
                    }
                }
                if (gen) {
                    fFreeMap[sp].push_back(newset);
                    fFreeItemsets.insert(newset);
                }
                else {
                    for (const auto &a : newset) {
                        auto x = subset(newset, a);
                        if (fFreeItemsets.find(x) == fFreeItemsets.end()) {
                            c = intersection(c, x);
                        }
                        if (!c.size()) break;
                    }
                }
                if (c.size()) {
                    suffix.emplace_back(tmpSuffix[e].fItem, exps[e]);
                    suffix.back().fCands = c;
                }
            }
        }
        // Sort suffix and recurse
        if (suffix.size()) {
            std::sort(suffix.begin(), suffix.end(), [this](const MinerNode<PartitionTidList>& a, const MinerNode<PartitionTidList>& b) {
                int aAtt = fDb.getAttrIndex(a.fItem);
                int bAtt = fDb.getAttrIndex(b.fItem);
                return fAttOrder[aAtt] < fAttOrder[bAtt] || (aAtt == bAtt && ((a.fItem > 0 && b.fItem < 0) || (a.fItem > 0 && b.fItem > 0 && support(a.fTids) > support(b.fTids))));
            });
            integratedDFS(iset, suffix);
        }
    }
}

void CFDDiscovery::integratedDFS(int minsup, int maxSize, double minconf) {
    fMinSup = minsup;
    fMinConf = minconf;
    fMaxSize = maxSize;
    PartitionTable::fDbSize = fDb.size();
    std::vector<MinerNode<PartitionTidList> > items = getAllSingletons(fMinSup);
    fStore[Itemset()] = convert(iota(fDb.size()));
    fAttOrder = Itemset(fDb.nrAttrs());
    auto attrs = iota(fDb.nrAttrs());
    std::sort(attrs.begin(), attrs.end(), [this](int a, int b) {
        return fDb.getDomain(a).size() > fDb.getDomain(b).size();
    });
    for (int ai = 0; ai < attrs.size(); ai++) {
        fAttOrder[attrs[ai]] = ai;
    }
    for (const auto& cp : items) {
        fAllAttrs.push_back(cp.fItem);
    }
    std::sort(fAllAttrs.begin(), fAllAttrs.end());
    fCandStore.insert(Itemset(), fAllAttrs);
    for (auto& a : items) {
        Itemset at;
        for (int cat : fAllAttrs) {
            if (a.fItem < 0 && cat >= 0) continue;
            if (a.fItem != cat && fDb.getAttrIndex(a.fItem) == fDb.getAttrIndex(cat)) continue;
            at.push_back(cat);
        }
        a.fCands = at;
        fFreeMap[std::make_pair(support(a.fTids),a.fTids.fNrSets)].push_back(itemset(a.fItem));
        fFreeItemsets.insert(itemset(a.fItem));
    }
    std::sort(items.begin(), items.end(), [this](const MinerNode<PartitionTidList>& a, const MinerNode<PartitionTidList>& b) {
        int aAtt = fDb.getAttrIndex(a.fItem);
        int bAtt = fDb.getAttrIndex(b.fItem);
        return fAttOrder[aAtt] > fAttOrder[bAtt] || (aAtt == bAtt && ((a.fItem > 0 && b.fItem < 0) || (a.fItem > 0 && b.fItem > 0 && support(a.fTids) > support(b.fTids))));
    });
    integratedDFS(Itemset(), items);
}

void CFDDiscovery::itemsetsFirstDFS(int minsup, int maxSize, SUBSTRATEGY ss, double minconf) {
    fMinSup = minsup;
    fMinConf = minconf;
    fMaxSize = maxSize;
    std::vector<MinerNode<SimpleTidList> > items = getSingletons(fMinSup);
    std::set<int> atts;
    for (auto& a : items) {
        atts.insert(-1-fDb.getAttrIndex(a.fItem));
    }
    Itemset attsV(atts.begin(), atts.end());
    fAllAttrs = setdiff(range(-1-fDb.nrAttrs(),0), attsV);
    for (auto& a : items) {
        fAllAttrs.push_back(a.fItem);
    }
    std::unordered_map<Itemset, SimpleTidList> tids(items.size());
    tids[Itemset()] = iota(fDb.size());
    for (auto& a : items) {
        a.fCands = where(fAllAttrs, [this,a](int i){ return (i == a.fItem) || fDb.getAttrIndex(a.fItem) != fDb.getAttrIndex(i); });
        fFreeMap[std::make_pair(support(a.fTids),1)].push_back(itemset(a.fItem));
        fFreeItemsets.insert(itemset(a.fItem));
    }
    if (ss == SUBSTRATEGY::DFS)
        mineFDsDFS(iota(fDb.size()), Itemset(), range(-fDb.nrAttrs(), 0));
    else if (ss == SUBSTRATEGY::BFS)
        mineFDs(iota(fDb.size()), Itemset(), range(-fDb.nrAttrs(), 0));
    itemsetsFirstDFS(Itemset(), items, ss);
}

void CFDDiscovery::itemsetsFirstDFS(const Itemset &prefix, std::vector<MinerNode<SimpleTidList> > &items,
                                    SUBSTRATEGY ss) {
    for (int ix = items.size()-1; ix >= 0; ix--) {
        MinerNode<SimpleTidList>& inode = items[ix];
        const Itemset iset = join(prefix, inode.fItem);
        std::vector<MinerNode<SimpleTidList> > suffix;
        if (ix < items.size() - 1) suffix.reserve(items.size()-ix-1);
        auto insect = intersection(iset, inode.fCands);
        if (iset.size() > 1) {
            for (int out : insect) {
                Itemset sub = subset(iset, out);
                auto storedSub = fTidStore.find(sub);
                if (storedSub == fTidStore.end()) {
                    continue;
                }
                if (support(storedSub->second) < fMinSup) continue;
                bool lhsGen = true;
                if (fRules.find(out) != fRules.end()) {
                    for (const auto &subRule : fRules[out]) {
                        if (precedes(subRule, sub)) {
                            lhsGen = false;
                            inode.fCands = intersection(inode.fCands, subRule);
                        }
                    }
                }
                if (fFreeItemsets.find(sub) == fFreeItemsets.end()) {
                    lhsGen = false;
                }
                if (!lhsGen) continue;
                double e = setdiff(storedSub->second, inode.fTids).size();
                double conf = 1 - (e / support(storedSub->second));

                if (conf >= fMinConf) {
                    fCFDs.emplace_back(sub, out);
                }
                if (conf >= 1) {
                    fRules[out].push_back(sub);
                    fRules[-1 - fDb.getAttrIndex(out)].push_back(sub);
                    inode.fCands = intersection(inode.fCands, sub);
                    pruneCands(items, sub, out);
                }
            }
        }
        if (iset.size() == fMaxSize) continue;

        if (fGens.addMinGen(iset, inode.fSupp, inode.fHash) && iset.size() < fMaxSize-1) {
            auto nodeAttrs = fDb.getAttrVectorItems(iset);
            auto cAs = fDb.getAttrVectorItems(inode.fCands);
            Itemset allAttrs = intersection(cAs, setdiff(range(-fDb.nrAttrs(), 0), nodeAttrs));
            if (inode.fSupp >= fMinSup && allAttrs.size()) {
                if (ss == SUBSTRATEGY::DFS)
                    mineFDsDFS(inode.fTids, iset, allAttrs);
                else if (ss == SUBSTRATEGY::BFS)
                    mineFDs(inode.fTids, iset, allAttrs);
            }
            fTidStore[iset] = inode.fTids;
        }
        fCandStore.insert(iset, inode.fCands);

        auto nodeAttrs = fDb.getAttrVector(iset);
        for (int jx = items.size()-1; jx > ix; jx--) {
            const auto &jnode = items[jx];
            if (jnode.fPrefix != inode.fPrefix) continue;
            if (std::binary_search(nodeAttrs.begin(), nodeAttrs.end(), -1-fDb.getAttrIndex(jnode.fItem))) continue;
            auto c = intersection(inode.fCands, jnode.fCands);
            SimpleTidList ijtids = intersection(inode.fTids, jnode.fTids);
            int ijsupp = ijtids.size();
            if (ijsupp >= fMinSup * fMinConf) {
                Itemset newset = join(iset, jnode.fItem);
                for (int zz : newset) {
                    if (!c.size()) break;
                    auto zsub = subset(newset, zz);
                    auto subCandsPtr = fCandStore.find(zsub);
                    if (subCandsPtr)
                        c = intersection(c, *subCandsPtr);
                    else {
                        c.clear();
                    }
                }
                if (c.size()) {
                    bool gen = true;
                    auto sp = std::make_pair(ijsupp, 1);
                    auto fm2 = fFreeMap.find(sp);
                    if (fm2 != fFreeMap.end()) {
                        auto freeCands = fm2->second;
                        for (const auto &subCand : freeCands) {
                            if (isSubsetOf(subCand, newset)) {
                                gen = false;
                                break;
                            }
                        }
                    }
                    if (gen) {
                        fFreeMap[sp].push_back(newset);
                        fFreeItemsets.insert(newset);
                    }

                    suffix.emplace_back(jnode.fItem, std::move(ijtids), ijsupp);
                    suffix.back().fCands = c;
                }
            }
        }
        // Sort suffix and recurse
        if (suffix.size()) {
            std::sort(suffix.begin(), suffix.end(), [this](const MinerNode<SimpleTidList>& a, const MinerNode<SimpleTidList>& b) {
                return support(a.fTids) < support(b.fTids);
            });
            itemsetsFirstDFS(iset, suffix, ss);
        }
    }
}

void CFDDiscovery::mineFDsDFS(const SimpleTidList &tids, const Itemset &tp, const Itemset &allAttrs) {
    PartitionTable::fDbSize = fDb.size();
    Itemset allAtt;
    for (int a : allAttrs) {
        allAtt.push_back(-1 - a);
    }
    std::vector<MinerNode<PartitionTidList> > items = getPartitionSingletons(tids, allAtt);
    for (auto &a : items) {
        a.fCands = allAttrs;
        if (a.fTids.fNrSets > 1) {
            fFreeMap[std::make_pair(support(a.fTids), a.fTids.fNrSets)].push_back(join(tp, a.fItem));
            fFreeItemsets.insert(join(tp, a.fItem));
        }
    }
    PrefixTree<Itemset, Itemset> candStore;
    fStore.clear();
    fStore[Itemset()] = convert(tids);
    candStore.insert(Itemset(), allAttrs);
    mineFDsDFS(items, Itemset(), tp);
}

void CFDDiscovery::mineFDsDFS(std::vector<MinerNode<PartitionTidList> > &items, const Itemset &prefix,
                              const Itemset &tp) {
    for (int ix = items.size() - 1; ix >= 0; ix--) {
        MinerNode<PartitionTidList> &inode = items[ix];
        const Itemset iset = join(prefix, inode.fItem);
        auto insect = intersection(iset, inode.fCands);
        if (iset.size() > 1) {
            for (int out : insect) {
                Itemset sub = subset(iset, out);
                if (!std::binary_search(inode.fCands.begin(), inode.fCands.end(), out)) continue;
                if (inode.fTids.fNrSets == 1 || isConstRule(inode.fTids, -1 - out)) continue;
                auto storedSub = fStore.find(sub);
                if (storedSub == fStore.end()) {
                    continue;
                }
                Itemset cSub = join(sub, tp);

                bool lhsGen = true;
                if (fFreeItemsets.find(cSub) == fFreeItemsets.end()) {
                    lhsGen = false;
                }
                if (fRules.find(out) != fRules.end()) {
                    for (const auto &subRule : fRules[out]) {
                        if (precedes(subRule, cSub)) {
                            auto subAtts = fDb.getAttrVector(subRule);
                            inode.fCands = intersection(inode.fCands, subRule);
                            lhsGen = false;
                        }
                    }
                }

                if (!lhsGen) continue;
                double e = PartitionTable::partitionError(storedSub->second, inode.fTids);
                double conf = 1 - (e / support(storedSub->second));
                if (conf >= fMinConf) {
                    fCFDs.emplace_back(cSub, out);
                }
                if (conf >= 1) {
                    fRules[out].push_back(cSub);
                    inode.fCands = intersection(inode.fCands, cSub);
                    pruneCands(items, sub, out);
                }
            }
        }
        if (inode.fCands.empty()) continue;
        if (iset.size() + tp.size() >= fMaxSize) continue;
        fStore[iset] = inode.fTids;
        fCandStore.insert(iset, inode.fCands);

        std::vector<MinerNode<PartitionTidList> > suffix;
        std::vector<const PartitionTidList*> expands;
        std::vector<MinerNode<PartitionTidList> > tmpSuffix;
        for (int jx = items.size()-1; jx > ix; jx--) {
            const auto &jnode = items[jx];
            auto nodeAttrs = fDb.getAttrVector(iset);
            if (jnode.fPrefix != inode.fPrefix) continue;
            if (std::binary_search(nodeAttrs.begin(), nodeAttrs.end(), fDb.getAttrIndex(jnode.fItem))) continue;
            Itemset newset = join(iset, items[jx].fItem);
            auto c = inode.fCands;
            for (int zz : newset) {
                if (!c.size()) break;
                auto subCandsPtr = fCandStore.find(subset(newset, zz));
                if (subCandsPtr)
                    c = intersection(c, *subCandsPtr);
                else
                    c.clear();
            }
            if (c.size()) {
                expands.push_back(&items[jx].fTids);
                tmpSuffix.emplace_back(items[jx].fItem);
                tmpSuffix.back().fCands = c;
                tmpSuffix.back().fPrefix = subset(newset, tmpSuffix.back().fItem);
            }
        }
        const auto exps = PartitionTable::intersection(items[ix].fTids, expands);
        for (int e = 0; e < exps.size(); e++) {
            bool gen = true;
            auto newset = join(tmpSuffix[e].fPrefix, tmpSuffix[e].fItem);
            newset = join(tp, newset);
            auto sp = std::make_pair(support(exps[e]), exps[e].fNrSets);
            auto fm2 = fFreeMap.find(sp);
            if (fm2 != fFreeMap.end()) {
                auto freeCands = fm2->second;
                for (const auto &subCand : freeCands) {
                    if (isSubsetOf(subCand, newset)) {
                        gen = false;
                        break;
                    }
                }
            }
            if (gen) {
                fFreeMap[sp].push_back(newset);
                fFreeItemsets.insert(newset);
            }
            suffix.emplace_back(tmpSuffix[e].fItem, exps[e]);
            suffix.back().fCands = tmpSuffix[e].fCands;
            suffix.back().fPrefix = tmpSuffix[e].fPrefix;
        }
        if (suffix.size()) {
            std::sort(suffix.begin(), suffix.end(), [this](const MinerNode<PartitionTidList>& a, const MinerNode<PartitionTidList>& b) {
                return a.fTids.fNrSets < b.fTids.fNrSets;
            });
            mineFDsDFS(suffix, iset, tp);
        }
    }
}

void CFDDiscovery::itemsetsFirstBFS(int minsup, int maxSize, SUBSTRATEGY ss, double minconf) {
    fMinSup = minsup;
    fMinConf = minconf;
    fMaxSize = maxSize;
    std::vector<MinerNode<SimpleTidList> > items = getSingletons(fMinSup);
    std::set<int> atts;
    for (auto& a : items) {
        atts.insert(-1-fDb.getAttrIndex(a.fItem));
    }
    Itemset attsV(atts.begin(), atts.end());
    fAllAttrs = setdiff(range(-1-fDb.nrAttrs(),0), attsV);
    for (auto& a : items) {
        fAllAttrs.push_back(a.fItem);
    }
    std::unordered_map<Itemset, SimpleTidList> tids(items.size());
    tids[Itemset()] = iota(fDb.size());
    for (auto& a : items) {
        a.fCands = where(fAllAttrs, [this,a](int i){ return (i == a.fItem) || fDb.getAttrIndex(a.fItem) != fDb.getAttrIndex(i); });
        fFreeMap[std::make_pair(support(a.fTids),1)].push_back(itemset(a.fItem));
        fFreeItemsets.insert(itemset(a.fItem));
    }
    if (ss == SUBSTRATEGY::DFS)
        mineFDsDFS(iota(fDb.size()), Itemset(), range(-fDb.nrAttrs(), 0));
    else if (ss == SUBSTRATEGY::BFS)
        mineFDs(iota(fDb.size()), Itemset(), range(-fDb.nrAttrs(), 0));
    while (items.size()) {
        std::vector<MinerNode<SimpleTidList> > suffix;
        std::unordered_map<Itemset, SimpleTidList> newTids(items.size());
        PrefixTree<Itemset, Itemset> candStore;
        for (int i = 0; i < items.size(); i++) {
            auto &inode = items[i];
            Itemset iset = inode.fPrefix;
            iset.push_back(inode.fItem);
            auto insect = intersection(iset, inode.fCands);
            if (iset.size() > 1) {
                for (int out : insect) {
                    Itemset sub = subset(iset, out);
                    auto storedSub = tids.find(sub);
                    if (storedSub == tids.end()) {
                        continue;
                    }
                    if (support(storedSub->second) < fMinSup) continue;
                    Itemset cSub = sub;
                    bool lhsGen = true;
                    if (fRules.find(out) != fRules.end()) {
                        for (const auto &subRule : fRules[out]) {
                            if (precedes(subRule, cSub)) {
                                lhsGen = false;
                                inode.fCands = intersection(inode.fCands, subRule);
                                //break;
                            }
                        }
                    }
                    if (fFreeItemsets.find(cSub) == fFreeItemsets.end()) {
                        lhsGen = false;
                    }
                    if (!lhsGen) continue;
                    double e = setdiff(storedSub->second, inode.fTids).size();
                    double conf = 1 - (e / support(storedSub->second));

                    if (conf >= fMinConf) {
                        fCFDs.emplace_back(cSub, out);
                    }
                    if (conf >= 1) {
                        fRules[out].push_back(sub);
                        fRules[-1 - fDb.getAttrIndex(out)].push_back(cSub);
                        inode.fCands = intersection(inode.fCands, cSub);
                        pruneCands(items, sub, out);
                    }
                }
            }
            if (fGens.addMinGen(iset, inode.fSupp, inode.fHash)) {
                auto nodeAttrs = fDb.getAttrVectorItems(iset);
                auto cAs = fDb.getAttrVectorItems(inode.fCands);
                Itemset allAttrs = intersection(cAs, setdiff(range(-fDb.nrAttrs(), 0), nodeAttrs));
                if (inode.fSupp >= fMinSup && allAttrs.size()) {
                    if (ss == SUBSTRATEGY::DFS)
                        mineFDsDFS(inode.fTids, iset, allAttrs);
                    else if (ss == SUBSTRATEGY::BFS)
                        mineFDs(inode.fTids, iset, allAttrs);
                }
                newTids[iset] = inode.fTids;
            }
        }
        for (int i = 0; i < items.size(); i++) {
            MinerNode<SimpleTidList> &inode = items[i];
            const Itemset iset = join(inode.fPrefix, inode.fItem);
            candStore.insert(iset, inode.fCands);
        }
        for (int i = 0; i < items.size(); i++) {
            auto &inode = items[i];
            Itemset iset = inode.fPrefix;
            iset.push_back(inode.fItem);
            auto nodeAttrs = fDb.getAttrVectorItems(iset);
            if (iset.size() == fMaxSize) continue;
            for (int j = i + 1; j < items.size(); j++) {
                const auto &jnode = items[j];
                if (jnode.fPrefix != inode.fPrefix) continue;
                if (std::binary_search(nodeAttrs.begin(), nodeAttrs.end(), -1-fDb.getAttrIndex(jnode.fItem))) continue;
                auto c = intersection(inode.fCands, jnode.fCands);
                Itemset jset = jnode.fPrefix;
                jset.push_back(jnode.fItem);
                Itemset newset = join(jset, inode.fItem);
                for (int zz : newset) {
                    if (!c.size()) break;
                    auto zsub = subset(newset, zz);
                    auto subCandsPtr = candStore.find(zsub);
                    if (subCandsPtr)
                        c = intersection(c, *subCandsPtr);
                    else
                        c.clear();
                }
                if (c.size()) {
                    SimpleTidList ijtids = intersection(inode.fTids, jnode.fTids);
                    int ijsupp = ijtids.size();
                    if (ijsupp >= fMinSup * fMinConf) {
                        bool gen = true;
                        auto sp = std::make_pair(ijsupp, 1);
                        auto fm2 = fFreeMap.find(sp);
                        if (fm2 != fFreeMap.end()) {
                            auto freeCands = fm2->second;
                            for (const auto &subCand : freeCands) {
                                if (isSubsetOf(subCand, newset)) {
                                    gen = false;
                                    break;
                                }
                            }
                        }
                        if (gen) {
                            fFreeMap[sp].push_back(newset);
                            fFreeItemsets.insert(newset);
                        }

                        int jtem = newset.back();
                        newset.pop_back();
                        suffix.emplace_back(jtem, std::move(ijtids), ijsupp, newset);
                        suffix.back().fCands = c;
                    }
                }
            }
        }
        tids.swap(newTids);
        items.swap(suffix);
    }
}

void CFDDiscovery::mineFDs(const SimpleTidList& tids, const Itemset& tp, const Itemset& allAttrs) {
    PartitionTable::fDbSize = fDb.size();
    Itemset allAtt;
    for (int a : allAttrs) {
        allAtt.push_back(-1-a);
    }
    std::vector<MinerNode<PartitionTidList> > items = getPartitionSingletons(tids, allAtt);
    for (auto& a : items) {
        a.fCands = allAttrs;
        if (a.fTids.fNrSets > 1) {
            fFreeMap[std::make_pair(support(a.fTids),a.fTids.fNrSets)].push_back(join(tp, a.fItem));
            fFreeItemsets.insert(join(tp, a.fItem));
        }
    }
    PrefixTree<Itemset, Itemset> candStore;
    fStore.clear();
    fStore[Itemset()] = convert(tids);
    candStore.insert(Itemset(), allAttrs);


    while (!items.empty()) {
        std::map<Itemset,PartitionTidList> newStore;
        PrefixTree<Itemset, Itemset> newCandStore;
        std::vector<std::vector<std::pair<Itemset,int> > > next(items.size());
        for (int i = 0; i < items.size(); i++) {
            MinerNode<PartitionTidList>& inode = items[i];
            const Itemset iset = join(inode.fPrefix, inode.fItem);
            auto insect = intersection(iset, inode.fCands);
            if (iset.size() > 1) {
                for (int out : insect) {
                    Itemset sub = subset(iset, out);
                    if (!std::binary_search(inode.fCands.begin(), inode.fCands.end(), out)) continue;
                    if (inode.fTids.fNrSets == 1 || isConstRule(inode.fTids, -1-out)) continue;
                    auto storedSub = fStore.find(sub);
                    if (storedSub == fStore.end()) {
                        continue;
                    }
                    Itemset cSub = join(sub, tp);

                    bool lhsGen = true;
                    if (fFreeItemsets.find(cSub) == fFreeItemsets.end()) {
                        lhsGen = false;
                    }
                    if (fRules.find(out) != fRules.end()) {
                        for (const auto &subRule : fRules[out]) {
                            if (precedes(subRule, cSub)) {
                                auto subAtts = fDb.getAttrVector(subRule);
                                inode.fCands = intersection(inode.fCands, subRule);
                                lhsGen = false;
                            }
                        }
                    }

                    if (!lhsGen) continue;
                    double e = PartitionTable::partitionError(storedSub->second, inode.fTids);
                    double conf = 1 - (e / support(storedSub->second));
                    if (conf >= fMinConf) {
                        fCFDs.emplace_back(cSub, out);
                    }
                    if (conf >= 1) {
                        fRules[out].push_back(cSub);
                        inode.fCands = intersection(inode.fCands, cSub);
                        pruneCands(items, sub, out);
                    }
                }
            }
            if (inode.fCands.empty()) continue;
            if (iset.size() + tp.size() >= fMaxSize) continue;
            newStore[iset] = inode.fTids;
            auto nodeAttrs = fDb.getAttrVector(iset);
            for (int j = i+1; j < items.size(); j++) {
                if (j == i) continue;
                const auto& jnode = items[j];
                if (jnode.fPrefix != inode.fPrefix) continue;
                if (std::binary_search(nodeAttrs.begin(), nodeAttrs.end(), fDb.getAttrIndex(jnode.fItem))) continue;
                next[i].emplace_back(iset, j);
            }
        }
        for (int i = 0; i < items.size(); i++) {
            MinerNode<PartitionTidList> &inode = items[i];
            const Itemset iset = join(inode.fPrefix, inode.fItem);
            newCandStore.insert(iset, inode.fCands);
        }

        std::vector<MinerNode<PartitionTidList> > suffix;
        for (int i = 0; i < items.size(); i++) {
            std::vector<PartitionTidList*> expands;
            std::vector<MinerNode<PartitionTidList> > tmpSuffix;
            for (auto& newsetTup : next[i]) {
                int j = newsetTup.second;
                Itemset newset = join(newsetTup.first, items[j].fItem);
                auto c = items[i].fCands;
                for (int zz : newset) {
                    if (!c.size()) break;
                    auto subCandsPtr = newCandStore.find(subset(newset, zz));
                    if (subCandsPtr)
                        c = intersection(c, *subCandsPtr);
                    else
                        c.clear();
                }
                if (c.size()) {
                    expands.push_back(&items[j].fTids);
                    tmpSuffix.emplace_back(items[j].fItem);
                    tmpSuffix.back().fCands = c;
                    tmpSuffix.back().fPrefix = newsetTup.first;
                }
            }
            const auto exps = PartitionTable::intersection(items[i].fTids, expands);
            for (int e = 0; e < exps.size(); e++) {
                bool gen = true;
                auto newset = join(tmpSuffix[e].fPrefix, tmpSuffix[e].fItem);
                newset = join(tp, newset);
                auto sp = std::make_pair(support(exps[e]), exps[e].fNrSets);
                auto fm2 = fFreeMap.find(sp);
                if (fm2 != fFreeMap.end()) {
                    auto freeCands = fm2->second;
                    for (const auto &subCand : freeCands) {
                        if (isSubsetOf(subCand, newset)) {
                            gen = false;
                            break;
                        }
                    }
                }
                if (gen) {
                    fFreeMap[sp].push_back(newset);
                    fFreeItemsets.insert(newset);
                }
                suffix.emplace_back(tmpSuffix[e].fItem, exps[e]);
                suffix.back().fCands = tmpSuffix[e].fCands;
                suffix.back().fPrefix = tmpSuffix[e].fPrefix;
            }
        }
        fStore.swap(newStore);
        items.swap(suffix);
    }
}

void CFDDiscovery::fdsFirstDFS(int minsup, int maxSize, SUBSTRATEGY ss, double minconf) {
    fMinSup = minsup;
    fMinConf = minconf;
    fMaxSize = maxSize;
    fAllAttrs = range(-fDb.nrAttrs(), 0);
    PartitionTable::fDbSize = fDb.size();
    std::vector<MinerNode<PartitionTidList> > items = getPartitionSingletons();
    for (auto& a : items) {
        a.fCands = fAllAttrs;
        fFreeMap[std::make_pair(support(a.fTids),a.fTids.fNrSets)].push_back(itemset(a.fItem));
        fFreeItemsets.insert(itemset(a.fItem));
    }
    fCandStore = PrefixTree<Itemset, Itemset>();
    fStore[Itemset()] = convert(iota(fDb.size()));
    fCandStore.insert(Itemset(), fAllAttrs);
    fdsFirstDFS(Itemset(), items, ss);
}

void CFDDiscovery::fdsFirstDFS(const Itemset &prefix, std::vector<MinerNode<PartitionTidList> > &items, SUBSTRATEGY ss) {
    for (int ix = items.size()-1; ix >= 0; ix--) {
        MinerNode<PartitionTidList>& inode = items[ix];
        const Itemset iset = join(prefix, inode.fItem);
        auto insect = intersection(iset, inode.fCands);
        for (int out : insect) {
            Itemset sub = subset(iset, out);
            if (inode.fTids.fNrSets == 1 || isConstRule(inode.fTids, -1-out)) continue;
            auto storedSub = fStore.find(sub);
            if (storedSub == fStore.end()) {
                continue;
            }
            bool lhsGen = true;
            if (fFreeItemsets.find(sub) == fFreeItemsets.end()) {
                lhsGen = false;
            }
            if (fRules.find(out) != fRules.end()) {
                for (const auto& subRule : fRules[out]) {
                    if (!has(subRule, [](int si) -> bool { return si < 0; })) continue;
                    if (precedes(subRule, sub)) {
                        lhsGen = false;
                    }
                }
            }
            if (lhsGen) {
                double e = PartitionTable::partitionError(storedSub->second, inode.fTids);
                double conf = 1 - (e / support(storedSub->second));
                if (conf >= fMinConf) {
                    fCFDs.emplace_back(sub, out);
                }
                if (conf >= 1) {
                    fRules[out].push_back(sub);
                }
            }
            if (ss == SUBSTRATEGY::DFS)
                minePatternsDFS(sub, out, inode.fTids);
            else if (ss == SUBSTRATEGY::BFS)
                minePatternsBFS(sub, out, inode.fTids);
        }
        if (inode.fCands.empty()) continue;
        if (iset.size() == fMaxSize) continue;
        fStore[iset] = inode.fTids;
        fCandStore.insert(iset, inode.fCands);

        auto nodeAttrs = fDb.getAttrVector(iset);
        std::vector<const PartitionTidList*> expands;
        std::vector<MinerNode<PartitionTidList> > tmpSuffix;
        for (int jx = items.size()-1; jx > ix; jx--) {
            const auto &jnode = items[jx];
            if (std::binary_search(nodeAttrs.begin(), nodeAttrs.end(), fDb.getAttrIndex(jnode.fItem))) continue;
            Itemset newset = join(iset, jnode.fItem);
            auto c = intersection(inode.fCands, jnode.fCands);
            for (int zz : newset) {
                auto zsub = subset(newset, zz);
                auto storedSub = fStore.find(zsub);
                if (storedSub == fStore.end()) {
                    c.clear();
                    break;
                }
                const Itemset &subCands = *fCandStore.find(zsub);
                c = intersection(c, subCands);
            }
            if (c.size()) {
                expands.push_back(&items[jx].fTids);
                tmpSuffix.emplace_back(items[jx].fItem);
                tmpSuffix.back().fCands = c;
                tmpSuffix.back().fPrefix = subset(newset, tmpSuffix.back().fItem);
            }
        }
        const auto exps = PartitionTable::intersection(items[ix].fTids, expands);
        std::vector<MinerNode<PartitionTidList> > suffix;

        for (int e = 0; e < exps.size(); e++) {
            bool gen = true;
            auto newset = join(tmpSuffix[e].fPrefix, tmpSuffix[e].fItem);
            auto sp = std::make_pair(support(exps[e]), exps[e].fNrSets);
            auto fm2 = fFreeMap.find(sp);
            if (fm2 != fFreeMap.end()) {
                auto freeCands = fm2->second;
                for (const auto &subCand : freeCands) {
                    if (isSubsetOf(subCand, newset)) {
                        gen = false;
                        break;
                    }
                }
            }
            if (gen) {
                fFreeMap[sp].push_back(newset);
                fFreeItemsets.insert(newset);
            }
            suffix.emplace_back(tmpSuffix[e].fItem, exps[e]);
            suffix.back().fCands = tmpSuffix[e].fCands;
            suffix.back().fPrefix = tmpSuffix[e].fPrefix;
        }
        if (suffix.size()) {
            std::sort(suffix.begin(), suffix.end(), [this](const MinerNode<PartitionTidList>& a, const MinerNode<PartitionTidList>& b) {
                return a.fTids.fNrSets < b.fTids.fNrSets;
            });
            fdsFirstDFS(iset, suffix, ss);
        }
    }
}

void CFDDiscovery::fdsFirstBFS(int minsup, int maxSize, SUBSTRATEGY ss, double minconf) {
    fMinSup = minsup;
    fMinConf = minconf;
    fMaxSize = maxSize;
    fAllAttrs = range(-fDb.nrAttrs(), 0);
    PartitionTable::fDbSize = fDb.size();
    std::vector<MinerNode<PartitionTidList> > items = getPartitionSingletons();
    for (auto& a : items) {
        a.fCands = fAllAttrs;
        fFreeMap[std::make_pair(support(a.fTids),a.fTids.fNrSets)].push_back(itemset(a.fItem));
        fFreeItemsets.insert(itemset(a.fItem));
    }
    fCandStore = PrefixTree<Itemset, Itemset>();
    fStore[Itemset()] = convert(iota(fDb.size()));
    fCandStore.insert(Itemset(), fAllAttrs);

    while (!items.empty()) {
        std::map<Itemset,PartitionTidList> newStore;
        PrefixTree<Itemset, Itemset> newCandStore;
        std::vector<std::vector<std::pair<Itemset,int> > > next(items.size());
        for (int i = 0; i < items.size(); i++) {
            MinerNode<PartitionTidList>& inode = items[i];
            const Itemset iset = join(inode.fPrefix, inode.fItem);
            auto insect = intersection(iset, inode.fCands);
            for (int out : insect) {
                Itemset sub = subset(iset, out);
                if (inode.fTids.fNrSets == 1 || isConstRule(inode.fTids, -1-out)) continue;
                auto storedSub = fStore.find(sub);
                if (storedSub == fStore.end()) {
                    continue;
                }
                bool lhsGen = true;
                if (fFreeItemsets.find(sub) == fFreeItemsets.end()) {
                    lhsGen = false;
                }
                if (fRules.find(out) != fRules.end()) {
                    for (const auto& subRule : fRules[out]) {
                        if (!has(subRule, [](int si) -> bool { return si < 0; })) continue;
                        if (precedes(subRule, sub)) {
                            lhsGen = false;
                        }
                    }
                }
                if (lhsGen) {
                    double e = PartitionTable::partitionError(storedSub->second, inode.fTids);
                    double conf = 1 - (e / support(storedSub->second));
                    if (conf >= fMinConf) {
                        fCFDs.emplace_back(sub, out);
                    }
                    if (conf >= 1) {
                        fRules[out].push_back(sub);
                    }
                }

                if (ss == SUBSTRATEGY::DFS)
                    minePatternsDFS(sub, out, inode.fTids);
                else if (ss == SUBSTRATEGY::BFS)
                    minePatternsBFS(sub, out, inode.fTids);
            }
            if (inode.fCands.empty()) continue;
            if (iset.size() == fMaxSize) continue;

            newStore[iset] = inode.fTids;
            auto nodeAttrs = fDb.getAttrVector(iset);
            for (int j = i+1; j < items.size(); j++) {
                if (j == i) continue;
                const auto& jnode = items[j];
                if (jnode.fPrefix != inode.fPrefix) continue;
                if (std::binary_search(nodeAttrs.begin(), nodeAttrs.end(), fDb.getAttrIndex(jnode.fItem))) continue;
                next[i].emplace_back(iset, j);
            }
        }
        for (int i = 0; i < items.size(); i++) {
            MinerNode<PartitionTidList> &inode = items[i];
            const Itemset iset = join(inode.fPrefix, inode.fItem);
            newCandStore.insert(iset, inode.fCands);
        }

        std::vector<MinerNode<PartitionTidList> > suffix;
        for (int i = 0; i < items.size(); i++) {
            std::vector<PartitionTidList*> expands;
            std::vector<MinerNode<PartitionTidList> > tmpSuffix;
            for (auto& newsetTup : next[i]) {
                int j = newsetTup.second;
                Itemset newset = join(newsetTup.first, items[j].fItem);
                auto c = items[i].fCands;
                for (int zz : newset) {

                    auto zsub = subset(newset, zz);
                    auto storedSub = newStore.find(zsub);
                    if (storedSub == newStore.end()) {
                        c.clear();
                        break;
                    }
                    const Itemset &subCands = *newCandStore.find(zsub);
                    c = intersection(c, subCands);
                }
                if (c.size()) {
                    expands.push_back(&items[j].fTids);
                    tmpSuffix.emplace_back(items[j].fItem);
                    tmpSuffix.back().fCands = c;
                    tmpSuffix.back().fPrefix = newsetTup.first;
                }
            }
            const auto exps = PartitionTable::intersection(items[i].fTids, expands);
            for (int e = 0; e < exps.size(); e++) {
                bool gen = true;
                auto newset = join(tmpSuffix[e].fPrefix, tmpSuffix[e].fItem);
                auto sp = std::make_pair(support(exps[e]), exps[e].fNrSets);
                auto fm2 = fFreeMap.find(sp);
                if (fm2 != fFreeMap.end()) {
                    auto freeCands = fm2->second;
                    for (const auto &subCand : freeCands) {
                        if (isSubsetOf(subCand, newset)) {
                            gen = false;
                            break;
                        }
                    }
                }
                if (gen) {
                    fFreeMap[sp].push_back(newset);
                    fFreeItemsets.insert(newset);
                }

                suffix.emplace_back(tmpSuffix[e].fItem, exps[e]);
                suffix.back().fCands = tmpSuffix[e].fCands;
                suffix.back().fPrefix = tmpSuffix[e].fPrefix;
            }
        }
        fStore.swap(newStore);
        fCandStore = newCandStore;
        items.swap(suffix);
    }
}

int CFDDiscovery::getPartitionSupport(const SimpleTidList& pids, const std::vector<int>& psupps) {
    int res = 0;
    for (int p : pids) {
        res += psupps[p];
    }
    return res;
}

int CFDDiscovery::getPartitionError(const SimpleTidList& pids, const std::vector<std::pair<Itemset, std::vector<int> > >& partitions) {
    int res = 0;
    for (int p : pids) {
        int max = 0;
        int total = 0;
        for (const auto& rhs : partitions[p].second) {
            total += rhs;
            if (rhs > max) {
                max = rhs;
            }
        }
        res += total - max;
    }
    return res;
}

bool CFDDiscovery::isConstRulePartition(const SimpleTidList &items,
                                        const std::vector<std::vector<std::pair<int, int> > > &rhses) {
    int rhsVal;
    bool first = true;
    for (int pi : items) {
        if (first) {
            first = false;
            rhsVal = rhses[pi][0].first;
        }
        for (auto rp : rhses[pi]) {
            int r = rp.first;
            if (r != rhsVal) return false;
        }
    }
    return true;
}

void CFDDiscovery::minePatternsBFS(const Itemset &lhs, int rhs, const PartitionTidList &allTids) {
    std::map<int, SimpleTidList> pidlists;
    auto lhsAttrs = fDb.getAttrVector(lhs);
    std::vector<std::pair<Itemset, std::vector<int> > > partitions;
    std::vector<std::vector<std::pair<int,int> > > rhses2;
    std::unordered_map<Itemset, int> ruleIxs;
    std::vector<int> rhses;
    int count = 0;
    for (int pi = 0; pi <= allTids.fTids.size(); pi++) {
        if (pi == allTids.fTids.size() || allTids.fTids[pi] == PartitionTidList::SEP ) {
            const Transaction& trans = fDb.getRow(allTids.fTids[pi-1]);
            Itemset lhsConstants = projection(trans, lhsAttrs);
            auto ruleI = ruleIxs.find(lhsConstants);
            if (ruleI == ruleIxs.end()) {
                rhses.push_back(trans[-1-rhs]);
                rhses2.emplace_back();
                rhses2.back().emplace_back(trans[-1-rhs],count);
                ruleIxs[lhsConstants] = partitions.size();
                partitions.emplace_back(lhsConstants, std::vector<int>());
                partitions.back().second.push_back(count);
            }
            else {
                partitions[ruleI->second].second.push_back(count);
                rhses2[ruleI->second].emplace_back(trans[-1-rhs], count);
            }
            count = 0;
        }
        else {
            count++;
        }
    }
    std::vector<MinerNode<SimpleTidList> > items;
    std::vector<int> psupps(partitions.size());
    int ri = 0;
    for (const auto& rule : partitions) {
        for (int i : rule.first) {
            pidlists[i].push_back(ri);
        }
        int s = 0;
        for (const auto& rhs : rule.second) {
            s += rhs;
        }
        psupps[ri] = s;
        ri++;
    }
    for (const auto& item : pidlists) {
        int psupp = getPartitionSupport(item.second, psupps);
        if (psupp >= fMinSup) {
            Itemset ns = join(itemset(item.first), subset(lhs, -1-fDb.getAttrIndex(item.first)));
            std::set<Itemset> nrParts;
            for (int pid : item.second) {
                nrParts.insert(partitions[pid].first);
            }
            bool gen = true;
            auto sp = std::make_pair(psupp,nrParts.size());
            auto fm2 = fFreeMap.find(sp);
            if (fm2 != fFreeMap.end()) {
                auto freeCands = fm2->second;
                for (const auto &subCand : freeCands) {
                    if (isSubsetOf(subCand, ns)) {
                        gen = false;
                        break;
                    }
                }
            }
            if (gen) {
                fFreeMap[sp].push_back(ns);
                fFreeItemsets.insert(ns);
            }
            items.emplace_back(item.first, item.second, psupp);
        }
    }
    while (items.size()) {
        std::vector<MinerNode<SimpleTidList> > suffix;
        for (int i = 0; i < items.size(); i++) {
            const auto &inode = items[i];
            Itemset iset = inode.fPrefix;
            iset.push_back(inode.fItem);
            auto nodeAttrs = fDb.getAttrVectorItems(iset);
            bool lhsGen = true;
            int out = (iset.size() == lhs.size()) ? getMaxElem(rhses2[inode.fTids[0]]) : rhs;
            auto sub = join(iset, setdiff(lhs, nodeAttrs));

            if (out > 0 || !isConstRulePartition(inode.fTids, rhses2)) {
                if (fRules.find(out) != fRules.end()) {
                    for (const auto &subRule : fRules[out]) {
                        if (out < 0 && !has(subRule, [](int si) -> bool { return si < 0; })) continue;
                        if (precedes(subRule, sub)) {
                            lhsGen = false;
                        }
                    }
                }
                if (lhsGen) {
                    int e = getPartitionError(inode.fTids, partitions);
                    double conf = 1.0 - ((double) e / (double) inode.fSupp);
                    if (conf >= fMinConf) {
                        fCFDs.emplace_back(sub, out);
                    }
                    if (conf >= 1) {
                        fRules[out].push_back(sub);
                        if (out > 0) fRules[-1 - fDb.getAttrIndex(out)].push_back(sub);
                    }
                }
            }
            for (int j = i + 1; j < items.size(); j++) {
                const auto &jnode = items[j];
                if (jnode.fPrefix != inode.fPrefix) continue;
                if (std::binary_search(nodeAttrs.begin(), nodeAttrs.end(), -1-fDb.getAttrIndex(jnode.fItem))) continue;
                Itemset jset = jnode.fPrefix;
                jset.push_back(jnode.fItem);
                Itemset newset = join(jset, inode.fItem);
                SimpleTidList ijtids = intersection(inode.fTids, jnode.fTids);
                int ijsupp = getPartitionSupport(ijtids, psupps);
                if (ijsupp >= fMinSup){
                    bool gen = true;
                    auto nas = fDb.getAttrVectorItems(newset);
                    Itemset ns = join(newset, setdiff(lhs, nas));
                    std::set<Itemset> nrParts;
                    for (int pid : ijtids) {
                        nrParts.insert(partitions[pid].first);
                    }
                    auto sp = std::make_pair(ijsupp, nrParts.size());
                    auto fm2 = fFreeMap.find(sp);
                    if (fm2 != fFreeMap.end()) {
                        auto freeCands = fm2->second;
                        for (const auto &subCand : freeCands) {
                            if (isSubsetOf(subCand, ns)) {
                                gen = false;
                                break;
                            }
                        }
                    }
                    if (gen) {
                        fFreeMap[sp].push_back(ns);
                        fFreeItemsets.insert(ns);
                    }
                    int jtem = newset.back();
                    newset.pop_back();
                    suffix.emplace_back(jtem, std::move(ijtids), ijsupp, newset);
                }
            }
        }
        items.swap(suffix);
    }
}

void CFDDiscovery::minePatternsDFS(const Itemset &lhs, int rhs, const PartitionTidList &allTids) {
    std::map<int, SimpleTidList> pidlists;
    auto lhsAttrs = fDb.getAttrVector(lhs);
    std::vector<std::pair<Itemset, std::vector<int> > > partitions;
    std::vector<std::vector<std::pair<int,int> > > rhses2;
    std::unordered_map<Itemset, int> ruleIxs;
    std::vector<int> rhses;
    int count = 0;
    for (int pi = 0; pi <= allTids.fTids.size(); pi++) {
        if (pi == allTids.fTids.size() || allTids.fTids[pi] == PartitionTidList::SEP ) {
            const Transaction& trans = fDb.getRow(allTids.fTids[pi-1]);
            Itemset lhsConstants = projection(trans, lhsAttrs);
            auto ruleI = ruleIxs.find(lhsConstants);
            if (ruleI == ruleIxs.end()) {
                rhses.push_back(trans[-1-rhs]);
                rhses2.emplace_back();
                rhses2.back().emplace_back(trans[-1-rhs],count);
                ruleIxs[lhsConstants] = partitions.size();
                partitions.emplace_back(lhsConstants, std::vector<int>());
                partitions.back().second.push_back(count);
            }
            else {
                partitions[ruleI->second].second.push_back(count);
                rhses2[ruleI->second].emplace_back(trans[-1-rhs], count);
            }
            count = 0;
        }
        else {
            count++;
        }
    }
    std::vector<MinerNode<SimpleTidList> > items;
    std::vector<int> psupps(partitions.size());
    int ri = 0;
    for (const auto& rule : partitions) {
        for (int i : rule.first) {
            pidlists[i].push_back(ri);
        }
        int s = 0;
        for (const auto& rhs : rule.second) {
            s += rhs;
        }
        psupps[ri] = s;
        ri++;
    }
    for (const auto& item : pidlists) {
        int psupp = getPartitionSupport(item.second, psupps);
        if (psupp >= fMinSup) {
            Itemset ns = join(itemset(item.first), subset(lhs, -1-fDb.getAttrIndex(item.first)));
            std::set<Itemset> nrParts;
            for (int pid : item.second) {
                nrParts.insert(partitions[pid].first);
            }
            bool gen = true;
            auto sp = std::make_pair(psupp,nrParts.size());
            auto fm2 = fFreeMap.find(sp);
            if (fm2 != fFreeMap.end()) {
                auto freeCands = fm2->second;
                for (const auto &subCand : freeCands) {
                    if (isSubsetOf(subCand, ns)) {
                        gen = false;
                        break;
                    }
                }
            }
            if (gen) {
                fFreeMap[sp].push_back(ns);
                fFreeItemsets.insert(ns);
            }
            items.emplace_back(item.first, item.second, psupp);
        }
    }
    minePatternsDFS(Itemset(), items, lhs, rhs, rhses2, partitions, psupps);
}

void CFDDiscovery::minePatternsDFS(const Itemset &prefix, std::vector<MinerNode<SimpleTidList> > &items,
                                   const Itemset &lhs, int rhs,
                                   std::vector<std::vector<std::pair<int, int> > > &rhses2,
                                   std::vector<std::pair<Itemset, std::vector<int> > > &partitions,
                                   std::vector<int> &psupps) {
    for (int ix = items.size()-1; ix >= 0; ix--) {
        const auto &inode = items[ix];
        Itemset iset = join(prefix, inode.fItem);
        auto nodeAttrs = fDb.getAttrVectorItems(iset);
        bool lhsGen = true;
        int out = (iset.size() == lhs.size()) ? getMaxElem(rhses2[inode.fTids[0]]) : rhs;
        auto sub = join(iset, setdiff(lhs, nodeAttrs));

        if (out > 0 || !isConstRulePartition(inode.fTids, rhses2)) {
            if (fRules.find(out) != fRules.end()) {
                for (const auto &subRule : fRules[out]) {
                    if (out < 0 && !has(subRule, [](int si) -> bool { return si < 0; })) continue;
                    if (precedes(subRule, sub)) {
                        lhsGen = false;
                    }
                }
            }
            if (lhsGen) {
                int e = getPartitionError(inode.fTids, partitions);
                double conf = 1.0 - ((double) e / (double) inode.fSupp);
                if (conf >= fMinConf) {
                    fCFDs.emplace_back(sub, out);
                }
                if (conf >= 1) {
                    fRules[out].push_back(sub);
                    if (out > 0) fRules[-1 - fDb.getAttrIndex(out)].push_back(sub);
                }
            }
        }
        std::vector<MinerNode<SimpleTidList> > suffix;
        for (int j = ix + 1; j < items.size(); j++) {
            const auto &jnode = items[j];
            if (std::binary_search(nodeAttrs.begin(), nodeAttrs.end(), -1-fDb.getAttrIndex(jnode.fItem))) continue;
            Itemset newset = join(iset, jnode.fItem);
            SimpleTidList ijtids = intersection(inode.fTids, jnode.fTids);
            int ijsupp = getPartitionSupport(ijtids, psupps);
            if (ijsupp >= fMinSup){
                bool gen = true;
                auto nas = fDb.getAttrVectorItems(newset);
                Itemset ns = join(newset, setdiff(lhs, nas));
                std::set<Itemset> nrParts;
                for (int pid : ijtids) {
                    nrParts.insert(partitions[pid].first);
                }
                auto sp = std::make_pair(ijsupp, nrParts.size());
                auto fm2 = fFreeMap.find(sp);
                if (fm2 != fFreeMap.end()) {
                    auto freeCands = fm2->second;
                    for (const auto &subCand : freeCands) {
                        if (isSubsetOf(subCand, ns)) {
                            gen = false;
                            break;
                        }
                    }
                }
                if (gen) {
                    fFreeMap[sp].push_back(ns);
                    fFreeItemsets.insert(ns);
                }
                suffix.emplace_back(jnode.fItem, std::move(ijtids), ijsupp);
            }
        }
        if (suffix.size()) {
            std::sort(suffix.begin(), suffix.end(), [this](const MinerNode<SimpleTidList>& a, const MinerNode<SimpleTidList>& b) {
                return support(a.fTids) < support(b.fTids);
            });
            minePatternsDFS(iset, suffix, lhs, rhs, rhses2, partitions, psupps);
        }
    }
}

std::vector<MinerNode<PartitionTidList> > CFDDiscovery::getPartitionSingletons() {
    std::vector<std::vector<SimpleTidList> > partitions(fDb.nrAttrs());
    std::unordered_map<int, std::pair<int,int> > attrIndices;
    for (int a = 0; a < fDb.nrAttrs(); a++) {
        const auto& dom = fDb.getDomain(a);
        partitions[a] = std::vector<SimpleTidList>(dom.size());
        for (int i = 0; i < dom.size(); i++) {
            partitions[a][i].reserve(fDb.frequency(dom[i]));
            attrIndices[dom[i]] = std::make_pair(a, i);
        }
    }
    for (unsigned row = 0; row < fDb.size(); row++) {
        const auto& tup = fDb.getRow(row);
        for (int item : tup) {
            const auto& attrNodeIx = attrIndices.at(item);
            partitions[attrNodeIx.first][attrNodeIx.second].push_back(row);
        }
    }
    std::vector<MinerNode<PartitionTidList> > singletons;
    for (int a = 0; a < fDb.nrAttrs(); a++) {
        int attrItem = -1 - a;
        singletons.push_back(MinerNode<PartitionTidList>(attrItem));
        const auto& dom = fDb.getDomain(a);
        singletons.back().fTids.fTids.reserve(fDb.size()+dom.size()-1);
        singletons.back().fTids.fNrSets = dom.size();
        for (int i = 0; i < dom.size(); i++) {
            auto& ts = singletons.back().fTids.fTids;
            ts.insert(ts.end(), partitions[a][i].begin(), partitions[a][i].end());
            if (i < dom.size() - 1) {
                ts.push_back(PartitionTidList::SEP);
            }
        }
        singletons.back().hashTids();
    }
    return singletons;
}

std::vector<MinerNode<PartitionTidList> > CFDDiscovery::getPartitionSingletons(const SimpleTidList& tids, const Itemset& attrs) {
    std::vector<std::vector<SimpleTidList> > partitions(fDb.nrAttrs());
    std::unordered_map<int, std::pair<int,int> > attrIndices;
    for (int a : attrs) {
        const auto& dom = fDb.getDomain(a);
        partitions[a] = std::vector<SimpleTidList>(dom.size());
        for (int i = 0; i < dom.size(); i++) {
            //partitions[a][i].reserve(fDb.frequency(dom[i]));
            attrIndices[dom[i]] = std::make_pair(a, i);
        }
    }
    for (int row: tids) {
        const auto& tup = fDb.getRow(row);
        for (int a : attrs) {
            int item = tup[a];
            const auto& attrNodeIx = attrIndices.at(item);
            partitions[attrNodeIx.first][attrNodeIx.second].push_back(row);
        }
    }
    std::vector<MinerNode<PartitionTidList> > singletons;
    for (int a : attrs) {
        int attrItem = -1 - a;
        singletons.push_back(MinerNode<PartitionTidList>(attrItem));
        const auto& dom = fDb.getDomain(a);
        singletons.back().fTids.fTids.reserve(tids.size()+dom.size()-1);
        singletons.back().fTids.fNrSets = 0;//dom.size();
        for (int i = 0; i < dom.size(); i++) {
            if (partitions[a][i].size()) {
                singletons.back().fTids.fNrSets++;
                auto &ts = singletons.back().fTids.fTids;
                ts.insert(ts.end(), partitions[a][i].begin(), partitions[a][i].end());
                ts.push_back(PartitionTidList::SEP);
            }
        }
        singletons.back().fTids.fTids.pop_back();
        singletons.back().hashTids();
    }
    return singletons;
}

std::vector<MinerNode<PartitionTidList> > CFDDiscovery::getAllSingletons(int minsup) {
    auto parts = getPartitionSingletons();
    auto singles = getSingletons(minsup);
    for (const auto& sing : singles) {
        parts.push_back(MinerNode<PartitionTidList>(sing.fItem));
        parts.back().fTids = {sing.fTids, 1};
    }
    return parts;
}

std::vector<MinerNode<SimpleTidList> > CFDDiscovery::getSingletons(int minsup) {
    std::vector<MinerNode<SimpleTidList> > singletons;
    std::unordered_map<int, int> nodeIndices;
    for (unsigned item = 1; item <= fDb.nrItems(); item++) {
        if (fDb.frequency(item) >= minsup) {
            singletons.push_back(MinerNode<SimpleTidList>(item, fDb.frequency(item)));
            nodeIndices[item] = singletons.size() - 1;
        }
    }

    for (unsigned row = 0; row < fDb.size(); row++) {
        const auto& tup = fDb.getRow(row);
        for (int item : tup) {
            const auto& it = nodeIndices.find(item);
            if (it != nodeIndices.end()) {
                singletons[it->second].fTids.push_back(row);
            }
        }
    }

    //std::sort(singletons.begin(), singletons.end());
    for (int i = singletons.size() - 1; i >= 0; i--) {
        auto& node = singletons[i];
        node.hashTids();
    }
    return singletons;
}
