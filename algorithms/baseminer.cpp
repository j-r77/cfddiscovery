#include "baseminer.h"
#include "../util/setutil.h"
#include "../util/output.h"
#include <set>

BaseMiner::BaseMiner(Database& db)
    :fDb(db) {
}

std::vector<int> BaseMiner::filterSameAttr(const Itemset& newset, const std::vector<int>& items) const {
    std::vector<int> attrs = fDb.getAttrVector(newset);
    std::vector<int> rightItems;
    for (int i : items) {
        if (i > 0 && std::find(attrs.begin(), attrs.end(), fDb.getAttrIndex(i)) == attrs.end()) {
            rightItems.push_back(i);
        }
        else if (i < 0 && std::find(attrs.begin(), attrs.end(), -1 - i) == attrs.end()) {
            rightItems.push_back(i);
        }
    }
    return rightItems;
}

std::vector<MinerNode<PartitionTidList> > BaseMiner::getPartitionSingletons() {
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

std::vector<MinerNode<PartitionTidList> > BaseMiner::getPartitionSingletons(const SimpleTidList& tids, const Itemset& attrs) {
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

std::vector<MinerNode<PartitionTidList> > BaseMiner::getAllSingletons(int minsup) {
    auto parts = getPartitionSingletons();
    auto singles = getSingletons(minsup);
    for (const auto& sing : singles) {
        parts.push_back(MinerNode<PartitionTidList>(sing.fItem));
        parts.back().fTids = {sing.fTids, 1};
    }
    return parts;
}

PartitionTidList BaseMiner::getTids(const Itemset& items, const Database& db) {
    auto attrs = fDb.getAttrVector(items);
    auto tp = where(items, [](int i){return i >= 0;});
    std::map<Itemset, SimpleTidList> values;
    for (int i = 0; i < db.size(); i++) {
        auto tup = db.getRow(i);
        if (isSubsetOf(tp, sorted(tup))) {
            auto vals = projection(tup, attrs);
            values[vals].push_back(i);
        }
    }
    PartitionTidList res { SimpleTidList(), 0 };
    for (const auto& kvp : values) {
        res.fTids.insert(res.fTids.end(), kvp.second.begin(), kvp.second.end());
        res.fTids.push_back(PartitionTidList::SEP);
        res.fNrSets++;
    }
    if (res.fTids.size()) res.fTids.pop_back();
    return res;
}

/*std::vector<MinerNode<PartitionTidList> > BaseMiner::getPartitionSingletons(int minsup, const Database& db) {
    std::vector<MinerNode<PartitionTidList> > singletons;
    std::unordered_map<int, int> nodeIndices;
    std::unordered_map<int, std::pair<int,int> > attrIndices;

    for (int a = 0; a < db.nrAttrs(); a++) {
        int attrItem = -1 - a;
        singletons.push_back(MinerNode<PartitionTidList>(attrItem));
        const auto& dom = db.getDomain(a);
        singletons.back().tidmap.resize(dom.size());
        for (int i = 0; i < dom.size(); i++) {
            attrIndices[dom[i]] = std::make_pair(singletons.size() - 1, i);
            singletons.back().tidmap[i].reserve(db.frequency(dom[i]));
        }
    }

    for (unsigned item = 1; item <= db.nrItems(); item++) {
        if (db.frequency(item) >= minsup) {
            singletons.push_back(MinerNode<PartitionTidList>(item, db.frequency(item)));
            singletons.back().tidmap.resize(1);
            nodeIndices[item] = singletons.size() - 1;
        }
    }

    for (unsigned row = 0; row < db.size(); row++) {
        const auto& tup = db.getRow(row);
        for (int item : tup) {
            const auto& attrNodeIx = attrIndices.at(item);
            singletons[attrNodeIx.first].tidmap[attrNodeIx.second].push_back(row);
            const auto& it = nodeIndices.find(item);
            if (it != nodeIndices.end()) {
                singletons[it->second].tidmap[0].push_back(row);
            }
        }
    }

    std::sort(singletons.begin(), singletons.end());
    //fItemFrequencies.reserve(singletons.size());
    for (int i = singletons.size() - 1; i >= 0; i--) {
        auto& node = singletons[i];
        node.hashTids();
        //fItemFrequencies.push_back(std::make_pair(node.fItem, node.fSupp));
    }
    return singletons;
}*/

std::vector<MinerNode<SimpleTidList> > BaseMiner::getSingletons(int minsup) {
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

bool BaseMiner::getDiffsets(std::vector<Diffset>& dss) {
    return getDiffsets(dss, fDb);
}

bool BaseMiner::getDiffsets(std::vector<Diffset>& dss, const Database& db) {
    std::set<Diffset> diffs2;
    for (int i = 0; i < db.size() - 1; i++) {
        const Transaction& ti = db.getRow(i);
        for (int j = i + 1; j < db.size(); j++) {
            const Transaction& tj = db.getRow(j);
            Diffset d;
            for (int ix = ti.size()-1; ix >= 0; ix--) {
                if (ti[ix] != tj[ix]) {
                    d.push_back(-1 -ix);
                }
            }
            if (d.size()) {
                diffs2.insert(d);
            }
        }
    }
    dss = std::vector<Diffset>(diffs2.begin(), diffs2.end());
    return true;
}

/*bool BaseMiner::getDiffsets(std::vector<Diffset>& dss, const Database& db) {
    std::set<Diffset> diffs;
    for (int i = 0; i < db.size() - 1; i++) {
        std::cout << "DIFFSET: " << i << std::endl;
        const Transaction& ti = db.getRow(i);
        for (int j = i + 1; j < db.size(); j++) {
            const Transaction& tj = db.getRow(j);
            Diffset d;
            for (int ix = ti.size()-1; ix >= 0; ix--) {
                if (ti[ix] != tj[ix]) {
                    d.push_back(-1 -ix);
                }
            }
            if (d.size()) {
                diffs.insert(d);
            }
        }
    }
    dss = std::vector<Diffset>(diffs.begin(), diffs.end());
    return true;
}*/

bool BaseMiner::getDiffsets(const SimpleTidList& tids, std::vector<Diffset>& dss) {
    std::set<Diffset> diffs;
    for (int i = 0; i < tids.size(); i++) {
        const Transaction& ti = fDb.getRow(tids[i]);
        for (int j = i + 1; j < tids.size(); j++) {
            const Transaction& tj = fDb.getRow(tids[j]);
            Diffset d;
            for (int ix = ti.size()-1; ix >= 0; ix--) {
                if (ti[ix] != tj[ix]) {
                    d.push_back(-1 -ix);
                }
            }
            if (d.size()) {
                diffs.insert(d);
            }
        }
    }
    dss = std::vector<Diffset>(diffs.begin(), diffs.end());
    return true;
}

std::vector<Diffset> BaseMiner::projectDiffsets(const std::vector<Diffset> diffs, int attr, std::vector<std::pair<int,int> >& attrCounts) {
    std::vector<Diffset> filterDiffs;
    attrCounts.reserve(fDb.nrAttrs());
    for (int a = 0; a < fDb.nrAttrs(); a++) {
        attrCounts.emplace_back(-1-a, 0);
    }
    //attrCounts = std::vector<std::pair<int,int> >(fDb.nrAttrs());
    for (const Diffset& diff : diffs) {
        if (std::find(diff.begin(), diff.end(), (-1-attr)) != diff.end()) {
            const auto ds = subset(diff, (-1-attr));
            if (!containsSubsetOf(filterDiffs,ds)) {
                filterDiffs.push_back(ds);
                for (int a : ds) {
                    if (a) attrCounts[-1-a].second++;
                }
            }
        }
    }
    std::sort(attrCounts.begin(), attrCounts.end(), 
        [](const std::pair<int,int>& a, const std::pair<int,int>& b) -> bool
    { 
        return a.second > b.second; 
    });
    std::vector<Diffset> res;
    for (const Diffset& ds : filterDiffs) {
        if (!containsStrictSubsetOf(filterDiffs, ds)) {
            res.push_back(ds);
        }
    }
    return filterDiffs;//std::vector<Diffset>(filterDiffs.begin(), filterDiffs.end());
}

std::vector<Diffset> BaseMiner::filterDiffsets(const std::vector<Diffset> diffs, int attr, std::vector<std::pair<int,int> >& attrCounts) {
    std::set<Diffset> filterDiffs;
    attrCounts.reserve(fDb.nrAttrs());
    for (int a = 0; a < fDb.nrAttrs(); a++) {
        attrCounts.emplace_back(-1-a, 0);
    }

    for (const Diffset& diff : diffs) {
        if (std::find(diff.begin(), diff.end(), attr) == diff.end()) {
            filterDiffs.insert(diff);
            for (int a : diff) {
                if (a) attrCounts[-1-a].second++;
            }
        }
    }
    std::sort(attrCounts.begin(), attrCounts.end(), 
        [](const std::pair<int,int>& a, const std::pair<int,int>& b) -> bool
    { 
        return a.second > b.second; 
    });
    return std::vector<Diffset>(filterDiffs.begin(), filterDiffs.end());
}

std::vector<SimpleTidList> BaseMiner::bucketTids(const std::vector<int>& items, const SimpleTidList& nodeTids) const {
    std::vector<SimpleTidList> ijtids(items.size());
    for (int ix = 0; ix < items.size(); ix++) {
        for (int t : nodeTids) {
            const Transaction& tup = fDb.getRow(t);
            if (std::find(tup.begin(), tup.end(), items[ix]) != tup.end()) {
                ijtids[ix].push_back(t);
            }
        }
    }
    return ijtids;
}

int BaseMiner::project(const SimpleTidList& tids, int attr) {
    int proj = -1;
    for (int t : tids) {
        if (proj < 0) {
            proj = fDb.getRow(t)[attr];
        }
        else if (proj != fDb.getRow(t)[attr]) {
            return -1;
        }
    }
    return proj;
}


int BaseMiner::projectConf(const SimpleTidList& tids, int attr, double& conf) {
    std::map<int, int> counts;
    for (int t : tids) {
        counts[fDb.getRow(t)[attr]]++;
    }
    auto max = std::max_element(counts.begin(), counts.end(),
                              [](const std::pair<int, int>& p1, const std::pair<int, int>& p2) {
                                  return p1.second < p2.second;
                              });
    conf = (max->second/(double)tids.size());
    if (max->first > fDb.nrItems()) {
        std::cout << max->first << " " << attr << " " << tids.size() << std::endl;
    }
    return max->first;
}


Itemset BaseMiner::subsetItems(const Itemset& items, const Itemset& attrs) {
    Itemset res(attrs.size());
    for (int ix = 0; ix < attrs.size(); ix++) {
        res[attrs.size() - 1 - ix] = items[-1-attrs[ix]];
    }
    std::sort(res.begin(), res.end());
    return res;
}

int BaseMiner::subsetItems(const Itemset& items, int attr) {
    return items[-1-attr];
}