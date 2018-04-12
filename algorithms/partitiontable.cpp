#include "partitiontable.h"
#include "../util/output.h"

int PartitionTable::fDbSize;

std::vector<PartitionTidList> PartitionTable::intersection(const PartitionTidList &lhs, const std::vector<PartitionTidList*> rhses) {
    std::unordered_map<int, int> eqIndices(support(lhs.fTids));
    std::vector<std::vector<int> > eqClasses(lhs.fNrSets);
    // Construct a lookup from tid to equivalence class
    int eix = 0;
    int count = 0;
    for (int ix = 0; ix <= lhs.fTids.size(); ix++) {
        count++;
        if (ix == lhs.fTids.size() || lhs.fTids[ix] == PartitionTidList::SEP) {
            eqClasses[eix++].reserve(count);
            count = 0;
        }
        else {
            eqIndices[lhs.fTids[ix]] = eix+1;
        }
    }
    std::vector<PartitionTidList> res;
    for (const PartitionTidList* rhs : rhses) {
        res.emplace_back();
        res.back().fNrSets = 0;
        res.back().fTids.reserve(lhs.fTids.size());
        for (int ix = 0; ix <= rhs->fTids.size(); ix++) {
            if (ix == rhs->fTids.size() || rhs->fTids[ix] == PartitionTidList::SEP) {
                for (auto& eqcl : eqClasses) {
                    if (eqcl.size()) {
                        res.back().fTids.insert(res.back().fTids.end(),eqcl.begin(),eqcl.end());
                        res.back().fTids.push_back(PartitionTidList::SEP);
                        res.back().fNrSets++;
                        eqcl.clear();
                    }
                }
            }
            else {
                const int jt = rhs->fTids[ix];
                if (eqIndices.count(jt)) {
                    eqClasses[eqIndices[jt] - 1].push_back(jt);
                }
            }
        }
        if (res.back().fTids.size() && res.back().fTids.back() == PartitionTidList::SEP) {
            res.back().fTids.pop_back();
        }
    }
    return res;
}

std::vector<PartitionTidList> PartitionTable::intersection(const PartitionTidList &lhs, const std::vector<const PartitionTidList*> rhses) {
    std::unordered_map<int, int> eqIndices(support(lhs.fTids));
    std::vector<std::vector<int> > eqClasses(lhs.fNrSets);
    // Construct a lookup from tid to equivalence class
    int eix = 0;
    int count = 0;
    for (int ix = 0; ix <= lhs.fTids.size(); ix++) {
        count++;
        if (ix == lhs.fTids.size() || lhs.fTids[ix] == PartitionTidList::SEP) {
            eqClasses[eix++].reserve(count);
            count = 0;
        }
        else {
            eqIndices[lhs.fTids[ix]] = eix+1;
        }
    }
    std::vector<PartitionTidList> res;
    for (const PartitionTidList* rhs : rhses) {
        res.emplace_back();
        res.back().fNrSets = 0;
        res.back().fTids.reserve(lhs.fTids.size());
        for (int ix = 0; ix <= rhs->fTids.size(); ix++) {
            if (ix == rhs->fTids.size() || rhs->fTids[ix] == PartitionTidList::SEP) {
                for (auto& eqcl : eqClasses) {
                    if (eqcl.size()) {
                        res.back().fTids.insert(res.back().fTids.end(),eqcl.begin(),eqcl.end());
                        res.back().fTids.push_back(PartitionTidList::SEP);
                        res.back().fNrSets++;
                        eqcl.clear();
                    }
                }
            }
            else {
                const int jt = rhs->fTids[ix];
                if (eqIndices.count(jt)) {
                    eqClasses[eqIndices[jt] - 1].push_back(jt);
                }
            }
        }
        if (res.back().fTids.size() && res.back().fTids.back() == PartitionTidList::SEP) {
            res.back().fTids.pop_back();
        }
    }
    return res;
}

PartitionTidList PartitionTable::intersection(const PartitionTidList &lhs, const PartitionTidList &rhs) {
    std::unordered_map<int, int> eqIndices;
    eqIndices.reserve(lhs.fTids.size() + 1 - lhs.fNrSets);
    std::vector<std::vector<int> > eqClasses(lhs.fNrSets);

    // Construct a lookup from tid to equivalence class
    int eix = 0;
    int count = 0;
    for (int ix = 0; ix <= lhs.fTids.size(); ix++) {
        count++;
        if (ix == lhs.fTids.size() || lhs.fTids[ix] == PartitionTidList::SEP) {
            eqClasses[eix].reserve(count);
            count = 0;
            eix++;
        }
        else {
            eqIndices[lhs.fTids[ix]] = eix+1;
        }
    }
    // For each rhs partition, for each eq class: spread all tids over eq classes in lhs
    PartitionTidList res;
    res.fNrSets = 0;
    for (int ix = 0; ix <= rhs.fTids.size(); ix++) {
        if (ix == rhs.fTids.size() || rhs.fTids[ix] == PartitionTidList::SEP) {
            for (auto& eqcl : eqClasses) {
                if (eqcl.size()) {
                    res.fTids.insert(res.fTids.end(),eqcl.begin(),eqcl.end());
                    res.fTids.push_back(PartitionTidList::SEP);
                    res.fNrSets++;
                    eqcl.clear();
                }
            }
        }
        else {
            const int jt = rhs.fTids[ix];
            if (eqIndices[jt]) {
                eqClasses[eqIndices[jt] - 1].push_back(jt);
            }
        }
    }
    if (res.fTids.size() && res.fTids.back() == PartitionTidList::SEP) {
        res.fTids.pop_back();
    }
    return res;
}

std::vector<int> PartitionTable::partitionMap2(const PartitionTidList& x, const PartitionTidList& xa) {
    std::vector<int> res;

    // Map all partitions in xa to their size; use a random tid as identifier
    std::unordered_map<int,int> bigt;
    bigt.reserve(xa.fNrSets);
    int count = 0;
    for (int pi = 0; pi <= xa.fTids.size(); pi++) {
        if (xa.fTids[pi] == PartitionTidList::SEP || pi == xa.fTids.size()) {
            bigt[xa.fTids[pi-1]] = count;
            count = 0;
        }
        else {
            count++;
        }
    }

    int eix = 0;
    int rep = -1;
    std::vector<int> bins;
    for (int cix = 0; cix < x.fTids.size(); cix++) {
        if (x.fTids[cix] == PartitionTidList::SEP) {
            //res[eix] = std::make_pair(rep, bins);
            res.push_back(PartitionTidList::SEP);
            rep = -1;
            bins.clear();
            eix++;
        }
        else {
            int t = x.fTids[cix];
            if (bigt.count(t)) {
                if (rep < 0) {
                    rep = t;
                    res.push_back(t);
                }
                res.push_back(bigt[t]);
                //res[eix].first = t;
                //res[eix].second.push_back(bigt[t]);
            }
        }
    }
    return res;
}

std::vector<std::pair<int, std::vector<int> > > PartitionTable::partitionMap(const PartitionTidList& x, const PartitionTidList& xa) {
    std::vector<std::pair<int, std::vector<int> > > res(x.fNrSets);

    // Map all partitions in xa to their size; use a random tid as identifier
    std::unordered_map<int,int> bigt;
    bigt.reserve(xa.fNrSets);
    int count = 0;
    for (int pi = 0; pi <= xa.fTids.size(); pi++) {
        if (xa.fTids[pi] == PartitionTidList::SEP || pi == xa.fTids.size()) {
            bigt[xa.fTids[pi-1]] = count;
            count = 0;
        }
        else {
            count++;
        }
    }

    int eix = 0;
    int rep = -1;
    std::vector<int> bins;
    for (int cix = 0; cix < x.fTids.size(); cix++) {
        if (x.fTids[cix] == PartitionTidList::SEP) {
            res[eix] = std::make_pair(rep, bins);
            rep = -1;
            bins.clear();
            eix++;
        }
        else {
            int t = x.fTids[cix];
            if (bigt.count(t)) {
                rep = t;
                bins.push_back(bigt[t]);
                //res[eix].first = t;
                //res[eix].second.push_back(bigt[t]);
            }
        }
    }
    return std::vector<std::pair<int, std::vector<int> > >();
}

int PartitionTable::partitionError(const PartitionTidList& x, const PartitionTidList& xa) {
    int e = 0;

    std::map<int,int> bigt;
    //bigt.reserve(xa.fNrSets);
    int count = 0;
    for (int pi = 0; pi <= xa.fTids.size(); pi++) {
        if (pi == xa.fTids.size() || xa.fTids[pi] == PartitionTidList::SEP) {
            bigt[xa.fTids[pi-1]] = count;
            count = 0;
        }
        else {
            count++;
        }
    }

    int eix = 0;
    count = 0;
    int m = 0;
    for (int cix = 0; cix <= x.fTids.size(); cix++) {
        if (cix == x.fTids.size() || x.fTids[cix] == PartitionTidList::SEP) {
            e += count - m;
            eix++;
            m = 0;
            count = 0;
        }
        else {
            count++;
            int t = x.fTids[cix];
            if (bigt.count(t)) {
                if (bigt[t] > m) {
                    m = bigt[t];
                }
            }
        }
    }
    return e;
}

bool PartitionTable::violatedInCleaned(const PartitionTidList& x, const PartitionTidList& xa, int rep) {
    std::unordered_map<int,int> bigt;
    bigt.reserve(xa.fNrSets);
    bool first = true;
    for (int pi = 0; pi <= xa.fTids.size(); pi++) {
        if (pi == xa.fTids.size() || xa.fTids[pi] == PartitionTidList::SEP ) {
            first = true;
        }
        else if (first) {
            bigt[xa.fTids[pi]] = 1;
            first = false;
        }
    }

    bool inCleaned = false;
    for (int cix = 0; cix <= x.fTids.size(); cix++) {
        if (cix == x.fTids.size() || x.fTids[cix] == PartitionTidList::SEP  ) {
            inCleaned = false;
        }
        else {
            int t = x.fTids[cix];
            if (bigt.count(t)) {
                if (inCleaned && t < rep) {
                    return true;
                }
                else if (t < rep) {
                    inCleaned = true;
                }
            }
        }
    }
    return false;
}

SimpleTidList PartitionTable::violations(const PartitionTidList& x, const PartitionTidList& xa) {
    SimpleTidList res;

    std::set<int> bigt;
    for (int pi = 0; pi <= xa.fTids.size(); pi++) {
        if (pi == xa.fTids.size() || xa.fTids[pi] == PartitionTidList::SEP) {
            bigt.insert(xa.fTids[pi-1]);
        }
    }

    int refs = 0;
    SimpleTidList part;
    for (int cix = 0; cix <= x.fTids.size(); cix++) {
        if (cix == x.fTids.size() || x.fTids[cix] == PartitionTidList::SEP) {
            if (refs != 1) {
                res.insert(res.end(), part.begin(), part.end());
            }
            part.clear();
            refs = 0;
        }
        else {
            int t = x.fTids[cix];
            part.push_back(t);
            if (bigt.count(t)) {
                refs++;
            }
        }
    }
    std::sort(res.begin(), res.end());
    return res;
}

SimpleTidList PartitionTable::violations(const PartitionTidList& x, const PartitionTidList& xa, int& e) {
    SimpleTidList res;
    e = 0;

    std::map<int,int> bigt;
    int count = 0;
    for (int pi = 0; pi <= xa.fTids.size(); pi++) {
        if (pi == xa.fTids.size() || xa.fTids[pi] == PartitionTidList::SEP) {
            bigt[xa.fTids[pi-1]] = count;
            count = 0;
        }
        else {
            count++;
        }
    }

    int eix = 0;
    count = 0;
    int m = 0;
    int refs = 0;
    SimpleTidList part;
    for (int cix = 0; cix <= x.fTids.size(); cix++) {
        if (cix == x.fTids.size() || x.fTids[cix] == PartitionTidList::SEP) {
            if (refs > 1) {
                res.insert(res.end(), part.begin(), part.end());
            }
            part.clear();
            refs = 0;
            e += count - m;
            eix++;
            m = 0;
            count = 0;
        }
        else {
            count++;
            int t = x.fTids[cix];
            part.push_back(t);
            if (bigt.count(t)) {
                if (bigt[t] > m) {
                    m = bigt[t];
                }
                refs++;
            }
        }
    }
    std::sort(res.begin(), res.end());
    return res;
}

int PartitionTable::partitionError(const std::vector<std::vector<int> >& x, const std::vector<std::vector<int> >& xa) {
    int e = 0;
    std::map<int,int> bigt;
    for (int pi = 0; pi < xa.size(); pi++) {
        bigt[xa[pi][0]] = xa[pi].size();
    }

    for (const auto& c : x) {
        int m = 1;
        for (int t : c) {
            if (bigt.count(t)) {
                if (bigt[t] > m) {
                    m = bigt[t];
                }
            }
        }
        e += c.size() - m;
    }
    return e;
}