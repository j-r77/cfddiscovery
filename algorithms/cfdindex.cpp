#include "cfdindex.h"
#include "../util/output.h"

CFDIndex::CFDIndex(const Database& db)
    :fDb(db) {

}

bool CFDIndex::redundant(const LHSValue& super, const LHSValue& sub) {
    const auto superAttrs = fDb.getAttrVectorItems(super);
    const auto subAttrs = fDb.getAttrVectorItems(sub);
    if (!isSubsetOf(subAttrs, superAttrs)) return false;
    for (int item : super) {
        int attr = item >= 0 ? -1 - fDb.getAttrIndex(item) : item;
        if (std::find(subAttrs.begin(), subAttrs.end(), attr) == subAttrs.end()) continue;
        if (item < 0 && std::find(sub.begin(), sub.end(), attr) == sub.end()) return false;
        if (item >= 0 && std::find(sub.begin(), sub.end(), attr) == sub.end()
            && std::find(sub.begin(), sub.end(), item) == sub.end()) {
            return false;
        }
    }
    return true;
}

bool CFDIndex::isLeftRedundant(const LHSValue& lhs, std::map<LHSAttrs, LHSValues>& cands) {
    const auto lhsAttrs = fDb.getAttrVectorItems(lhs);
    for (const std::pair<LHSAttrs, LHSValues>& storedLHS : cands) {
        if (!storedLHS.second.size()) continue;
        if (isStrictSubsetOf(storedLHS.first, lhsAttrs)) return true;
        if (storedLHS.first == lhsAttrs) {
            for (const LHSValue& storedValue : storedLHS.second) {
                if (redundant(lhs, storedValue)) return true;
            }
        }
    }
    return false;
}

/*void CFDIndex::deleteRightRedundant(const LHSValue& lhs, std::map<LHSAttrs, LHSValues>& cands) {
    
}*/

void CFDIndex::deleteLeftRedundant(const LHSValue& lhs, std::map<LHSAttrs, LHSValues>& cands) {
    const auto lhsAttrs = fDb.getAttrVectorItems(lhs);
    std::vector<LHSAttrs> dels;
    for (const std::pair<LHSAttrs, LHSValues>& storedLHS : cands) {
        if (!storedLHS.second.size()) {
            dels.push_back(storedLHS.first);
            continue;
        }
        if (isStrictSubsetOf(lhsAttrs, storedLHS.first)) {
            dels.push_back(storedLHS.first);
            continue;
        }
        if (storedLHS.first == lhsAttrs) {
            std::vector<LHSValue> subdels;
            for (const LHSValue& storedValue : storedLHS.second) {
                if (redundant(storedValue, lhs)) {
                    subdels.push_back(storedValue);
                }
            }
            for (const LHSValue& delValue : subdels) {
                cands[storedLHS.first].erase(std::find(cands[storedLHS.first].begin(), cands[storedLHS.first].end(), delValue));
            }
        }
        if (!storedLHS.second.size()) {
            dels.push_back(storedLHS.first);
        }
    }
    for (const LHSAttrs& delAttrs : dels) {
        cands.erase(delAttrs);
    }
}

void CFDIndex::addCFD(const LHSValue& lhs, const int rhs) {
    const auto lhsAttrs = fDb.getAttrVectorItems(lhs);
    int rhsAttr = rhs >= 0 ? -1 - fDb.getAttrIndex(rhs) : rhs;

    // Process all CFDs with same RHS. If LHS superset of current -> delete
    //deleteLeftRedundant(lhs, fCFDs[rhsAttr][rhs]);
    //if (rhs != rhsAttr) deleteLeftRedundant(lhs, fCFDs[rhsAttr][rhsAttr]);
    //if (rhs != rhsAttr) deleteRightRedundant(lhs, fCFDs[rhsAttr][rhsAttr]);

   /*for (const auto& cfd : fSimple) {
        if (cfd.second == rhs) {
            if (isSubsetOf(cfd.first, lhs)) return;
            if (rhs >= 0) {
                continue;
            }

            auto lAttr = fDb.getAttrVectorItems(lhs);
            if (isSubsetOf(cfd.first, lAttr)) return;
            auto cAttr = fDb.getAttrVector(cfd.first);
            Itemset lConstPart;
            Itemset lVarPart;
            for (int i : lhs) {
                if (i >= 0) lConstPart.push_back(i);
                else lVarPart.push_back(i);
            }
            std::sort(lConstPart.begin(), lConstPart.end());
            std::sort(lVarPart.begin(), lVarPart.end());
            Itemset cConstPart;
            Itemset cVarPart;
            for (int i : cfd.first) {
                if (i >= 0) cConstPart.push_back(i);
                else cVarPart.push_back(i);
            }
            std::sort(cConstPart.begin(), cConstPart.end());
            std::sort(cVarPart.begin(), cVarPart.end());
            if (isSubsetOf(cVarPart, lVarPart) && isStrictSubsetOf(cConstPart, lConstPart)) return;
            Itemset lVars = fDb.getAttrVectorItems(lhs);
            //if (isStrictSubsetOf(cVarPart, lVars)) return;
        }
    }*/
    //if (isLeftRedundant(lhs, fCFDs[rhsAttr][rhs])) return;
    //if (rhs == rhsAttr) {
    //    if (isLeftRedundant(lhs, fCFDs[rhsAttr][rhs])) continue;
    //}
    //if (rhs == rhsAttr) {
    //    if (isRightRedundant(lhs, fCFDs[rhsAttr][rhs])) continue;
    //}
    //Output::printCFD(lhs, rhs, fDb);
    //fCFDs[rhsAttr][rhs][lhsAttrs].push_back(lhs);
    fSimple.emplace_back(lhs, rhs);
    /*
    for (auto it = fCFDs[rhsAttr][rhs].cbegin(); it != fCFDs[rhsAttr][rhs].cend(); ) {
        bool del = false;
        const auto& ll = it->first;
        if (isStrictSubsetOf(lhsAttrs, ll) && it->second.size()) {
            del = true;
        }
        else if (ll == lhsAttrs) {
            for (const auto& llvs : it->second) {
                if (redundant(llvs, lhs)) {
                    del = true;
                    break;
                }
            }
        }
        /*if (del) fCFDs[rhsAttr][rhs].erase(it++);
        else ++it;    
    }
    // Constant RHS: If a CFD with '_' already exists -> delete
    if (rhs > 0) {
        if (fCFDs[rhsAttr][rhsAttr].find(lhsAttrs) != fCFDs[rhsAttr][rhsAttr].end()) {
            bool del = false;
            for (auto llit = fCFDs[rhsAttr][rhsAttr][lhsAttrs].cbegin(); llit != fCFDs[rhsAttr][rhsAttr][lhsAttrs].cend(); ++llit) {
                if (lhs == *llit) {
                    del = true;
                    break;
                }
            }
            /*if (del) llit->second.erase(llit++);
            else ++llit; 
        }
    }
    // Process all CFDs with same RHS. If LHS subset of current -> redundant, return
    for (auto it = fCFDs[rhsAttr][rhs].cbegin(); it != fCFDs[rhsAttr][rhs].cend(); ++it) {
        
    }
    // Variable RHS: If a constant CFD already exists -> redundant, return
    if (rhs < 0) {
        for (auto it = fCFDs[rhsAttr].begin(); it != fCFDs[rhsAttr].end(); ++it) {
            if (it->first == rhs) continue;
            for (auto llit = it->second.cbegin(); llit != it->second.cend(); ++llit) {
                bool del = false;
                const auto& ll = llit->first;
                if (isStrictSubsetOf(lhsAttrs, ll) && it->second.size()) {
                    return;
                }
                else if (ll == lhsAttrs) {
                    for (const auto& llvs : llit->second) {
                        if (redundant(llvs, lhs)) {
                            return;
                        }
                    }
                }
            }
        }
    }
    /*if (rhs > 0) {
        for (const auto& llp : fCFDs[rhsAttr][rhsAttr]) {
            const auto& ll = llp.first;
            if (isSubsetOf(ll, lhsAttrs) || ll == lhsAttrs) {
                for (const auto& llvs : llp.second) {
                    if (redundant(lhs, llvs)) return;
                }
            }
        }
    }*/
    // Variable RHS: Erase all stored Constant RHSes
    /*else {
        for (auto it = fCFDs[rhsAttr].begin(); it != fCFDs[rhsAttr].end(); ++it) {
            if (it->first == rhs) continue;
            for (auto llit = it->second.cbegin(); llit != it->second.cend(); ) {
                bool del = false;
                const auto& ll = llit->first;
                if (ll == lhsAttrs || isSubsetOf(lhsAttrs, ll)) {
                    for (const auto& llvs : llit->second) {
                        if (redundant(llvs, lhs)) {
                            del = true;
                            break;
                        }
                    }
                }
                if (del) it->second.erase(llit++);
                else ++llit;    
            }
        }
    }*/
    //Output::printCFD(lhs, rhs, fDb);
    //std::cout << "Added! " << count() << std::endl;
    //if (!rhsAttr) std::cout << rhsAttr << " " << rhs << std::endl;
}

int CFDIndex::count() const {
    return fSimple.size();
    int cnt = 0;
    std::set<Itemset> attrsets;
    std::set<std::pair<Itemset, int> > test;
    std::map<Itemset,std::vector<std::pair<Itemset, int> > > cfds;
    for (const auto& outerit: fSimple) {
        auto lhsAttrs = fDb.getAttrVectorItems(outerit.first);
        cfds[lhsAttrs].push_back(std::make_pair(outerit.first, outerit.second));
    }
    for (const auto& it: cfds) {
        std::set<int> rhses;
        for (const auto& pat : it.second) {
            rhses.insert(pat.second > 0 ? -1-fDb.getAttrIndex(pat.second) : pat.second);
        }
        cnt += rhses.size();//it.second.size(); 
    }
    return cnt;
}

CFDList CFDIndex::getCFDs() {
    return fSimple;
}