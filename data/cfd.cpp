#include "cfd.h"
#include <sstream>
#include <string>

bool isValid(const Itemset& lhs, int rhs) {
    // Discard Variable -> Constant and Constant -> Variable CFDs
    if (rhs < 0) {
        if (lhs.size() && !has(lhs, [](int si) -> bool { return si < 0; })) return false;
    }
    else {
        if (has(lhs, [](int si) -> bool { return si < 0; })) false;
    }
    return true;
}

SimpleTidList getConstantVio(const Itemset& lhs, int rhs, const SimpleTidList& tids, const Database& db) {
    SimpleTidList res;
    int rAttr = db.getAttrIndex(rhs);
    for (int t : tids) {
        const Transaction& trans = db.getRow(t);
        if (trans[rAttr] != rhs) {
            res.push_back(t);
        }
    }
    return res;
}

SimpleTidList getVariableVio(const Itemset& lhs, int rhs, const SimpleTidList& tids, const Database& db) {
    SimpleTidList res;
    int rAttr = -1-rhs;
    std::map<Itemset, std::set<int> > rules;
    for (int t : tids) {
        const Transaction& trans = db.getRow(t);
        Itemset lhsConstants;
        for (int a : lhs) {
            lhsConstants.push_back(trans[-1-a]);
        }
        rules[lhsConstants].insert(trans[rAttr]);
    }
    for (int t : tids) {
        const Transaction& trans = db.getRow(t);
        Itemset lhsConstants;
        for (int a : lhs) {
            lhsConstants.push_back(trans[-1-a]);
        }
        if (rules[lhsConstants].size() > 1) {
            res.push_back(t);
        }
    }
    return res;
}

SimpleTidList getConstantVio(const Itemset& lhs, int rhs, const Database& db) {
    SimpleTidList res;
    int rAttr = db.getAttrIndex(rhs);
    for (int t = 0; t < db.size(); t++) {
        const Transaction& trans = db.getRow(t);
        Transaction sTrans = trans;
        std::sort(sTrans.begin(), sTrans.end());
        if (isSubsetOf(lhs, sTrans)) {
            if (trans[rAttr] != rhs) {
                res.push_back(t);
            }
        }
    }
    return res;
}

SimpleTidList getVariableVio(const Itemset& lhs, int rhs, const Database& db) {
    auto tp = where(lhs, [](int i) { return i >= 0; });
    auto attrs = where(lhs, [](int i) { return i < 0; });
    SimpleTidList res;
    int rAttr = -1-rhs;
    std::map<Itemset, std::set<int> > rules;
    SimpleTidList tids;
    for (int t = 0; t < db.size(); t++) {
        const Transaction& trans = db.getRow(t);
        Transaction sTrans = trans;
        std::sort(sTrans.begin(), sTrans.end());
        if (isSubsetOf(tp, sTrans)) {
            tids.push_back(t);
            Itemset lhsConstants;
            for (int a : attrs) {
                lhsConstants.push_back(trans[-1-a]);
            }
            rules[lhsConstants].insert(trans[rAttr]);
        }
    }
    for (int t : tids) {
        const Transaction& trans = db.getRow(t);
        Itemset lhsConstants;
        for (int a : attrs) {
            lhsConstants.push_back(trans[-1-a]);
        }
        if (rules[lhsConstants].size() > 1) {
            res.push_back(t);
        }
    }
    return res;
}

std::string printCFD(const std::string& relation, const char* chars, const Itemset& lhs, int rhs, Database& db) {
    if (rhs < 0) {
        std::stringstream res;
        res << relation << "(";
        for (int i = 0; i < db.nrAttrs(); i++) {
            res << db.getAttrName(i);
            res << ": $" << chars[i] << "1";
            if (i < db.nrAttrs()-1) {
                res << ", ";
            }
        }
        res << "), ";
        res << relation << "(";
        for (int i = 0; i < db.nrAttrs(); i++) {
            res << db.getAttrName(i);
            res << ": $" << chars[i] << "2";
            if (i < db.nrAttrs()-1) {
                res << ", ";
            }
        }
        res << "), ";
        for (int v : lhs) {
            int a = v > 0 ? db.getAttrIndex(v) : -1-v;
            res << "$" << chars[a] << "1";
            res << " == ";
            res << "$" << chars[a] << "2";
            res << ", ";
            if (v > 0) {
                res << "$" << chars[a] << "1";
                res << " == ";
                res << "\"" << db.getValue(v) << "\"";
                res << ", ";
            }
        }
        int a = -1-rhs;
        res << "$" << chars[a] << "1";
        res << " != ";
        res << "$" << chars[a] << "2";
        res << " -> #fail.";
        return res.str();
    }
    else {
        std::stringstream res;
        res << relation << "(";
        for (int i = 0; i < db.nrAttrs(); i++) {
            res << db.getAttrName(i);
            res << ": $" << chars[i] << "1";
            if (i < db.nrAttrs()-1) {
                res << ", ";
            }
        }
        res << "), ";
        for (int v : lhs) {
            int a = db.getAttrIndex(v);
            res << "$" << chars[a] << "1";
            res << " == ";
            res << "\"" << db.getValue(v) << "\"";
            res << ", ";
        }
        int a = db.getAttrIndex(rhs);
        res << "$" << chars[a] << "1";
        res << " != ";
        res << "\"" << db.getValue(rhs) << "\"";
        res << " -> #fail.";
        return res.str();
    }
}

CFD getCFD(std::vector<int> attrIxs, std::vector<std::string> vals, int rhsAttr, const std::string& rhsVal, Database& db) {
    Itemset lhs(attrIxs.size());
    int rhs = 0;
    for (int i = 0; i < attrIxs.size(); i++) {
        if (vals[i] == "-") {
            lhs[i] = -1 - attrIxs[i];
        }
        else {
            lhs[i] = db.translateToken(attrIxs[i], vals[i]);
        }
    }
    std::sort(lhs.begin(), lhs.end());
    if (rhsVal == "-") {
        rhs = -1 - rhsAttr;
    }
    else {
        rhs = db.translateToken(rhsAttr, rhsVal);
    }
    return std::make_pair(lhs, rhs);
}

CFDList getCFDsFromFile(const std::string& fileName, const Database& db) {
    CFDList res;
    std::ifstream infile(fileName);
    std::string line;
    std::getline(infile, line);
    bool eof = false;
    while (!eof && line.size())
    {
        if (line[0] == '/') {
            std::string lhs = line.substr(line.find(":")+3, line.find(" => ")-line.find(":")-3);
            trim(lhs);
            lhs.pop_back();
            std::string rhs = line.substr(line.find(" => ")+3);
            trim(rhs);
            int rhsItem = getItem(rhs, db);

            std::istringstream iss(lhs);
            std::string lhsPart;
            std::vector<int> lhsItems;
            int i = 0;
            while (std::getline(iss, lhsPart, ',')) {
                trim(lhsPart);
                lhsItems.push_back(getItem(lhsPart, db));
                i++;
            }
            if (i > 0) {
                res.emplace_back(lhsItems, rhsItem);
            }
            if (infile.eof()) {
                eof = true;
            }
        }
        std::getline(infile, line);
    }
    return res;
}

int getItem(const std::string& str, const Database& db) {
    if (str.find('=') == std::string::npos) {
        return -1 - db.getAttr(str);
    }
    else {
        std::string attrStr = str.substr(0, str.find('='));
        trim(attrStr);
        std::string valStr = str.substr(str.find('=')+1);
        trim(valStr);
        int attr = db.getAttr(attrStr);
        return db.getItem(attr, valStr);
    }
}