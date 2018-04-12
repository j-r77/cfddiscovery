#include "database.h"
#include <sstream>
#include <iostream>
#include <algorithm>
#include "../util/output.h"

Database::Database()
    :fNumTokens(1) {

}

unsigned Database::size() const {
    return fData.size();
}

unsigned Database::nrAttrs() const {
    return fAttributes.size();
}

unsigned Database::nrItems() const {
    return fItems.size();
}

const Transaction& Database::getRow(int row) const {
    return fData.at(row);
}

void Database::setRow(int r, const Transaction& row) {
    for (int i = 0; i < row.size(); i++) {
        if (fData[r][i] != row[i]) {
            fItems[fData[r][i] - 1].fFrequency--;
            fItems[row[i] - 1].fFrequency++;
        }
    }
    fData[r] = row;
}

void Database::addRow(const Transaction& t) {
    fData.push_back(t);
}

void Database::setAttributes(const std::vector<std::string>& attrs) {
    fAttributes.reserve(attrs.size());
    for (int i = 0; i < attrs.size(); i++) {
        fAttributes.emplace_back(attrs[i]);
    }
}

const std::vector<int>& Database::getDomainOfItem(int item) const {
    return fAttributes.at(fItems[item - 1].fAttribute).fValues;
}

const std::vector<int>& Database::getDomain(int attr) const {
    return fAttributes.at(attr).fValues;
}

int Database::translateToken(int attr, const std::string& strVal) {
    auto ptr = fItemDictionary.find(std::make_pair(attr, strVal));
    if (ptr != fItemDictionary.end()) return ptr->second;

    fItems.emplace_back(strVal, attr);
    fAttributes[attr].fValues.push_back(fNumTokens);
    fItemDictionary[std::make_pair(attr, strVal)] = fNumTokens;
    return fNumTokens++;
}

int Database::getAttr(const std::string& s) const {
    for (int i = 0; i < fAttributes.size(); i++) {
        const AttributeInfo& ai = fAttributes[i];
        if (ai.fName == s) {
            return i;
        }
    }
    return -1;
}

int Database::getItem(int attr, const std::string& strVal) const {
    return fItemDictionary.at(std::make_pair(attr, strVal));
}

void Database::sort() {
    std::sort(fData.begin(), fData.end(),
        [](const Transaction& a, const Transaction& b)
        {
            for (int i = 0; i < a.size(); i++) {
                if (a[i] < b[i]) return true;
                if (b[i] < a[i]) return false;
            }
            return false;
        });
}

void Database::toFront(const SimpleTidList& tids) {
    for (int i = 0; i < tids.size(); i++) {
        std::swap(fData[i], fData[tids[i]]);
    }
}


void Database::incFreq(int i) {
    fItems[i - 1].fFrequency++;
}

int Database::frequency(int i) const {
    return fItems[i - 1].fFrequency;
}

int Database::getAttrIndex(int i) const {
    if (i > 0) return fItems[i - 1].fAttribute;
    return -1 - i;
}

void Database::writeInt(std::ofstream& file, char delim) const {
    for (const auto& row : fData) {
        for (const int& i : row) {
            file << i << delim;
        }
        file << std::endl;
    }
}

void Database::write(std::ofstream& file, char delim) const {
    for (int ai = 0; ai < fAttributes.size(); ai++) {
        const auto& attr = fAttributes[ai];
        file << attr.fName;
        if (ai < fAttributes.size()-1)
            file << delim;
        else
            file << std::endl;
    }
    for (const auto& row : fData) {
        for (int ri = 0; ri < row.size(); ri++) {
            const auto& item = row[ri];
            file << fItems[item-1].fValue;
            if (ri < row.size()-1)
                file << delim;
            else
                file << std::endl;
        }
    }
}

void Database::write(std::ofstream& file, const SimpleTidList& subset, char delim) const {
    for (int ai = 0; ai < fAttributes.size(); ai++) {
        const auto& attr = fAttributes[ai];
        file << attr.fName;
        if (ai < fAttributes.size()-1)
            file << delim;
        else
            file << std::endl;
    }
    for (int i : subset) {
        auto& row = fData[i];
        for (int ri = 0; ri < row.size(); ri++) {
            const auto& item = row[ri];
            file << fItems[item-1].fValue;
            if (ri < row.size()-1)
                file << delim;
            else
                file << std::endl;
        }
    }
}

const std::string& Database::getAttrName(int i) const {
    return fAttributes[i].fName;
}

const std::string& Database::getValue(int i) const {
    return fItems[i - 1].fValue;
}

std::vector<int> Database::getAttrVector(const Itemset& items) const {
    std::vector<int> attrs;
    attrs.reserve(items.size());
    for (int i : items) {
        if (i > 0)
            attrs.push_back(getAttrIndex(i));
        else
            attrs.push_back(-1 - i);
    }
    std::sort(attrs.begin(), attrs.end());
    return attrs;
}

std::vector<int> Database::getAttrVectorItems(const Itemset& items) const {
    std::vector<int> attrs;
    attrs.reserve(items.size());
    for (int i : items) {
        if (i > 0)
            attrs.push_back(-1 - getAttrIndex(i));
        else
            attrs.push_back(i);
    }
    std::sort(attrs.begin(), attrs.end());
    return attrs;
}

std::vector<int> Database::getDiffs(const Database& lhs, const Database& rhs) {
    std::vector<int> diffs;
    for (int i = 0; i < lhs.size(); i++) {
        if (lhs.getRow(i) != rhs.getRow(i)) {
            diffs.push_back(i);
        }
    }
    return diffs;
}

std::unordered_map<std::pair<int, std::string>, int, pairhash> Database::getDictionary() const {
    return fItemDictionary;
}

void Database::setDictionary(const std::unordered_map<std::pair<int, std::string>, int, pairhash>& dict) {
    fItems.resize(dict.size());
    for (const auto& kvp : dict) {
        int attr = kvp.first.first;
        std::string strVal = kvp.first.second;
        int val = kvp.second;
        fItemDictionary[std::make_pair(attr, strVal)] = val;
        fItems[val-1] = ItemInfo(strVal, attr);
        fAttributes[attr].fValues.push_back(val);
        if (val >= fNumTokens) {
            fNumTokens = val + 1;
        }
    }
    for (auto& att : fAttributes) {
        std::sort(att.fValues.begin(), att.fValues.end());
    }
}