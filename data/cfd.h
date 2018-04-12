#ifndef DATA_CFD_H_
#define DATA_CFD_H_

#include "types.h"
#include "../data/database.h"
#include "../util/stringutil.h"

typedef std::pair<Itemset, int> CFD;
typedef std::vector<CFD> CFDList;

SimpleTidList getConstantVio(const Itemset&, int, const SimpleTidList&, const Database&);
SimpleTidList getVariableVio(const Itemset&, int, const SimpleTidList&, const Database&);
SimpleTidList getConstantVio(const Itemset& lhs, int rhs, const Database& db);
SimpleTidList getVariableVio(const Itemset& lhs, int rhs, const Database& db);
bool isValid(const Itemset& lhs, int rhs);

std::string printCFD(const std::string& relation, const char* chars, const Itemset& lhs, int rhs, Database& db);
CFD getCFD(std::vector<int> attrIxs, std::vector<std::string> vals, int rhsAttr, const std::string& rhsVal, Database& db);
CFDList getCFDsFromFile(const std::string& fileName, const Database& db);
int getItem(const std::string&, const Database&);

#endif //DATA_CFD_H_