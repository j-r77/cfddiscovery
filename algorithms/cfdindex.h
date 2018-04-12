#ifndef ALGORITHMS_CFDINDEX_H_
#define ALGORITHMS_CFDINDEX_H_

#include "../data/database.h"
#include "../util/setutil.h"
#include "../data/cfd.h"

class CFDIndex {
    public:
        typedef int RHSAttr;
        typedef int RHSValue;
        typedef Itemset LHSAttrs;
        typedef Itemset LHSValue;
        typedef std::vector<Itemset> LHSValues;

        CFDIndex(const Database&);
        void addCFD(const LHSValue& lhs, const int rhs);
        int count() const;
        bool redundant(const LHSValue&, const LHSValue&);
        bool isLeftRedundant(const LHSValue&, std::map<LHSAttrs, LHSValues>&);
        void deleteLeftRedundant(const LHSValue&, std::map<LHSAttrs, LHSValues>&);
        std::map<RHSAttr, std::map<RHSValue, std::map<LHSAttrs, LHSValues> > > fCFDs;
        CFDList getCFDs();

    private:
        std::vector<std::pair<Itemset, int> > fSimple;
        const Database& fDb;
};

#endif //ALGORITHMS_CFDINDEX_H_