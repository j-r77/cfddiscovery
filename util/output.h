#ifndef UTIL_OUTPUT_H_
#define UTIL_OUTPUT_H_

#include <iostream>
#include <ostream>
#include "../data/database.h"
#include "../algorithms/minernode.h"
#include "../data/cfd.h"

class Output {
public:
	template <typename T>
    static void printCollection(const T& coll, std::ostream& out=std::cout) {
		for (const typename T::value_type& i : coll) {
			out << i << " ";
		}
	}

    template<typename T>
    static void printCollection(const T& items, const std::string& join, std::ostream& out=std::cout) {
        bool comma = false;
        for (const typename T::value_type& item : items) {
            if (comma) {
                out << join;
            }
            else {
                comma = true;
            }
            out << item;
        }
    }

    static void printItemset(const Itemset& items, const Database& db, std::ostream& out=std::cout, bool endl=true) {
        out << "(";
        std::vector<std::string> parts;
        for (uint ix = 0; ix < items.size(); ix++) {
            int item = items[ix];
            if (item < 0) {
                parts.push_back(db.getAttrName(-1-item));
            }
            else if (item == 0) {
                parts.push_back(db.getAttrName(ix) + "=N/A");
            }
            else {
                parts.push_back(db.getAttrName(db.getAttrIndex(item)) + "=" + db.getValue(item));
            }
        }
        printCollection(parts, ", ", out);
        out << ")";
        if (endl) {
            out << std::endl;
        }
    }

    static void printCFDList(const CFDList& cs, const Database& db, std::ostream& out=std::cout, bool endl=true) {
        for (const auto& c : cs) {
            printCFD(c.first, c.second, db, out, endl);
        }
    }

    static void printCFD(const CFD& c, const Database& db, std::ostream& out=std::cout, bool endl=true) {
        printCFD(c.first, c.second, db, out, endl);
    }

    static void printCFD(const Itemset& lhs, const int rhs, const Database& db, std::ostream& out=std::cout, bool endl=true) {
        printItemset(lhs, db, out, false);
        out << " => ";
        if (rhs < 0) {
            out << db.getAttrName(-1-rhs);
        }
        else {
            out << (db.getAttrName(db.getAttrIndex(rhs)) + "=" + db.getValue(rhs));
        }
        if (endl) {
            out << std::endl;
        }
    }
};

#endif //UTIL_OUTPUT_H_
