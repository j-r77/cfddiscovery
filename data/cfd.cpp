#include "cfd.h"

bool isValid(const Itemset& lhs, int rhs) {
    // Discard Variable -> Constant and Constant -> Variable CFDs
    if (rhs < 0) {
        if (lhs.size() && !has(lhs, [](int si) -> bool { return si < 0; })) return false;
    }
    else {
        if (has(lhs, [](int si) -> bool { return si < 0; })) return false;
    }
    return true;
}