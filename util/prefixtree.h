#ifndef UTIL_PREFIXTREE_H_
#define UTIL_PREFIXTREE_H_

#include <map>
#include <unordered_set>
#include <queue>
#include <iostream>
#include <algorithm>

template<typename T>
size_t hashCollection(const T& xs) {
	size_t res = 0;
	for (const typename T::value_type& x : xs) {
		res ^= std::hash<typename T::value_type>()(x) + 0x9e3779b9 + (res << 6) + (res >> 2);
	}
	return res;
}

template <typename Key, typename Value>
class PrefixTree {
public:
	PrefixTree();
	void reserve(int);
	void insert(const Key&, const Value&);
	Value* find(const Key&) const;
	void erase(const Key&) const;
	bool hasSubset(const Key&, const Value&) const;
	bool hasStrictSubset(const Key&, const Value&) const;
	std::vector<Key> getSubsets(const Key&, const Value&) const;
	std::vector<const Value*> getSubsets(const Key&) const;
	std::vector<Key> getSets() const;
	std::map<Key, int> getSupports(const Key&) const;
private:
	struct PrefixNode {
		Value fValue;
		std::map<typename Key::value_type, PrefixNode> fSubTrees;
		PrefixNode* fParent;
		typename Key::value_type fKey;
		int fDepth;
	};
private:
	int size;
	PrefixNode fRoot;
	std::unordered_set<size_t> fHashes;
	std::unordered_map<size_t, std::vector<PrefixNode*> > fJumps;
};

template <typename Key, typename Value>
PrefixTree<Key,Value>::PrefixTree() {
	fRoot.fDepth = 0;
	fRoot.fValue = Value();
}

template <typename Key, typename Value>
void PrefixTree<Key,Value>::reserve(int r) {
	size = r;
	//fRoot.fSubTrees.reserve(r);
}

template <typename Key, typename Value>
void PrefixTree<Key,Value>::insert(const Key& k, const Value& v) {
	typename Key::const_iterator it = k.begin();
	PrefixNode* insertionPoint = &fRoot;
	while (it != k.end()) {
		PrefixNode* next = &insertionPoint->fSubTrees[*it];
		if (!next->fParent) {
			next->fKey = *it;
			next->fDepth = insertionPoint->fDepth + 1;
			next->fValue = Value();
			next->fParent = insertionPoint;
		}
		insertionPoint = next;
		it++;
	}
	fJumps[hashCollection(k)].push_back(insertionPoint);
	insertionPoint->fValue = v;
}

template <typename Key, typename Value>
Value* PrefixTree<Key,Value>::find(const Key& k) const {
	const auto& aaa = fJumps.find(hashCollection(k));
	if (aaa == fJumps.end()) return 0;

	for (PrefixNode* pn : aaa->second) {
		if (k.size() != pn->fDepth) continue;
		PrefixNode* src = pn;
		typename Key::const_reverse_iterator rit = k.rbegin();
		while (rit != k.rend() && pn && pn->fKey == *rit) {
			pn = pn->fParent;
			rit++;
		}
		if (rit == k.rend() && pn == &fRoot) {
			if (src->fValue != Value())
				return &src->fValue;
			else
				return 0;
		}
	}
	return 0;
}

template <typename Key, typename Value>
void PrefixTree<Key,Value>::erase(const Key& k) const {
	const auto& aaa = fJumps.find(hashCollection(k));
	if (aaa == fJumps.end()) return;

	for (PrefixNode* pn : aaa->second) {
		if (k.size() != pn->fDepth) continue;
		PrefixNode* src = pn;
		typename Key::const_reverse_iterator rit = k.rbegin();
		while (rit != k.rend() && pn && pn->fKey == *rit) {
			pn = pn->fParent;
			rit++;
		}
		if (rit == k.rend() && pn == &fRoot) {
			src->fValue = Value();
			return;
		}
	}
}

template <typename Key, typename Value>
bool PrefixTree<Key,Value>::hasSubset(const Key& k, const Value& v) const {
	std::queue<const PrefixNode*> fringe;
	fringe.push(&fRoot);
	while (!fringe.empty()) {
		const PrefixNode* elem = fringe.front();
		fringe.pop();
		if (elem->fValue == v) return true;
		for (const auto& sub : elem->fSubTrees) {
			if (std::binary_search(k.begin(), k.end(), sub.first)) {
				fringe.push(&sub.second);
			}
		}
	}
	return false;
}

template <typename Key, typename Value>
bool PrefixTree<Key,Value>::hasStrictSubset(const Key& k, const Value& v) const {
	std::queue<const PrefixNode*> fringe;
	fringe.push(&fRoot);
	int d = 0;
	while (!fringe.empty() && d < k.size()) {
		const PrefixNode* elem = fringe.front();
		fringe.pop();
		if (elem->fValue == v) return true;
		for (const auto& sub : elem->fSubTrees) {
			if (std::binary_search(k.begin(), k.end(), sub.first)) {
				fringe.push(&sub.second);
			}
		}
		d++;
	}
	return false;
}

template <typename Key, typename Value>
std::vector<Key> PrefixTree<Key,Value>::getSubsets(const Key& k, const Value& v) const {
    std::vector<Key> subs;
    std::queue<std::pair<const PrefixNode*, Key> > fringe;
    fringe.push(std::make_pair(&fRoot, Key()));
    while (!fringe.empty()) {
        const PrefixNode* elem = fringe.front().first;
        Key elemKey = fringe.front().second;
        fringe.pop();
        if (elem->fValue == v) {
            subs.push_back(elemKey);
        }
        for (const auto& sub : elem->fSubTrees) {
            if (std::binary_search(k.begin(), k.end(), sub.first)) {
                Key subKey(elemKey);
                subKey.push_back(sub.first);
                fringe.push(std::make_pair(&sub.second, subKey));
            }
        }
    }
    return subs;
}


template <typename Key, typename Value>
std::vector<const Value*> PrefixTree<Key,Value>::getSubsets(const Key& k) const {
    std::vector<const Value*> subs;
    std::queue<std::pair<const PrefixNode*, Key> > fringe;
    fringe.push(std::make_pair(&fRoot, Key()));
    while (!fringe.empty()) {
        const PrefixNode* elem = fringe.front().first;
        Key elemKey = fringe.front().second;
        fringe.pop();
        if (elem->fValue != Value()) {
            subs.push_back(&elem->fValue);
        }
        for (const auto& sub : elem->fSubTrees) {
            if (std::binary_search(k.begin(), k.end(), sub.first)) {
                Key subKey(elemKey);
                subKey.push_back(sub.first);
                fringe.push(std::make_pair(&sub.second, subKey));
            }
        }
    }
    return subs;
}

template <typename Key, typename Value>
std::vector<Key> PrefixTree<Key,Value>::getSets() const {
	std::vector<Key> subs;
	std::queue<std::pair<const PrefixNode*, Key> > fringe;
	fringe.push(std::make_pair(&fRoot, Key()));
	while (!fringe.empty()) {
		const PrefixNode* elem = fringe.front().first;
		Key elemKey = fringe.front().second;
		fringe.pop();
		if (elem->fValue != Value()) {
			subs.push_back(elemKey);
		}
		for (const auto& sub : elem->fSubTrees) {
			Key subKey(elemKey);
			subKey.push_back(sub.first);
			fringe.push(std::make_pair(&sub.second, subKey));
		}
	}
	return subs;
}

template <typename Key, typename Value>
std::map<Key, int> PrefixTree<Key,Value>::getSupports(const Key& k) const {
	std::map<Key, int> result;
	std::queue<std::pair<Key, const PrefixNode*> > fringe;
	fringe.push(std::make_pair(Key(), &fRoot));
	while (!fringe.empty()) {
		const PrefixNode* elem = fringe.front().second;
		const Key& pk = fringe.front().first;
		for (const auto& sub : elem->fSubTrees) {
			if (std::binary_search(k.begin(), k.end(), sub.first)) {
				Key newk(pk.begin(), pk.end());
				newk.push_back(sub.first);
				result[newk] = sub.second.fValue;
				fringe.push(std::make_pair(newk, &sub.second));
			}
		}
		fringe.pop();
	}
	return result;
}

#endif //UTIL_PREFIXTREE_H_