#ifndef DATA_DATABASE_H_
#define DATA_DATABASE_H_


#include <fstream>
#include "types.h"

class Database {
public:
	Database();
	unsigned size() const;
	unsigned nrAttrs() const;
    unsigned nrItems() const;

    const Transaction& getRow(int) const;
	void setRow(int i, const Transaction&);
	void write(std::ofstream& file, char delim=' ') const;
    void write(std::ofstream& file, const SimpleTidList& subset, char delim=' ') const;
	void writeInt(std::ofstream& file, char delim=' ') const;
	void sort();
	void toFront(const SimpleTidList&);

	int getAttrIndex(int t) const;
	int frequency(int i) const;
	const std::string& getValue(int i) const;
	const std::vector<int>& getDomainOfItem(int) const;
	const std::vector<int>& getDomain(int) const;
    const std::string& getAttrName(int) const;
	std::vector<int> getAttrVector(const Itemset&) const;
    std::vector<int> getAttrVectorItems(const Itemset&) const;
    int translateToken(int, const std::string&);
	int getAttr(const std::string&) const;
    int getItem(int, const std::string&) const;

    static std::vector<int> getDiffs(const Database&, const Database&);
	std::unordered_map<std::pair<int, std::string>, int, pairhash> getDictionary() const;
	void setDictionary(const std::unordered_map<std::pair<int, std::string>, int, pairhash>&);
protected:
    void setAttributes(const std::vector<std::string>&);
    void incFreq(int i);
    void addRow(const Transaction&);

private:
	struct AttributeInfo {
		AttributeInfo(const std::string& name):fName(name) { }
		std::string fName;
		std::vector<int> fValues;
	};
	struct ItemInfo {
		ItemInfo() { }
		ItemInfo(const std::string& val, int attr):fValue(val), fAttribute(attr), fFrequency(0) { }
		std::string fValue;
		int fAttribute;
		int fFrequency;
	};
	std::vector<Transaction> fData; //array of data represented as rows of integers
	std::unordered_map<std::pair<int, std::string>, int, pairhash> fItemDictionary; // maps a value string for the given attribute index to corresponding item id
	std::vector<AttributeInfo> fAttributes;
	std::vector<ItemInfo> fItems;
	int fNumTokens;
};

#endif // DATA_DATABASE_H_
