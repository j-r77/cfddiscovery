#ifndef DATA_DATABASEREADER_H_
#define DATA_DATABASEREADER_H_

#include "database.h"

class DatabaseReader : private Database {
public:
	static Database fromTable(const std::string&, char delim=' ', double rSample=1, double cSample=1, bool headers=true);
	static Database fromTable(const std::string&, const Database&, char delim=' ', bool headers=true);
    static Database fromTable(std::ifstream&, char delim=' ', double rSample=1, double cSample=1, bool headers=true);
	static Database fromTable(std::ifstream&, const Database&, char delim=' ', bool headers=true);
	static Database fromTable(std::ifstream&, int nrCols, char delim=' ', bool headers=true);
	static Database fromTablePart(std::ifstream&, int nrCols, int nrTuples, char delim=' ', bool headers=true);
	static Database fromItemsets(std::ifstream&, char delim=' ');
};

#endif // DATA_DATABASEREADER_H_