#include "databasereader.h"
#include "../util/output.h"

#include <sstream>
#include <random>

Database DatabaseReader::fromItemsets(std::ifstream& infile, char delim) {
	DatabaseReader db = DatabaseReader();
	std::vector<std::string> heads;
	heads.push_back("Item");
    db.setAttributes(heads);

	std::string line;
	std::getline(infile, line);
	bool eof = false;
	while (!eof)
	{
	    std::istringstream iss(line);
	    std::string val;
        std::vector<int> row;
    	int i = 0;
	    while (std::getline(iss, val, delim)) {
            row.push_back(db.translateToken(1,val));
    		i++;
	    }
	    if (i > 0) {
	    	db.addRow(row);
	    }
	    if (infile.eof()) {
	    	eof = true;
	    }
    	std::getline(infile, line);
	}
	return db;
}

Database DatabaseReader::fromTable(const std::string& infile, char delim, double rSample, double cSample, bool headers) {
	std::ifstream stream(infile);
	return DatabaseReader::fromTable(stream, delim, rSample, cSample, headers);
}

Database DatabaseReader::fromTable(const std::string& infile, const Database& orig, char delim, bool headers) {
	std::ifstream stream(infile);
	return DatabaseReader::fromTable(stream, orig, delim, headers);
}

Database DatabaseReader::fromTable(std::ifstream& infile, const Database& orig, char delim, bool headers) {
	DatabaseReader db = DatabaseReader();
	std::string line;
	std::vector<std::string> heads;
	std::getline(infile, line);
	std::istringstream iss(line);
	std::string val;
	int i = 1;
	while (std::getline(iss, val, delim)) {
		if (headers)
			heads.push_back(val);
		else
            heads.push_back("Attr" + std::to_string(i));
		i++;
	}
	db.setAttributes(heads);
	db.setDictionary(orig.getDictionary());
	if (headers)
		std::getline(infile, line);
	int size = i - 1;
	bool eof = false;
	while (!eof)
	{
	    std::istringstream iss(line);
	    std::string val;
        std::vector<int> row(size);
    	int i = 0;
	    while (std::getline(iss, val, delim)) {
	    	int it = db.translateToken(i, val);
	    	db.incFreq(it);
            row[i] = it;
    		i++;
	    }
	    if (i == heads.size()) {
            db.addRow(row);
	    }
	    if (infile.eof()) {
	    	eof = true;
	    }
    	std::getline(infile, line);

	}
	return db;
}

Database DatabaseReader::fromTable(std::ifstream& infile, int nrCols, char delim, bool headers) {
	DatabaseReader db = DatabaseReader();
	std::string line;
	std::vector<std::string> heads;
	std::getline(infile, line);
	std::istringstream iss(line);
	std::string val;
	int i = 1;
	while (std::getline(iss, val, delim) && heads.size() < nrCols) {
		if (headers)
			heads.push_back(val);
		else
			heads.push_back("Attr" + std::to_string(i));
		i++;
	}
	db.setAttributes(heads);
	if (headers)
		std::getline(infile, line);
	int size = i - 1;
	bool eof = false;
	while (!eof)
	{
		std::istringstream iss(line);
		std::string val;
		std::vector<int> row(size);
		int i = 0;
		while (std::getline(iss, val, delim) && i < nrCols) {
			int it = db.translateToken(i, val);
			db.incFreq(it);
			row[i] = it;
			i++;
		}
		if (i == heads.size()) {
			db.addRow(row);
		}
		if (infile.eof()) {
			eof = true;
		}
		std::getline(infile, line);

	}
	return db;
}

Database DatabaseReader::fromTablePart(std::ifstream& infile, int nrCols, int nrTuples, char delim, bool headers) {
	DatabaseReader db = DatabaseReader();
	std::string line;
	std::vector<std::string> heads;
	std::getline(infile, line);
	std::istringstream iss(line);
	std::string val;
	int i = 1;
	while (std::getline(iss, val, delim) && heads.size() < nrCols) {
		if (headers)
			heads.push_back(val);
		else
			heads.push_back("Attr" + std::to_string(i));
		i++;
	}
	db.setAttributes(heads);
	if (headers)
		std::getline(infile, line);
	int size = i - 1;
	bool eof = false;
	while (!eof && db.size() < nrTuples)
	{
		std::istringstream iss(line);
		std::string val;
		std::vector<int> row(size);
		int i = 0;
		while (std::getline(iss, val, delim) && i < nrCols) {
			int it = db.translateToken(i, val);
			db.incFreq(it);
			row[i] = it;
			i++;
		}
		if (i == heads.size()) {
			db.addRow(row);
		}
		if (infile.eof()) {
			eof = true;
		}
		std::getline(infile, line);

	}
	return db;
}


Database DatabaseReader::fromTable(std::ifstream& infile, char delim, double rSample, double cSample, bool headers) {
	DatabaseReader db = DatabaseReader();
	std::string line;
	std::vector<std::string> allHeads;
	std::getline(infile, line);
	std::istringstream iss(line);
	std::string val;
    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_real_distribution<double> uni(0.0,1.0); // guaranteed unbiased
    std::vector<int> sColumns;
	int i = 1;
	while (std::getline(iss, val, delim)) {
        if (headers)
            allHeads.push_back(val);
        else
            allHeads.push_back("Attr" + std::to_string(i));
        sColumns.push_back(i-1);
        i++;
	}
    int size = sColumns.size() * cSample;
    shuffle(sColumns);
    sColumns = std::vector<int>(sColumns.begin(), sColumns.begin() + size);
    std::sort(sColumns.begin(), sColumns.end());
    std::vector<std::string> heads;
    for (int cc : sColumns) {
        heads.push_back(allHeads[cc]);
    }
	db.setAttributes(heads);
	if (headers)
		std::getline(infile, line);
	bool eof = false;
	while (!eof)
	{
        if (uni(rng) < rSample) {
            std::istringstream iss(line);
            std::string val;
            std::vector<int> row(size);
            int i = 0;
            int j = 0;
            while (std::getline(iss, val, delim)) {
                if (std::binary_search(sColumns.begin(), sColumns.end(), i)) {
                    int it = db.translateToken(j, val);
                    db.incFreq(it);
                    row[j] = it;
                    j++;
                }
                i++;
            }
            if (j > 0) {
                db.addRow(row);
            }
            if (infile.eof()) {
                eof = true;
            }
        }
        std::getline(infile, line);
    }
	return db;
}
