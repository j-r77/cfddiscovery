#include <fstream>
#include <iostream>
#include <cstring>
#include <chrono>
#include "util/output.h"
#include "util/stringutil.h"
#include "data/databasereader.h"
#include "algorithms/cfddiscovery.h"

int main(int argc, char *argv[]) {
    if (argc == 5 || argc == 6) {
        std::ifstream dataFile(argv[1]);
        if (!dataFile.good()) {
            std::cout << "[ERROR] File not found: " << argv[1] << std::endl;
            return -1;
        }
        Database db = DatabaseReader::fromTable(dataFile, ',');

        int supp = atoi(argv[2]);
        if (supp < 1 || supp > db.size()) {
            std::cout << "[ERROR] Illegal Support value: " << argv[2] << " not in [1," << db.size() << "]" << std::endl;
            return -1;
        }
        
        double conf = atof(argv[3]);
        if (conf < 0 || conf >= 1) {
            std::cout << "[ERROR] Illegal Confidence value: " << argv[3] << " not in [0,1[" << std::endl;
            return -1;
        }
    
        int maxSize = atoi(argv[4])+1;
        if (maxSize < 2) {
            std::cout << "[ERROR] Illegal Max Size value: " << argv[4] << " less than 1" << std::endl;
            return -1;
        }
        
        auto t1 = std::chrono::high_resolution_clock::now();
        CFDDiscovery cfdd(db);
        if (argc == 5) {
            cfdd.fdsFirstDFS(supp, maxSize, SUBSTRATEGY::DFS, conf);
        }
        else if (strcmp(argv[5], "FD-First-DFS-dfs") == 0) {
            cfdd.fdsFirstDFS(supp, maxSize, SUBSTRATEGY::DFS, conf);
        }
        else if (strcmp(argv[5], "FD-First-DFS-bfs") == 0) {
            cfdd.fdsFirstDFS(supp, maxSize, SUBSTRATEGY::BFS, conf);
        }
        else if (strcmp(argv[5], "FD-First-BFS-dfs") == 0) {
            cfdd.fdsFirstBFS(supp, maxSize, SUBSTRATEGY::DFS, conf);
        }
        else if (strcmp(argv[5], "FD-First-BFS-bfs") == 0) {
            cfdd.fdsFirstBFS(supp, maxSize, SUBSTRATEGY::BFS, conf);
        }
        else if (strcmp(argv[5], "Itemset-First-DFS-dfs") == 0) {
            cfdd.itemsetsFirstDFS(supp, maxSize, SUBSTRATEGY::DFS, conf);
        }
        else if (strcmp(argv[5], "Itemset-First-DFS-bfs") == 0) {
            cfdd.itemsetsFirstDFS(supp, maxSize, SUBSTRATEGY::BFS, conf);
        }
        else if (strcmp(argv[5], "Itemset-First-BFS-dfs") == 0) {
            cfdd.itemsetsFirstBFS(supp, maxSize, SUBSTRATEGY::DFS, conf);
        }
        else if (strcmp(argv[5], "Itemset-First-BFS-bfs") == 0) {
            cfdd.itemsetsFirstBFS(supp, maxSize, SUBSTRATEGY::BFS, conf);
        }
        else if (strcmp(argv[5], "Integrated-DFS") == 0) {
            cfdd.integratedDFS(supp, maxSize, conf);
        }
        else if (strcmp(argv[5], "Integrated-BFS") == 0) {
            cfdd.ctane(supp, maxSize, conf);
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        CFDList cfds = cfdd.getCFDs();
        for (const auto& c : cfds) {
            Output::printCFD(c, db);
        }
        std::cout << "Mined " << cfds.size() << " cfds in " << time << " milliseconds" << std::endl;
    }
    else {
        std::cout << "[ERROR] Wrong number of Arguments: " << (argc-1) << " instead of 4 [5]" << std::endl;
        std::cout << "Usage: ./CFDD csvfile minsupport minconfidence maxsize [algorithm]" << std::endl;
        std::cout << "\t For Example: ./CFDD datasets/adult.csv 1000 0.95 5" << std::endl;
        return -1;
    }
    return 0;
}
