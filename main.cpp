#include <fstream>
#include <iostream>
#include <chrono>
#include "util/output.h"
#include "data/databasereader.h"
#include "algorithms/ctane.h"

int NR_DATASETS = 4;
const char *DATASETS[] = {
        "abalone",//0
        "sp500",//1
        "Soccer",//2
        "adult",//3
        "Tax",//4
        "Hospital",//5
        "Employees",//6
        "GenLarge",//7
        "cancer",//8
        "iris",//9
        "poker",//10
        "mushroom",//11
        "nursery",//12
        "ipums",//13
};

const int ABALONE = 0;
const int SP500 = 1;
const int SOCCER = 2;
const int ADULT = 3;
const int CANCER = 8;
const int IRIS = 9;
const int POKER = 10;
const int MUSHROOM = 11;
const int NURSERY = 12;
const int IPUMS = 13;

const char *EXE_PATH = "";
//const char *DATA_PATH = "/home/jr/phd/Code/CFDDiscovery/datasets/";
const char *DIRTY_PATH = "dirtysets/";
const char *DIRTY_PATH_FALCON = "falconsets/";

const char *DATA_PATH = "/Users/jr/Documents/PhD Git/Code/CFDDiscovery/datasets/";


void expBfsDfs() {
    int dbs[] = {ADULT};
    std::map<int, std::vector<double> > minsups;
    minsups[NURSERY] = {0.1, 0.5, 1, 5, 10, 15};
    //minsups[ADULT] = {0.1, 0.5, 1, 5, 10, 15};
    //minsups[ADULT] = {5, 10, 15};
    minsups[ADULT] = {0.1};
    minsups[MUSHROOM] = {0.1, 0.5, 1, 5, 10, 15};
    //minsups[NURSERY] = {5};
    //minsups[ADULT] = {5};
    //minsups[MUSHROOM] = {5};
    for (int dbIx : dbs) {
        std::string outpath = concat(3, EXE_PATH, DATASETS[dbIx], "_bfsdfs.csv");
        std::ofstream ofile(outpath.c_str());
        ofile << "algorithm,strategy,minsup,minsupPct,nrCfds,runtime" << std::endl;
        std::ifstream dbFile(concat(3, DATA_PATH, DATASETS[dbIx], ".csv"));
        Database db = DatabaseReader::fromTable(dbFile, 11, ',');
        int onePct = db.size() / 100.0;
        int maxSize = 7;
        for (double msP:minsups[dbIx]) {
            CFDList c1;
            CFDList c3;
            int ms = msP * onePct;
            std::cout << "DATASET " << DATASETS[dbIx] << " minsup " << ms << std::endl;
            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CTane ct(db);
                ct.mine(ms, maxSize);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                ofile << "CTane,BFS," << ms << "," << msP << "," << ct.nrCFDs() << "," << iftime << std::endl;
                //std::cout << "CTane - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }

            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CTane ct(db);
                ct.mineFreeDepth(ms, maxSize);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                ofile << "CTane,DFS," << ms << "," << msP << "," << ct.nrCFDs() << "," << iftime << std::endl;
                //std::cout << "CTane - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }
/*
            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CTane ct(db);
                ct.mineItemsetsFirst(ms, maxSize, SUBSTRATEGY::BFS);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                ofile << "ItemsetsFirst,BFSbfs," << ms << "," << msP << "," << ct.nrCFDs() << "," << iftime << std::endl;
                //std::cout << "ItemsetsFirst - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }

            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CTane ct(db);
                ct.mineItemsetsFirst(ms, maxSize, SUBSTRATEGY::DFS);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                ofile << "ItemsetsFirst,BFSdfs," << ms << "," << msP << "," << ct.nrCFDs() << "," << iftime << std::endl;
                //std::cout << "ItemsetsFirst - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }

            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CTane ct(db);
                ct.mineItemsetsFirstDepth(ms, maxSize, SUBSTRATEGY::BFS);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                ofile << "ItemsetsFirst,DFSbfs," << ms << "," << msP << "," << ct.nrCFDs() << "," << iftime << std::endl;
                //std::cout << "ItemsetsFirst - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }

            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CTane ct(db);
                ct.mineItemsetsFirstDepth(ms, maxSize, SUBSTRATEGY::DFS);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                ofile << "ItemsetsFirst,DFSdfs," << ms << "," << msP << "," << ct.nrCFDs() << "," << iftime << std::endl;
                //std::cout << "ItemsetsFirst - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }

            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CTane ct(db);
                ct.mineFDsFirst(ms, maxSize, SUBSTRATEGY::BFS);
                c1 = ct.getCFDs();
                auto t2 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                ofile << "FDsFirst,BFSbfs," << ms << "," << msP << "," << ct.nrCFDs() << "," << iftime << std::endl;
                //std::cout << "FDsFirst - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }

            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CTane ct(db);
                ct.mineFDsFirst(ms, maxSize, SUBSTRATEGY::DFS);
                c3 = ct.getCFDs();
                auto t2 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                ofile << "FDsFirst,BFSdfs," << ms << "," << msP << "," << ct.nrCFDs() << "," << iftime << std::endl;
                //std::cout << "FDsFirst - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }

            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CTane ct(db);
                ct.mineFDsFirstDepth(ms, maxSize, SUBSTRATEGY::BFS);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                ofile << "FDsFirst,DFSbfs," << ms << "," << msP << "," << ct.nrCFDs() << "," << iftime << std::endl;
                //std::cout << "FDsFirst - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }

            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CTane ct(db);
                ct.mineFDsFirstDepth(ms, maxSize, SUBSTRATEGY::DFS);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                ofile << "FDsFirst,DFSdfs," << ms << "," << msP << "," << ct.nrCFDs() << "," << iftime << std::endl;
                //std::cout << "FDsFirst - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }
            /*std::sort(c1.begin(), c1.end());
            std::sort(c3.begin(), c3.end());
            if (c1 != c3) {
                std::cout << "FAIL" << std::endl;
                auto diff = setdiff(c3, c1);
                std::cout << "MISSING: " << std::endl;
                for (const auto& cfd : diff) {
                    Output::printCFD(cfd, db);
                    Output::printCollection(cfd.first);
                    std::cout << " => ";
                    Output::printCollection(itemset(cfd.second));
                    std::cout << std::endl;
                }
                std::cout << "MISSING: " << std::endl;
                for (const auto& cfd : setdiff(c1, c3)) {
                    Output::printCFD(cfd, db);
                    Output::printCollection(cfd.first);
                    std::cout << " => ";
                    Output::printCollection(itemset(cfd.second));
                    std::cout << std::endl;
                    Output::printCollection(db.getAttrVectorItems(cfd.first));
                    std::cout << " => ";
                    Output::printCollection(db.getAttrVectorItems(itemset(cfd.second)));
                    std::cout << std::endl;
                }
            }*/
        }
        ofile.close();
    }
}

void expSupport() {
    int dbs[] = {MUSHROOM, NURSERY, ADULT};
    std::map<int, std::vector<double> > minsups;
    minsups[NURSERY] = {0.1, 0.5, 1, 5, 10, 15};
    minsups[ADULT] = {0.1, 0.5, 1, 5, 10, 15};
    minsups[MUSHROOM] = {0.1, 0.5, 1, 5, 10, 15};
    for (int dbIx : dbs) {
        std::string outpath = concat(3, EXE_PATH, DATASETS[dbIx], "_runtime.csv");
        std::ofstream ofile(outpath.c_str());
        ofile << "algorithm,minsup,minsupPct,nrCfds,runtime" << std::endl;
        std::ifstream dbFile(concat(3, DATA_PATH, DATASETS[dbIx], ".csv"));
        Database db = DatabaseReader::fromTable(dbFile, 11, ',');
        int onePct = db.size() / 100.0;
        int maxSize = 7;
        for (double msP:minsups[dbIx]) {
            int ms = msP * onePct;
            std::cout << "DATASET " << DATASETS[dbIx] << " minsup " << ms << std::endl;
            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CTane ct(db);
                ct.mineFreeDepth(ms, maxSize);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                ofile << "CTane," << ms << "," << msP << "," << ct.nrCFDs() << "," << time << std::endl;
                std::cout << "CTane - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << time << " ms" << std::endl;
            }

            {
                auto t3 = std::chrono::high_resolution_clock::now();
                CTane ct2(db);
                ct2.mineItemsetsFirst(ms, maxSize, SUBSTRATEGY::DFS);
                auto t4 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
                ofile << "ItemsetsFirst," << ms << "," << msP << "," << ct2.nrCFDs() << "," << iftime << std::endl;
                std::cout << "ItemsetsFirst - Minsup " << ms << ": " << ct2.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }

            {
                auto t5 = std::chrono::high_resolution_clock::now();
                CTane ct3(db);
                ct3.mineFDsFirstDepth(ms, maxSize, SUBSTRATEGY::DFS);
                auto t6 = std::chrono::high_resolution_clock::now();
                auto fftime = std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count();
                ofile << "FDsFirst," << ms << "," << msP << "," << ct3.nrCFDs() << "," << fftime << std::endl;
                std::cout << "FDsFirst - Minsup " << ms << ": " << ct3.nrCFDs() << " cfds in " << fftime << " ms"
                          << std::endl;
            }
        }
        ofile.close();
    }
}

void expAttributes() {
    int dbs[] = {MUSHROOM, NURSERY, ADULT};
    std::map<int, std::vector<int> > attSizes;
    attSizes[NURSERY] = {5,7,9};
    attSizes[ADULT] = {7,9,11};
    attSizes[MUSHROOM] = {7,9,11,13,15};
    for (int dbIx : dbs) {
        std::string outpath = concat(3, EXE_PATH, DATASETS[dbIx], "_attributes.csv");
        std::ofstream ofile(outpath.c_str());
        ofile << "algorithm,minsup,atts,nrCfds,runtime" << std::endl;
        auto atts = attSizes[dbIx];
        for (int as:atts) {
            std::ifstream dbFile(concat(3, DATA_PATH, DATASETS[dbIx], ".csv"));
            Database db = DatabaseReader::fromTable(dbFile, as, ',');
            int onePct = db.size() / 100.0;
            int maxSize = 7;
            int minsups[] = {10 * onePct};
            for (int ms:minsups) {
                std::cout << "DATASET " << DATASETS[dbIx] << " minsup " << ms << std::endl;

                {
                    auto t1 = std::chrono::high_resolution_clock::now();
                    CTane ct(db);
                    ct.mineFreeDepth(ms, maxSize);
                    auto t2 = std::chrono::high_resolution_clock::now();
                    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                    ofile << "CTane," << ms << "," << as << "," << ct.nrCFDs() << "," << time << std::endl;
                }

                {
                    auto t3 = std::chrono::high_resolution_clock::now();
                    CTane ct2(db);
                    ct2.mineItemsetsFirst(ms, maxSize, SUBSTRATEGY::DFS);
                    auto t4 = std::chrono::high_resolution_clock::now();
                    auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
                    ofile << "ItemsetsFirst," << ms << "," << as << "," << ct2.nrCFDs() << "," << iftime << std::endl;
                }

                {
                    auto t5 = std::chrono::high_resolution_clock::now();
                    CTane ct3(db);
                    ct3.mineFDsFirstDepth(ms, maxSize, SUBSTRATEGY::DFS);
                    auto t6 = std::chrono::high_resolution_clock::now();
                    auto fftime = std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count();
                    ofile << "FDsFirst," << ms << "," << as << "," << ct3.nrCFDs() << "," << fftime << std::endl;
                }
            }
        }
        ofile.close();
    }
}

void expTuples() {
    int dbs[] = {MUSHROOM};
    std::map<int,int> dbSizes;
    dbSizes[NURSERY] = 12960;
    dbSizes[ADULT] = 48842;
    dbSizes[MUSHROOM] = 8124;
    for (int dbIx : dbs) {
        std::string outpath = concat(3, EXE_PATH, DATASETS[dbIx], "_tuples.csv");
        std::ofstream ofile(outpath.c_str());
        ofile << "algorithm,minsup,tuples,nrCfds,runtime" << std::endl;
        int aPct = dbSizes[dbIx] / 100.0;
        int tupPcts[] = {10, 25, 50, 75, 100};
        for (int tupPct:tupPcts) {
            int nrTups = aPct * tupPct;
            int maxSize = 7;
            int nrAtts = 11;
            std::ifstream dbFile(concat(3, DATA_PATH, DATASETS[dbIx], ".csv"));
            Database db = DatabaseReader::fromTablePart(dbFile, nrAtts, nrTups, ',');
            int onePct = db.size() / 100.0;
            int minsups[] = {10 * onePct};
            for (int ms:minsups) {
                std::cout << "DATASET " << DATASETS[dbIx] << " minsup " << ms << " size " << db.size() << std::endl;

                {
                    auto t1 = std::chrono::high_resolution_clock::now();
                    CTane ct(db);
                    ct.mineFreeDepth(ms, maxSize);
                    auto t2 = std::chrono::high_resolution_clock::now();
                    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                    ofile << "CTane," << ms << "," << tupPct << "," << ct.nrCFDs() << "," << time << std::endl;
                }

                {
                    auto t3 = std::chrono::high_resolution_clock::now();
                    CTane ct2(db);
                    ct2.mineItemsetsFirst(ms, maxSize, SUBSTRATEGY::DFS);
                    auto t4 = std::chrono::high_resolution_clock::now();
                    auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
                    ofile << "ItemsetsFirst," << ms << "," << tupPct << "," << ct2.nrCFDs() << "," << iftime
                          << std::endl;
                }

                {
                    auto t5 = std::chrono::high_resolution_clock::now();
                    CTane ct3(db);
                    ct3.mineFDsFirstDepth(ms, maxSize, SUBSTRATEGY::DFS);
                    auto t6 = std::chrono::high_resolution_clock::now();
                    auto fftime = std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count();
                    ofile << "FDsFirst," << ms << "," << tupPct << "," << ct3.nrCFDs() << "," << fftime << std::endl;
                }
            }
        }
        ofile.close();
    }
}

void expMaxSize() {
    int dbs[] = {MUSHROOM, NURSERY, ADULT};
    std::map<int, std::vector<double> > minsups;
    minsups[NURSERY] = {10};
    minsups[ADULT] = {10};
    minsups[MUSHROOM] = {10};
    std::map<int, std::vector<int> > maxSizes;
    maxSizes[NURSERY] = {3,5,7,9};
    maxSizes[ADULT] = {3,5,7,9};
    maxSizes[MUSHROOM] = {5,7,9,11};
    for (int dbIx : dbs) {
        std::string outpath = concat(3, EXE_PATH, DATASETS[dbIx], "_maxSize.csv");
        std::ofstream ofile(outpath.c_str());
        ofile << "algorithm,minsup,minsupPct,maxSize,nrCfds,runtime" << std::endl;
        std::ifstream dbFile(concat(3, DATA_PATH, DATASETS[dbIx], ".csv"));
        Database db = DatabaseReader::fromTable(dbFile, 11, ',');
        int onePct = db.size() / 100.0;
        for (int maxSize : maxSizes[dbIx]) {
            for (double msP:minsups[dbIx]) {
                int ms = msP * onePct;
                std::cout << "DATASET " << DATASETS[dbIx] << " minsup " << ms << std::endl;

                {
                    auto t1 = std::chrono::high_resolution_clock::now();
                    CTane ct(db);
                    ct.mineFreeDepth(ms, maxSize);
                    auto t2 = std::chrono::high_resolution_clock::now();
                    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                    ofile << "CTane," << ms << "," << msP << "," << maxSize << "," << ct.nrCFDs() << "," << time
                          << std::endl;
                    std::cout << "CTane - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << time << " ms"
                              << std::endl;
                }

                {
                    auto t3 = std::chrono::high_resolution_clock::now();
                    CTane ct2(db);
                    ct2.mineItemsetsFirst(ms, maxSize, SUBSTRATEGY::DFS);
                    auto t4 = std::chrono::high_resolution_clock::now();
                    auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
                    ofile << "ItemsetsFirst," << ms << "," << msP << "," << maxSize << "," << ct2.nrCFDs() << ","
                          << iftime << std::endl;
                    std::cout << "ItemsetsFirst - Minsup " << ms << ": " << ct2.nrCFDs() << " cfds in " << iftime
                              << " ms" << std::endl;
                }

                {
                    auto t5 = std::chrono::high_resolution_clock::now();
                    CTane ct3(db);
                    ct3.mineFDsFirstDepth(ms, maxSize, SUBSTRATEGY::DFS);
                    auto t6 = std::chrono::high_resolution_clock::now();
                    auto fftime = std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count();
                    ofile << "FDsFirst," << ms << "," << msP << "," << maxSize << "," << ct3.nrCFDs() << "," << fftime
                          << std::endl;
                    std::cout << "FDsFirst - Minsup " << ms << ": " << ct3.nrCFDs() << " cfds in " << fftime << " ms"
                              << std::endl;
                }
            }
        }
        ofile.close();
    }
}

void expConfidence() {
    int dbs[] = {MUSHROOM, NURSERY, ADULT};
    std::map<int, std::vector<double> > minsups;
    minsups[NURSERY] = {10};
    minsups[ADULT] = {10};
    minsups[MUSHROOM] = {10};
    for (int dbIx : dbs) {
        std::string outpath = concat(3, EXE_PATH, DATASETS[dbIx], "_confidence.csv");
        std::ofstream ofile(outpath.c_str());
        ofile << "algorithm,minsup,minsupPct,conf,nrCfds,runtime" << std::endl;
        std::ifstream dbFile(concat(3, DATA_PATH, DATASETS[dbIx], ".csv"));
        Database db = DatabaseReader::fromTable(dbFile, 11, ',');
        int onePct = db.size() / 100.0;
        int maxSize = 7;
        double confs[] = {0.6, 0.75, 0.90, 0.95, 1};
        for (double conf : confs) {
            for (double msP:minsups[dbIx]) {
                int ms = msP * onePct;
                std::cout << "DATASET " << DATASETS[dbIx] << " minsup " << ms << std::endl;

                {
                    auto t1 = std::chrono::high_resolution_clock::now();
                    CTane ct(db);
                    ct.mineFreeDepth(ms, maxSize, conf);
                    auto t2 = std::chrono::high_resolution_clock::now();
                    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                    ofile << "CTane," << ms << "," << msP << "," << conf << "," << ct.nrCFDs() << "," << time
                          << std::endl;
                    std::cout << "CTane - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << time << " ms" << std::endl;
                }

                {
                    auto t3 = std::chrono::high_resolution_clock::now();
                    CTane ct2(db);
                    ct2.mineItemsetsFirst(ms, maxSize, SUBSTRATEGY::DFS, conf);
                    auto t4 = std::chrono::high_resolution_clock::now();
                    auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
                    ofile << "ItemsetsFirst," << ms << "," << msP << "," << conf << "," << ct2.nrCFDs() << "," << iftime
                          << std::endl;
                    std::cout << "ItemsetsFirst - Minsup " << ms << ": " << ct2.nrCFDs() << " cfds in " << iftime
                              << " ms" << std::endl;
                }

                {
                    auto t5 = std::chrono::high_resolution_clock::now();
                    CTane ct3(db);
                    ct3.mineFDsFirstDepth(ms, maxSize, SUBSTRATEGY::DFS, conf);
                    auto t6 = std::chrono::high_resolution_clock::now();
                    auto fftime = std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count();
                    ofile << "FDsFirst," << ms << "," << msP << "," << conf << "," << ct3.nrCFDs() << "," << fftime
                          << std::endl;
                    std::cout << "FDsFirst - Minsup " << ms << ": " << ct3.nrCFDs() << " cfds in " << fftime << " ms"
                              << std::endl;
                }
            }
        }
        ofile.close();
    }
}

void expNrCFD() {
    int dbs[] = {MUSHROOM, NURSERY, ADULT};
    std::map<int, std::vector<double> > minsups;
    minsups[NURSERY] = {5,10,15};
    minsups[ADULT] = {5,10,15};
    minsups[MUSHROOM] = {5,10,15};
    for (int dbIx : dbs) {
        std::string outpath = concat(3, EXE_PATH, DATASETS[dbIx], "_nrcfd.csv");
        std::ofstream ofile(outpath.c_str());
        ofile << "algorithm,minsup,minsupPct,conf,nrCfds,runtime" << std::endl;
        std::ifstream dbFile(concat(3, DATA_PATH, DATASETS[dbIx], ".csv"));
        Database db = DatabaseReader::fromTable(dbFile, 15, ',');
        int onePct = db.size() / 100.0;
        int maxSize = 7;
        double confs[] = {0.95, 0.99, 1};
        for (double conf : confs) {
            for (double msP:minsups[dbIx]) {
                int ms = msP * onePct;
                std::cout << "DATASET " << DATASETS[dbIx] << " minsup " << ms << std::endl;

                {
                    auto t5 = std::chrono::high_resolution_clock::now();
                    CTane ct3(db);
                    ct3.mineFDsFirstDepth(ms, maxSize, SUBSTRATEGY::DFS, conf);
                    auto t6 = std::chrono::high_resolution_clock::now();
                    auto fftime = std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count();
                    ofile << "FDsFirst," << ms << "," << msP << "," << conf << "," << ct3.nrCFDs() << "," << fftime
                          << std::endl;
                    std::cout << "FDsFirst - Minsup " << ms << ": " << ct3.nrCFDs() << " cfds in " << fftime << " ms"
                              << std::endl;
                }
            }
        }
        ofile.close();
    }
}

/*void expPruning() {
    int dbs[] = {MUSHROOM, NURSERY, SP500};

    std::map<int, std::vector<int> > attSizes;
    attSizes[NURSERY] = {9};
    attSizes[SP500] = {7};
    attSizes[MUSHROOM] = {12};//,14,15};
    for (int dbIx : dbs) {
        std::string outpath = concat(3, EXE_PATH, DATASETS[dbIx], "_pruning.csv");
        std::ofstream ofile(outpath.c_str());
        ofile << "algorithm,minsup,atts,nrCfds,runtime,visited,validated" << std::endl;

        auto atts = attSizes[dbIx];
        for (int as:atts) {
            std::ifstream dbFile(concat(3, DATA_PATH, DATASETS[dbIx], ".csv"));
            Database db = DatabaseReader::fromTable(dbFile, as, ',');
            int onePct = db.size() / 100.0;
            int minsups[] = {1 * onePct, 5 * onePct, 10 * onePct, 25 * onePct};
            //int minsups[] = {1 * onePct};
            for (int ms:minsups) {
                std::cout << "DATASET " << DATASETS[dbIx] << " minsup " << ms << std::endl;
                auto t1 = std::chrono::high_resolution_clock::now();
                CTane ct(db);
                ct.mine(ms);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                ofile << "bfs," << ms << "," << as << "," << ct.nrCFDs() << "," << time << "," << ct.fVisited << "," << ct.fValidated << std::endl;
                //std::cout << "BFS " << ms << ": " << ct.nrCFDs() << " cfds in " << time << " ms" << std::endl;
                //std::cout << "BFS " << ms << ": " << ct.fVisited << " visited and " << ct.fValidated << " validated" << std::endl;

                auto t3 = std::chrono::high_resolution_clock::now();
                CTane ct2(db);
                ct2.mineFreeDepth(ms);
                auto t4 = std::chrono::high_resolution_clock::now();
                auto ftime = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
                ofile << "dfs," << ms << "," << as << "," << ct2.nrCFDs() << "," << ftime << "," << ct2.fVisited << "," << ct2.fValidated << std::endl;
                //std::cout << "DFS " << ms << ": " << ct2.nrCFDs() << " cfds in " << ftime << " ms" << std::endl;
                //std::cout << "DFS " << ms << ": " << ct2.fVisited << " visited and " << ct2.fValidated << " validated" << std::endl;

                auto c1 = ct.getCFDs();
                std::sort(c1.begin(), c1.end());
                auto c3 = ct2.getCFDs();
                std::sort(c3.begin(), c3.end());
                if (c1 != c3) {
                    std::cout << "FAIL" << std::endl;
                    auto diff = setdiff(c3, c1);
                    Output::printCFDList(c1, db);
                    std::cout << "----" << std::endl;
                    Output::printCFDList(c3, db);
                    std::cout << "MISSING: " << std::endl;
                    for (const auto& cfd : diff) {
                        Output::printCFD(cfd, db);
                        Output::printCollection(join(cfd.first,cfd.second));
                        std::cout << std::endl;
                    }
                    std::cout << "MISSING: " << std::endl;
                    for (const auto& cfd : setdiff(c1, c3)) {
                        Output::printCFD(cfd, db);
                        Output::printCollection(join(cfd.first,cfd.second));
                        std::cout << std::endl;
                    }
                }
            }
        }
        ofile.close();
    }
}*/

int main(int argc, char *argv[]) {
    //expSupport();
    //expAttributes();
    //expMaxSize();
    //expConfidence();
    // expPruning();
    //expTuples();
    expBfsDfs();
    //expNrCFD();
    return 0;
}
