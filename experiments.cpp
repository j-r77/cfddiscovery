#include <fstream>
#include <iostream>
#include <chrono>
#include "util/output.h"
#include "util/stringutil.h"
#include "data/databasereader.h"
#include "algorithms/cfddiscovery.h"

const char *DATASETS[] = {
        "adult",//0
        "mushroom",//1
        "nursery",//2
        "contraceptive",//3
};

const int ADULT = 0;
const int MUSHROOM = 1;
const int NURSERY = 2;
const int CONTRACEPTIVE = 3;

const int MAX_ATTRS = 15;

const char *EXE_PATH = "";
const char *DATA_PATH = "/datasets/";

void expBfsDfs() {
    int dbs[] = {ADULT, MUSHROOM, NURSERY};
    std::map<int, std::vector<double> > minsups;
    minsups[NURSERY] = {0.1, 0.5, 1, 5, 10, 15};
    minsups[ADULT] = {0.1, 0.5, 1, 5, 10, 15};
    minsups[MUSHROOM] = {0.1, 0.5, 1, 5, 10, 15};

    for (int dbIx : dbs) {
        std::ifstream dbFile(concat(3, DATA_PATH, DATASETS[dbIx], ".csv"));
        Database db = DatabaseReader::fromTable(dbFile, MAX_ATTRS, ',');
        int onePct = db.size() / 100.0;
        int maxSize = 7;
        for (double msP:minsups[dbIx]) {
            int ms = msP * onePct;
            std::cout << "DATASET " << DATASETS[dbIx] << " minsup " << ms << std::endl;
            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CFDDiscovery ct(db);
                ct.ctane(ms, maxSize);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                std::cout << "Integrated BFS (Ctane) - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }

            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CFDDiscovery ct(db);
                ct.integratedDFS(ms, maxSize);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                std::cout << "Integrated DFS - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }

            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CFDDiscovery ct(db);
                ct.itemsetsFirstBFS(ms, maxSize, SUBSTRATEGY::BFS);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                std::cout << "ItemsetsFirst BFSbfs - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }

            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CFDDiscovery ct(db);
                ct.itemsetsFirstBFS(ms, maxSize, SUBSTRATEGY::DFS);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                std::cout << "ItemsetsFirst BFSdfs - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }

            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CFDDiscovery ct(db);
                ct.itemsetsFirstDFS(ms, maxSize, SUBSTRATEGY::BFS);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                std::cout << "ItemsetsFirst DFSbfs - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }

            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CFDDiscovery ct(db);
                ct.itemsetsFirstDFS(ms, maxSize, SUBSTRATEGY::DFS);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                std::cout << "ItemsetsFirst DFSdfs - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }

            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CFDDiscovery ct(db);
                ct.fdsFirstBFS(ms, maxSize, SUBSTRATEGY::BFS);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                std::cout << "FDsFirst BFSbfs - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }

            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CFDDiscovery ct(db);
                ct.fdsFirstBFS(ms, maxSize, SUBSTRATEGY::DFS);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                std::cout << "FDsFirst BFSdfs - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }

            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CFDDiscovery ct(db);
                ct.fdsFirstDFS(ms, maxSize, SUBSTRATEGY::BFS);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                std::cout << "FDsFirst DFSbfs - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }

            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CFDDiscovery ct(db);
                ct.fdsFirstDFS(ms, maxSize, SUBSTRATEGY::DFS);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                std::cout << "FDsFirst DFSdfs - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }
        }
    }
}

void expSupport() {
    int dbs[] = {ADULT, MUSHROOM, NURSERY};
    std::map<int, std::vector<double> > minsups;
    minsups[NURSERY] = {0.1, 0.5, 1, 5, 10, 15};
    minsups[ADULT] = {0.1, 0.5, 1, 5, 10, 15};
    minsups[MUSHROOM] = {0.1, 0.5, 1, 5, 10, 15};

    for (int dbIx : dbs) {
        std::ifstream dbFile(concat(3, DATA_PATH, DATASETS[dbIx], ".csv"));
        Database db = DatabaseReader::fromTable(dbFile, MAX_ATTRS, ',');
        int onePct = db.size() / 100.0;
        int maxSize = 7;
        for (double msP:minsups[dbIx]) {
            int ms = msP * onePct;
            std::cout << "DATASET " << DATASETS[dbIx] << " minsup " << ms << std::endl;
            {
                auto t1 = std::chrono::high_resolution_clock::now();
                CFDDiscovery ct(db);
                ct.integratedDFS(ms, maxSize);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                std::cout << "Integrated DFS - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << time << " ms" << std::endl;
            }

            {
                auto t3 = std::chrono::high_resolution_clock::now();
                CFDDiscovery ct2(db);
                ct2.itemsetsFirstBFS(ms, maxSize, SUBSTRATEGY::DFS);
                auto t4 = std::chrono::high_resolution_clock::now();
                auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
                std::cout << "ItemsetsFirst BFSdfs - Minsup " << ms << ": " << ct2.nrCFDs() << " cfds in " << iftime << " ms" << std::endl;
            }

            {
                auto t5 = std::chrono::high_resolution_clock::now();
                CFDDiscovery ct3(db);
                ct3.fdsFirstDFS(ms, maxSize, SUBSTRATEGY::DFS);
                auto t6 = std::chrono::high_resolution_clock::now();
                auto fftime = std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count();
                std::cout << "FDsFirst DFSdfs - Minsup " << ms << ": " << ct3.nrCFDs() << " cfds in " << fftime << " ms" << std::endl;
            }
        }
    }
}

void expAttributes() {
    int dbs[] = {ADULT, MUSHROOM, NURSERY};
    std::map<int, std::vector<int> > attSizes;
    attSizes[NURSERY] = {5,7,9};
    attSizes[ADULT] = {7,9,11};
    attSizes[MUSHROOM] = {7,9,11,13,15};

    for (int dbIx : dbs) {
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
                    CFDDiscovery ct(db);
                    ct.integratedDFS(ms, maxSize);
                    auto t2 = std::chrono::high_resolution_clock::now();
                    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                    std::cout << "Integrated DFS - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << time << " ms" << std::endl;
                }

                {
                    auto t3 = std::chrono::high_resolution_clock::now();
                    CFDDiscovery ct2(db);
                    ct2.itemsetsFirstBFS(ms, maxSize, SUBSTRATEGY::DFS);
                    auto t4 = std::chrono::high_resolution_clock::now();
                    auto iftime = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
                    std::cout << "ItemsetsFirst BFSdfs - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << time << " ms" << std::endl;
                }

                {
                    auto t5 = std::chrono::high_resolution_clock::now();
                    CFDDiscovery ct3(db);
                    ct3.fdsFirstDFS(ms, maxSize, SUBSTRATEGY::DFS);
                    auto t6 = std::chrono::high_resolution_clock::now();
                    auto fftime = std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count();
                    std::cout << "FDsFirst DFSdfs - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << time << " ms" << std::endl;
                }
            }
        }
    }
}

void expTuples() {
    int dbs[] = {ADULT, MUSHROOM, NURSERY};
    std::map<int,int> dbSizes;
    dbSizes[NURSERY] = 12960;
    dbSizes[ADULT] = 48842;
    dbSizes[MUSHROOM] = 8124;
    for (int dbIx : dbs) {
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
                    CFDDiscovery ct(db);
                    ct.integratedDFS(ms, maxSize);
                    auto t2 = std::chrono::high_resolution_clock::now();
                    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                    std::cout << "Integrated DFS - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << time << " ms" << std::endl;
                }

                {
                    auto t3 = std::chrono::high_resolution_clock::now();
                    CFDDiscovery ct(db);
                    ct.itemsetsFirstBFS(ms, maxSize, SUBSTRATEGY::DFS);
                    auto t4 = std::chrono::high_resolution_clock::now();
                    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
                    std::cout << "ItemsetsFirst BFSdfs - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << time << " ms" << std::endl;
                }

                {
                    auto t5 = std::chrono::high_resolution_clock::now();
                    CFDDiscovery ct(db);
                    ct.fdsFirstDFS(ms, maxSize, SUBSTRATEGY::DFS);
                    auto t6 = std::chrono::high_resolution_clock::now();
                    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count();
                    std::cout << "FDsFirst DFSdfs - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << time << " ms" << std::endl;
                }
            }
        }
    }
}

void expMaxSize() {
    int dbs[] = {ADULT, MUSHROOM, NURSERY};
    std::map<int, std::vector<double> > minsups;
    minsups[NURSERY] = {10};
    minsups[ADULT] = {10};
    minsups[MUSHROOM] = {10};
    std::map<int, std::vector<int> > maxSizes;
    maxSizes[NURSERY] = {3,5,7,9};
    maxSizes[ADULT] = {3,5,7,9};
    maxSizes[MUSHROOM] = {5,7,9,11};
    for (int dbIx : dbs) {
        std::ifstream dbFile(concat(3, DATA_PATH, DATASETS[dbIx], ".csv"));
        Database db = DatabaseReader::fromTable(dbFile, 11, ',');
        int onePct = db.size() / 100.0;
        for (int maxSize : maxSizes[dbIx]) {
            for (double msP:minsups[dbIx]) {
                int ms = msP * onePct;
                std::cout << "DATASET " << DATASETS[dbIx] << " minsup " << ms << std::endl;

                {
                    auto t1 = std::chrono::high_resolution_clock::now();
                    CFDDiscovery ct(db);
                    ct.integratedDFS(ms, maxSize);
                    auto t2 = std::chrono::high_resolution_clock::now();
                    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                    std::cout << "Integrated DFS - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << time << " ms" << std::endl;
                }

                {
                    auto t3 = std::chrono::high_resolution_clock::now();
                    CFDDiscovery ct(db);
                    ct.itemsetsFirstBFS(ms, maxSize, SUBSTRATEGY::DFS);
                    auto t4 = std::chrono::high_resolution_clock::now();
                    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
                    std::cout << "ItemsetsFirst BFSdfs - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << time << " ms" << std::endl;
                }

                {
                    auto t5 = std::chrono::high_resolution_clock::now();
                    CFDDiscovery ct(db);
                    ct.fdsFirstDFS(ms, maxSize, SUBSTRATEGY::DFS);
                    auto t6 = std::chrono::high_resolution_clock::now();
                    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count();
                    std::cout << "FDsFirst DFSdfs - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << time << " ms" << std::endl;
                }
            }
        }
    }
}

void expConfidence() {
    int dbs[] = {ADULT, MUSHROOM, NURSERY};
    std::map<int, std::vector<double> > minsups;
    minsups[NURSERY] = {10};
    minsups[ADULT] = {10};
    minsups[MUSHROOM] = {10};
    for (int dbIx : dbs) {
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
                    CFDDiscovery ct(db);
                    ct.integratedDFS(ms, maxSize, conf);
                    auto t2 = std::chrono::high_resolution_clock::now();
                    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                    std::cout << "Integrated DFS - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << time << " ms" << std::endl;
                }

                {
                    auto t3 = std::chrono::high_resolution_clock::now();
                    CFDDiscovery ct(db);
                    ct.itemsetsFirstBFS(ms, maxSize, SUBSTRATEGY::DFS, conf);
                    auto t4 = std::chrono::high_resolution_clock::now();
                    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
                    std::cout << "ItemsetsFirst BFSdfs - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << time << " ms" << std::endl;
                }

                {
                    auto t5 = std::chrono::high_resolution_clock::now();
                    CFDDiscovery ct(db);
                    ct.fdsFirstDFS(ms, maxSize, SUBSTRATEGY::DFS, conf);
                    auto t6 = std::chrono::high_resolution_clock::now();
                    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count();
                    std::cout << "FDsFirst DFSdfs - Minsup " << ms << ": " << ct.nrCFDs() << " cfds in " << time << " ms" << std::endl;
                }
            }
        }
    }
}

void expNrCFD() {
    int dbs[] = {ADULT, MUSHROOM, NURSERY};
    std::map<int, std::vector<double> > minsups;
    minsups[NURSERY] = {5,10,15};
    minsups[ADULT] = {5,10,15};
    minsups[MUSHROOM] = {5,10,15};
    for (int dbIx : dbs) {
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
                    CFDDiscovery ct3(db);
                    ct3.fdsFirstDFS(ms, maxSize, SUBSTRATEGY::DFS, conf);
                    auto t6 = std::chrono::high_resolution_clock::now();
                    auto fftime = std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count();
                    std::cout << "FDsFirst DFSdfs - Minsup " << ms << ": " << ct3.nrCFDs() << " cfds in " << fftime << " ms" << std::endl;
                }
            }
        }
    }
}

int main(int argc, char *argv[]) {
    expSupport();
    std::cout << "*----*" << std::endl;
    std::cout << "*----*" << std::endl;
    expAttributes();
    std::cout << "*----*" << std::endl;
    std::cout << "*----*" << std::endl;
    expMaxSize();
    std::cout << "*----*" << std::endl;
    std::cout << "*----*" << std::endl;
    expConfidence();
    std::cout << "*----*" << std::endl;
    std::cout << "*----*" << std::endl;
    expTuples();
    std::cout << "*----*" << std::endl;
    std::cout << "*----*" << std::endl;
    expBfsDfs();
    std::cout << "*----*" << std::endl;
    std::cout << "*----*" << std::endl;
    expNrCFD();
    return 0;
}
