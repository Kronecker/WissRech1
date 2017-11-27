//
// Created by grabiger on 27.11.2017.
//
#include <string>
#ifndef WISSRECH1_EXECUTIONDATA_H
#define WISSRECH1_EXECUTIONDATA_H
#include <vector>

struct executionData {
    std::string name;
    double executionTimeMs;
    int iterations;
    int samplePoints;
    double h;
    double veloZero;
    std::string matrixMethod;
    std::string diffMethod;
};


#endif //WISSRECH1_EXECUTIONDATA_H
