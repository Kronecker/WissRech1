//
// Created by Looky on 25/11/2017.
//

#include <math.h>
#include "jacobiIterationSolver.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <chrono>
#include "executionData.h"

using namespace std;

double* calculateDiffusion(bool centralDiff, bool jacobiIteration, int n, double veloZero, double f, double valBoundary, bool saveToFile, executionData *execData) ;
void buildAndSaveGnuPlotSkript(std::vector<executionData> execData,string fileName);

void solveTask5() {


    double h, veloZero, veloX, veloY, hSquare, f, valBoundary;
    vector <double *> solutionArrays;
    bool centralDiff, jacobiIteration,saveToFile;

    std::vector<executionData> overviewExecData;

    f = 1;
    valBoundary = 0;
    vector<int> numPointsVec{25, 50, 100};
    vector<double> vZeroVec = {50, 100,200};
    executionData execData;
    centralDiff = true;
    saveToFile=true;

    for (int i = 0; i < numPointsVec.size(); i++) {

        for (int k = 0;k < vZeroVec.size(); k++) {
            jacobiIteration = true;
            solutionArrays.push_back(calculateDiffusion(centralDiff, jacobiIteration, numPointsVec[i], vZeroVec[k], f, valBoundary,saveToFile,&execData));
            overviewExecData.push_back(execData);
            jacobiIteration = false;
            solutionArrays.push_back(calculateDiffusion(centralDiff, jacobiIteration, numPointsVec[i], vZeroVec[k], f, valBoundary,saveToFile,&execData));
            overviewExecData.push_back(execData);
        }
    }



    centralDiff = false;
    saveToFile=true;
    for (int i = 0; i < numPointsVec.size(); i++) {

        for (int k = 0;k < vZeroVec.size(); k++) {
            jacobiIteration = true;
            solutionArrays.push_back(calculateDiffusion(centralDiff, jacobiIteration, numPointsVec[i], vZeroVec[k], f, valBoundary,saveToFile,&execData));
            overviewExecData.push_back(execData);
            jacobiIteration = false;
            solutionArrays.push_back(calculateDiffusion(centralDiff, jacobiIteration, numPointsVec[i], vZeroVec[k], f, valBoundary,saveToFile,&execData));
            overviewExecData.push_back(execData);
        }
    }

    for(int i=0;i<solutionArrays.size();i++) {
        delete(solutionArrays[i]);
    }

    std::cout<<"n\tdd\tv0\titera\ttime:ms\tAlgo\tDiff"<<std::endl;
    for(std::vector<executionData>::size_type i = 0; i != overviewExecData.size(); i++) {
        if(overviewExecData[i].executionTimeMs<1) {
            std::cout << overviewExecData[i].samplePoints << "\t" << (overviewExecData[i].diagonalDomiant ? "y" : "n")
                      << "\t" << overviewExecData[i].veloZero << "\t" << overviewExecData[i].iterations << "\t"
                      << "<1" << "\t" << overviewExecData[i].matrixMethod << "\t"
                      << overviewExecData[i].diffMethod << std::endl;
        } else {
            std::cout << overviewExecData[i].samplePoints << "\t" << (overviewExecData[i].diagonalDomiant ? "y" : "n")
                      << "\t" << overviewExecData[i].veloZero << "\t" << overviewExecData[i].iterations << "\t"
                      << overviewExecData[i].executionTimeMs << "\t" << overviewExecData[i].matrixMethod << "\t"
                      << overviewExecData[i].diffMethod << std::endl;
        }
    }
    buildAndSaveGnuPlotSkript(overviewExecData,"gplotScript");
}

double* calculateDiffusion(bool centralDiff, bool jacobiIteration, int n, double veloZero, double f, double valBoundary, bool saveToFile, executionData *execData) {

    double *solutionVector;
    double h, hSquare, veloX, veloY, valLowBlockDiag, valLowMinDiag, valMainDiag, valUpDiag, valUpBlockDiag;
    int numberOfIterations=0;
    bool diagonalDominant;


    h = 1 / (double(n -1));
    hSquare = h * h;
    veloX = 1 * veloZero / sqrt(5);
    veloY = 2 * veloZero / sqrt(5);

    if (centralDiff) {
        valLowBlockDiag = -1 / hSquare - veloY / 2 / h;   // central
        valLowMinDiag = -1 / hSquare - veloX / 2 / h;
        valUpDiag = -1 / hSquare + veloX / 2 / h;
        valUpBlockDiag = -1 / hSquare + veloY / 2 / h;
        valMainDiag = 4 / hSquare;
    } else {
        valLowBlockDiag = -1 / hSquare - veloY / h;  // upwind
        valLowMinDiag = -1 / hSquare - veloX / h;
        valUpDiag = -1 / hSquare;
        valUpBlockDiag = -1 / hSquare;
        valMainDiag = 4 / hSquare + veloX / h + veloY / h;
    }
    std::cout<<"Starting calculation for n="<<n<<"(h="<<h<<")"<<" and v0="<<veloZero<<" with "<<(jacobiIteration?"Jacobi":"Gauss-Seidel")<<" and "<<(centralDiff?"central diff.":"upwind diff.")<<std::endl;


    // std::cout<<valLowBlockDiag<<" "<< valLowMinDiag<<" "<< valUpDiag<<" "<<valUpBlockDiag <<" "<< valMainDiag<<endl;
    auto start = std::chrono::high_resolution_clock::now();
    if (jacobiIteration) {
        solutionVector = jacobiIterOfBlockMatrixFourDiags(valLowBlockDiag, valLowMinDiag, valMainDiag, valUpDiag,
                                                          valUpBlockDiag, n, f, valBoundary,&numberOfIterations,&diagonalDominant);
    } else {
        solutionVector = gaussSeidelIterOfBlockMatrixFourDiags(valLowBlockDiag, valLowMinDiag, valMainDiag, valUpDiag,
                                                               valUpBlockDiag, n, f, valBoundary,&numberOfIterations,&diagonalDominant);
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start);

    //std::chrono::duration<double> elapsed = finish - start;

    std::cout<<"Execution time was "<<elapsed.count()*1000<<"ms."<<std::endl;

    std::ostringstream fileName;
    fileName << "T5" << (jacobiIteration ? "Jacobi" : "Gauss") << (centralDiff ? "Central" : "Upwind");
    fileName << "_n=" << std::setfill('0') << std::setw(4) << n;
    fileName << "_vo=" << std::setfill('0') << std::setw(3) << veloZero;


    *execData={fileName.str(),round(elapsed.count()*1000),numberOfIterations,n,h,veloZero,diagonalDominant,(jacobiIteration ? "Jacobi" : "Gauss"),(centralDiff ? "Central" : "Upwind")};

    if(saveToFile) {
        int index;
        ofstream myfile;
        fileName << ".dat";

      //  std::cout << fileName.str() << endl;
        myfile.open(fileName.str());

        for (int k = 0; k < n; k++) {
            for (int i = 0; i < n; i++) {
                index = k * n + i;
                myfile<<k*h<<" "<<i*h<<" "<<solutionVector[index];
                myfile << "\n";
            }
            myfile << "\n";
        }
        myfile.close();
      std::cout<<"Results saved to "<<fileName.str()<<" ."<<std::endl;
    }
    std::cout<<endl;



    return solutionVector;
}


void buildAndSaveGnuPlotSkript(std::vector<executionData> execData,string fileName) {

    ofstream myfile;

    myfile.open(fileName+"win");
    myfile<<"set pm3d"<<std::endl;


    for(std::vector<executionData>::size_type i = 0; i != execData.size(); i++) {
        myfile<<"set term wxt "<<(i+1)<<std::endl;
        myfile<<"set title \"T5 n="<<execData[i].samplePoints<<" v0="<<execData[i].veloZero<<" "<<execData[i].matrixMethod<<" "<< execData[i].diffMethod<<(execData[i].diagonalDomiant?" diag dom ":" not diag dom ")<<"\""<<std::endl;
        myfile<<"splot \""<<execData[i].name<<".dat\""<<std::endl;
    }


}
























