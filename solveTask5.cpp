//
// Created by Looky on 25/11/2017.
//

#include <math.h>
#include "jacobiIterationSolver.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <iomanip>

using namespace std;

double* calculateDiffusion(bool centralDiff, bool jacobiIteration, int n, double veloZero, double f, double valBoundary, bool saveToFile) ;

void solveTask5() {


    double h, veloZero, veloX, veloY, hSquare, f, valBoundary;
    vector <double *> solutionArrays;
    bool centralDiff, jacobiIteration,saveToFile;


    f = 1;
    valBoundary = 0;
    vector<int> numPointsVec{25, 50, 100};
    vector<double> vZeroVec = {50, 100};

    centralDiff = true;
    jacobiIteration = true;
    saveToFile=true;
    for (int i = 0; i < numPointsVec.size(); i++) {
        for (int k = 0;k < vZeroVec.size(); k++) {
            solutionArrays.push_back(calculateDiffusion(centralDiff, jacobiIteration, numPointsVec[i], vZeroVec[k], f, valBoundary,saveToFile));
        }
    }






    for(int i=0;i<solutionArrays.size();i++) {
        delete(solutionArrays[i]);
    }
}

double* calculateDiffusion(bool centralDiff, bool jacobiIteration, int n, double veloZero, double f, double valBoundary, bool saveToFile) {

    double *solutionVector;
    double h, hSquare, veloX, veloY, valLowBlockDiag, valLowMinDiag, valMainDiag, valUpDiag, valUpBlockDiag;




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



    // std::cout<<valLowBlockDiag<<" "<< valLowMinDiag<<" "<< valUpDiag<<" "<<valUpBlockDiag <<" "<< valMainDiag<<endl;
    if (jacobiIteration) {
        solutionVector = jacobiIterOfBlockMatrixFourDiags(valLowBlockDiag, valLowMinDiag, valMainDiag, valUpDiag,
                                                          valUpBlockDiag, n, f, valBoundary);
    } else {
        solutionVector = gaussSeidelIterOfBlockMatrixFourDiags(valLowBlockDiag, valLowMinDiag, valMainDiag, valUpDiag,
                                                               valUpBlockDiag, n, f, valBoundary);
    }

    if(saveToFile) {
        int index;
        ofstream myfile;
        std::ostringstream  fileName;
        char* formattedNumber=new char[20];

        fileName<<"T5"<<(jacobiIteration?"Jacobi":"Gauss")<<(centralDiff?"Central":"Upwind");
        fileName<<"_n="<< std::setfill('0') << std::setw(4)<<n;
        fileName<<"_vo="<< std::setfill('0') << std::setw(3)<<veloZero;
        fileName<<".dat";

        std::cout<<fileName.str()<<endl;
        myfile.open(fileName.str());

        for (int k = 0; k < n; k++) {
            for (int i = 0; i < n; i++) {
                index = k * n + i;

                sprintf(formattedNumber, "%.10lf", k*h);
                myfile << formattedNumber << " ";

                sprintf(formattedNumber, "%.10lf", i*h);
                myfile << formattedNumber << " ";



                sprintf(formattedNumber, "%.10lf", solutionVector[index]);
                myfile << formattedNumber<< " ";
                myfile << "\n\n";
            }

        }
        myfile.close();

    }
    return solutionVector;
}



























