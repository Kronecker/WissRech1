//
// Created by looky on 22.11.17.
//

#include <iostream>
#include <iomanip>
#include <random>
#include "jacobiIterationSolver.h"
#include "mUtil.h"


double* inversePowerIterationForTriDiagonalMatrix(double valueOfLowerDiag, double valueOfUpperDiag, double* valuesOfMainDiag, int n, double h) {

    std::random_device r;
    std::mt19937 gen(r());
    std::uniform_real_distribution<> dis(0.0,1.0);

    double sum=0;
    double* eigenvec;
    double* eigenvecGuess=new double[n];
    double eigenValue, eigenValueGuess;
    int iteration=0;
    int maxIter=10000000;
    double tol=0.00001;


    for(int i=0;i<n;i++) {
        eigenvecGuess[i]=dis(gen);
        sum+=eigenvecGuess[i];
    }
    for(int i=0;i<n;i++) {
        eigenvecGuess[i]=eigenvecGuess[i]/sum;
    }
    double resi=tol+1;
    eigenvec = solveLaeWJacobiIterOfTridiagonalMatrixVariMainDiag(valueOfLowerDiag, valueOfUpperDiag,
                                                                  valuesOfMainDiag, h, n, eigenvecGuess);
    eigenValueGuess=vectorTranspTimesTridiagTimesvector(eigenvec,valueOfLowerDiag, valueOfUpperDiag,valuesOfMainDiag, n)/vectorTimesVectorScalar(eigenvec,eigenvec,n)/vectorTimesVectorScalar(eigenvec,eigenvec,n);

    while(iteration<maxIter && resi>tol) {
        eigenvec = solveLaeWJacobiIterOfTridiagonalMatrixVariMainDiag(valueOfLowerDiag, valueOfUpperDiag,
                                                                      valuesOfMainDiag, h, n, eigenvecGuess);
        sum = 0;
        for (int i = 0; i < n; i++) {
            sum += eigenvec[i];

        }
        for (int i = 0; i < n; i++) {
            eigenvec[i] = eigenvec[i] / sum;
            //std::cout<<i<<": "<<eigenvecGuess[i]<<" << "<<eigenvec[i]<<std::endl;
        }
        eigenValue=vectorTranspTimesTridiagTimesvector(eigenvec,valueOfLowerDiag, valueOfUpperDiag,valuesOfMainDiag, n)/vectorTimesVectorScalar(eigenvec,eigenvec,n);
        resi=fabs(eigenValue-eigenValueGuess);
        eigenValueGuess=eigenValue;

        delete(eigenvecGuess);  // very bad memory managment...i'm sorry
        eigenvecGuess = eigenvec;
        std::cout<<"Iteration: "<<iteration<<": "<<"Residual EW: "<<resi<<"  Ew: "<< eigenValue<<std::endl;
        iteration++;
    }

    std::cout<<"Residual EW: "<<resi<<"  Ew: "<< eigenValue<<std::endl;

    return eigenvec;
}