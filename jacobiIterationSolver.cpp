//
// Created by grabiger on 15.11.2017.
//

#include "jacobiIterationSolver.h"
#include <iostream>

double* solveLaeWJacobiIterOfTridiagonalMatrix(double valueOFLowMinDiag, double valOfUpMinDiag, double valueOfMainDiag, double h, int n, double* f) {

    double* actualIteration=new double[n];
    double* lastIterSol=new double[n]();
    int maxIter=1000;
    double tol=0.0001;
    int iteration=0;

    double* fTimesHsquare=new double[n];
    double hsquare=h*h;

    for(int i=0;i<n;i++) {
        fTimesHsquare[i]=f[i]*hsquare;
        lastIterSol[i]=fTimesHsquare[i]/valueOfMainDiag;
    }





    int end=n-1;

    while(iteration<maxIter){
        actualIteration[0]=1/valueOfMainDiag*(valOfUpMinDiag*lastIterSol[1]+fTimesHsquare[0]);
        for(int i=1;i<end;i++) {
            actualIteration[i]=1/valueOfMainDiag*(valOfUpMinDiag*lastIterSol[i+1]+valueOFLowMinDiag*lastIterSol[i-1]+fTimesHsquare[i]);
        }
        actualIteration[end]=1/valueOfMainDiag*(valueOFLowMinDiag*lastIterSol[end-1]+fTimesHsquare[end]);

        for(int i=0;i<n;i++) {
            lastIterSol[i]=actualIteration[i];
        }
        iteration++;
    }





    return actualIteration;
}





