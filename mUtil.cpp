//
// Created by looky on 22.11.17.
//

#include "mUtil.h"
#include <iostream>


double *linspace(double start, double end, int n) {
    double range = end - start;
    double h = range / (n - 1);

    if (range == 0) {
        return nullptr;
    }

    double *linSpacedVector = new double[n]();

    for (int i = 0; i < n; i++) {
        linSpacedVector[i] = start + h * i;
        //     std::cout<<i<<": "<<linSpacedVector[i]<<std::endl;
    }

    return linSpacedVector;


}


double * vectorTimesTridiag(double *vector, double valueOfLowMinDiag, double valOfUpMinDiag, double *valuesOfMainDiag, int n) {

    double *resultVec = new double[n];
    int end = n - 1;

    resultVec[0] = valuesOfMainDiag[0] * vector[0] + valueOfLowMinDiag * vector[1];
    for (int i = 1; i < end; i++) {
        resultVec[i] =
                valuesOfMainDiag[i] * vector[i] + valueOfLowMinDiag * vector[i + 1] + valOfUpMinDiag * vector[i - 1];
    }

    resultVec[end] = valuesOfMainDiag[end] * vector[end] + valOfUpMinDiag * vector[end - 1];

return resultVec;
}

double * tridiagTimesVector(double *vector, double valueOfLowMinDiag, double valOfUpMinDiag, double *valuesOfMainDiag, int n) {
    return vectorTimesTridiag(vector, valOfUpMinDiag,valueOfLowMinDiag,valuesOfMainDiag,n);
}

double vectorTimesVectorScalar(double* vecLeft, double*vecRight, int n) {
    double result;
    for(int i=0;i<n;i++) {
        result+=vecLeft[i]*vecRight[i];
    }
    return result;
}

double vectorTranspTimesTridiagTimesvector(double* vector, double valueOfLowMinDiag, double valOfUpMinDiag, double *valuesOfMainDiag, int n) {
    double result;

    int end = n - 1;

    result+= (valuesOfMainDiag[0] * vector[0] + valueOfLowMinDiag * vector[1])*vector[0];
    for (int i = 1; i < end; i++) {
        result+=(valuesOfMainDiag[i] * vector[i] + valueOfLowMinDiag * vector[i + 1] + valOfUpMinDiag * vector[i - 1])*vector[i];
    }

    result+=(valuesOfMainDiag[end] * vector[end] + valOfUpMinDiag * vector[end - 1])*vector[end];


return result;
}



