//
// Created by looky on 22.11.17.
//

#ifndef WISSRECH1_MUTIL_H
#define WISSRECH1_MUTIL_H

double *linspace(double start, double end, int n);

double *
vectorTimesTridiag(double *vector, double valueOfLowMinDiag, double valOfUpMinDiag, double *valuesOfMainDiag, int n);

double *
tridiagTimesVector(double *vector, double valueOfLowMinDiag, double valOfUpMinDiag, double *valuesOfMainDiag, int n);

double vectorTimesVectorScalar(double* vecLeft, double*vecRight, int n);
double vectorTranspTimesTridiagTimesvector(double* vector, double valueOfLowMinDiag, double valOfUpMinDiag, double *valuesOfMainDiag, int n);

#endif //WISSRECH1_MUTIL_H
