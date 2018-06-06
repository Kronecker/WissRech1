//
// Created by grabiger on 15.11.2017.
//

#ifndef WISSRECH1_JACOBIITERATIONSOLVER_H
#define WISSRECH1_JACOBIITERATIONSOLVER_H

double* solveLaeWJacobiIterOfTridiagonalMatrix(double valueOFLowMinDiag, double valOfUpMinDiag, double valueOfMainDiag,double solBoundLow, double solBoundHigh, double h, int N, double* f);
double* solveLaeWJacobiIterOfTridiagonalMatrixVariMainDiag(double valueOFLowMinDiag, double valOfUpMinDiag, double* valuesOfMainDiag, double h,int n, double* f);
double* solveLaeWJacobiIterOfBlockMatrix(double valueOfLowMinDiag, double valOfUpMinDiag, double valueOfMainDiag,double valueOfLowBlockDia, double valueOfUpBlockDia, double h, int N, double** f);
double* inversePowerIterationForTriDiagonalMatrix(double valueOfLowerDiag, double valueOfUpperDiag, double* valuesOfMainDiag, int n, double h);
double* jacobiIterOfBlockMatrixFourDiags(double valLowBlockDiag,double valLowMinDiag,double valMainDiag, double valUpDiag,double valUpBlockDiag, int n, double f, double valBoundary, int* numberOfIterations, bool* diagonalDominant);
float* jacobiIterOfBlockMatrixFourDiags(float valLowBlockDiag,float valLowMinDiag,float valMainDiag, float valUpDiag,float valUpBlockDiag, int n, float f, float valBoundary, int* numberOfIterations, bool* diagonalDominant);
double* gaussSeidelIterOfBlockMatrixFourDiags(double valLowBlockDiag,double valLowMinDiag,double valMainDiag, double valUpDiag,double valUpBlockDiag, int n, double f, double valBoundary, int* numberOfIterations, bool* diagonalDominant);
#endif //WISSRECH1_JACOBIITERATIONSOLVER_H
