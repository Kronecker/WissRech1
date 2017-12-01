//
// Created by grabiger on 15.11.2017.
//

#ifndef WISSRECH1_JACOBIITERATIONSOLVER_H
#define WISSRECH1_JACOBIITERATIONSOLVER_H

double* jacobiIterOfBlockMatrixFourDiags(double valLowBlockDiag,double valLowMinDiag,double valMainDiag, double valUpDiag,double valUpBlockDiag, int n, double f, double valBoundary, int* numberOfIterations, bool* diagonalDominant);
double* gaussSeidelIterOfBlockMatrixFourDiags(double valLowBlockDiag,double valLowMinDiag,double valMainDiag, double valUpDiag,double valUpBlockDiag, int n, double f, double valBoundary, int* numberOfIterations, bool* diagonalDominant);
#endif //WISSRECH1_JACOBIITERATIONSOLVER_H
