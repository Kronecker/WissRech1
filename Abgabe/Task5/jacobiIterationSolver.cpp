//
// Created by grabiger on 15.11.2017.
//

#include <iostream>
#include <math.h>

double* jacobiIterOfBlockMatrixFourDiags(double valLowBlockDiag,double valLowMinDiag,double valMainDiag, double valUpDiag,double valUpBlockDiag, int n, double f, double valBoundary, int* numberOfIterations, bool* diagonalDominant) {

    double* actualIteration=new double[n*n]();
    double* lastIterSol=new double[n*n]();
    int maxIter=9999;
    double tol=0.0001;
    int iteration=0;
    double resi=tol+1;
    int step=(maxIter / 100);
    double relativeTolerance=10000;

    step =2;
    *diagonalDominant=fabs(fabs(valMainDiag)-(fabs(valLowBlockDiag)+fabs(valLowMinDiag)+fabs(valUpDiag)+fabs(valUpBlockDiag)))<=fabs(valMainDiag)/relativeTolerance;

    if(*diagonalDominant) {
        std::cout<<"Matrix is (weak) diagonal dominant"<<std::endl;
    } else {
        std::cout<<"Matrix is not diagonal dominant"<<std::endl;
    }
    // boundary values init (outer)
    for(int i=0;i<n;i++) {
        actualIteration[i]=valBoundary;
        lastIterSol[i]=valBoundary;
        actualIteration[n*(n-1)+i]=valBoundary;
        lastIterSol[n*(n-1)+i]=valBoundary;
    }
    for(int k=1;k<n-1;k++) { // iterate through blocks
        actualIteration[k*n]=valBoundary;
        lastIterSol[k*n]=valBoundary;
        actualIteration[(k+1)*n-1]=valBoundary;
        lastIterSol[(k+1)*n-1]=valBoundary;
    }

    int nm1=n-1;
    int index;
    while(iteration<maxIter&&resi>tol) {

        // first block = boundary with u=0
            // already done above
        // last block = boundary with u=max
            // already done above


        // consecutive blocks

        for(int k=1;k<nm1;k++) { // iterate through blocks

            for(int i=1;i<nm1;i++) {  // iterate in block
                index=k*n+i;
                actualIteration[index]=1/valMainDiag*(f-valLowBlockDiag*lastIterSol[index-n]-valLowMinDiag*lastIterSol[index-1]-valUpDiag*lastIterSol[index+1]-valUpBlockDiag*lastIterSol[index+n]);
            }

        }


        if (!(iteration % step)) {
            resi=0;
            for(int i=0;i<n*n;i++) {
                resi+=fabs(actualIteration[i]- lastIterSol[i]);
            }
         //   std::cout << iteration <<": "<< resi<< std::endl;
        }


        for (int i = 0; i < n*n; i++) {
            lastIterSol[i] = actualIteration[i];
        }
        iteration++;


    }
    std::cout << "Calculation finished after "<<iteration<<" Iterations.(%"<<step<<")"<<std::endl;
    *numberOfIterations=iteration;

    return actualIteration;
}
double* gaussSeidelIterOfBlockMatrixFourDiags(double valLowBlockDiag,double valLowMinDiag,double valMainDiag, double valUpDiag,double valUpBlockDiag, int n, double f, double valBoundary, int* numberOfIterations, bool* diagonalDominant) {

    double* actualIteration=new double[n*n]();
    double* lastIterSol=new double[n*n]();
    int maxIter=9999;
    double tol=0.0001;
    int iteration=0;
    double resi=tol+1;
    int step=(maxIter / 100);
    double relativeTolerance=10000;

    step =2;
    *diagonalDominant=fabs(fabs(valMainDiag)-(fabs(valLowBlockDiag)+fabs(valLowMinDiag)+fabs(valUpDiag)+fabs(valUpBlockDiag)))<=fabs(valMainDiag)/relativeTolerance;

    if(*diagonalDominant) {
        std::cout<<"Matrix is (weak) diagonal dominant"<<std::endl;
    } else {
        std::cout<<"Matrix is not diagonal dominant"<<std::endl;
    }

    // boundary values init (outer)
    for(int i=0;i<n;i++) {
        actualIteration[i]=valBoundary;
        lastIterSol[i]=valBoundary;
        actualIteration[n*(n-1)+i]=valBoundary;
        lastIterSol[n*(n-1)+i]=valBoundary;
    }
    for(int k=1;k<n-1;k++) { // iterate through blocks
        actualIteration[k*n]=valBoundary;
        lastIterSol[k*n]=valBoundary;
        actualIteration[(k+1)*n-1]=valBoundary;
        lastIterSol[(k+1)*n-1]=valBoundary;
    }

    int nm1=n-1;
    int index;
    while(iteration<maxIter&&resi>tol) {

        // first block = boundary with u=0
        // already done above
        // last block = boundary with u=max
        // already done above


        // consecutive blocks

        for(int k=1;k<nm1;k++) { // iterate through blocks

            for(int i=1;i<nm1;i++) {  // iterate in block
                index=k*n+i;
                actualIteration[index]=1/valMainDiag*(f-valLowBlockDiag*actualIteration[index-n]-valLowMinDiag*actualIteration[index-1]-valUpDiag*lastIterSol[index+1]-valUpBlockDiag*lastIterSol[index+n]);
            }

        }


        if (!(iteration % step)) {
            resi=0;
            for(int i=0;i<n*n;i++) {
                resi+=fabs(actualIteration[i]- lastIterSol[i]);
            }
          //     std::cout << iteration <<": "<< resi<< std::endl;
        }


        for (int i = 0; i < n*n; i++) {
            lastIterSol[i] = actualIteration[i];
        }
        iteration++;


    }
    std::cout << "Calculation finished after "<<iteration<<" Iterations.(%"<<step<<")"<<std::endl;

    *numberOfIterations=iteration;

    return actualIteration;
}
