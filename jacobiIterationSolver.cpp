//
// Created by grabiger on 15.11.2017.
//

#include "jacobiIterationSolver.h"
#include <iostream>
#include <math.h>


double* solveLaeWJacobiIterOfTridiagonalMatrix(double valueOFLowMinDiag, double valOfUpMinDiag, double valueOfMainDiag,double solBoundLow, double solBoundHigh, double h, int N, double* f) {

    double* actualIteration=new double[N];
    actualIteration[0]=solBoundLow;
    actualIteration[N-1]=solBoundHigh;
    double* lastIterSol=new double[N]();
    int maxIter=10;
    double tol=0.0001;
    int iteration=0;

    int n=N-2;

    double* fTimesHsquare=new double[n];
    double hsquare=h*h;



    for(int i=0;i<n;i++) {
        fTimesHsquare[i]=f[i+1]*hsquare;
        lastIterSol[i]=fTimesHsquare[i]/valueOfMainDiag;
    }



    double resi;


    while(iteration<maxIter) {
        actualIteration[1] = 1 / valueOfMainDiag * (-valOfUpMinDiag * lastIterSol[2] + fTimesHsquare[0]);
        for (int i = 2; i < n; i++) {
            actualIteration[i] = 1 / valueOfMainDiag *
                                 (-valOfUpMinDiag * lastIterSol[i + 1] - valueOFLowMinDiag * lastIterSol[i - 1] +
                                  fTimesHsquare[i-1]);
        }
        actualIteration[n] = 1 / valueOfMainDiag * (-valueOFLowMinDiag * lastIterSol[n] + fTimesHsquare[n-1]);


        if (!(iteration % (maxIter / 20))) {
            resi=0;
            for(int i=1;i<n;i++) {
                resi+=fabs(actualIteration[i]- lastIterSol[i]);
            }
//            std::cout << iteration <<": "<< resi<< std::endl;
        }


        for (int i = 0; i < n; i++) {
            lastIterSol[i] = actualIteration[i];
        }
        iteration++;


    }





    return actualIteration;
}

double* solveLaeWJacobiIterOfTridiagonalMatrixVariMainDiag(double valueOFLowMinDiag, double valOfUpMinDiag, double* valuesOfMainDiag, double h, int n, double* f) {

    double* actualIteration=new double[n];
    double* lastIterSol=new double[n]();
    int maxIter=5000000;
    double tol=0.0001;
    int iteration=0;


    double resi;
    int end=n-1;



    while(iteration<maxIter) {
        actualIteration[0] = 1 / valuesOfMainDiag[0] * (-valOfUpMinDiag * lastIterSol[1] + f[0]);
        for (int i = 1; i < end; i++) {
            actualIteration[i] = 1 / valuesOfMainDiag[i] *
                                 (-valOfUpMinDiag * lastIterSol[i + 1] - valueOFLowMinDiag * lastIterSol[i - 1] +
                                  f[i]);
        }
        actualIteration[end] = 1 / valuesOfMainDiag[end] * (-valueOFLowMinDiag * lastIterSol[end - 1] + f[end]);


        if (!(iteration % (maxIter / 20))) {
            resi=0;
            for(int i=0;i<n;i++) {
                resi+=fabs(actualIteration[i]- lastIterSol[i]);
            }
//            std::cout << iteration <<": "<< resi<< std::endl;
            if(resi<tol) {
//                std::cout << "Residual below tolerance"<< std::endl;
                break;
            }
        }


        for (int i = 0; i < n; i++) {
            lastIterSol[i] = actualIteration[i];
        }
        iteration++;


    }





    return actualIteration;
}


double* solveLaeWJacobiIterOfBlockMatrix(double valueOfLowMinDiag, double valOfUpMinDiag, double valueOfMainDiag,double valueOfLowBlockDia, double valueOfUpBlockDia, double h, int N, double** f) {

    int n=N-2;
    double* actualIteration=new double[n*n]();
    double* lastIterSol=new double[n*n]();
    int maxIter=5000;
    double tol=0.0001;
    int iteration=0;



    double* fTimesHsquare=new double[n*n];
    double hsquare=h*h;


    int index=0;
    for(int i=1;i<=n;i++) {
        for(int j=1;j<=n;j++) {
            index=(i-1)*n+j-1;  // row major
            fTimesHsquare[index] = f[i][j] * hsquare;
            lastIterSol[index] = fTimesHsquare[index] / valueOfMainDiag;
        }
    }



    double resi;


    while(iteration<maxIter) {

        // 1st block
        actualIteration[0] = 1 / valueOfMainDiag * (-valOfUpMinDiag * lastIterSol[1]
                                                    - valueOfUpBlockDia*lastIterSol[n] + fTimesHsquare[0]);
        for(int i=1;i<n-1;i++) {
            actualIteration[i] = 1 / valueOfMainDiag *
                                 (-valOfUpMinDiag * lastIterSol[i + 1] - valueOfLowMinDiag * lastIterSol[i - 1]
                                  -valueOfUpBlockDia*lastIterSol[n+i]+
                                  fTimesHsquare[i-1]);
        }
        actualIteration[n-1] = 1 / valueOfMainDiag * (-valueOfLowMinDiag * lastIterSol[n-2]- valueOfUpBlockDia*lastIterSol[2*n-1] + fTimesHsquare[n-1]);


        // consecutive blocks
        for(int k=1;k<n-1;k++) {

            actualIteration[k*n] = 1 / valueOfMainDiag *
                                 (-valOfUpMinDiag * lastIterSol[k*n+1] - valueOfUpBlockDia * lastIterSol[(k+1)*n]
                                  -valueOfLowBlockDia*lastIterSol[(k-1)*n]+
                                  fTimesHsquare[k*n]);
            for (int i = 1; i < n - 1; i++) {
                actualIteration[k*n+i] = 1 / valueOfMainDiag *
                                       (-valOfUpMinDiag * lastIterSol[k*n+1+i]-valueOfLowMinDiag*lastIterSol[k*n-1+i]
                                        - valueOfUpBlockDia * lastIterSol[(k+1)*n+i] -valueOfLowBlockDia*lastIterSol[(k-1)*n+i]+
                                        fTimesHsquare[k*n+i]);
            }
            actualIteration[(k+1)*n-1] = 1 / valueOfMainDiag *
                                     (-valueOfLowMinDiag*lastIterSol[(k+1)*n-2] - valueOfUpBlockDia * lastIterSol[(k+2)*n-1]
                                      -valueOfLowBlockDia*lastIterSol[(k)*n-1]+
                                      fTimesHsquare[(k+1)*n-1]);
        }


        // last block


        actualIteration[(n-1)*n] = 1 / valueOfMainDiag *
                               (-valOfUpMinDiag * lastIterSol[(n-1)*n+1]
                                -valueOfLowBlockDia*lastIterSol[(n-2)*n]+
                                fTimesHsquare[(n-1)*n]);
        for (int i = 1; i < n - 1; i++) {
            actualIteration[(n-1)*n+i] = 1 / valueOfMainDiag *
                                     (-valOfUpMinDiag * lastIterSol[(n-1)*n+1+i]-valueOfLowMinDiag*lastIterSol[(n-1)*n-1+i]
                                      -valueOfLowBlockDia*lastIterSol[(n-2)*n+i]+
                                      fTimesHsquare[(n-1)*n+i]);
        }
        actualIteration[n*n-1] = 1 / valueOfMainDiag *
                                     (-valueOfLowMinDiag*lastIterSol[n*n-2]
                                      -valueOfLowBlockDia*lastIterSol[(n-1)*n-1]+
                                      fTimesHsquare[n*n-1]);






        if (!(iteration % (maxIter / 10))) {
            resi=0;
            for(int i=0;i<n*n;i++) {
                resi+=fabs(actualIteration[i]- lastIterSol[i]);
            }
//            std::cout << iteration <<": "<< resi<< std::endl;
        }


        for (int i = 0; i < n*n; i++) {
            lastIterSol[i] = actualIteration[i];
        }
        iteration++;


    }




    return actualIteration;
}

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

float* jacobiIterOfBlockMatrixFourDiagsFloatSIMD(float valLowBlockDiag,float valLowMinDiag,float valMainDiag, float valUpDiag,float valUpBlockDiag, int n, float f, float valBoundary, int* numberOfIterations, bool* diagonalDominant) {

    float* actualIteration=new float[n*n]();
    float* lastIterSol=new float[n*n]();
    int maxIter=9999;
    float tol=0.0001;
    int iteration=0;
    float resi=tol+1;
    int step=(maxIter / 100);
    float relativeTolerance=10000;

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


float* jacobiIterOfBlockMatrixFourDiags(float valLowBlockDiag,float valLowMinDiag,float valMainDiag, float valUpDiag,float valUpBlockDiag, int n, float f, float valBoundary, int* numberOfIterations, bool* diagonalDominant) {

    float* actualIteration=new float[n*n]();
    float* lastIterSol=new float[n*n]();
    int maxIter=9999;
    float tol=0.0001;
    int iteration=0;
    float resi=tol+1;
    int step=(maxIter / 100);
    float relativeTolerance=10000;

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
















