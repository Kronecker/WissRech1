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
            std::cout << iteration <<": "<< resi<< std::endl;
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

    double* fTimesHsquare=new double[n];
    double hsquare=h*h;

    for(int i=0;i<n;i++) {
        fTimesHsquare[i]=f[i]*hsquare;
        lastIterSol[i]=fTimesHsquare[i]/valuesOfMainDiag[i];
    }




    double resi;
    int end=n-1;

    while(iteration<maxIter) {
        actualIteration[0] = 1 / valuesOfMainDiag[0] * (-valOfUpMinDiag * lastIterSol[1] + fTimesHsquare[0]);
        for (int i = 1; i < end; i++) {
            actualIteration[i] = 1 / valuesOfMainDiag[i] *
                                 (-valOfUpMinDiag * lastIterSol[i + 1] - valueOFLowMinDiag * lastIterSol[i - 1] +
                                  fTimesHsquare[i]);
        }
        actualIteration[end] = 1 / valuesOfMainDiag[end] * (-valueOFLowMinDiag * lastIterSol[end - 1] + fTimesHsquare[end]);


        if (!(iteration % (maxIter / 20))) {
            resi=0;
            for(int i=0;i<n;i++) {
                resi+=fabs(actualIteration[i]- lastIterSol[i]);
            }
            std::cout << iteration <<": "<< resi<< std::endl;
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
            std::cout << iteration <<": "<< resi<< std::endl;
        }


        for (int i = 0; i < n*n; i++) {
            lastIterSol[i] = actualIteration[i];
        }
        iteration++;


    }




    return actualIteration;
}

