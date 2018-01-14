//
// Created by looky on 14.01.18.
//
#include <chrono>
#include <iostream>
#include <math.h>
#include <pthread.h>
#include <fstream>
#include <sstream>


using namespace std;



double* jacobiIterOfBlockMatrixFourDiags2(double valLowBlockDiag,double valLowMinDiag,double valMainDiag, double valUpDiag,double valUpBlockDiag, int n, double f, double valBoundary);
double* jacobiIterOfBlockMatrixFourDiagsPThread(double valLowBlockDiag,double valLowMinDiag,double valMainDiag, double valUpDiag,double valUpBlockDiag, int n, double f, double valBoundary, int procs);

void solveTask12();




void solveTask12() {

    int n=1024;
    double f=0;
    double valBoundary=0;
    double veloZero=200;

    double *solutionVector;
    double h, hSquare, veloX, veloY, valLowBlockDiag, valLowMinDiag, valMainDiag, valUpDiag, valUpBlockDiag;



    h = 1 / (double(n -1));
    hSquare = h * h;
    veloX = 1 * veloZero / sqrt(5);
    veloY = 2 * veloZero / sqrt(5);


    valLowBlockDiag = -1 / hSquare - veloY / h;  // upwind
    valLowMinDiag = -1 / hSquare - veloX / h;
    valUpDiag = -1 / hSquare;
    valUpBlockDiag = -1 / hSquare;
    valMainDiag = 4 / hSquare + veloX / h + veloY / h;

    auto start = std::chrono::high_resolution_clock::now();

        solutionVector = jacobiIterOfBlockMatrixFourDiags2(valLowBlockDiag, valLowMinDiag, valMainDiag, valUpDiag,
                                                          valUpBlockDiag, n, f, valBoundary);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start);


    bool saveToFile=true;


    if(saveToFile) {

        int index;
        ofstream myfile;
        myfile.open("single.dat");

        for (int k = 0; k < n; k++) {
            for (int i = 0; i < n; i++) {
                index = k * n + i;
                myfile<<k*h<<" "<<i*h<<" "<<solutionVector[index];
                myfile << "\n";
            }
            myfile << "\n";
        }
        myfile.close();
        std::cout<<"Results saved."<<std::endl;
    }
    std::cout<<endl;




    int procs=4;
    start = std::chrono::high_resolution_clock::now();

    solutionVector = jacobiIterOfBlockMatrixFourDiagsPThread(valLowBlockDiag, valLowMinDiag, valMainDiag, valUpDiag,
                                                       valUpBlockDiag, n, f, valBoundary, procs);

    finish = std::chrono::high_resolution_clock::now();
    elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start);


    if(saveToFile) {

        int index;
        ofstream myfile;
        myfile.open("threads.dat");

        for (int k = 0; k < n; k++) {
            for (int i = 0; i < n; i++) {
                index = k * n + i;
                myfile<<k*h<<" "<<i*h<<" "<<solutionVector[index];
                myfile << "\n";
            }
            myfile << "\n";
        }
        myfile.close();
        std::cout<<"Results saved."<<std::endl;
    }
    std::cout<<endl;



























}









double* jacobiIterOfBlockMatrixFourDiags2(double valLowBlockDiag,double valLowMinDiag,double valMainDiag, double valUpDiag,double valUpBlockDiag, int n, double f, double valBoundary) {

    double* actualIteration=new double[n*n]();
    double* lastIterSol=new double[n*n]();
    double* temp;
    int maxIter=2000;
//    double tol=0.0001;
    int iteration=0;
  //  double resi=tol+1;

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
    while(iteration<maxIter) {  // removed &&resi>tol to compare by iteration speed

        for(int k=1;k<nm1;k++) { // iterate through blocks

            for(int i=1;i<nm1;i++) {  // iterate in block
                index=k*n+i;
                actualIteration[index]=1/valMainDiag*(f-valLowBlockDiag*lastIterSol[index-n]-valLowMinDiag*lastIterSol[index-1]-valUpDiag*lastIterSol[index+1]-valUpBlockDiag*lastIterSol[index+n]);
            }

        }

        temp=actualIteration;
        actualIteration=lastIterSol;
        lastIterSol=actualIteration;


    }

    return lastIterSol;
}
double* jacobiIterOfBlockMatrixFourDiagsPThread(double valLowBlockDiag,double valLowMinDiag,double valMainDiag, double valUpDiag,double valUpBlockDiag, int n, double f, double valBoundary, int procs) {

    double* actualIteration=new double[n*n]();
    double* lastIterSol=new double[n*n]();
    double* temp;
    int maxIter=2000;
//    double tol=0.0001;
    int iteration=0;
    //  double resi=tol+1;

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
    while(iteration<maxIter) {  // removed &&resi>tol to compare by iteration speed

        for(int k=1;k<nm1;k++) { // iterate through blocks

            for(int i=1;i<nm1;i++) {  // iterate in block
                index=k*n+i;
                actualIteration[index]=1/valMainDiag*(f-valLowBlockDiag*lastIterSol[index-n]-valLowMinDiag*lastIterSol[index-1]-valUpDiag*lastIterSol[index+1]-valUpBlockDiag*lastIterSol[index+n]);
            }

        }

        temp=actualIteration;
        actualIteration=lastIterSol;
        lastIterSol=actualIteration;


    }

    return lastIterSol;
}
