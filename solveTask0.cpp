//
// Created by looky on 22.11.17.
//
#include <iostream>
#include "jacobiIterationSolver.h"

using namespace std;

void solveTask0() {



    double* sol;

    int n=50;
    double h=1/(n+1.);
    double umax,umin,ustart,uend;



    double* diagonal=new double[n]();

    for(int i=0;i<n;i++) {
        diagonal[i]=(i>20)&&(i<30)?1:0;
        std::cout<<i<<": " << diagonal[i] << std::endl;
    }


    cout<<"hi"<<std::endl;
    sol=solveLaeWJacobiIterOfTridiagonalMatrix(-1,-1,2,0,0,h,n,diagonal);
    cout<<"bye"<<std::endl;
    ustart=sol[0];
    uend=sol[n-1];
    umax=sol[0];
    umin=umax;
    for(int i=0;i<n;i++) {
        if(sol[i]>umax) {
            umax=sol[i];
        }
        if(sol[i]<umin) {
            umin=sol[i];
        }
       // std::cout << sol[i] << std::endl;
    }
    std::cout <<std::endl<<ustart<<" "<<uend<<" "<<umin<<" "<<umax<<" "<<std::endl;







}