#include <iostream>
#include "jacobiIterationSolver.h"



int main() {
    std::cout << "Hello, World!" << std::endl;


    double* sol;

    int n=1000;
    double h=1/(n+1.);



    double* f=new double[n]();
    for(int i=0;i<n;i++) {
        f[i]=i*h*(i*h-1);
      //  std::cout << f[i] << std::endl;

    }



    sol=solveLaeWJacobiIterOfTridiagonalMatrix(-1,-1,2,h,n,f);


    for(int i=0;i<n;i++) {
        std::cout << sol[i] << std::endl;
    }


    return 0;
}