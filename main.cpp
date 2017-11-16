#include <iostream>
#include "jacobiIterationSolver.h"
#include <fstream>
using namespace std;

int main() {
    std::cout << "Hello, World!" << std::endl;

    ofstream myfile;
    myfile.open ("M:\\Uni\\WissRech1\\a.txt");





    double* sol;

    int n=50;
    double h=1/(n+1.);
    double umax,umin,ustart,uend;

/*

    double* f=new double[n]();
    for(int i=0;i<n;i++) {
        f[i]=i*h*(i*h-1);
        std::cout << f[i] << std::endl;

    }



    sol=solveLaeWJacobiIterOfTridiagonalMatrix(-1,-1,2,0,0,h,n,f);

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


*/




    double** f2Dim=new double*[n];
    for(int i=0;i<n;i++) {
        f2Dim[i]=new double[n]();
    }

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            f2Dim[i][j]=i*h*(i*h-1)*j*h*(j*h-1);
           // std::cout <<f2Dim[i][j] << std::endl;
        }
    }
    std::cout<<"starte 2dim"<<std::endl;
    sol=solveLaeWJacobiIterOfBlockMatrix(-1,-1,4,-1,-1,h,n,f2Dim);

    char* formattedNumber=new char[50];

    umax=sol[0];
    umin=umax;
    int index;
    for(int i=0;i<(n-2);i++) {
        for(int j=0;j<(n-2);j++) {
            index=i*(n-2)+j;
            if (sol[index] > umax) {
                umax = sol[index];
            }
            if (sol[index] < umin) {
                umin = sol[index];
            }
            sprintf(formattedNumber,"%.10lf",sol[index]);
            std::cout <<formattedNumber << " ";
            myfile<<formattedNumber << "\t";

        }
        std::cout<<std::endl;
        myfile<<formattedNumber << "\r\n";
    }
    std::cout <<std::endl<<umin<<" "<<umax<<" "<<std::endl;



    myfile.close();
    return 0;
}