//
// Created by looky on 23.11.17.
//

#include "jacobiIterationSolver.h"

#include "solveTasks.h"
#include "mUtil.h"
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

void solveTask4() {

    ofstream myfile;
   // myfile.open ("M:\\Uni\\WissRech1\\a.txt");

    myfile.open ("./T4a.txt");

    int n = 1000;

    double *r = linspace(-10, 10, n);
    double h = fabs(r[1] - r[0]);
    double *diagonal = new double[n]();
    double* eigenvec;
    double cachedDiagonal = 2 / h / h;
   // std::cout << h << std::endl;


    std::cout<<"Task 4 a)"<<endl;
    for (int i = 0; i < n; i++) {
        //diagonal[i] = ((i > ((int) (double(n) / 3.))) && (i < (int) (double(n) * 2. / 3.)) ? 1 : 0) + cachedDiagonal;
        diagonal[i]=cachedDiagonal+(fabs(r[i])<1?1:0);
        //std::cout << i << ": " << r[i] << std::endl;
    }


    eigenvec=inversePowerIterationForTriDiagonalMatrix(-1 / h / h, -1 / h / h, diagonal, n, h);
    for(int i=0;i<n;i++) {
        myfile << eigenvec[i] << " ";
    }
    myfile.close();
    std::cout<<"Eigenvector saved to T4a.txt"<<std::endl;

    std::cout<<"Task 4 b)"<<endl;
    myfile.open ("./T4b.txt");
    for (int i = 0; i < n; i++) {
        //diagonal[i] = ((i > ((int) (double(n) / 3.))) && (i < (int) (double(n) * 2. / 3.)) ? 1 : 0) + cachedDiagonal;
        diagonal[i]=cachedDiagonal+(fabs(r[i])<1?0:1);
        //std::cout << i << ": " << r[i] << std::endl;
    }


    eigenvec=inversePowerIterationForTriDiagonalMatrix(-1 / h / h, -1 / h / h, diagonal, n, h);
    for(int i=0;i<n;i++) {
        myfile << eigenvec[i] << " ";
    }
    myfile.close();
    std::cout<<"Eigenvector saved to T4b.txt"<<std::endl;

    std::cout<<"Task 4 c)"<<endl;
    myfile.open ("./T4c.txt");
    for (int i = 0; i < n; i++) {
        //diagonal[i] = ((i > ((int) (double(n) / 3.))) && (i < (int) (double(n) * 2. / 3.)) ? 1 : 0) + cachedDiagonal;
        diagonal[i]=cachedDiagonal+r[i]*r[i];
        //std::cout << i << ": " << r[i] << std::endl;
    }


    eigenvec=inversePowerIterationForTriDiagonalMatrix(-1 / h / h, -1 / h / h, diagonal, n, h);
    for(int i=0;i<n;i++) {
        myfile << eigenvec[i] << " ";
    }
    myfile.close();
    std::cout<<"Eigenvector saved to T4c.txt"<<std::endl;

}
















