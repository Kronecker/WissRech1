//
// Created by looky on 22.11.17.
//

#include "jacobiIterationSolver.h"
#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

void solveTask3() {

    ofstream myfile;

    myfile.open("./T3");

    double *sol;
    double **f2Dim;
    int n;
    double h;
    double umax, umin, ustart, uend;
    char *formattedNumber = new char[50];
    int index;
    int lengthNs = 5;
    double *nValues = new double[lengthNs];
    nValues[0] = 50;
    nValues[1] = 200;
    nValues[2] = 300;
    nValues[3] = 400;
    nValues[4] = 600;


    for (int i = 0; i < lengthNs; i++) {

        n = nValues[i];
        myfile<<"# n= "<<n<<std::endl;
        h = 1 / (n + 1.);

        f2Dim = new double *[n];
        for (int i = 0; i < n; i++) {
            f2Dim[i] = new double[n]();
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                f2Dim[i][j] = i * h * (i * h - 1) * j * h * (j * h - 1);
                // std::cout <<f2Dim[i][j] << std::endl;
            }
        }
        std::cout << "starte 2dim" << std::endl;
        sol = solveLaeWJacobiIterOfBlockMatrix(-1, -1, 4, -1, -1, h, n, f2Dim);

        umax = sol[0];
        umin = umax;

        for (int i = 0; i < (n - 2); i++) {
            for (int j = 0; j < (n - 2); j++) {
                index = i * (n - 2) + j;
                if (sol[index] > umax) {
                    umax = sol[index];
                }
                if (sol[index] < umin) {
                    umin = sol[index];
                }
                sprintf(formattedNumber, "%.10lf", sol[index]);
              //  std::cout << formattedNumber << " ";
                myfile << formattedNumber << "\t";

            }
            //std::cout << std::endl;
            myfile << formattedNumber << "\r\n";
        }
        std::cout << std::endl << umin << " " << umax << " " << std::endl;


        f2Dim = new double *[n];
        for (int i = 0; i < n; i++) {
            delete(f2Dim[i]);
        }
        delete(f2Dim);
        delete(sol);
    }
    myfile.close();
}