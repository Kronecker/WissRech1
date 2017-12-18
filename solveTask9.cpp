//
// Created by grabiger on 05.12.2017.
//

#include "memAlignOS.h"
#include "mUtil.h"
#include <iostream>

using namespace std;

void solveTask9() {

    float *a,*b;
    memAlignOS((void**) &a,16,16);
    memAlignOS((void**) &b,16,16);

    a[0]=1;a[1]=2;a[2]=-1;a[3]=-2;
    b[0]=4;b[1]=4;b[2]=-3;b[3]=-4;

    std::cout<<subtractAndSumVectorsAligned16(a,b,4)<<std::endl;

};


// changes to jacobi

// temporary vector for values that are lowMinDiag and upMinDiag to align to 16bit
// seperate first iteration and last iteration for whole matrix and for each block to account for boundary values