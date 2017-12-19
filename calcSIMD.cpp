//
// Created by looky on 19.12.17.
//

#include "calcSIMD.h"
#include <math.h>
#include <iostream
#include <xmmintrin.h>
#include "memAlignOS.h"


using namespace std;

float* jacobiIterOfBlockMatrixFourDiagsFloatSIMD(float valLowBlockDiag,float valLowMinDiag,float valMainDiag, float valUpDiag,float valUpBlockDiag, int n, float f, float valBoundary, int* numberOfIterations, bool* diagonalDominant) {

    n=1024;
    float *actualIteration, *lastIterSol, *lowItaDy, *upItaDy, *temp;
    memAlignOS((void**) &actualIteration, 16,4*n*n);
    memAlignOS((void**) &lastIterSol, 16,4*n*n);
    memAlignOS((void**) &lowItaDy, 16,16);
    memAlignOS((void**) &upItaDy, 16,16);

    __m128 centralUSIMD,lowDiagUSIMD,upDiagUSIMD, lowBlockUSIMD, upBlockUSIMD,valLowBlockDiagSIMD, valLowMinDiagSIMD, valMainDiagInvertSIMD, valUpDiagSIMD, valUpBlockDiagSIMD, fSIMD, valBoundarySMD;

    int maxIter=9999;
    float tol=0.0001;
    int iteration=0;
    float resi=tol+1;
    int step=(maxIter / 100);step =2;
    float relativeTolerance=10000;

    *diagonalDominant=fabs(fabs(valMainDiag)-(fabs(valLowBlockDiag)+fabs(valLowMinDiag)+fabs(valUpDiag)+fabs(valUpBlockDiag)))<=fabs(valMainDiag)/relativeTolerance;

    // Store maatrix into SIMD vectors, lowItaDy is used as dummy variable
    lowItaDy[0]=lowItaDy[1]=lowItaDy[2]=lowItaDy[3]=valLowBlockDiag;
    valLowBlockDiagSIMD=_mm_load_ps(lowItaDy);
    lowItaDy[0]=lowItaDy[1]=lowItaDy[2]=lowItaDy[3]=valLowMinDiag;
    valLowMinDiagSIMD=_mm_load_ps(lowItaDy);
    lowItaDy[0]=lowItaDy[1]=lowItaDy[2]=lowItaDy[3]=1/valMainDiag;
    valMainDiagInvertSIMD=_mm_load_ps(lowItaDy);
    lowItaDy[0]=lowItaDy[1]=lowItaDy[2]=lowItaDy[3]=valUpDiag;
    valUpDiagSIMD=_mm_load_ps(lowItaDy);
    lowItaDy[0]=lowItaDy[1]=lowItaDy[2]=lowItaDy[3]=valUpBlockDiag;
    valUpBlockDiagSIMD=_mm_load_ps(lowItaDy);
    lowItaDy[0]=lowItaDy[1]=lowItaDy[2]=lowItaDy[3]=f;
    fSIMD=_mm_load_ps(lowItaDy);


    if(*diagonalDominant) {
        std::cout<<"Matrix is (weak) diagonal dominant"<<std::endl;
    } else {
        std::cout<<"Matrix is not diagonal dominant"<<std::endl;
    }


    int nm1=n-1;
    int nm3=n-3;
    int index;
    while(iteration<maxIter&&resi>tol) {

        // first block = boundary with u=0
        for(int i=4;i<nm3;i+=4) {  // iterate in block
            actualIteration[i]=1/valMainDiag*(f-valLowBlockDiag*lastIterSol[i-n]-valLowMinDiag*lastIterSol[index-1]-valUpDiag*lastIterSol[index+1]-valUpBlockDiag*lastIterSol[index+n]);

        }
        // consecutive blocks
        for(int k=2;k<nm1;k++) { // iterate through blocks

            for(int i=1;i<nm1;i+=4) {  // iterate in block
                index=k*n+i;
                actualIteration[index]=1/valMainDiag*(f-valLowBlockDiag*lastIterSol[index-n]-valLowMinDiag*lastIterSol[index-1]-valUpDiag*lastIterSol[index+1]-valUpBlockDiag*lastIterSol[index+n]);

            }

        }

        // last block = boundary with u=max
        for(int i=1;i<nm1;i+=4) {  // iterate in block
            index=i;
            actualIteration[index]=1/valMainDiag*(f-valLowBlockDiag*lastIterSol[index-n]-valLowMinDiag*lastIterSol[index-1]-valUpDiag*lastIterSol[index+1]-valUpBlockDiag*lastIterSol[index+n]);

        }

        if (!(iteration % step)) {
            resi=0;
            for(int i=0;i<n*n;i++) {
                resi+=fabs(actualIteration[i]- lastIterSol[i]);
            }
            //   std::cout << iteration <<": "<< resi<< std::endl;
        }


        temp=lastIterSol;  // #swap actualIteration and lastIterSol
        lastIterSol=actualIteration;
        actualIteration=temp;
        iteration++;


    }
    std::cout << "Calculation finished after "<<iteration<<" Iterations.(%"<<step<<")"<<std::endl;
    *numberOfIterations=iteration;

    freeAlignedMemOS(actualIteration);
    freeAlignedMemOS(lastIterSol);
    freeAlignedMemOS(lowItaDy);
    freeAlignedMemOS(upItaDy);

    return lastIterSol;   // lastIterSol isntead of actual is returned, to account for swapping of Vectors at the end of while loop (#swap)
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


















