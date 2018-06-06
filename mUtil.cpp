//
// Created by looky on 22.11.17.
//

#include "mUtil.h"
#include <iostream>
#include <xmmintrin.h>
#include <pmmintrin.h>
#include "memAlignOS.h"


float subtractAndSumVectorsAligned16(float *vec_a,float *vec_b, int n) {

    int blocks;
    float  *result;
    memAlignOS((void**) &result,16,(size_t) 16);
    float resultDirect;
    __m128 a, b, c, diff, zero, signMask;
    if ((n % 4)) {
        n = (n / 4) * 4; //floor to multiple of 4
    }
    std::cout << "n = " << n << std::endl;

    c = _mm_setzero_ps();
    zero=_mm_setzero_ps();

    blocks = n / 4;
    for (int i = 0; i < blocks; i++) {
        a = _mm_load_ps(&(vec_a[i * 4]));
        b = _mm_load_ps(&(vec_b[i * 4]));
        diff=_mm_sub_ps(a, b);
        signMask=_mm_cmplt_ps(diff, zero);
        diff=_mm_sub_ps(_mm_andnot_ps(signMask,diff),_mm_and_ps(signMask,diff));
        c = _mm_add_ps(c,diff);
    }

    c = _mm_hadd_ps(_mm_hadd_ps(c, c), c);
    _mm_store_ps1(result, c);
    resultDirect=result[0];
    freeAlignedMemOS(result);
    return resultDirect;
}





double *linspace(double start, double end, int n) {
    double range = end - start;
    double h = range / (n - 1);

    if (range == 0) {
        return nullptr;
    }

    double *linSpacedVector = new double[n]();

    for (int i = 0; i < n; i++) {
        linSpacedVector[i] = start + h * i;
        //     std::cout<<i<<": "<<linSpacedVector[i]<<std::endl;
    }

    return linSpacedVector;


}


double * vectorTimesTridiag(double *vector, double valueOfLowMinDiag, double valOfUpMinDiag, double *valuesOfMainDiag, int n) {

    double *resultVec = new double[n];
    int end = n - 1;

    resultVec[0] = valuesOfMainDiag[0] * vector[0] + valueOfLowMinDiag * vector[1];
    for (int i = 1; i < end; i++) {
        resultVec[i] =
                valuesOfMainDiag[i] * vector[i] + valueOfLowMinDiag * vector[i + 1] + valOfUpMinDiag * vector[i - 1];
    }

    resultVec[end] = valuesOfMainDiag[end] * vector[end] + valOfUpMinDiag * vector[end - 1];

return resultVec;
}

double * tridiagTimesVector(double *vector, double valueOfLowMinDiag, double valOfUpMinDiag, double *valuesOfMainDiag, int n) {
    return vectorTimesTridiag(vector, valOfUpMinDiag,valueOfLowMinDiag,valuesOfMainDiag,n);
}

double vectorTimesVectorScalar(double* vecLeft, double*vecRight, int n) {
    double result;
    for(int i=0;i<n;i++) {
        result+=vecLeft[i]*vecRight[i];
    }
    return result;
}

double vectorTranspTimesTridiagTimesvector(double* vector, double valueOfLowMinDiag, double valOfUpMinDiag, double *valuesOfMainDiag, int n) {
    double result;

    int end = n - 1;

    result+= (valuesOfMainDiag[0] * vector[0] + valueOfLowMinDiag * vector[1])*vector[0];
    for (int i = 1; i < end; i++) {
        result+=(valuesOfMainDiag[i] * vector[i] + valueOfLowMinDiag * vector[i + 1] + valOfUpMinDiag * vector[i - 1])*vector[i];
    }

    result+=(valuesOfMainDiag[end] * vector[end] + valOfUpMinDiag * vector[end - 1])*vector[end];


return result;
}



