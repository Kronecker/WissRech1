//
// Created by grabiger on 05.12.2017.
//

#include <iostream>
#include "memAlignOS.h"
#include <xmmintrin.h>
#include "math.h"
#include <chrono>

// Abbruchkriterium exp : |x| <= 1+1/2N   ==> N=(|x|-1)*2 für x>1w ?

int factorialInt(int n);
float* calcExpWithDivisionFloat(int numTerms, float *xValues);
void compareResultsWithExp(float *x, float *results);
float* calcExpWithExternalFactorial(int numTerms, float *xValues);
float* calcExpWithExternalFactorialSpreadIntoMem(int numTerms, float *xValues);

void solveTask8() {

    int numTerms=100, repeat=10000000;
    float *results, *x;
    x=new float[4];
    x[0]=0.5;
    x[1]=0.25;
    x[2]=10;
    x[3]=15;


    auto start = std::chrono::high_resolution_clock::now();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start);
    //std::cout<<(elapsed.count()*1000)/repeat<<std::endl;

    start = std::chrono::high_resolution_clock::now();
    for(int i=0;i<repeat;i++) {
        results = calcExpWithDivisionFloat(numTerms, x);
        //compareResultsWithExp(x, results);
        freeAlignedMemOS(results);
    }
    finish = std::chrono::high_resolution_clock::now();
    elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start);
    std::cout<<(elapsed.count()*1000)/repeat<<std::endl;
    start = std::chrono::high_resolution_clock::now();
    for(int i=0;i<repeat;i++) {
        results = calcExpWithExternalFactorial(numTerms, x);
      //  compareResultsWithExp(x, results);
        freeAlignedMemOS(results);
    }
    finish = std::chrono::high_resolution_clock::now();
    elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start);
    std::cout<<(elapsed.count()*1000)/repeat<<std::endl;
    start = std::chrono::high_resolution_clock::now();
    for(int i=0;i<repeat;i++) {
        results=calcExpWithExternalFactorialSpreadIntoMem(numTerms,x);
        //compareResultsWithExp(x,results);
        freeAlignedMemOS(results);
    }
    finish = std::chrono::high_resolution_clock::now();
    elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start);
    std::cout<<(elapsed.count()*1000)/repeat<<std::endl;

};

float* calcExpWithDivisionFloat(int numTerms, float *xValues) {
    // Pure SIMD
    float *numTInverseFactorial, *numTermsFloat, *x, *ones, *resultExp;
    __m128 coeffCurrentOrder, xSMD, currentOrder, oneSMD, resultExpSMD;

    memAlignOS((void**) &numTInverseFactorial,16,16);
    memAlignOS((void**) &numTermsFloat,16,16);
    memAlignOS((void**) &x,16,16);
    memAlignOS((void**) &ones,16,16);
    memAlignOS((void**) &resultExp,16,16);
    numTInverseFactorial[0]=numTInverseFactorial[1]=numTInverseFactorial[2]=numTInverseFactorial[3]=float(1)/factorialInt(numTerms);
    numTermsFloat[0]=numTermsFloat[1]=numTermsFloat[2]=numTermsFloat[3]=(float)numTerms;

    x[0]=xValues[0];
    x[1]=xValues[1];
    x[2]=xValues[2];
    x[3]=xValues[3];

    ones[0]=ones[1]=ones[2]=ones[3]=1;
    currentOrder=_mm_load_ps(numTermsFloat);
    coeffCurrentOrder=_mm_load_ps(numTInverseFactorial);
    oneSMD=_mm_load_ps(ones);
    xSMD=_mm_load_ps(x);

    resultExpSMD=_mm_mul_ps(coeffCurrentOrder,xSMD);
    for(int i=numTerms;i>1;i--) {
        coeffCurrentOrder=_mm_mul_ps(coeffCurrentOrder,currentOrder);
        resultExpSMD=_mm_add_ps(resultExpSMD,coeffCurrentOrder);
        currentOrder=_mm_sub_ps(currentOrder,oneSMD);
        resultExpSMD=_mm_mul_ps(resultExpSMD,xSMD);
    }
    resultExpSMD=_mm_add_ps(resultExpSMD,oneSMD);
    _mm_store_ps(resultExp,resultExpSMD);

    freeAlignedMemOS(numTInverseFactorial);freeAlignedMemOS(numTermsFloat);freeAlignedMemOS(x);freeAlignedMemOS(ones);

    return resultExp;
}

float* calcExpWithExternalFactorialSpreadIntoMem(int numTerms, float *xValues) {

    float *coeffCurrentOrder,*x, *ones, *resultExp;
    __m128 coeffCurrentOrderSMD, xSMD, currentOrder, oneSMD, resultExpSMD, dummyFactorial, currentN;

    memAlignOS((void**) &coeffCurrentOrder,16,16*(numTerms-1));
    memAlignOS((void**) &x,16,16);
    memAlignOS((void**) &ones,16,16);
    memAlignOS((void**) &resultExp,16,16);
    ones[0]=ones[1]=ones[2]=ones[3]=1;
    oneSMD=_mm_load_ps(ones);


    dummyFactorial=_mm_load_ps(ones);
    currentN=_mm_load_ps(ones);
    int index;


    for(int i=2;i<=numTerms;i++) {
        currentN=_mm_add_ps(oneSMD,currentN);
        dummyFactorial=_mm_mul_ps(dummyFactorial,currentN);
        index=(i-2)*4;
        _mm_store_ps(&coeffCurrentOrder[index],_mm_div_ps(oneSMD,dummyFactorial));
    }

    x[0]=xValues[0];
    x[1]=xValues[1];
    x[2]=xValues[2];
    x[3]=xValues[3];



    coeffCurrentOrderSMD=_mm_load_ps(&coeffCurrentOrder[4*(numTerms-2)]);

    xSMD=_mm_load_ps(x);

    resultExpSMD=_mm_mul_ps(coeffCurrentOrderSMD,xSMD);
    for(int i=numTerms-3;i>=0;i--) {
        coeffCurrentOrderSMD=_mm_load_ps(&coeffCurrentOrder[i*4]);
        resultExpSMD=_mm_add_ps(resultExpSMD,coeffCurrentOrderSMD);
        resultExpSMD=_mm_mul_ps(resultExpSMD,xSMD);
    }
    resultExpSMD=_mm_add_ps(resultExpSMD,oneSMD);
    resultExpSMD=_mm_mul_ps(resultExpSMD,xSMD);
    resultExpSMD=_mm_add_ps(resultExpSMD,oneSMD);
    _mm_store_ps(resultExp,resultExpSMD);
    freeAlignedMemOS(coeffCurrentOrder);freeAlignedMemOS(x);freeAlignedMemOS(ones);
    return resultExp;
}

float* calcExpWithExternalFactorial(int numTerms, float *xValues) {

    float *numFactorial, *coeffCurrentOrder,*x, *ones, *resultExp;
    __m128 coeffCurrentOrderSMD, xSMD, currentOrder, oneSMD, resultExpSMD;
    float dummyFactorial;

    memAlignOS((void**) &coeffCurrentOrder,16,16);
    memAlignOS((void**) &x,16,16);
    memAlignOS((void**) &ones,16,16);
    memAlignOS((void**) &resultExp,16,16);


    numFactorial=new float[numTerms-1];
    dummyFactorial=1;

    for(int i=2;i<=numTerms;i++) {
        dummyFactorial*=i;
        numFactorial[i-2]=float(1)/dummyFactorial;
    }

    x[0]=xValues[0];
    x[1]=xValues[1];
    x[2]=xValues[2];
    x[3]=xValues[3];

    ones[0]=ones[1]=ones[2]=ones[3]=1;
    coeffCurrentOrder[0]=coeffCurrentOrder[1]=coeffCurrentOrder[2]=coeffCurrentOrder[3]=numFactorial[numTerms-2];

    coeffCurrentOrderSMD=_mm_load_ps(coeffCurrentOrder);
    oneSMD=_mm_load_ps(ones);
    xSMD=_mm_load_ps(x);

    resultExpSMD=_mm_mul_ps(coeffCurrentOrderSMD,xSMD);
    for(int i=numTerms-1;i>1;i--) {
        coeffCurrentOrder[0]=coeffCurrentOrder[1]=coeffCurrentOrder[2]=coeffCurrentOrder[3]=numFactorial[i-2];
        coeffCurrentOrderSMD=_mm_load_ps(coeffCurrentOrder);
        resultExpSMD=_mm_add_ps(resultExpSMD,coeffCurrentOrderSMD);
        resultExpSMD=_mm_mul_ps(resultExpSMD,xSMD);
    }
    resultExpSMD=_mm_add_ps(resultExpSMD,oneSMD);
    resultExpSMD=_mm_mul_ps(resultExpSMD,xSMD);
    resultExpSMD=_mm_add_ps(resultExpSMD,oneSMD);
    _mm_store_ps(resultExp,resultExpSMD);

    delete(numFactorial);
    freeAlignedMemOS(coeffCurrentOrder); freeAlignedMemOS(x); freeAlignedMemOS(ones);

    return resultExp;
}




void compareResultsWithExp(float *x, float *results) {
    std::cout<<results[0]<<" "<<results[1]<<" "<<results[2]<<" "<<results[3]<<std::endl;
    std::cout<<results[0]-exp(x[0])<<" "<<results[1]-exp(x[1])<<" "<<results[2]-exp(x[2])<<" "<<results[3]-exp(x[3])<<std::endl;
}


int factorialInt(int n) {
    int factorial=1;
    for(;n>1;n--)
        factorial*=n;
    return factorial;
}