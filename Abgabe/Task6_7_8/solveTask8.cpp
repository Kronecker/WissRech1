//
// Created by grabiger on 05.12.2017.
//

#include <iostream>
#include "../../memAlignOS.h"
#include <xmmintrin.h>
#include "math.h"
#include <chrono>

// Abbruchkriterium exp : |x| <= 1+1/2N   ==> N=(|x|-1)*2 für x>1w ?

int factorialInt(int n);
void compareResultsWithExp(float *x, float *results);
void calcExpWithExternalFactorialSpreadIntoMemExternResults(int numTerms, float *xValues, float *resultsExp);
int estimateIterations(float *x);
// some additional tries
float* calcExpWithExternalFactorial(int numTerms, float *xValues);
float* calcExpWithExternalFactorialSpreadIntoMem(int numTerms, float *xValues);
float* calcExpWithDivisionFloat(int numTerms, float *xValues);


void solveTask8() {

    int numTerms=100, repeat=1000000;
    float *results, *x,*resultExp;
    x=new float[4];
    x[0]=1;
    x[1]=4;
    x[2]=-1;
    x[3]=-4;

    numTerms=estimateIterations(x);

    std::cout<<"Doing "<<numTerms<<" orders."<<std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed;
    //std::cout<<(elapsed.count()*1000)/repeat<<std::endl;
    start = std::chrono::high_resolution_clock::now();
    memAlignOS((void**) &resultExp,16,16);
    for(int i=0;i<repeat;i++) {
        resultExp[0]=exp(x[0]);
        resultExp[1]=exp(x[1]);
        resultExp[2]=exp(x[2]);
        resultExp[3]=exp(x[3]);
    }
    finish = std::chrono::high_resolution_clock::now();
    elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start);
    std::cout<<"Std:  "<<(elapsed.count()*1000)/repeat<<"ms/iteration"<<std::endl;

    start = std::chrono::high_resolution_clock::now();
    memAlignOS((void**) &resultExp,16,16);
    for(int i=0;i<repeat;i++) {
        calcExpWithExternalFactorialSpreadIntoMemExternResults(numTerms,x,resultExp);

    }
    finish = std::chrono::high_resolution_clock::now();
    elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start);
    std::cout<<"SIMD: "<<(elapsed.count()*1000)/repeat<<"ms/iteration"<<std::endl;
    compareResultsWithExp(x,resultExp);
    freeAlignedMemOS(resultExp);
};

void calcExpWithExternalFactorialSpreadIntoMemExternResults(int numTerms, float *xValues, float *resultExp) {

    float *coeffCurrentOrder,*x, *ones;
    __m128 coeffCurrentOrderSMD, xSMD, currentOrder, oneSMD, resultExpSMD, dummyFactorial, currentN, zero, signMask;

    memAlignOS((void**) &coeffCurrentOrder,16,16*(numTerms-1));
    memAlignOS((void**) &x,16,16);
    memAlignOS((void**) &ones,16,16);

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
    zero=_mm_setzero_ps();
    signMask=_mm_cmplt_ps(xSMD,zero);

    xSMD=_mm_sub_ps(_mm_andnot_ps(signMask,xSMD),_mm_and_ps(signMask,xSMD));  // absolut values

    resultExpSMD=_mm_mul_ps(coeffCurrentOrderSMD,xSMD);
    for(int i=numTerms-3;i>=0;i--) {
        coeffCurrentOrderSMD=_mm_load_ps(&coeffCurrentOrder[i*4]);
        resultExpSMD=_mm_add_ps(resultExpSMD,coeffCurrentOrderSMD);
        resultExpSMD=_mm_mul_ps(resultExpSMD,xSMD);
    }
    resultExpSMD=_mm_add_ps(resultExpSMD,oneSMD);
    resultExpSMD=_mm_mul_ps(resultExpSMD,xSMD);
    resultExpSMD=_mm_add_ps(resultExpSMD,oneSMD);
    resultExpSMD=_mm_add_ps(_mm_andnot_ps(signMask,resultExpSMD),_mm_and_ps(signMask,_mm_div_ps(oneSMD,resultExpSMD)));

    _mm_store_ps(resultExp,resultExpSMD);
    freeAlignedMemOS(coeffCurrentOrder);freeAlignedMemOS(x);freeAlignedMemOS(ones);
}



void compareResultsWithExp(float *x, float *results) {
    std::cout<<"x: \t"<<x[0]<<"\t"<<x[1]<<"\t"<<x[2]<<"\t"<<x[3]<<std::endl;
    std::cout<<"exp(x): \t"<<results[0]<<"\t"<<results[1]<<"\t"<<results[2]<<"\t"<<results[3]<<std::endl;
    std::cout<<"resi: \t"<<results[0]-exp(x[0])<<"\t"<<results[1]-exp(x[1])<<"\t"<<results[2]-exp(x[2])<<"\t"<<results[3]-exp(x[3])<<std::endl;
}

int factorialInt(int n) {
    int factorial=1;
    for(;n>1;n--)
        factorial*=n;
    return factorial;
}

int estimateIterations(float *x) {
    // Abbruchkriterium exp : |x| <= 1+1/2N   ==> N=(|x|-1)*2 für x>1w ?
    // aber mind. 4 Iterationen
    float max=fabs(x[0]);
    if(fabs(x[1]))
        max=fabs(x[1]);
    if(fabs(x[2]))
        max=fabs(x[2]);
    if(fabs(x[3]))
        max=fabs(x[3]);
    if(max>3) {
        return (max-1)*2*4;  //twicefor better accuracy
    }else {
        return 4;
    }
}
// some additional tries, x<0 is not handled correctly
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

float* calcExpWithDivisionFloat(int numTerms, float *xValues) {
    // Pure SIMD
    float *numTInverseFactorial, *numTermsFloat, *x, *ones, *resultExp;
    __m128 coeffCurrentOrder, xSMD, currentOrder, oneSMD, resultExpSMD, signMask, zero;

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
