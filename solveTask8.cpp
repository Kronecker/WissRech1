//
// Created by grabiger on 05.12.2017.
//

#include <iostream>
#include "memAlignOS.h"
#include <xmmintrin.h>
#include "math.h"

int factorialInt(int n);

void solveTask8() {

    int numTerms=12;

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

    x[0]=0.5;
    x[1]=1.;
    x[2]=2.;
    x[3]=3.;

    ones[0]=ones[1]=ones[2]=ones[3]=1;
    //std::cout<<float(1)/factorialInt(numTerms)<<std::endl;
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

    std::cout<<resultExp[0]<<" "<<resultExp[1]<<" "<<resultExp[2]<<" "<<resultExp[3]<<" "<<std::endl;
    std::cout<<(resultExp[0]-exp(x[0]))<<" "<<(resultExp[1]-exp(x[1]))<<" "<<(resultExp[2]-exp(x[2]))<<" "<<(resultExp[3]-exp(x[3]))<<" "<<std::endl;


    std::cout<<std::endl;

    resultExp[0]=resultExp[1]=resultExp[2]=resultExp[3]=0;
    for(int i=numTerms;i>=0;i--) {
        for(int k=0;k<4;k++) {
            resultExp[k] = resultExp[k]*x[k]+1./float(factorialInt(i));

        }
    //    std::cout<<resultExp[2]<<std::endl;
     //   std::cout<<"    >"<<1./float(factorialInt(i))<<std::endl;
    }

    std::cout<<resultExp[0]<<" "<<resultExp[1]<<" "<<resultExp[2]<<" "<<resultExp[3]<<" "<<std::endl;
    std::cout<<(resultExp[0]-exp(x[0]))<<" "<<(resultExp[1]-exp(x[1]))<<" "<<(resultExp[2]-exp(x[2]))<<" "<<(resultExp[3]-exp(x[3]))<<" "<<std::endl;


    resultExp[0]=resultExp[1]=resultExp[2]=resultExp[3]=0;
    for(int i=numTerms;i>0;i--) {
        for(int k=0;k<4;k++) {
            resultExp[k] += pow(x[k],i)/float(factorialInt(i));

        }
       // std::cout<<resultExp[2]<<std::endl;
       // std::cout<<"    >"<<1./float(factorialInt(i))<<std::endl;
    }
    for(int k=0;k<4;k++) {
        resultExp[k] = resultExp[k]+1;

    }
    std::cout<<std::endl<<resultExp[0]<<" "<<resultExp[1]<<" "<<resultExp[2]<<" "<<resultExp[3]<<" "<<std::endl;
    std::cout<<(resultExp[0]-exp(x[0]))<<" "<<(resultExp[1]-exp(x[1]))<<" "<<(resultExp[2]-exp(x[2]))<<" "<<(resultExp[3]-exp(x[3]))<<" "<<std::endl;




};



int factorialInt(int n) {
    int factorial=1;
    for(;n>1;n--)
        factorial*=n;
    return factorial;
}