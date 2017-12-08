#include <xmmintrin.h>
#include "memAlignOS.h"
#include <iostream>

using namespace std;
//
// Created by grabiger on 05.12.2017.
//
float* smdNewtonReciprocalFourVals(float *a);
float* smdNewtonSqrtRcFourVals(float *a);
void solveTask7() {




    float *a=new float[4];
    a[0]=1;
    a[1]=0.5;
    a[2]=2;
    a[3]=4;

    smdNewtonReciprocalFourVals(a);





    delete(a);
};



float* smdNewtonReciprocalFourVals(float *a) {

    __m128 xvalSMD, yvalsSMD,tempTwoSMD;
    int iterations=101;
    float *xvals, *result, *tempTwo;
    memAlignOS((void**) &result,16,16);
    memAlignOS((void**) &xvals,16,16);
    memAlignOS((void**) &tempTwo,16,16);
    xvals[0]=a[0];xvals[1]=a[1];xvals[2]=a[2];xvals[3]=a[3];
    tempTwo[0]=2;tempTwo[1]=2;tempTwo[2]=2;tempTwo[3]=2;


    //c = _mm_add_ps(c, _mm_mul_ps(a, b));
    xvalSMD = _mm_load_ps(xvals);
    tempTwoSMD=_mm_load_ps(tempTwo);
    yvalsSMD= _mm_rcp_ps (xvalSMD);
    _mm_store_ps(result, yvalsSMD);

    std::cout<<result[0]<<" "<<result[1]<<" "<<result[2]<<" "<<result[3]<<" "<<std::endl;
    std::cout<<(result[0]-1/a[0])<<" "<<(result[1]-1/a[1])<<" "<<(result[2]-1/a[2])<<" "<<(result[3]-1/a[3])<<" "<<std::endl;

    for(int i=0;i<iterations;i++) {
        yvalsSMD = _mm_mul_ps(yvalsSMD, _mm_sub_ps(tempTwoSMD, _mm_mul_ps(xvalSMD, yvalsSMD)));
    }

    _mm_store_ps(result, yvalsSMD);
    std::cout<<result[0]<<" "<<result[1]<<" "<<result[2]<<" "<<result[3]<<" "<<std::endl;
    std::cout<<(result[0]-1/a[0])<<" "<<(result[1]-1/a[1])<<" "<<(result[2]-1/a[2])<<" "<<(result[3]-1/a[3])<<" "<<std::endl;

}

float* smdNewtonSqrtRcFourVals(float *a) {

    __m128 xvalSMD, yvalsSMD,tempTwoSMD;
    int iterations=101;
    float *xvals, *result, *tempTwo;
    memAlignOS((void**) &result,16,16);
    memAlignOS((void**) &xvals,16,16);
    memAlignOS((void**) &tempTwo,16,16);
    xvals[0]=a[0];xvals[1]=a[1];xvals[2]=a[2];xvals[3]=a[3];
    tempTwo[0]=2;tempTwo[1]=2;tempTwo[2]=2;tempTwo[3]=2;


    //c = _mm_add_ps(c, _mm_mul_ps(a, b));
    xvalSMD = _mm_load_ps(xvals);
    tempTwoSMD=_mm_load_ps(tempTwo);
    yvalsSMD= _mm_rcp_ps (xvalSMD);
    _mm_store_ps(result, yvalsSMD);

    std::cout<<result[0]<<" "<<result[1]<<" "<<result[2]<<" "<<result[3]<<" "<<std::endl;
    std::cout<<(result[0]-1/a[0])<<" "<<(result[1]-1/a[1])<<" "<<(result[2]-1/a[2])<<" "<<(result[3]-1/a[3])<<" "<<std::endl;

    for(int i=0;i<iterations;i++) {
        yvalsSMD = _mm_mul_ps(yvalsSMD, _mm_sub_ps(tempTwoSMD, _mm_mul_ps(xvalSMD, yvalsSMD)));
    }

    _mm_store_ps(result, yvalsSMD);
    std::cout<<result[0]<<" "<<result[1]<<" "<<result[2]<<" "<<result[3]<<" "<<std::endl;
    std::cout<<(result[0]-1/a[0])<<" "<<(result[1]-1/a[1])<<" "<<(result[2]-1/a[2])<<" "<<(result[3]-1/a[3])<<" "<<std::endl;






}