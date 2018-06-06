#include <xmmintrin.h>
#include "memAlignOS.h"
#include <iostream>
#include <math.h>
#include <chrono>

using namespace std;
//
// Created by grabiger on 05.12.2017.
//
float* smdNewtonReciprocalFourVals(float *a);
float* smdNewtonSqrtRcFourVals(float *a);
void solveTask7() {
std::cout<<"Task 7"<<std::endl;
float* result;
    float *a=new float[4];
    a[0]=1;
    a[1]=0.5;
    a[2]=2;
    a[3]=4;

    std::cout<<"x: "<<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<std::endl<<std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed;

    start = std::chrono::high_resolution_clock::now();
    result=smdNewtonReciprocalFourVals(a);
    finish = std::chrono::high_resolution_clock::now();
    elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start);
    std::cout<<"SIMD: "<<(elapsed.count()*1000)<<"ms"<<std::endl;

    std::cout<<"1/x: "<<result[0]<<" "<<result[1]<<" "<<result[2]<<" "<<result[3]<<" "<<std::endl;
    std::cout<<"resi: "<<(result[0]-1/a[0])<<" "<<(result[1]-1/a[1])<<" "<<(result[2]-1/a[2])<<" "<<(result[3]-1/a[3])<<" "<<std::endl;

    start = std::chrono::high_resolution_clock::now();
    for(int i=0;i<4;i++) {
        result[i]=1./a[i];
    }
    finish = std::chrono::high_resolution_clock::now();
    elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start);
    std::cout<<"Std: "<<(elapsed.count()*1000)<<"ms"<<std::endl;
    std::cout<<"1/x: "<<result[0]<<" "<<result[1]<<" "<<result[2]<<" "<<result[3]<<" "<<std::endl;

    freeAlignedMemOS(result);
    std::cout<<std::endl;


    start = std::chrono::high_resolution_clock::now();
    result=smdNewtonSqrtRcFourVals(a);
    finish = std::chrono::high_resolution_clock::now();
    elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start);
    std::cout<<"SIMD: "<<(elapsed.count()*1000)<<"ms"<<std::endl;
    std::cout<<"1/sqrt(x): "<<result[0]<<" "<<result[1]<<" "<<result[2]<<" "<<result[3]<<" "<<std::endl;
    std::cout<<"resi: "<<(result[0]-1/sqrt(a[0]))<<" "<<(result[1]-1/sqrt(a[1]))<<" "<<(result[2]-1/sqrt(a[2]))<<" "<<(result[3]-1/sqrt(a[3]))<<" "<<std::endl;

    start = std::chrono::high_resolution_clock::now();
    for(int i=0;i<4;i++) {
        result[i]=1./sqrt(a[i]);
    }
    finish = std::chrono::high_resolution_clock::now();
    elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start);
    std::cout<<"Std: "<<(elapsed.count()*1000)<<"ms"<<std::endl;
    std::cout<<"1/sqrt(x): "<<result[0]<<" "<<result[1]<<" "<<result[2]<<" "<<result[3]<<" "<<std::endl;
    freeAlignedMemOS(result);

    delete(a);
};



float* smdNewtonReciprocalFourVals(float *a) {

    __m128 xvalSMD, yvalsSMD,tempTwoSMD;

    float *xvals, *result, *tempTwo;
    memAlignOS((void**) &result,16,16);
    memAlignOS((void**) &xvals,16,16);
    memAlignOS((void**) &tempTwo,16,16);
    xvals[0]=a[0];xvals[1]=a[1];xvals[2]=a[2];xvals[3]=a[3];
    tempTwo[0]=2;

    xvalSMD = _mm_load_ps(xvals);
    tempTwoSMD=_mm_load_ps1(tempTwo);
    yvalsSMD= _mm_rcp_ps (xvalSMD);
    _mm_store_ps(result, yvalsSMD);


        yvalsSMD = _mm_mul_ps(yvalsSMD, _mm_sub_ps(tempTwoSMD, _mm_mul_ps(xvalSMD, yvalsSMD)));


    _mm_store_ps(result, yvalsSMD);

    freeAlignedMemOS(xvals);freeAlignedMemOS(tempTwo);
    return result;
}

float* smdNewtonSqrtRcFourVals(float *a) {

    __m128 xvalSMD, yvalsSMD,tempOneHalfSMD, tempThreeHalfSMD;

    float *xvals, *result, *tempOneHalf, *tempThreeHalf;
    memAlignOS((void**) &result,16,16);
    memAlignOS((void**) &xvals,16,16);
    memAlignOS((void**) &tempOneHalf,16,16);
    memAlignOS((void**) &tempThreeHalf,16,16);
    xvals[0]=a[0];xvals[1]=a[1];xvals[2]=a[2];xvals[3]=a[3];
    tempOneHalf[0]=1./2;
    tempThreeHalf[0]=3./2;

    xvalSMD = _mm_load_ps(xvals);
    tempOneHalfSMD=_mm_load_ps1(tempOneHalf);
    tempThreeHalfSMD=_mm_load_ps1(tempThreeHalf);
    yvalsSMD= _mm_rsqrt_ps (xvalSMD);
    _mm_store_ps(result, yvalsSMD);


        yvalsSMD = _mm_mul_ps(yvalsSMD,_mm_sub_ps(tempThreeHalfSMD, _mm_mul_ps(_mm_mul_ps(tempOneHalfSMD,xvalSMD),_mm_mul_ps(yvalsSMD,yvalsSMD))));


    _mm_store_ps(result, yvalsSMD);

    freeAlignedMemOS(xvals);freeAlignedMemOS(tempOneHalf);freeAlignedMemOS(tempThreeHalf);
    return result;
}


























