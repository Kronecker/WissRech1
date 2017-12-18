//
// Created by grabiger on 05.12.2017.
//

#include <stdlib.h>
#include <iostream>
#include "xmmintrin.h"
#include "memAlignOS.h"
#include "pmmintrin.h"
#include <chrono>
#include <math.h>
using namespace std;



void solveTask6() {
std::cout<<"Task 6"<<std::endl;
    int n=5000, blocks;
    int numOfIta=10000;
    float *vec_a, *vec_b, *result;
    float resultDirect;
    __m128 a,b,c;
    if((n%4)) {
       n=(n/4)*4; //floor to multiple of 4
    }
    std::cout<<"n = "<<n<<std::endl;

    memAlignOS((void**) &vec_a,16,(size_t) 4*n);

    memAlignOS((void**)&vec_b,16,(size_t) 4*n);

    memAlignOS((void**) &result,16,(size_t) 16);


    for(int i=0;i<n;i++) {
        vec_a[i]=(i+1);
        vec_b[i]=(n-i);
    }

    auto start = std::chrono::high_resolution_clock::now();
    for(int i=0;i<numOfIta;i++) {
        c = _mm_setzero_ps();
        blocks = n / 4;
        for (int i = 0; i < blocks; i++) {
            a = _mm_load_ps(&(vec_a[i * 4]));
            b = _mm_load_ps(&(vec_b[i * 4]));
            c = _mm_add_ps(c, _mm_mul_ps(a, b));
        }

        c = _mm_hadd_ps(_mm_hadd_ps(c, c),c);
        _mm_store_ps1(result, c);
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start);
    std::cout<<"SIMD: "<<(elapsed.count()*1000)/numOfIta<<"ms"<<std::endl;



    start = std::chrono::high_resolution_clock::now();
    for(int i=0;i<numOfIta;i++) {
        resultDirect=0;
        for(int i=0;i<n;i++) {
            resultDirect+=(vec_a[i]) * (vec_b[i]);
        }
    }
    finish = std::chrono::high_resolution_clock::now();
    elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start);
    std::cout<<"Std: "<<(elapsed.count()*1000)/numOfIta<<"ms"<<std::endl;

    freeAlignedMemOS(vec_a);
    freeAlignedMemOS(vec_b);
    freeAlignedMemOS(result);
};