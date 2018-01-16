//
// Created by looky on 29.12.17.
//

#include <stdio.h>
#include <chrono>
#include <iostream>
#include <stdlib.h>
#include <pthread.h>
#define NUM_THREADS 8

using namespace std;


double scalarprodSingleThread(double *a,double *b, int n);
double scalarprodPThread(double *a, double *b, int n, int procs);
double scalarprodPThreadGlobalNoMutex(double *a, double *b, int n, int procs) ;
double scalarprodPThreadGlobalSummationMutex(double *a, double *b, int n, int procs);
double scalarprodPThreadGlobalResultMutex(double *a, double *b, int n, int procs);

void* subrScalarprodPThread(void* param);
void* subrScalarprodPThreadGlobalNoMutex(void* param);
void* subrScalarprodPThreadGlobalSummationMutex(void* param);
void* subrScalarprodPThreadGlobalResultMutex(void* param);
void* subrScalarprodPThreadNotCached(void* param);

double globalTestSum=0;




typedef struct messagetype{
   double *a,*b;
    int n;
    double sum;
};





void solveTask10() {

    int n=500000;
    double *a,*b, sum=0;
    a=new double[n];
    b=new double[n];

    for(int i=0;i<n;i++) {
        a[i]=i+1;
        b[i]=n-1;
    }

    std::cout<<"Direct       : "<<((double)n*(double)n-1)*(double)n/2<<std::endl;

    // single-thread

    auto start = std::chrono::high_resolution_clock::now();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed;



    int runs=30;
    double procList[]={1,2,4,8,16,32,64,128,256};
    int lengthProcList=9;


    start = std::chrono::high_resolution_clock::now();
    for(int run=0;run<runs;run++) {
        sum = scalarprodSingleThread(a, b, n);
    }
    finish = std::chrono::high_resolution_clock::now();
    elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start);

    std::cout<<"Single Thread: "<<(elapsed.count()*1000)<<"ms"<<"\t"<<sum<<std::endl;




    std::cout << "Multi Thread"<<std::endl;
    for(int k=0;k<lengthProcList;k++) {
        start = std::chrono::high_resolution_clock::now();
        for (int run = 0; run < runs; run++) {
            sum = scalarprodPThread(a, b, n, procList[k]);
        }
        finish = std::chrono::high_resolution_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start);

        std::cout << "Procs: "<<procList[k]<<" \t" << (elapsed.count() * 1000) << "ms  " << "\t" << sum << std::endl;
    }



    std::cout << "Multi Thread Global Variable No Mutex"<<std::endl;
    for(int k=0;k<lengthProcList;k++) {
        start = std::chrono::high_resolution_clock::now();
        for (int run = 0; run < runs; run++) {
            sum = scalarprodPThreadGlobalNoMutex(a, b, n, procList[k]);
        }
        finish = std::chrono::high_resolution_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start);

        std::cout << "Procs: "<<procList[k]<<" \t" << (elapsed.count() * 1000) << "ms  " << "\t" << sum << std::endl;
    }


    std::cout << "Multi Thread Global Variable Summation Mutex"<<std::endl;
    for(int k=0;k<lengthProcList;k++) {
        start = std::chrono::high_resolution_clock::now();
        for (int run = 0; run <runs; run++) {
            sum = scalarprodPThreadGlobalSummationMutex(a, b, n, procList[k]);
        }
        finish = std::chrono::high_resolution_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start);

        std::cout << "Procs: "<<procList[k]<<" \t" << (elapsed.count() * 1000) << "ms  " << "\t" << sum << std::endl;
    }

    std::cout << "Multi Thread Global Variable Result Mutex"<<std::endl;
    for(int k=0;k<lengthProcList;k++) {
        start = std::chrono::high_resolution_clock::now();
        for (int run = 0; run < runs; run++) {
            sum = scalarprodPThreadGlobalResultMutex(a, b, n, procList[k]);
        }
        finish = std::chrono::high_resolution_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start);

        std::cout << "Procs: "<<procList[k]<<" \t" << (elapsed.count() * 1000) << "ms  " << "\t" << sum << std::endl;
    }

}

double scalarprodSingleThread(double *a,double *b, int n) {
    double sum=0;
    for(int i=0;i<n;i++) {
        sum+=a[i]*b[i];
    }
    return sum;
}

double scalarprodPThread(double *a, double *b, int n, int procs) {

    int r,numOfVals=n/procs;    // evenly distribute, add rest to the last proc

    messagetype *messArray=new messagetype[procs];
    pthread_t tid, *thread_id;
    thread_id=new pthread_t[procs];

    for(int i=0;i<procs;i++) {

        messArray[i].a=&a[i*numOfVals];
        messArray[i].b=&b[i*numOfVals];
        messArray[i].n=(i<procs-1)?numOfVals:n-numOfVals*(procs-1);
        messArray[i].sum=0;

        pthread_create(&tid, NULL, subrScalarprodPThreadNotCached, (void*) &messArray[i]);
        thread_id[i]=tid;

    }

    double sum=0;
    for(int i=0;i<procs;i++) {

        pthread_join(thread_id[i],NULL);
        sum+=messArray[i].sum;

    }




    delete(messArray);
    delete(thread_id);
    return sum;
}

double scalarprodPThreadGlobalNoMutex(double *a, double *b, int n, int procs) {

    int r,numOfVals=n/procs;    // evenly distribute, add rest to the last proc

    messagetype *messArray=new messagetype[procs];
    pthread_t tid, *thread_id;
    thread_id=new pthread_t[procs];
    globalTestSum=0;

    for(int i=0;i<procs;i++) {

        messArray[i].a=&a[i*numOfVals];
        messArray[i].b=&b[i*numOfVals];
        messArray[i].n=(i<procs-1)?numOfVals:n-numOfVals*(procs-1);
        messArray[i].sum=0;

        pthread_create(&tid, NULL, subrScalarprodPThreadGlobalNoMutex, (void*) &messArray[i]);
        thread_id[i]=tid;

    }
    for(int i=0;i<procs;i++) {
        pthread_join(thread_id[i],NULL);
    }

    delete(messArray);
    delete(thread_id);
    return globalTestSum;
}

pthread_mutex_t mtx;

double scalarprodPThreadGlobalSummationMutex(double *a, double *b, int n, int procs) {
// Copy of scalarprodPThreadGlobalNoMutex
    int r,numOfVals=n/procs;    // evenly distribute, add rest to the last proc

    messagetype *messArray=new messagetype[procs];
    pthread_t tid, *thread_id;
    thread_id=new pthread_t[procs];
    globalTestSum=0;

    pthread_mutex_init(&mtx,NULL);

    for(int i=0;i<procs;i++) {

        messArray[i].a=&a[i*numOfVals];
        messArray[i].b=&b[i*numOfVals];
        messArray[i].n=(i<procs-1)?numOfVals:n-numOfVals*(procs-1);
        messArray[i].sum=0;

        pthread_create(&tid, NULL, subrScalarprodPThreadGlobalSummationMutex, (void*) &messArray[i]);
        thread_id[i]=tid;

    }
    for(int i=0;i<procs;i++) {
        pthread_join(thread_id[i],NULL);
    }

    delete(messArray);
    delete(thread_id);
    return globalTestSum;
}

double scalarprodPThreadGlobalResultMutex(double *a, double *b, int n, int procs) {
// Copy of scalarprodPThreadGlobalNoMutex
    int r,numOfVals=n/procs;    // evenly distribute, add rest to the last proc

    messagetype *messArray=new messagetype[procs];
    pthread_t tid, *thread_id;
    thread_id=new pthread_t[procs];
    globalTestSum=0;

    pthread_mutex_init(&mtx,NULL);

    for(int i=0;i<procs;i++) {

        messArray[i].a=&a[i*numOfVals];
        messArray[i].b=&b[i*numOfVals];
        messArray[i].n=(i<procs-1)?numOfVals:n-numOfVals*(procs-1);
        messArray[i].sum=0;

        pthread_create(&tid, NULL, subrScalarprodPThreadGlobalResultMutex, (void*) &messArray[i]);
        thread_id[i]=tid;

    }
    for(int i=0;i<procs;i++) {
        pthread_join(thread_id[i],NULL);
    }

    delete(messArray);
    delete(thread_id);
    return globalTestSum;
}

void* subrScalarprodPThread(void* param) {
    messagetype *mess=(messagetype *) param;
    double sum=0;
    double *a=mess->a,*b=mess->b;
    int n=mess->n;
    for(int i=0;i<n;i++) {
        sum+=a[i]*b[i];
    }
    mess->sum=sum;
}
void* subrScalarprodPThreadNotCached(void* param) {
    messagetype *mess=(messagetype *) param;
    mess->sum=0;
    for(int i=0;i<mess->n;i++) {
        mess->sum+=mess->a[i]*mess->b[i];
    }

}

void* subrScalarprodPThreadGlobalNoMutex(void* param) {
    messagetype *mess=(messagetype *) param;
    for(int i=0;i<mess->n;i++) {
        globalTestSum+=mess->a[i]*mess->b[i];
    }
}

void* subrScalarprodPThreadGlobalSummationMutex(void* param) {
    messagetype *mess=(messagetype *) param;
    for(int i=0;i<mess->n;i++) {
        pthread_mutex_lock(&mtx);
        globalTestSum+=mess->a[i]*mess->b[i];
        pthread_mutex_unlock(&mtx);
    }
}

void* subrScalarprodPThreadGlobalResultMutex(void* param) {
    messagetype *mess=(messagetype *) param;
    double sum=0;
    for(int i=0;i<mess->n;i++) {

        sum+=mess->a[i]*mess->b[i];

    }
    pthread_mutex_lock(&mtx);
    globalTestSum+=sum;
    pthread_mutex_unlock(&mtx);
}

