//
// Created by looky on 14.01.18.
//
#include <chrono>
#include <iostream>
#include <math.h>
#include <pthread.h>
#include <fstream>
#include <sstream>


using namespace std;

double* jacobiIterOfBlockMatrixFourDiags2(double valLowBlockDiag,double valLowMinDiag,double valMainDiag, double valUpDiag,double valUpBlockDiag, int n, double f, double valBoundary);
double* jacobiIterOfBlockMatrixFourDiagsPThread(double valLowBlockDiag,double valLowMinDiag,double valMainDiag, double valUpDiag,double valUpBlockDiag, int n, double f, double valBoundary, int procs);
double* jacobiIterOfBlockMatrixFourDiagsPThreadBarrier(double valLowBlockDiag,double valLowMinDiag,double valMainDiag, double valUpDiag,double valUpBlockDiag, int n, double f, double valBoundary, int procs);
void* subrJacobiThreads(void* param);
void* subrJacobiThreadsBarrierCached(void* param);
void solveTask12();

int globalMaxIteration=2000;


int main() {
    solveTask12();
}

void solveTask12() {

    int n=1024;
    double f=1;
    double valBoundary=0;
    double veloZero=200;

    double *solutionVector,*solutionVector2;
    double h, hSquare, veloX, veloY, valLowBlockDiag, valLowMinDiag, valMainDiag, valUpDiag, valUpBlockDiag;



    h = 1 / (double(n -1));
    hSquare = h * h;
    veloX = 1 * veloZero / sqrt(5);
    veloY = 2 * veloZero / sqrt(5);


    valLowBlockDiag = -1 / hSquare - veloY / h;  // upwind
    valLowMinDiag = -1 / hSquare - veloX / h;
    valUpDiag = -1 / hSquare;
    valUpBlockDiag = -1 / hSquare;
    valMainDiag = 4 / hSquare + veloX / h + veloY / h;

    auto start = std::chrono::high_resolution_clock::now();

        solutionVector = jacobiIterOfBlockMatrixFourDiags2(valLowBlockDiag, valLowMinDiag, valMainDiag, valUpDiag,
                                                          valUpBlockDiag, n, f, valBoundary);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start);
    std::cout<<"Threadless Jacobi finished after "<<elapsed.count()*1000<<"ms."<<std::endl;


    std::cout<<std::endl<<"Jacobi with PThreads started"<<std::endl;
    int procs[]={1,2,4,8,16,32};
    int procNum=6;
    for(int i=0;i<procNum;i++) {
        start = std::chrono::high_resolution_clock::now();

        solutionVector2 = jacobiIterOfBlockMatrixFourDiagsPThread(valLowBlockDiag, valLowMinDiag, valMainDiag, valUpDiag,
                                                                 valUpBlockDiag, n, f, valBoundary, procs[i]);

        finish = std::chrono::high_resolution_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start);
        std::cout<<(procs[i]>10?"Calculation finished with ":"Calculation finished with  ")<<procs[i]<<" threads after "<<elapsed.count()*1000<<"ms."<<std::endl;

      delete (solutionVector2);

    }

    std::cout<<std::endl<<"Jacobi with PThreads Sync started"<<std::endl;

    for(int i=0;i<procNum;i++) {
        start = std::chrono::high_resolution_clock::now();

        solutionVector2 = jacobiIterOfBlockMatrixFourDiagsPThreadBarrier(valLowBlockDiag, valLowMinDiag, valMainDiag, valUpDiag,
                                                                         valUpBlockDiag, n, f, valBoundary, procs[i]);

        finish = std::chrono::high_resolution_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start);
        std::cout<<(procs[i]>10?"Calculation finished with ":"Calculation finished with  ")<<procs[i]<<" threads after "<<elapsed.count()*1000<<"ms."<<std::endl;
    }
    delete(solutionVector);
}

double* jacobiIterOfBlockMatrixFourDiags2(double valLowBlockDiag,double valLowMinDiag,double valMainDiag, double valUpDiag,double valUpBlockDiag, int n, double f, double valBoundary) {

    double* actualIteration=new double[n*n]();
    double* lastIterSol=new double[n*n]();
    double* temp;
    int maxIter=globalMaxIteration;
    int iteration=0;

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
    while(iteration<maxIter) {  // removed &&resi>tol to compare by iteration speed

        for(int k=1;k<nm1;k++) { // iterate through blocks
            for(int i=1;i<nm1;i++) {  // iterate in block
                index=k*n+i;
                actualIteration[index]=1/valMainDiag*(f-valLowBlockDiag*lastIterSol[index-n]-valLowMinDiag*lastIterSol[index-1]-valUpDiag*lastIterSol[index+1]-valUpBlockDiag*lastIterSol[index+n]);
            }
        }
        temp=actualIteration;
        actualIteration=lastIterSol;
        lastIterSol=temp;
        iteration++;
    }
    delete(actualIteration);
    return lastIterSol;
}


typedef struct messJacobiShared{
    double *actualIteration;
    double *lastIterSol;
    double valMainDiag;
    double valLowBlockDiag;
    double valLowMinDiag;
    double valUpDiag;
    double valUpBlockDiag;
    double f;
    int n;
    int iteration;
    int maxIteration;
};

typedef struct messJacobiPrivate{
    int threadId;
    int startBlockIndex;
    int endBlockIndex;
    messJacobiShared* mJShared;
};


double* jacobiIterOfBlockMatrixFourDiagsPThread(double valLowBlockDiag,double valLowMinDiag,double valMainDiag, double valUpDiag,double valUpBlockDiag, int n, double f, double valBoundary, int procs) {

    double* actualIteration=new double[n*n]();
    double* lastIterSol=new double[n*n]();
    double* temp;
    int maxIter=globalMaxIteration;
    int iteration=0;

    // prepare data with boundary values
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

    // first and last block are boundary value only and do not need any computation
    // -> number ofblocks to iterate thorugh is n-2

    int numBlocksPerThread=(n-2)/procs;  // integer division gives number of blocks per thread to distribute, we also need the rest to
    int restBlocks=(n-2)-procs*numBlocksPerThread; //these blocks must be divided as well (1 per thread til none are left)

    messJacobiShared mJShared={actualIteration, lastIterSol, valMainDiag, valLowBlockDiag, valLowMinDiag,valUpDiag, valUpBlockDiag,f,n,0,maxIter};

    messJacobiPrivate* mJPrivates=new messJacobiPrivate[procs];
    mJPrivates[0].threadId=0;
    mJPrivates[0].mJShared=&mJShared;
    mJPrivates[0].startBlockIndex=1;
    mJPrivates[0].endBlockIndex=mJPrivates[0].startBlockIndex+numBlocksPerThread+(restBlocks>0?1:0)-1;
    restBlocks--;

    for(int l=1;l<procs;l++) {
        mJPrivates[l].threadId=l;
        mJPrivates[l].mJShared=&mJShared;
        mJPrivates[l].startBlockIndex=1+mJPrivates[l-1].endBlockIndex;
        mJPrivates[l].endBlockIndex=mJPrivates[l].startBlockIndex+numBlocksPerThread+(restBlocks>0?1:0)-1;
        restBlocks--;
    }

    pthread_t tid, *thread_id;
    thread_id=new pthread_t[procs];

    while(iteration<maxIter) {  // removed &&resi>tol to compare by iteration speed

        for(int i=0;i<procs;i++) {

            pthread_create(&tid, NULL, subrJacobiThreads, (void*) &mJPrivates[i]);
            thread_id[i]=tid;

        }

        for(int i=0;i<procs;i++) {

            pthread_join(thread_id[i],NULL);

        }

        temp=mJShared.actualIteration;
        mJShared.actualIteration=mJShared.lastIterSol;
        mJShared.lastIterSol=temp;

        iteration++;
    }
    delete(mJShared.actualIteration);
    return mJShared.lastIterSol;
}

pthread_barrier_t iterationBarrier;
double* jacobiIterOfBlockMatrixFourDiagsPThreadBarrier(double valLowBlockDiag,double valLowMinDiag,double valMainDiag, double valUpDiag,double valUpBlockDiag, int n, double f, double valBoundary, int procs) {

    double* actualIteration=new double[n*n]();
    double* lastIterSol=new double[n*n]();
    int maxIter=globalMaxIteration;

    // prepare data with boundary values
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
    // first and last block are boundary value only and do not need any computation
    // -> number of blocks to iterate thorugh is n-2

    int numBlocksPerThread=(n-2)/procs;  // integer division gives number of blocks per thread to distribute, we also need the rest to
    int restBlocks=(n-2)-procs*numBlocksPerThread; //these blocks must be divided as well (1 per thread til none are left)

    messJacobiShared mJShared={actualIteration, lastIterSol, valMainDiag, valLowBlockDiag, valLowMinDiag,valUpDiag, valUpBlockDiag,f,n,0,maxIter};

    messJacobiPrivate* mJPrivates=new messJacobiPrivate[procs];
    mJPrivates[0].threadId=0;
    mJPrivates[0].mJShared=&mJShared;
    mJPrivates[0].startBlockIndex=1;
    mJPrivates[0].endBlockIndex=mJPrivates[0].startBlockIndex+numBlocksPerThread+(restBlocks>0?1:0)-1;
    restBlocks--;

    for(int l=1;l<procs;l++) {
        mJPrivates[l].threadId=l;
        mJPrivates[l].mJShared=&mJShared;
        mJPrivates[l].startBlockIndex=1+mJPrivates[l-1].endBlockIndex;
        mJPrivates[l].endBlockIndex=mJPrivates[l].startBlockIndex+numBlocksPerThread+(restBlocks>0?1:0)-1;
        restBlocks--;
    }

    pthread_t tid, *thread_id;
    thread_id=new pthread_t[procs];
    pthread_barrier_init(&iterationBarrier,NULL,procs);

        for(int i=0;i<procs;i++) {

            pthread_create(&tid, NULL, subrJacobiThreadsBarrierCached, (void*) &mJPrivates[i]);
            thread_id[i]=tid;

        }

        for(int i=0;i<procs;i++) {

            pthread_join(thread_id[i],NULL);

        }

    delete(mJShared.actualIteration);
    return mJShared.lastIterSol;
}
void* subrJacobiThreads(void* param) {
    messJacobiPrivate *m=(messJacobiPrivate *) param;
   // Sharing is caring, but caching is smashing.
    double *actualIteration=m->mJShared->actualIteration, *lastIterSol=m->mJShared->lastIterSol;
    double valLowMinDiag=m->mJShared->valLowMinDiag, valLowBlockDiag=m->mJShared->valLowBlockDiag, valMainDiag=m->mJShared->valMainDiag, valUpDiag=m->mJShared->valUpDiag,
            valUpBlockDiag=m->mJShared->valUpBlockDiag, f=m->mJShared->f;
    double *temp;
    int nm1=m->mJShared->n-1, n=m->mJShared->n;
    int startBlockIndex=m->startBlockIndex, endBlockIndex=m->endBlockIndex;
    int index;
        for (int k =startBlockIndex; k <= endBlockIndex; k++) { // iterate through blocks -> divide into small for loops per thread
            for (int i = 1; i < nm1; i++) {  // iterate in block
                index = k * n + i;
                actualIteration[index] = 1 / valMainDiag *
                                                      (f - valLowBlockDiag * lastIterSol[index-n] -
                                                       valLowMinDiag * lastIterSol[index - 1] -
                                                       valUpDiag * lastIterSol[index + 1] -
                                                       valUpBlockDiag * lastIterSol[index + n]);
            }
        }
}

void* subrJacobiThreadsBarrierCached(void* param) {
    messJacobiPrivate *m=(messJacobiPrivate *) param;
    // Sharing is caring, but caching is smashing.
    double *actualIteration=m->mJShared->actualIteration, *lastIterSol=m->mJShared->lastIterSol;
    double valLowMinDiag=m->mJShared->valLowMinDiag, valLowBlockDiag=m->mJShared->valLowBlockDiag, valMainDiag=m->mJShared->valMainDiag, valUpDiag=m->mJShared->valUpDiag,
            valUpBlockDiag=m->mJShared->valUpBlockDiag, f=m->mJShared->f;
    double *temp;
    int nm1=m->mJShared->n-1, n=m->mJShared->n;
    int startBlockIndex=m->startBlockIndex, endBlockIndex=m->endBlockIndex;
    int index;
    while(m->mJShared->iteration<m->mJShared->maxIteration) {
        for (int k = startBlockIndex;
             k <= endBlockIndex; k++) { // iterate through blocks -> divide into small for loops per thread
            for (int i = 1; i < nm1; i++) {  // iterate in block
                index = k * n + i;
                actualIteration[index] = 1 / valMainDiag *
                                                      (f - valLowBlockDiag *  lastIterSol[index - n] -
                                                       valLowMinDiag * lastIterSol[index - 1] -
                                                       valUpDiag * lastIterSol[index + 1] -
                                                       valUpBlockDiag *    lastIterSol[index + n]);
            }
        }
        pthread_barrier_wait(&iterationBarrier);
        if(m->threadId==0) {  // I'm the Master.
            m->mJShared->iteration++;
        }
        temp=actualIteration;
        actualIteration=lastIterSol;
        lastIterSol=temp;
        pthread_barrier_wait(&iterationBarrier);
    }
    if(m->threadId==0) {  // I'm the Master.
        m->mJShared->lastIterSol=lastIterSol;
        m->mJShared->actualIteration=actualIteration;
    }

}
// Dump, pls ignore everything below
void* subrJacobiThreadsBarrier(void* param) {  // deprecated
    messJacobiPrivate *m=(messJacobiPrivate *) param;
    int nm1=m->mJShared->n-1;
    int index;
    while(m->mJShared->iteration<m->mJShared->maxIteration) {
        for (int k = m->startBlockIndex;
             k <= m->endBlockIndex; k++) { // iterate through blocks -> divide into small for loops per thread
            for (int i = 1; i < nm1; i++) {  // iterate in block
                index = k * m->mJShared->n + i;
                m->mJShared->actualIteration[index] = 1 / m->mJShared->valMainDiag *
                                                      (m->mJShared->f - m->mJShared->valLowBlockDiag *
                                                                        m->mJShared->lastIterSol[index -
                                                                                                 m->mJShared->n] -
                                                       m->mJShared->valLowMinDiag *
                                                       m->mJShared->lastIterSol[index - 1] -
                                                       m->mJShared->valUpDiag * m->mJShared->lastIterSol[index + 1] -
                                                       m->mJShared->valUpBlockDiag *
                                                       m->mJShared->lastIterSol[index + m->mJShared->n]);
            }
        }
        pthread_barrier_wait(&iterationBarrier);
        if(m->threadId==0) {  // I'm the Master.
            m->mJShared->iteration++;
            double *temp=m->mJShared->actualIteration;
            m->mJShared->actualIteration=m->mJShared->lastIterSol;
            m->mJShared->lastIterSol=temp;
        }
        pthread_barrier_wait(&iterationBarrier);
    }
}
