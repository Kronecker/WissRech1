//
// Created by looky on 29.12.17.
//

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#define NUM_THREADS 8


void *testfunc(void *param){
    double *message=(double*)param;
    printf("Parameter: %f\n", *message);
}


void solveTask10() {
    int i, r;
    pthread_t tid, thread_id[NUM_THREADS];
    double message[NUM_THREADS];
    printf("Threads starten:\n");
    for(i=0; i<NUM_THREADS; i++){
        message[i] = 1.0/(i+1);
        r = pthread_create(&tid, NULL, testfunc, (void*) &message[i]);
        thread_id[i] = tid;
    }
    printf("Threads einsammeln:\n");
    for(i=0; i<NUM_THREADS; i++){
        r = pthread_join(thread_id[i], NULL);
    }
}















