//
// Created by grabiger on 16.01.2018.
//

int globalMaxIteration=100;


void solveTask13() {












}


double* gradientForBlockDiag(double valLowBlockDiag,double valLowMinDiag,double valMainDiag, double valUpDiag,double valUpBlockDiag, int n, double f, double valBoundary, int procs) {

    double* r=new double[n*n]();
    double* a=new double[n*n]();

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








}










