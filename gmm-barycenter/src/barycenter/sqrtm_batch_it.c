#include "mex.h"
#include "matrix.h"
#include "cblas.h"
#include "lapacke.h"
#include <math.h>

/*
 S = sqrtm_batch_it(V, w, m);
 
 V is a d x d x n matrices
 w is a n length vector
 m is a d x d symmetric matrix
 S = \sum_{i=1}^n w(i) * sqrtm(m * V(:,:,i) * m)
 */

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    mwSize *dimensions;
    double *V, *w, *m, *S;
    
    dimensions = mxGetDimensions(prhs[0]);
    V = mxGetPr(prhs[0]);
    w = mxGetPr(prhs[1]);
    m = mxGetPr(prhs[2]);
    
    mwSize d, n;
    d=dimensions[0];
    n=dimensions[2];
    
    plhs[0] = mxCreateDoubleMatrix(d, d, mxREAL);
    S = mxGetPr(plhs[0]);
    
    double *Vtmp1, *Vtmp2, *Q, *L, *work;
    Vtmp1 = (double*) malloc(d*d*n*sizeof(double));
    Vtmp2 = (double*) malloc(d*d*sizeof(double));
    Q = (double*) malloc(d*d*sizeof(double));
    lapack_int lwork = 26*d, liwork=10*d, *iwork, suppz; 
    int info;
    L = (double*) malloc(d*3*sizeof(double));
    work = (double*) malloc(lwork*sizeof(double));
    iwork = (lapack_int*) malloc(liwork*sizeof(lapack_int));
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d*n, d, 1.0, m, d, V, d, 0.0, Vtmp1, d);
    mwSize i,j,k;
    lapack_int M;
    char cS='S';
    for (i=0; i<n; ++i) {
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, 1.0, Vtmp1+d*d*i, d, m, d, 0.0, Vtmp2, d);
        info=LAPACKE_dsyevr_work(LAPACK_COL_MAJOR, 'V', 'A', 
                                 'L', d, Vtmp2, 
                                 d, 0, 0, 
                                 0, 0, LAPACK_dlamch(&cS), 
                                 &M, L, Q, 
                                 d, &suppz, 
                                 work, lwork, 
                                 iwork, liwork);        
        for (j=0; j<d; ++j) {
            double lambda = sqrt(fabs(L[j]));
            for (k=0; k<d; ++k) {
                Vtmp2[j*d+k] = Q[j*d+k] * lambda;
            }
        }
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, d, d, d, w[i], Vtmp2, d, Q, d, 1.0, S, d);
    }
    
    free(Vtmp1); free(Vtmp2); free(Q); free(L); free(work); free(iwork); 
}