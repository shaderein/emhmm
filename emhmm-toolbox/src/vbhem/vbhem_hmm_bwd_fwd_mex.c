/*==========================================================
 * this is a MEX version of vbhem_hmm_bwd_fwd.m
 *
 * for diagonal covariances:
 *   [LL_elbo ...
 *    sum_nu_1 ...
 *    update_emit_pr ...
 *    update_emit_mu ...
 *    update_emit_Mu  ...
 *    sum_xi] = vbhem_hmm_bwd_fwd(h3m_b,h3m_r,T,smooth,maxN,maxN2)
 *
 * for full covariances:
 *   ... = = vbhem_hmm_bwd_fwd(h3m_b,h3m_r,T,smooth,maxN,maxN2,logdetCovPlusDdivlamR, invCovR)
 *
 * ---
 * Eye-Movement analysis with HMMs (emhmm-toolbox)
 * Copyright (c) 2021-07-29
 * Antoni B. Chan, Janet H. Hsiao, Lan Hui
 * City University of Hong Kong, University of Hong Kong
 *========================================================*/

/* Version info
 *  2017/08/07 : ABC - initial version
 *  2017/08/18 : ABC - added full covariance matrices
 *  2018/05/15 : v0.72 - full C compatability for Linux
 *  2021/07/29 : v0.80 - initial version (adapted from hem_hmm_bwd_fwd)
 */

/* compile commands
 *
 * for debugging
 *   mex -g vbhem_hmm_bwd_fwd_mex.c
 *
 * optimized
 *   mex -g vbhem_hmm_bwd_fwd_mex.c
 */

/* Assumptions:
 *    ncentres = 1 for all g3m emissions 
 *   covar_type = 'diag' 
 */

#include "mex.h"
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef enum {COV_DIAG, COV_FULL} covar_type;

/*#define DEBUG     */
/*#define DEBUGMORE */

/* use pointer arithmetic for faster matrix indexing
   (compiler optimizations give similar speed) */
#define USEPTRS

/* safe functions for indexing */
int IX(int i, int j, int M, int N) {
  mxAssert( (i>=0) && (i<M) , "IX: i out of bounds");
  mxAssert( (j>=0) && (j<N) , "IX: j out of bounds");    
  return i+j*M;
}

int IX3(int i, int j, int k, int M, int N, int D) {
  mxAssert( (i>=0) && (i<M) , "IX3: i out of bounds");
  mxAssert( (j>=0) && (j<N) , "IX3: j out of bounds");    
  mxAssert( (k>=0) && (k<D) , "IX3: k out of bounds");
  return i+j*M+k*M*N;
}

int IX4(int i, int j, int k, int l, int M, int N, int K, int L) {
  mxAssert( (i>=0) && (i<M) , "IX4: i out of bounds");
  mxAssert( (j>=0) && (j<N) , "IX4: j out of bounds");    
  mxAssert( (k>=0) && (k<K) , "IX4: k out of bounds");    
  mxAssert( (l>=0) && (l<L) , "IX4: l out of bounds");
  return i + j*M + k*M*N + l*M*N*K;
}

double parseScalar(const mxArray *mx) {
  if( !mxIsDouble(mx) || mxIsComplex(mx) ) {
    mexErrMsgIdAndTxt("vbhmm_fb_mex:invalidinput", "arg must be scalar.");
  }
  return mxGetScalar(mx);
}

/* get MxN matrix */
/* pass M=0 or N=0 to ignore that dimension when checking dimensions */
double *parseMatrix(const mxArray *mx, int M, int N) {
  mwSize ndims;  
  size_t nrows;
  size_t ncols;

  
  if( !mxIsDouble(mx) || mxIsComplex(mx) )
    mexErrMsgIdAndTxt("vbhmm_fb_mex:invalidinput", "arg must be double.");
 
  ndims = mxGetNumberOfDimensions(mx);
  
#ifdef DEBUG
  printf("ndims=%d; ", ndims);
#endif
  
  if (ndims != 2)
    mexErrMsgIdAndTxt("vbhmm_fb_mex:invalidinput", "invalid num dimensions.");

  nrows = mxGetM(mx);
  ncols = mxGetN(mx);

#ifdef DEBUG
  printf("nr=%d; nc=%d;\n", nrows, ncols);
#endif

  if (((nrows != M) && (M!=0)) || ((ncols != N) && (N!=0)))
    mexErrMsgIdAndTxt("vbhmm_fb_mex:invalidinput", "invalid size.");

  return mxGetPr(mx);
}


/* get MxNxD matrix */
double *parseMatrix3(const mxArray *mx, int M, int N, int D) {
  mwSize ndims;
  const mwSize *tmpdims;
  
  if( !mxIsDouble(mx) || mxIsComplex(mx) )
    mexErrMsgIdAndTxt("vbhmm_fb_mex:invalidinput", "arg must be double.");

  ndims = mxGetNumberOfDimensions(mx);
  
#ifdef DEBUG
  printf("ndims=%d; ", ndims);
#endif
  
  if (ndims != 3)
    mexErrMsgIdAndTxt("vbhmm_fb_mex:invalidinput", "invalid num dimensions.");
  
  tmpdims = mxGetDimensions(mx);
  
#ifdef DEBUG
  printf("nr=%d; nc=%d; np=%d;\n", tmpdims[0], tmpdims[1], tmpdims[2]);
#endif

  if ((tmpdims[0] != M) || (tmpdims[1] != N) || (tmpdims[2] != D))
    mexErrMsgIdAndTxt("vbhmm_fb_mex:invalidinput", "invalid size.");

  return mxGetPr(mx);
}

/* create 3D matrix */
mxArray *create3Dmatrix(int D1, int D2, int D3) {
  mwSize tmp_dims[3] = {D1, D2, D3};
  return mxCreateNumericArray(3, tmp_dims,  mxDOUBLE_CLASS, mxREAL);
}
    
/* debugging output */
void dumpMatrix(const char *name, const double *p, int M, int N) {
  int i, j;
  printf("-- %s --\n", name);
  for(i=0; i<M; i++) {
    for (j=0; j<N; j++) {
      printf("%0.5f ", p[IX(i,j,M,N)]);
    }
    printf("\n");
  }
}

void dumpCovMatrix(const char *name, const double *p, int dim, covar_type covmode) {
  switch(covmode) {
    case COV_DIAG:
      dumpMatrix(name, p, 1, dim);
      break;
    case COV_FULL:
      dumpMatrix(name, p, dim, dim);
      break;
  }
}

void dumpMatrix3(const char *name, const double *p, int M, int N, int D) {
  int i, j, d;
  for (d=0; d<D; d++) {
    printf("-- %s(:,:,%d) --\n", name, d);
    for(i=0; i<M; i++) {
      for (j=0; j<N; j++) {
        printf("%0.5f ", p[IX3(i,j,d,M,N,D)]);
      }
      printf("\n");
    }
  }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s] = logtrick(lA)
% logtrick - "log sum trick" - calculate log(sum(A)) using only log(A) 
%
%   s = logtrick(lA)
%
%   lA = column vector of log values
%
%   if lA is a matrix, then the log sum is calculated over each column
% 

[mv, mi] = max(lA, [], 1);
temp = lA - repmat(mv, size(lA,1), 1);
cterm = sum(exp(temp),1);
s = mv + log(cterm);
end
*/
void logtrick(double *lA, int M, int N, double *s) {
  int i, j;
  double mv, tmps;
  
#ifndef USEPTRS
  /* for each column */
  for(j=0; j<N; j++) {
    /* find max of column*/
    /* [mv, mi] = max(lA, [], 1); */
    mv = lA[IX(0,j,M,N)];
    for(i=1; i<M; i++)
      if (lA[IX(i,j,M,N)] > mv)
        mv = lA[IX(i,j,M,N)];
    
    /*
     * temp = lA - repmat(mv, size(lA,1), 1);
     * cterm = sum(exp(temp),1);
     * s = mv + log(cterm);
     */
    tmps = 0.0;
    for(i=0; i<M; i++) {
      tmps += exp(lA[IX(i,j,M,N)] - mv);
    }
    s[j] = mv+log(tmps);
  }
  
#else    
  
  double *ptr1, *ptr2, *ptr3;
  
  /* for each column */
  ptr2 = s;
  ptr3 = lA;
  for(j=0; j<N; j++, ++ptr2, ptr3+=M) {
    
    /* find max of column*/
    /* [mv, mi] = max(lA, [], 1); */
    ptr1 = ptr3;
    mv = *ptr1;
    ++ptr1;
    for(i=1; i<M; i++, ++ptr1)
      if (*ptr1 > mv)
        mv = *ptr1;
    
    /*
     * temp = lA - repmat(mv, size(lA,1), 1);
     * cterm = sum(exp(temp),1);
     * s = mv + log(cterm);
     */
    tmps = 0.0;
    ptr1 = ptr3;
    for(i=0; i<M; i++, ++ptr1) {
      tmps += exp(*ptr1 - mv);
    }
    *ptr2 = mv+log(tmps);
  }
#endif
   
}



/*
function [LL_elbo ...
                sum_nu_1 ...
                update_emit_pr ...
                update_emit_mu ...
                update_emit_Mu  ...
                sum_xi] = vbhem_hmm_bwd_fwd(h3m_b,h3m_r,T,smooth,maxN,maxN2)
 
function [LL_elbo ...
                sum_nu_1 ...
                update_emit_pr ...
                update_emit_mu ...
                update_emit_Mu  ...
                sum_xi] = vbhem_hmm_bwd_fwd(h3m_b,h3m_r,T,smooth,maxN,maxN2,logdetCovPlusDdivlamR, invCovR)
*/
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  covar_type covmode;
  
  /* size constants */
  int T, Kr, Kb, N, N2, dim, maxN, maxN2, dim2;
  
  
  /* input matrices and values */
  const mxArray *h3m_b, *h3m_r;
  const mxArray *logdetCovPlusDdivlamR, *invCovR;
  
  
  /* output matrices */
  double *LL_elbo;
  mxArray *out_sum_nu_1, *out_update_emit_pr, *out_update_emit_mu, *out_update_emit_Mu, *out_sum_xi;
  double *sum_nu_1, *sum_t_sum_g_xi, *update_emit_pr, *update_emit_mu, *update_emit_Mu; /* pointers to output */
  double *sum_t_nu;

  /* working variables */
  const mxArray *hmm_r, *hmm_b, *mxtmp, *mxtmp2;
  const double **h3m_r_A, **h3m_r_prior, **h3m_r_emit_m,  **h3m_r_emit_W, **h3m_r_emit_v, **h3m_r_loglambdaTildePDdivlam;
  const double **h3m_b_A, **h3m_b_prior, **h3m_b_emit_centres, **h3m_b_emit_covars;
  const double **h3m_r_invCov,  **h3m_r_logdetCovPlusDdivlam; 
  const double *Ar, *Ab;
  int *h3m_r_N2, *h3m_b_N;
  double *LLG_elbo, *sum_w_pr, *sum_w_mu, *sum_w_Mu; 
  double ELLs, tmpd, tmpd2;
  double *Theta, *LL_old, *LL_new, *logtheta, *logsumtheta, *Theta_1;
  double *nu, *foo, *xi_foo;
  double *tmpx;
  mxArray *tmpmx;
  
  
  /* working pointers */
#ifdef USEPTRS
  double *ptr1, *ptr2, *ptr3, *ptr4;
  const double *cptra,*cptrb, *cptrc, *cptrd, *cptre, *cptrf, *cptrx,*cptry;  
  const double **cpptra, **cpptrb, **cpptrc, **cpptrd, **cpptre, **cpptrf, **cpptrx, **cpptry;
#endif 
  
  /* counters */
  int n, k, t, ti, tj, tk;
  int i_h3m, j_h3m, rho, beta, sigma;
  
/* check for proper number of arguments */
  if((nrhs!=5) && (nrhs!=7)) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","5 or 7 inputs required.");
  }
  if(nlhs!=6) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","6 output required.");
  }

  if (nrhs == 7)
    covmode = COV_FULL;
  else
    covmode = COV_DIAG;
  
  /*
   * Check and Read Inputs:
   *  h3m_b - {1xKb}
   *  h3m_r - {1xKr}
   *    T    - scalar
   *  smooth - scalar
   *  maxN   - max N (states) in h3m_b
   *  maxN2  - max N2 (states) in h3m_r
   *  logdetCovPlusDdivlamR - {1xKr} [1 x N2]           (full cov only)
   *  invCovR    - {1xKr} [dim x dim x N2]   (full cov only)
   */
  
  if( !mxIsCell(prhs[0]) ) {
    mexErrMsgIdAndTxt("vbhem_mex:invalidinput", "1st arg must be cell");
  }
  h3m_b = prhs[0];
  
  if( !mxIsCell(prhs[1]) ) {
    mexErrMsgIdAndTxt("vbhem_mex:invalidinput", "2nd arg must be cell");
  }
  h3m_r = prhs[1];
  
  Kr = (int) mxGetNumberOfElements(h3m_r);
  Kb = (int) mxGetNumberOfElements(h3m_b);
  T       = (int) parseScalar(prhs[2]);
  maxN    = (int) parseScalar(prhs[3]);
  maxN2   = (int) parseScalar(prhs[4]);
  if (covmode == COV_FULL) {
    logdetCovPlusDdivlamR = prhs[5];
    invCovR    = prhs[6];
  }

  
#ifdef DEBUG
  printf("Kb=%d; Kr=%d; T=%d;  maxN=%d; maxN2=%d; covmode=%d\n", 
          Kb, Kr, T,  maxN, maxN2, covmode);
#endif
  
  /*
   * Create Outputs:
   * LL_elbo  [Kb x Kr]
   * sum_nu_1  {Kb x Kr}
   * update_emit_pr {Kb x Kr}
   * update_emit_mu {Kb x Kr}
   * update_emit_Mu {Kb x Kr}
   * sum_xi  {Kb x Kr}
   */

  /* create the output matrix */
  plhs[0] = mxCreateDoubleMatrix(Kb, Kr, mxREAL);
  plhs[1] = mxCreateCellMatrix(Kb, Kr);
  plhs[2] = mxCreateCellMatrix(Kb, Kr);
  plhs[3] = mxCreateCellMatrix(Kb, Kr);
  plhs[4] = mxCreateCellMatrix(Kb, Kr);
  plhs[5] = mxCreateCellMatrix(Kb, Kr);
 
  
  LL_elbo = mxGetPr(plhs[0]);
  out_sum_nu_1 = plhs[1];
  out_update_emit_pr  = plhs[2];
  out_update_emit_mu  = plhs[3];
  out_update_emit_Mu  = plhs[4];
  out_sum_xi          = plhs[5];
  

  
  /** cache the matrices into an array of pointers for easier processing */
  h3m_r_N2 = mxMalloc(sizeof(int)*Kr);                      /* [1 x Kr] */
  h3m_r_A  = mxMalloc(sizeof(double *)*Kr);                 /* [1 x Kr] */
  h3m_r_prior = mxMalloc(sizeof(double *)*Kr);              /* [1 x Kr] */
  h3m_r_loglambdaTildePDdivlam = mxMalloc(sizeof(double *)*maxN2*Kr);  /* [maxN2 x Kr] */
  h3m_r_emit_m = mxMalloc(sizeof(double *)*maxN2*Kr);  /* [maxN2 x Kr] */
  h3m_r_emit_W = mxMalloc(sizeof(double *)*maxN2*Kr);  /* [maxN2 x Kr] */
  h3m_r_emit_v = mxMalloc(sizeof(double *)*maxN2*Kr);  /* [maxN2 x Kr] */
  h3m_b_N  = mxMalloc(sizeof(int)*Kb);                      /* [1 x Kb] */
  h3m_b_A  = mxMalloc(sizeof(double *)*Kb);                 /* [1 x Kb] */
  h3m_b_prior = mxMalloc(sizeof(double *)*Kb);              /* [1 x Kb] */
  h3m_b_emit_centres = mxMalloc(sizeof(double *)*maxN*Kb);  /* [maxN x Kb] */
  h3m_b_emit_covars  = mxMalloc(sizeof(double *)*maxN*Kb);  /* [maxN x Kb] */
  if (covmode == COV_FULL) {
    h3m_r_invCov    = mxMalloc(sizeof(double *)*Kr);  /* [1 x Kr] */
    h3m_r_logdetCovPlusDdivlam = mxMalloc(sizeof(double *)*Kr);  /* [1 x Kr] */
  }

  
  /* reduced HMMs */
  for(j_h3m=0; j_h3m<Kr; j_h3m++) {
    hmm_r = mxGetCell(h3m_r, j_h3m);
    mxtmp = mxGetField(hmm_r, 0, "logATilde");
    N2                 = (int) mxGetM(mxtmp);  /* num states in hmm_r */
    h3m_r_N2[j_h3m]    = N2;
    h3m_r_A[j_h3m]     = mxGetPr(mxtmp);
    h3m_r_prior[j_h3m] = mxGetPr(mxGetField(hmm_r, 0, "logPiTilde"));
    mxtmp = mxGetField(hmm_r, 0, "emit");
    for(k=0; k<N2; k++) {
      mxtmp2 = mxGetCell(mxtmp, k);
      h3m_r_emit_m[IX(k,j_h3m,maxN2,Kr)] = mxGetPr(mxGetField(mxtmp2, 0, "m"));
      h3m_r_emit_W[IX(k,j_h3m,maxN2,Kr)]  = mxGetPr(mxGetField(mxtmp2, 0, "W")); 
      h3m_r_emit_v[IX(k,j_h3m,maxN2,Kr)]  = mxGetPr(mxGetField(mxtmp2, 0, "v")); 
      h3m_r_loglambdaTildePDdivlam[IX(k,j_h3m,maxN2,Kr)] = mxGetPr(mxGetField(mxtmp2, 0, "logLambdaTildePlusDdivlamda"));
      dim = (int) mxGetN(mxGetField(mxtmp2, 0, "m"));
    }
    
    dim2 = dim*dim;
    
    if (covmode == COV_FULL) {
       h3m_r_logdetCovPlusDdivlam[j_h3m] = mxGetPr(mxGetCell(logdetCovPlusDdivlamR, j_h3m));
       h3m_r_invCov[j_h3m]    = mxGetPr(mxGetCell(invCovR, j_h3m));
       dim2 = dim*dim;
    }
  }
  /* base HMMs */
  for(i_h3m=0; i_h3m<Kb; i_h3m++) {
    hmm_b = mxGetCell(h3m_b, i_h3m);    
    mxtmp = mxGetField(hmm_b, 0, "A");
    N     = (int) mxGetM(mxtmp);  /* num states in hmm_b */
    h3m_b_N[i_h3m] = N;    
    h3m_b_A[i_h3m]     = mxGetPr(mxtmp);
    h3m_b_prior[i_h3m] = mxGetPr(mxGetField(hmm_b, 0, "prior"));
    mxtmp = mxGetField(hmm_b, 0, "emit");
    for(k=0; k<N; k++) {
      mxtmp2 = mxGetCell(mxtmp, k);
      h3m_b_emit_centres[IX(k,i_h3m,maxN,Kb)] = mxGetPr(mxGetField(mxtmp2, 0, "centres"));
      h3m_b_emit_covars[IX(k,i_h3m,maxN,Kb)]  = mxGetPr(mxGetField(mxtmp2, 0, "covars"));
      dim = (int) mxGetScalar(mxGetField(mxtmp2, 0, "nin"));
    }
  }
  
  
  /* working variables */
  /* changes with the hmm_b and hmm_r -- the used sizes will be different */
  LLG_elbo = mxMalloc(sizeof(double)*maxN*maxN2);         /* [maxN x maxN2] --> [N x N2] */
  sum_w_pr = mxMalloc(sizeof(double)*maxN*maxN2);         /* [maxN x maxN2] --> [N x N2] */  
  sum_w_mu = mxMalloc(sizeof(double)*maxN*dim*maxN2);     /* [maxN x dim x maxN2] --> [N x dim x N2] */  
  Theta    = mxMalloc(sizeof(double)*maxN2*maxN2*maxN*T); /* [maxN2 x maxN2 x maxN x T] --> [N2 x N2 x N x T] */
  LL_new   = mxMalloc(sizeof(double)*maxN*maxN2);         /* [maxN x maxN2] --> [N x N2] */
  LL_old   = mxMalloc(sizeof(double)*maxN*maxN2);         /* [maxN x maxN2] --> [N x N2] */
  logtheta = mxMalloc(sizeof(double)*maxN2*maxN);         /* [maxN2 x maxN] --> [N2 x N] */
  logsumtheta = mxMalloc(sizeof(double)*1*maxN);          /* [1 x maxN]    --> [1 x N] */
  Theta_1  = mxMalloc(sizeof(double)*maxN2*maxN);         /* [maxN2 x maxN] --> [N2 x N] */
  nu       = mxMalloc(sizeof(double)*maxN2*maxN);         /* [maxN2 x maxN] --> [N2 x N] */
  foo      = mxMalloc(sizeof(double)*maxN2*maxN);         /* [maxN2 x maxN] --> [N2 x N] */
  xi_foo   = mxMalloc(sizeof(double)*maxN2*maxN);         /* [maxN2 x maxN] --> [N2 x N] */
  sum_t_nu = mxMalloc(sizeof(double)*maxN2*maxN);         /* [maxN2 x maxN] --> [N2 x N] */
  
  tmpx     = mxMalloc(sizeof(double)*dim);                 /* [dim x 1] */
  
  switch(covmode) {
    case COV_DIAG:
      sum_w_Mu = mxMalloc(sizeof(double)*maxN*dim*maxN2);     /* [maxN x dim x maxN2] --> [N x dim x N2] */
      break;
    case COV_FULL:
      sum_w_Mu = mxMalloc(sizeof(double)*maxN*dim*dim*maxN2);  /* [maxN x dim x dim x maxN2] --> [N x dim x dim x N2] */
      tmpx     = mxMalloc(sizeof(double)*dim);                 /* [dim x 1] */
      break;
  }
      
  /** loop over Kr and Kb **/
  for(j_h3m=0; j_h3m<Kr; j_h3m++) {
    N2 = h3m_r_N2[j_h3m];
    
#ifdef DEBUGMORE
    printf("== j_h3m=%d ========================================\n", j_h3m);
    printf("N2=%d; dim=%d\n", N2, dim);
    dumpMatrix("hmm_r_A", h3m_r_A[j_h3m], N2, N2);
    dumpMatrix("hmm_r_prior", h3m_r_prior[j_h3m], 1, N2);
    for (k=0; k<N2; k++) {
      dumpMatrix("hmm_r_emit_centres", h3m_r_emit_m[IX(k,j_h3m,maxN2,Kr)], 1, dim);
      dumpCovMatrix("hmm_r_emit_covars",  h3m_r_emit_W[IX(k,j_h3m,maxN2,Kr)], dim, covmode);
    }
#endif
    
    /* loop over base mixture components */
    for(i_h3m=0; i_h3m<Kb; i_h3m++) {
      N = h3m_b_N[i_h3m];
          
#ifdef DEBUGMORE
      printf("== j_h3m=%d, i_h3m=%d ========================================\n", j_h3m, i_h3m);
      printf("N=%d; dim=%d\n", N, dim);
      dumpMatrix("hmm_b_A", h3m_b_A[i_h3m], N, N);
      dumpMatrix("hmm_b_prior", h3m_b_prior[i_h3m], 1, N);
      for (k=0; k<N; k++) {
        dumpMatrix("hmm_b_emit_centres", h3m_b_emit_centres[IX(k,i_h3m,maxN,Kb)], 1, dim);
        dumpCovMatrix("hmm_b_emit_covars",  h3m_b_emit_covars[IX(k,i_h3m,maxN,Kb)], dim, covmode);
      }
#endif
      
      /*
       * % first, get the elbo of the E_gauss(b) [log p (y | gauss (r))]  (different b on different rows)
       * % and the sufficient statistics for the later updates
       * % sum_w_pr is a N2 by 1 cell of N by M
       * % sum_w_mu is a N2 by 1 cell of N by dim by M
       * % sum_w_Mu is a N2 by 1 cell of N by dim by M (for the diagonal case)
       * %               N2 by 1 cell of N by dim by dim by M (for the full case)
       *
       *[LLG_elbo sum_w_pr sum_w_mu sum_w_Mu] = g3m_stats(hmm_b.emit,hmm_r.emit, h3m_b_hmm_emit_cache);
       */
      
      /*** BEGIN g3m_stats ***********************************************************************************/
      /* Our variables are different (since we assume M=1):
       * LLG_elbo = [N x N2]
       * sum_w_pr = [N x N2]
       * sum_w_mu = [N x dim x N2]
       * sum_w_Mu = [N x dim x N2] (diagonal)
       *            [N x dim x dim x N2] (full)
       */
      
#ifndef USEPTRS
      /*
       * for rho = 1 : N2
       */
      for(rho=0; rho<N2; rho++) {
        /*
         * [UNNEEDED]
         * foo_sum_w_pr = zeros(N,M);
         * foo_sum_w_mu = zeros(N,dim,M);
         * foo_sum_w_Mu = zeros(N,dim,M);
         *
         * gmmR = g3m_r{rho};
         *
         * for beta = 1 : N       
         */
        for(beta=0; beta<N; beta++) {
          /* 
           * gmmB = g3m_b{beta};
           *
           * % compute the expected log-likelihood between the Gaussian components
           * % i.e., E_M(b),beta,m [log p(y | M(r),rho,l)], for m and l 1 ...M
           * [ELLs] =  compute_exp_lls(gmmB,gmmR);
           */
          
          switch(covmode) {
            case COV_DIAG:
              /* (note here: 1st arg gmmA=gmmB, 2nd arg gmmB=gmmR)
               * ELLs(a,b) = -.5 * ( ...
               *      dim*log(2*pi) + sum( log(gmmB.covars(b,:)),2) ...
               *     +  sum(gmmA.covars(a,:)./gmmB.covars(b,:),2)  ...
               *     + sum(((gmmA.centres(a,:) - gmmB.centres(b,:)).^2) ./ gmmB.covars(b,:),2) ...
               *     );
               */
              ELLs = dim*log(2.0*M_PI);
              for (ti=0; ti<dim; ti++) {
                ELLs += log(h3m_r_emit_W[IX(rho,j_h3m,maxN2,Kr)][ti]); /* [this can be precomputed] */
                ELLs += h3m_b_emit_covars[IX(beta,i_h3m,maxN,Kb)][ti] / h3m_r_emit_W[IX(rho,j_h3m,maxN2,Kr)][ti];
                tmpd = h3m_b_emit_centres[IX(beta,i_h3m,maxN,Kb)][ti] - h3m_r_emit_m[IX(rho,j_h3m,maxN2,Kr)][ti];
                ELLs += tmpd*tmpd / h3m_r_emit_W[IX(rho,j_h3m,maxN2,Kr)][ti];
              }
              ELLs *= -0.5;
              break;
              
            case COV_FULL:
              /*
               *  inv_covR = inv(gmmR.covars(:,:,b));
               * 
               *  ELLs(a,b) = -.5 * ( ...
               *     dim*log(2*pi) + log(det(gmmR.covars(:,:,b))) ...
               *      + trace( inv_covR * gmmB.covars(:,:,a) ) +  ...
               *      + (gmmB.centres(a,:) - gmmR.centres(b,:)) * inv_covR * (gmmB.centres(a,:) - gmmR.centres(b,:))' ...
               *    );
               */
              ELLs = dim*log(2.0*M_PI);
              ELLs += h3m_r_logdetCovPlusDdivlam[j_h3m][rho];
              
              /* trace term (use the fact that covars are symmetric) */
              for (ti=0; ti<dim; ti++)
                for (tj=0; tj<dim; tj++)
                  ELLs += h3m_r_invCov[j_h3m][IX3(ti,tj,rho,dim,dim,N2)]*
                          h3m_b_emit_covars[IX(beta,i_h3m,maxN,Kb)][IX(ti,tj,dim,dim)];             
              
              /* mahal term */
              for (ti=0; ti<dim; ti++)
                tmpx[ti] = h3m_b_emit_centres[IX(beta,i_h3m,maxN,Kb)][ti] - h3m_r_emit_m[IX(rho,j_h3m,maxN2,Kr)][ti];
              for (tj=0; tj<dim; tj++) {
                tmpd = 0.0;
                for (ti=0; ti<dim; ti++) {
                  tmpd += tmpx[ti]*h3m_r_invCov[j_h3m][IX3(ti,tj,rho,dim,dim,N2)];
                }
                ELLs += tmpd*tmpx[tj];
              }
              ELLs *= -0.5;

              break;
              
          }
                    
          /* [UNNEEDED since M=1, prior=1]
           * % compute log(omega_r) + E_M(b),beta,m [log p(y | M(r),rho,l)]
           * log_theta = ELLs + ones(M,1) * log(gmmR.priors);
      
           * % compute log Sum_b omega_b exp(-D(fa,gb))
           * log_sum_theta = logtrick(log_theta')';
           */
          
          /*
           * % compute L_variational(M(b)_i,M(r)_j) = Sum_a pi_a [  log (Sum_b omega_b exp(-D(fa,gb)))]
           * LLG_elbo(beta,rho) = gmmB.priors * log_sum_theta;
           */
          LLG_elbo[IX(beta,rho,N,N2)] = ELLs;
      
          /* [UNNEEDED - since M=1, priors=1 ==> theta=1.
           * % cache theta
           * theta = exp(log_theta -  log_sum_theta * ones(1,M));
           *
           * % aggregate in the output ss
           *
           * foo_sum_w_pr(beta,:) = gmmB.priors * theta;
           *
           * %          for l = 1 :  M
           * %
           * %              foo_sum_w_mu(beta,:,l) = foo_sum_w_mu(beta,:,l) + (gmmB.priors .* theta(:,l)') * gmmB.centres;
           * %
           * %          end
           * A
           * foo_sum_w_mu(beta,:,:) =  ( ((ones(M,1) * gmmB.priors) .* theta') * gmmB.centres )';
           *
           * foo_sum_w_Mu(beta,:,:) =  ( ((ones(M,1) * gmmB.priors) .* theta') * (gmmB.centres.^2 + gmmB.covars) )';
           */
          
          /* 
           * sum_w_pr{rho} = foo_sum_w_pr;
           * sum_w_mu{rho} = foo_sum_w_mu;
           * sum_w_Mu{rho} = foo_sum_w_Mu;
           */

          /* [TODO: these are unnecessary, and can be used directly later]; */
          /* Since M=1, prior=1, theta=1
           * foo_sum_w_pr(beta,:) = 1
           */
          sum_w_pr[IX(beta,rho,N,N2)] = 1;
          
          switch(covmode) {
            case COV_DIAG:
              /*
               * foo_sum_w_mu(beta,:) = gmmB.centres
               * foo_sum_w_Mu(beta,:) = gmmB.centres.^2 + gmmB.covars
               */
              /* [TODO: these are unnecessary, and can be used directly later]; */
              for(ti=0; ti<dim; ti++) {
                tmpd = h3m_b_emit_centres[IX(beta,i_h3m,maxN,Kb)][ti];
                sum_w_mu[IX3(beta,ti,rho,N,dim,N2)] = tmpd;
                sum_w_Mu[IX3(beta,ti,rho,N,dim,N2)] = tmpd*tmpd+h3m_b_emit_covars[IX(beta,i_h3m,maxN,Kb)][ti];
                /* [above can be precomputed] */
              }
              break;
              
            case COV_FULL:
              /*
               * foo_sum_w_mu(beta,:) = gmmB.centres
               *
               * tmp = gmmB.covars;
               * for m=1:M
               *   tmp(:,:,m) = tmp(:,:,m) + gmmB.centres(m,:)'*gmmB.centres(m,:);
               * end
               */
              for(ti=0; ti<dim; ti++) {
                tmpd = h3m_b_emit_centres[IX(beta,i_h3m,maxN,Kb)][ti];
                sum_w_mu[IX3(beta,ti,rho,N,dim,N2)] = tmpd;
                for(tj=0; tj<dim; tj++) {
                  tmpd2 = h3m_b_emit_centres[IX(beta,i_h3m,maxN,Kb)][tj];
                  sum_w_Mu[IX4(beta,ti,tj,rho,N,dim,dim,N2)] = tmpd*tmpd2+h3m_b_emit_covars[IX(beta,i_h3m,maxN,Kb)][IX(ti,tj,dim,dim)];
                }
              }
              
          }
          
        } /* end for beta */
      } /* end for rho    */
      
#else
      switch(covmode) {
        /** DIAGONAL ****************************************************/
        case COV_DIAG:
            
          
          cpptra = h3m_r_emit_W + j_h3m*maxN2;
          cpptrb = h3m_r_emit_m + j_h3m*maxN2;

          cpptrx = h3m_r_loglambdaTildePDdivlam + j_h3m*maxN2;
          cpptry = h3m_r_emit_v + j_h3m*maxN2;

          
          cpptre = h3m_b_emit_covars + i_h3m*maxN;
          cpptrf = h3m_b_emit_centres + i_h3m*maxN;

          ptr1 = LLG_elbo;
          ptr2 = sum_w_pr;

          
          for(rho=0; rho<N2; rho++, ++cpptra, ++cpptrb, ++cpptrx, ++cpptry) {
            cpptrc = cpptre;
            cpptrd = cpptrf;
            for(beta=0; beta<N; beta++, ++cpptrc, ++cpptrd, ++ptr1, ++ptr2) {
                            
          
              cptrx = *cpptrx;

                            
              ELLs = dim*log(2.0*M_PI) + *cptrx;

              cptra = *cpptra;
              cptrb = *cpptrb;
              cptrc = *cpptrc;
              cptrd = *cpptrd;
              cptry = *cpptry;

              
              for (ti=0; ti<dim; ti++, ++cptra, ++cptrb, ++cptrc, ++cptrd) {
               
                ELLs += *cptry * (*cptrc) * (*cptra);
                tmpd = *cptrd - *cptrb;
                ELLs += *cptry * (tmpd*tmpd) * (*cptra);
              }
              ELLs *= -0.5;                            
              

              *ptr1 = ELLs;
              *ptr2 = 1;


              
              cptrc = *cpptrc;
              cptrd = *cpptrd;
              ptr3 = sum_w_mu + beta + N*dim*rho;
              ptr4 = sum_w_Mu + beta + N*dim*rho;
              for(ti=0; ti<dim; ti++, ++cptrd, ++cptrc, ptr3+=N, ptr4+=N) {
                tmpd = *cptrd;
                *ptr3 = tmpd;
                *ptr4 = tmpd*tmpd + *cptrc;
              }
              
            } /* end for beta */
          } /* end for rho */
          break;
          
        case COV_FULL:
          /** FULL ****************************************************/
          cptra  = *(h3m_r_invCov + j_h3m);
          cpptrb = h3m_r_emit_m + j_h3m*maxN2;
          cpptre = h3m_b_emit_covars  + i_h3m*maxN;
          cpptrf = h3m_b_emit_centres + i_h3m*maxN;
          cptrf  = *(h3m_r_logdetCovPlusDdivlam + j_h3m);
          ptr1 = LLG_elbo;
          ptr2 = sum_w_pr;
          
          for(rho=0; rho<N2; rho++, cptra+=dim2, ++cptrf, ++cpptrb) {
            cpptrc = cpptre; /* h3m_b_emit_covars[0,i_h3m] */
            cpptrd = cpptrf; /* h3m_b_emit_centres[0,i_h3m] */
            for(beta=0; beta<N; beta++, ++cpptrc, ++cpptrd, ++ptr1, ++ptr2) {
              
              ELLs = dim*log(2.0*M_PI);
              
              /* ELLs += h3m_r_logdetCovPlusDdivlam[j_h3m][rho]; */
              ELLs += *cptrf;
                            
              /* trace term (use the fact that covars are symmetric) */
              /* for (ti=0; ti<dim; ti++)
               *   for (tj=0; tj<dim; tj++)
               *     ELLs += h3m_r_invCov[j_h3m][IX3(ti,tj,rho,dim,dim,N2)]*
               *       h3m_b_emit_covars[IX(beta,i_h3m,maxN,Kb)][IX(ti,tj,dim,dim)];
              */              
              cptrb = cptra;    /* h3m_r_invCov[:,:,rho] */
              cptrc = *cpptrc;  /* h3m_b_emit_covars[rho,i_h3m] */
              for (ti=0; ti<dim; ti++)
                for (tj=0; tj<dim; tj++, ++cptrb, ++cptrc)
                  ELLs += (*cptrb)*(*cptrc);             
              
              /* mahal term */
              /* for (ti=0; ti<dim; ti++)
               *   tmpx[ti] = h3m_b_emit_centres[IX(beta,i_h3m,maxN,Kb)][ti] - h3m_r_emit_m[IX(rho,j_h3m,maxN2,Kr)][ti];
              */
              ptr3 = tmpx;
              cptrc = *cpptrd;
              cptrd = *cpptrb;
              for (ti=0; ti<dim; ti++, ++ptr3, ++cptrc, ++cptrd)
                *ptr3 = *cptrc - *cptrd;
              
              /* for (tj=0; tj<dim; tj++) {
               *   tmpd = 0.0;
               *   for (ti=0; ti<dim; ti++) {
               *     tmpd += tmpx[ti]*h3m_r_invCov[j_h3m][IX3(ti,tj,rho,dim,dim,N2)];
               *   }
               *   ELLs += tmpd*tmpx[tj];
               * }
               */
              
              cptrb = tmpx;
              cptrd = cptra;                
              for (tj=0; tj<dim; tj++, ++cptrb) {
                cptrc = tmpx;
                tmpd = 0.0;
                for (ti=0; ti<dim; ti++, ++cptrc, ++cptrd) {
                  tmpd += (*cptrc) * (*cptrd);
                }
                ELLs += tmpd * (*cptrb);
              }
              
              ELLs *= -0.5;
              
              *ptr1 = ELLs;
              *ptr2 = 1;
              
              cptrc = *cpptrc;  /* h3m_b_emit_covars[rho,i_h3m] */
              cptrd = *cpptrd;  /* h3m_b_emit_centres[rho,i_h3m] */
              ptr3 = sum_w_mu + beta + N*dim*rho;
              ptr4 = sum_w_Mu + beta + N*dim*dim*rho;
              for(ti=0; ti<dim; ti++, ++cptrd, ptr3+=N) {
                tmpd = *cptrd;
                *ptr3 = tmpd;
                
                cptre = *cpptrd;
                for(tj=0; tj<dim; tj++, ++cptre, ++cptrc, ptr4+=N) {
                  *ptr4 = tmpd * (*cptre) + (*cptrc);
                }
              }             
              
            } /* end for beta */
          } /* end for rho */
          break;
      }
       

      
#endif
      
      /*** END g3m_stats ***********************************************************************************/
      
#ifdef DEBUGMORE
      dumpMatrix("LLG_elbo", LLG_elbo, N, N2); /* OK */
#endif
      
 
      
      /*
       * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       * %%%% do the backward recursion %%%%
       * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       *
       * Ab = hmm_b.A;
       * Ar = hmm_r.A;
       */
      Ab = h3m_b_A[i_h3m];
      Ar = h3m_r_A[j_h3m];      
      
      /* % allocate Theta, i.e., the VA parameters in the form of a CMC
       * Theta = zeros(N2,N2,N,T); % the dimensions refer to rho sigma (states of M(r)) gamma (state of M(b)) and time
       * % rho sigma gamma t
       *
       * % allocate the log-likelihood
       * LL_old = zeros(N,N2);       % the dimensions refer to gamma (state of M(b)) and sigma (states of M(r))
       */
#ifndef USEPTRS
      for(rho=0; rho<N2; rho++)
        for(beta=0; beta<N; beta++)
          LL_old[IX(beta,rho,N,N2)] = 0.0;
#else
      ptr1 = LL_old;
      for(rho=0; rho<N2; rho++)
        for(beta=0; beta<N; beta++, ++ptr1)
          *ptr1 = 0.0;
#endif
      
      /* % the log-likelihood can be intialized by all zeros ...
       *
       * for t = T : -1 : 2
       *
       * % [N,N2]
       * LL_new = zeros(size(LL_old));
       */
      for(t=T-1;t>=1; t--) {
        /*
         * %% old code
         * for rho = 1 : N2
         */
        for(rho=0; rho<N2; rho++) {
          /* % [N2xN] =          [1xN2]' x [1xN]      + [N,N2]'   +  [N,N2]'
           * logtheta = (log(Ar(rho,:))' * ones(1,N)) + LLG_elbo' + LL_old';
           */
#ifndef USEPTRS
          for(int tj=0; tj<N2; tj++) {
            for(int ti=0; ti<N; ti++) {
              logtheta[IX(tj,ti,N2,N)] = log(Ar[IX(rho,tj,N2,N2)])
                           + LLG_elbo[IX(ti,tj,N,N2)] + LL_old[IX(ti,tj,N,N2)];
            }
          }
#else
          cptra = Ar + rho;
          ptr2 = LLG_elbo;
          ptr3 = LL_old;
          for(tj=0; tj<N2; tj++, cptra+=N2) {
            ptr1 = logtheta + tj;
            for(ti=0; ti<N; ti++, ++ptr2, ++ptr3, ptr1+=N2) {
              *ptr1 = *cptra + *ptr2 + *ptr3;
            }
          }
#endif

          /* % [1xN]
           * logsumtheta = logtrick(logtheta);
           */
          logtrick(logtheta, N2, N, logsumtheta);
#ifdef DEBUGMORE
          dumpMatrix("logtheta", logtheta, N2, N);
          dumpMatrix("logsumtheta", logsumtheta, 1, N);
#endif
          
          /* % [Nx1]       =  [NxN] * [1xN]'
           * LL_new(:,rho) =  Ab * logsumtheta';
           */
#ifndef USEPTRS
          for(ti=0; ti<N; ti++) {
            tmpd = 0.0;
            for(tj=0; tj<N; tj++) 
              tmpd += Ab[IX(ti,tj,N,N)]*logsumtheta[tj];
            LL_new[IX(ti,rho,N,N2)] = tmpd;
          }
#else
          ptr1 = LL_new+rho*N;
          for(ti=0; ti<N; ti++, ++ptr1) {
            cptra = Ab+ti;
            ptr3 = logsumtheta;
            tmpd = 0.0;
            for(tj=0; tj<N; tj++, cptra+=N, ++ptr3) 
              tmpd += *cptra * *ptr3;
            *ptr1 = tmpd;
          }
#endif
          
          /* % normalize so that each clmn sums to 1 (may be not necessary ...)
           * % [N2xN] = [N2xN]   -       [N2x1]*[1xN]
           * theta = exp(logtheta - ones(N2,1) * logsumtheta);
           *
           * % and store for later
           * % [1xN2xNx1]
           * Theta(rho,:,:,t) = theta;
           */
#ifndef USEPTRS
          for(ti=0; ti<N; ti++) {
            for(tj=0; tj<N2; tj++) {
              Theta[IX4(rho,tj,ti,t,N2,N2,N,T)] = exp(logtheta[IX(tj,ti,N2,N)] - logsumtheta[ti]);
            }
          }
#else
          ptr1 = Theta + rho + N2*N2*N*t;
          ptr2 = logtheta;
          ptr3 = logsumtheta;
          for(ti=0; ti<N; ti++, ++ptr3) {
            for(tj=0; tj<N2; tj++, ptr1+=N2, ++ptr2) {
              *ptr1 = exp(*ptr2 - *ptr3);
            }
          }
#endif
            
        } /* end for rho */
        
        /*
         * LL_old = LL_new;
         */
#ifndef USEPTRS
        for(tj=0; tj<N2; tj++)
          for(ti=0; ti<N; ti++)
            LL_old[IX(ti,tj,N,N2)] = LL_new[IX(ti,tj,N,N2)];
#else
        ptr1 = LL_old;
        ptr2 = LL_new;
        for(tj=0; tj<N2; tj++)
          for(ti=0; ti<N; ti++, ++ptr1, ++ptr2)
            *ptr1 = *ptr2;
#endif
      } /* end for t */
      
      /*
       * % terminate the recursion
       *
       * logtheta = (log(hmm_r.prior) * ones(1,N)) + LLG_elbo' + LL_old';
       */
#ifndef USEPTRS
      for(int tj=0; tj<N2; tj++) {
        for(int ti=0; ti<N; ti++) {
          logtheta[IX(tj,ti,N2,N)] = log(h3m_r_prior[j_h3m][tj]) + LLG_elbo[IX(ti,tj,N,N2)] + LL_old[IX(ti,tj,N,N2)];
        }
      }
#else
      ptr2 = LLG_elbo;
      ptr3 = LL_old;
      cptra = h3m_r_prior[j_h3m];
      for(tj=0; tj<N2; tj++, ++cptra) {
        ptr1 = logtheta + tj;
        for(ti=0; ti<N; ti++, ptr1+=N2, ++ptr2, ++ptr3) {
          *ptr1 = *cptra + *ptr2 + *ptr3;
        }
      }
#endif
      
      /*
       * logsumtheta = logtrick(logtheta);
       */
      logtrick(logtheta, N2, N, logsumtheta);
      
      /*
       * LL_elbo =  hmm_b.prior' * logsumtheta';
       */
#ifndef USEPTRS
      tmpd = 0.0;
      for(ti=0; ti<N; ti++)
        tmpd += h3m_b_prior[i_h3m][ti] * logsumtheta[ti];
      LL_elbo[IX(i_h3m, j_h3m,Kb,Kr)] = tmpd;
#else
      cptra = h3m_b_prior[i_h3m];
      ptr2 = logsumtheta;
      tmpd = 0.0;
      for(ti=0; ti<N; ti++, ++cptra, ++ptr2)
        tmpd += *cptra * *ptr2;
      LL_elbo[IX(i_h3m, j_h3m,Kb,Kr)] = tmpd;
#endif
      
      /*
       * % normalize so that each clmn sums to 1 (may be not necessary ...)
       * theta = exp(logtheta - ones(N2,1) * logsumtheta);
       *
       * % and store for later
       * Theta_1 = theta;
       */
#ifndef USEPTRS
      for(ti=0; ti<N; ti++)
        for(tj=0; tj<N2; tj++) 
          Theta_1[IX(tj,ti,N2,N)] = exp(logtheta[IX(tj,ti,N2,N)] - logsumtheta[ti]);
#else
      ptr1 = Theta_1;
      ptr2 = logsumtheta;
      ptr3 = logtheta;
      for(ti=0; ti<N; ti++, ++ptr2)
        for(tj=0; tj<N2; tj++, ++ptr1, ++ptr3) 
          *ptr1 = exp(*ptr3 - *ptr2);
#endif
      
      /*
       * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       * %%%% do the forward recursion  %%%%
       * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       *
       * % rather then saving all intermediate values, just cache the cumulative
       * % values that are needed for the updates (it saves a lot of memory)
       *
       * % nu is N by N, first dimension indexed by sigma (M(r)), second by gamma (M(b))
       *
       * % initialize [N2 x N]
       * nu = (ones(N2,1) * hmm_b.prior') .* Theta_1; % N by N (indexed by sigma and gamma)
       */
#ifndef USEPTRS
      for(ti=0; ti<N; ti++)
        for(tj=0; tj<N2; tj++)
          nu[IX(tj,ti,N2,N)] = h3m_b_prior[i_h3m][ti] * Theta_1[IX(tj,ti,N2,N)];
#else
      ptr1 = nu;
      cptra = h3m_b_prior[i_h3m];
      ptr2 = Theta_1;
      for(ti=0; ti<N; ti++, ++cptra)
        for(tj=0; tj<N2; tj++, ++ptr1, ++ptr2)
          *ptr1 = *cptra * *ptr2;
#endif
      
      /* setup output cell entry */
      tmpmx = mxCreateDoubleMatrix(1,N2,mxREAL);
      mxSetCell(out_sum_nu_1, IX(i_h3m,j_h3m,Kb,Kr), tmpmx);
      sum_nu_1 = mxGetPr(tmpmx);
      
      /* [OUTPUT NOT NEEDED...ALREADY DECLARED ABOVE] */
      /*tmpmx = mxCreateDoubleMatrix(N2,N,mxREAL);
      mxSetCell(out_sum_t_nu, IX(i_h3m,j_h3m,Kb,Kr), tmpmx);
      sum_t_nu = mxGetPr(tmpmx);*/
      
      tmpmx = mxCreateDoubleMatrix(N2,N2,mxREAL);
      mxSetCell(out_sum_xi, IX(i_h3m,j_h3m,Kb,Kr), tmpmx);
      sum_t_sum_g_xi = mxGetPr(tmpmx);
      
      /*
       * % CACHE: sum_gamma nu_1(sigma,gamma) (this is one of the outputs ...)
       * sum_nu_1 = sum(nu,2)';
       */
#ifndef USEPTRS
      for(tj=0; tj<N2; tj++) {
        tmpd = 0.0;
        for(ti=0; ti<N; ti++)
          tmpd += nu[IX(tj,ti,N2,N)];
        sum_nu_1[tj] = tmpd;
      }
#else
      ptr1 = sum_nu_1;
      for(tj=0; tj<N2; tj++, ++ptr1) {
        ptr2 = nu + tj;      
        tmpd = 0.0;
        for(ti=0; ti<N; ti++, ptr2+=N2)
          tmpd += *ptr2;
        *ptr1 = tmpd;
      }
#endif
      /*
       * % CACHE: sum_t nu_t(sigma, gamma)
       * % sum_t_nu = zeros(N,N);
       * sum_t_nu = nu;
       */
#ifndef USEPTRS
      for(ti=0; ti<N; ti++)
        for(tj=0; tj<N2; tj++)
          sum_t_nu[IX(tj,ti,N2,N)] = nu[IX(tj,ti,N2,N)];
#else
      ptr1 = sum_t_nu;
      ptr2 = nu;
      for(ti=0; ti<N; ti++)
        for(tj=0; tj<N2; tj++, ++ptr1, ++ptr2)
          *ptr1 = *ptr2;
#endif
      
      /*
       * % CACHE: sum_t sum_gamma xi(rho,sigma,gamma,t)
       * sum_t_sum_g_xi = zeros(N2,N2); % N by N (indexed by rho and sigma)
       */
#ifndef USEPTRS
      for(ti=0; ti<N2; ti++)
        for(tj=0; tj<N2; tj++)
          sum_t_sum_g_xi[IX(tj,ti,N2,N2)] = 0.0;
#else
      ptr1 = sum_t_sum_g_xi;
      for(ti=0; ti<N2; ti++)
        for(tj=0; tj<N2; tj++, ++ptr1)
          *ptr1 = 0.0;
#endif
      
      /*
       * for t = 2 : T
       */
      for(t=1; t<T; t++) {        
        /*
         * % compute the inner part of the update of xi (does not depend on sigma)
         * foo = nu * Ab; % indexed by rho gamma [N2xN]
         */
#ifndef USEPTRS
        for(ti=0; ti<N; ti++) {
          for(tj=0; tj<N2; tj++) {
            tmpd = 0.0;
            for(tk=0; tk<N; tk++) {
              tmpd += nu[IX(tj,tk,N2,N)] * Ab[IX(tk,ti,N,N)];
            }
            foo[IX(tj,ti,N2,N)] = tmpd;
          }
        }
#else
        ptr1 = foo;
        cptra = Ab;
        for(ti=0; ti<N; ti++, cptra+=N) {
          for(tj=0; tj<N2; tj++, ++ptr1) {
            ptr2 = nu+tj;
            cptrb = cptra;
            tmpd = 0.0;
            for(tk=0; tk<N; tk++, ptr2+=N2, ++cptrb) {
              tmpd += *ptr2 * *cptrb;
            }
            *ptr1 = tmpd;
          }
        }
#endif
        
                
        /*
         * %% old code using loop
         * for sigma = 1 : N2
         */
        for(sigma=0; sigma<N2; sigma++) {
          /*
           * % new xi
           * % xi(:,sigma,:,t) = foo .* squeeze(Theta(:,sigma,:,t));
           * %xi_foo = foo .* squeeze(Theta(:,sigma,:,t)); % (indexed by rho gamma);
           *
           * % ABC: bug fix when another dim is 1
           * xi_foo = foo .* reshape(Theta(:,sigma,:,t), [size(Theta,1), size(Theta,3)]); % (indexed by rho gamma);
           */
#ifndef USEPTRS
          for(ti=0; ti<N; ti++)
            for(tj=0; tj<N2; tj++)
              xi_foo[IX(tj,ti,N2,N)] = foo[IX(tj,ti,N2,N)] * Theta[IX4(tj,sigma,ti,t,N2,N2,N,T)];
#else
          ptr1 = xi_foo;
          ptr2 = foo;
          ptr3 = Theta + N2*sigma + N2*N2*N*t;
          for(ti=0; ti<N; ti++, ptr3+=N2*N2) {
            ptr4 = ptr3;
            for(tj=0; tj<N2; tj++, ++ptr1, ++ptr2, ++ptr4)
              *ptr1 = *ptr2 * *ptr4;
          }
#endif
          /*
           * % CACHE:
           * sum_t_sum_g_xi(:,sigma) = sum_t_sum_g_xi(:,sigma) + sum(xi_foo,2);
           */
#ifndef USEPTRS
          for(tj=0; tj<N2; tj++) {
            tmpd = 0.0;
            for(ti=0; ti<N; ti++)
              tmpd += xi_foo[IX(tj,ti,N2,N)];
            sum_t_sum_g_xi[IX(tj,sigma,N2,N2)] += tmpd;
          }
#else
          ptr1 = sum_t_sum_g_xi+sigma*N2;
          for(tj=0; tj<N2; tj++, ++ptr1) {
            ptr2 = xi_foo+tj;
            tmpd = 0.0;
            for(ti=0; ti<N; ti++, ptr2+=N2)
              tmpd += *ptr2;
            *ptr1 += tmpd;
          }
#endif
          
          /*
           * % new nu
           * % nu(sigma,:) = ones(1,N) * squeeze(xi(:,sigma,:,t));
           * nu(sigma,:) = ones(1,N2) * xi_foo;
           */
#ifndef USEPTRS
          for(ti=0; ti<N; ti++) {
            tmpd = 0.0;
            for(tj=0; tj<N2; tj++)
              tmpd += xi_foo[IX(tj,ti,N2,N)];
            nu[IX(sigma,ti,N2,N)] = tmpd;
          }
#else
          ptr1 = nu+sigma;
          ptr2 = xi_foo;
          for(ti=0; ti<N; ti++, ptr1+=N2) {
            tmpd = 0.0;
            for(tj=0; tj<N2; tj++, ++ptr2)
              tmpd += *ptr2;
            *ptr1 = tmpd;
          }
#endif
        }
        
        /*
         * % CACHE: in the sum_t nu_t(sigma, gamma)
         * sum_t_nu = sum_t_nu + nu;
         */
#ifndef USEPTRS
        for(ti=0; ti<N; ti++)
          for(tj=0; tj<N2; tj++)
            sum_t_nu[IX(tj,ti,N2,N)] += nu[IX(tj,ti,N2,N)];
#else
        ptr1 = sum_t_nu;
        ptr2 = nu;
        for(ti=0; ti<N; ti++)
          for(tj=0; tj<N2; tj++, ++ptr1, ++ptr2)
            *ptr1 += *ptr2;
#endif
      }
      
      /* [STORED ALREADY]
       * % this is one of the outputs ...
       * sum_xi = sum_t_sum_g_xi;
       */

      /*
       * %%%% now prepare the cumulative sufficient statistics for the reestimation
       * %%%% of the emission distributions
       */
      
      /*
       * update_emit_pr = zeros(N2,M);
       */
      /* [ASSUME THAT M=1] */
      tmpmx = mxCreateDoubleMatrix(N2,1,mxREAL);
      mxSetCell(out_update_emit_pr, IX(i_h3m,j_h3m,Kb,Kr), tmpmx);
      update_emit_pr = mxGetPr(tmpmx);
            
      /*
       * update_emit_mu = zeros(N2,dim,M);
       */
      tmpmx = mxCreateDoubleMatrix(N2,dim,mxREAL);
      mxSetCell(out_update_emit_mu, IX(i_h3m,j_h3m,Kb,Kr), tmpmx);
      update_emit_mu = mxGetPr(tmpmx);
      
      /*
       * switch hmm_b.emit{1}.covar_type
       * case 'diag'
       *   update_emit_Mu = zeros(N2,dim,M);
       * case 'full'
       *   update_emit_Mu = zeros(N2,dim,dim,M);
       * end
       */
      switch(covmode) {
        case COV_DIAG:     
          tmpmx = mxCreateDoubleMatrix(N2,dim,mxREAL);
          break;
        case COV_FULL:
          tmpmx = create3Dmatrix(N2,dim,dim);
          break;
      }
      mxSetCell(out_update_emit_Mu, IX(i_h3m,j_h3m,Kb,Kr), tmpmx);
      update_emit_Mu = mxGetPr(tmpmx);            
      
      /*
       * % loop all the emission GMM of each state
       * for sigma = 1 : N2
       */
      for(sigma=0; sigma<N2; sigma++) {
        /*
         * update_emit_pr(sigma,:) = sum_t_nu(sigma,:) * sum_w_pr{sigma};
         */
#ifndef USEPTRS
        tmpd = 0.0;
        for(tj=0; tj<N; tj++)
          tmpd += sum_t_nu[IX(sigma,tj,N2,N)] * sum_w_pr[IX(tj,sigma,N,N2)];
        update_emit_pr[IX(sigma,0,N2,1)] = tmpd;
#else
        ptr1 = sum_t_nu+sigma;
        ptr2 = sum_w_pr+sigma*N;
        tmpd = 0.0;
        for(tj=0; tj<N; tj++, ptr1+=N2, ++ptr2)
          tmpd += *ptr1 * *ptr2;
        update_emit_pr[sigma] = tmpd;
#endif
        
        /*
         * foo_sum_w_mu = sum_w_mu{sigma};
         * foo_sum_w_Mu = sum_w_Mu{sigma};
         */
        
        /* [NO LOOP SINCE M=1]
         * for l = 1 : M
         */
        
        /*        
         * update_emit_mu(sigma,:,l) = sum_t_nu(sigma,:) * foo_sum_w_mu(:,:,l);
         */
#ifndef USEPTRS
        for(tk=0; tk<dim; tk++) {
          tmpd = 0.0;
          for(ti=0; ti<N; ti++)
            tmpd += sum_t_nu[IX(sigma,ti,N2,N)] * sum_w_mu[IX3(ti,tk,sigma,N,dim,N2)];
          update_emit_mu[IX(sigma,tk,N2,dim)] = tmpd;
        }
#else
        ptr1 = update_emit_mu+sigma;
        ptr3 = sum_w_mu + sigma*N*dim;
        for(tk=0; tk<dim; tk++, ptr1+=N2) {
          ptr2 = sum_t_nu+sigma;
          tmpd = 0.0;
          for(ti=0; ti<N; ti++, ++ptr3, ptr2+=N2)
            tmpd += *ptr2 * *ptr3;
          *ptr1 = tmpd;
        }
#endif
                
        /*
         * switch hmm_b.emit{1}.covar_type
         * case 'diag'
         *   update_emit_Mu(sigma,:,l) = sum_t_nu(sigma,:) * foo_sum_w_Mu(:,:,l);
         * case 'full'
         *   update_emit_Mu(sigma,:,:,l) = sum(bsxfun(@times,  ...
         *      reshape(sum_t_nu(sigma,:), [N 1 1]), foo_sum_w_Mu(:,:,:,l)), 1);
         * end
         */
#ifndef USEPTRS
        switch(covmode) {
          case COV_DIAG:
            for(tk=0; tk<dim; tk++) {
              tmpd = 0.0;
              for(ti=0; ti<N; ti++)
                tmpd += sum_t_nu[IX(sigma,ti,N2,N)] * sum_w_Mu[IX3(ti,tk,sigma,N,dim,N2)];
              update_emit_Mu[IX(sigma,tk,N2,dim)] = tmpd;              
            }
            break;
            
          case COV_FULL:
            for (tk=0; tk<dim; tk++) {
              for (tj=0; tj<dim; tj++) {
                tmpd = 0.0;
                for (ti=0; ti<N; ti++)
                  tmpd += sum_t_nu[IX(sigma,ti,N2,N)]* sum_w_Mu[IX4(ti,tk,tj,sigma,N,dim,dim,N2)];
                update_emit_Mu[IX3(sigma,tk,tj,N2,dim,dim)] = tmpd;
              }
            }
            break;
        }
#else
        switch(covmode) {
          case COV_DIAG:
            ptr1 = update_emit_Mu + sigma;
            ptr3 = sum_w_Mu + sigma*N*dim;
            for(tk=0; tk<dim; tk++, ptr1+=N2) {
              ptr2 = sum_t_nu + sigma;
              tmpd = 0.0;
              for(ti=0; ti<N; ti++, ptr2+=N2, ++ptr3)
                tmpd += *ptr2 * *ptr3;
              *ptr1 = tmpd;
            }
            break;
            
          case COV_FULL:
            /* for (tk=0; tk<dim; tk++) {
             *   for (tj=0; tj<dim; tj++) {
             *     tmpd = 0.0;
             *     for (ti=0; ti<N; ti++)
             *       tmpd += sum_t_nu[IX(sigma,ti,N2,N)]* sum_w_Mu[IX4(ti,tk,tj,sigma,N,dim,dim,N2)];
             *     update_emit_Mu[IX3(sigma,tk,tj,N2,dim,dim)] = tmpd;
             *   }
             * }
             */
            ptr1 = update_emit_Mu + sigma;
            ptr3 = sum_w_Mu + sigma*N*dim2;
            
            for (tk=0; tk<dim; tk++) {
              for (tj=0; tj<dim; tj++, ptr1+=N2) {
                ptr2 = sum_t_nu + sigma;                
                tmpd = 0.0;
                for (ti=0; ti<N; ti++, ptr2+=N2, ++ptr3)
                  tmpd += (*ptr2) * (*ptr3);
                *ptr1 = tmpd;
              }
            }
            break;
        }
#endif
        
        
      } /* end for sigma */
      
    } /* end for i_h3m */
  } /* end for j_hem, */
  
  /** clean up **/
  mxFree(h3m_r_N2);
  mxFree(h3m_r_A);
  mxFree(h3m_r_prior);
  mxFree(h3m_r_emit_m);
  mxFree(h3m_r_emit_W);
  mxFree(h3m_b_N);
  mxFree(h3m_b_A);
  mxFree(h3m_b_prior);
  mxFree(h3m_b_emit_centres);
  mxFree(h3m_b_emit_covars);
  mxFree(LLG_elbo);
  mxFree(sum_w_pr);
  mxFree(sum_w_mu);
  mxFree(sum_w_Mu);
  mxFree(Theta);
  mxFree(LL_new);
  mxFree(LL_old);
  mxFree(logtheta);
  mxFree(logsumtheta);
  mxFree(Theta_1);
  mxFree(nu);
  mxFree(foo);
  mxFree(xi_foo);
  mxFree(sum_t_nu);
  if (covmode == COV_FULL) {
    mxFree(h3m_r_logdetCovPlusDdivlam);
    mxFree(h3m_r_invCov);
    mxFree(tmpx);
  }
}
