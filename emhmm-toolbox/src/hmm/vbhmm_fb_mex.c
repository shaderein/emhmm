/*==========================================================
 * vbhmm_fb_mex - MEX file for running forward-backward algorithm
 *
 * this is an internal function called by vbhmm_fb
 *
 * Inputs
 *   data = {Nx1} - data{i} = [dim x T]
 *   K, N, dim, maxT = scalars [1x1]
 *   m = [dim x K]
 *   W = [dim x dim x K]
 *   v = [K x 1]
 *   beta = [K x 1]
 *   logLambdaTilde = [1 x K]
 *   const_denominator = scalar [1 x 1]
 *   t_pz1 = [1 x K]
 *   t_tpz1zt1 = [K x K]
 * 
 * Outputs:
 *   logrho_Saved = [K x N x maxT]
 *   gamma_all    = [K x N x maxT]
 *   xi_sum       = [K x K x N]
 *   phi_norm     = [1 x N]
 *   xi_Saved (unsupported, since currently unused elsewhere)
 * [logrho_Saved, gamma_all, xi_sum, phi_norm] = ...
 *       vbhmm_fb_mex(data, K, N, dim, maxT, m, W, v, beta, logLambdaTilde, const_denominator, t_pz1, t_tpztzt1);
 *
 *
 * ---
 * Eye-Movement analysis with HMMs (emhmm-toolbox)
 * Copyright (c) 2017-08-02
 * Antoni B. Chan, Janet H. Hsiao
 * City University of Hong Kong, University of Hong Kong
 *========================================================*/

/* Version info
 *  2017/08/02 : ABC - initial version
 *  2017/08/03 : ABC - update compatability with Windows
 *  2017/08/25 : ABC - scaling to handle numerical precision problems
 *  2018/05/15 : v0.72 - full C compatability for Linux
 */

/* compile commands
 *
 * for debugging
 *   mex -g vbhmm_fb_mex.c
 *
 * optimized
 *   mex vbhmm_fb_mex.c
 */

#include "mex.h"
#include <math.h>
#include <stdio.h>

/* #define DEBUG */

/* use pointer arithmetic for faster matrix indexing
   (compiler optimizations give similar speed) */
#define USEPTRS

/* Macros for indexing matrices */
/* (the last argument is not used, but there so that you can see the matrix size) */
/* #define IX(i,j,M,N)      (i+j*M)       */
/* #define IX3(i,j,k,M,N,D) (i+j*M+k*M*N) */

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
    mexErrMsgIdAndTxt("vbhmm_fb_mex:invalidinput", "parseMatrix: invalid num dimensions.");

  nrows = mxGetM(mx);
  ncols = mxGetN(mx);

#ifdef DEBUG
  printf("nr=%d; nc=%d;\n", nrows, ncols);
#endif

  if (((nrows != M) && (M!=0)) || ((ncols != N) && (N!=0)))
    mexErrMsgIdAndTxt("vbhmm_fb_mex:invalidinput", "parseMatrix: invalid size.");

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
  
  if (! (  ((ndims == 3) && (D>1)) || 
           ((ndims == 2) && (D==1))) )
    mexErrMsgIdAndTxt("vbhmm_fb_mex:invalidinput", "parseMatrix3: invalid num dimensions.");
  
  tmpdims = mxGetDimensions(mx);
  
#ifdef DEBUG
  printf("nr=%d; nc=%d; np=%d;\n", tmpdims[0], tmpdims[1], tmpdims[2]);
#endif

  if ((tmpdims[0] != M) || (tmpdims[1] != N) || ((D>1) && (tmpdims[2] != D)))
    mexErrMsgIdAndTxt("vbhmm_fb_mex:invalidinput", "parseMatrix3: invalid size.");

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



/* The gateway function */  
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  /*double multiplier;              input scalar */
  /*double *inMatrix;                1xN input matrix */
  /*size_t ncols;                    size of matrix */
  /*double *outMatrix;               output matrix */
  
  /* size constants */
  int K, N, dim, maxT;
  int dim2, NK, NK1, K2;  /* products of constants */
  
  /* input matrices and values */
  const double *tdata_t;
  const double *m, *W, *v, *beta, *logLambdaTilde, *t_pz1, *t_tpztzt1;
  double const_denominator;
  
  /* output matrices */
  double *logrho_Saved, *gamma_all, *xi_sum, *phi_norm;
  
  /* working variables */
  double *delta, *diff, *diff2;
  double *t_alpha, *t_Delta, *t_beta, *t_Eta, *t_c;
  double *t_pxtzt, *t_gamma, *t_sumxi, *bpi;
  int tT;
  double tmp;
  double *max_fb_logrho;
  
  /* working pointers */
#ifdef USEPTRS
  double *ptr1, *ptr2, *ptr3, *ptr4;
  const double *cptr4, *cptr5;
#endif  
  
  /* counters */
  int n, k, t, ti, tj, tk;
  
  /* check for proper number of arguments */
  if(nrhs!=13) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","13 inputs required.");
  }
  if(nlhs!=4) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
  }

  
  /*
   * Check and Read Inputs:
   *   data = {Nx1} - data{i} = [dim x T]
   *   K, N, dim, maxT = scalars [1x1]
   *   m = [dim x K]
   *   W = [dim x dim x K]
   *   v = [K x 1]
   *   beta = [K x 1]
   *   logLambdaTilde = [1 x K]
   *   const_denominator = scalar [1 x 1]
   *   t_pz1 = [1 x K]
   *   t_tpztzt1 = [K x K]
   */
  
  if( !mxIsCell(prhs[0]) ) {
    mexErrMsgIdAndTxt("vbhmm_fb_mex:invalidinput", "1st arg must be cell");
  }
  
  K    = (int) parseScalar(prhs[1]);
  N    = (int) parseScalar(prhs[2]);
  dim  = (int) parseScalar(prhs[3]);
  maxT = (int) parseScalar(prhs[4]);
  
  dim2 = dim*dim;
  NK   = N*K;
  NK1  = NK-K;
  K2   = K*K;
  
#ifdef DEBUG
  printf("K=%d; N=%d; dim=%d; maxT=%d\n", K, N, dim, maxT);
#endif
  
  m                 = parseMatrix(prhs[5], dim, K);
  W                 = parseMatrix3(prhs[6], dim, dim, K);
  v                 = parseMatrix(prhs[7], K, 1);
  beta              = parseMatrix(prhs[8], K, 1);
  logLambdaTilde    = parseMatrix(prhs[9], 1, K);  
  const_denominator = parseScalar(prhs[10]);
  t_pz1             = parseMatrix(prhs[11], 1, K);
  t_tpztzt1         = parseMatrix(prhs[12], K, K);
  
#ifdef DEBUG
  dumpMatrix("m", m, dim, K);
  dumpMatrix3("W", W, dim, dim, K);
  dumpMatrix("v", v, K, 1);
  dumpMatrix("pz1", t_pz1, 1, K);
#endif
  
  /*
   * Create Outputs:
   *   logrho_Saved = [K x N x maxT]
   *   gamma_all    = [K x N x maxT]
   *   xi_sum       = [K x K x N]
   *   phi_norm     = [1 x N]
   */

  /* create the output matrix */
  plhs[0] = create3Dmatrix(K, N, maxT);
  plhs[1] = create3Dmatrix(K, N, maxT);
  plhs[2] = create3Dmatrix(K, K, N);
  plhs[3] = mxCreateDoubleMatrix(1, (mwSize)N, mxREAL);

  logrho_Saved = mxGetPr(plhs[0]);
  gamma_all    = mxGetPr(plhs[1]);
  xi_sum       = mxGetPr(plhs[2]);
  phi_norm     = mxGetPr(plhs[3]);
  
  /* initialize storage */  
  delta = mxMalloc(K*maxT*sizeof(double));    /* [K x maxT]   */
  diff  = mxMalloc(dim*maxT*sizeof(double));  /* [dim x maxT] */
  diff2 = mxMalloc(dim*maxT*sizeof(double));  /* [dim x maxT] */
  t_alpha = mxMalloc(maxT*K*sizeof(double));  /* [maxT x K]   */
  t_Delta = mxMalloc(maxT*K*sizeof(double));  /* [maxT x K]   */
  t_beta  = mxMalloc(maxT*K*sizeof(double));  /* [maxT x K]   */
  t_Eta   = mxMalloc(maxT*K*sizeof(double));  /* [maxT x K]   */
  t_c     = mxMalloc(maxT*sizeof(double));    /* [maxT]       */
  t_pxtzt = mxMalloc(maxT*K*sizeof(double));  /* [maxT x K]   */
  t_gamma = mxMalloc(K*maxT*sizeof(double));  /* [K x maxT]   */
  t_sumxi = mxMalloc(K*K*sizeof(double));     /* [K x K]      */
  bpi     = mxMalloc(K*sizeof(double));       /* [1 x K]      */
  max_fb_logrho = mxMalloc(maxT*sizeof(double)); /* [1 x maxT] */
     
  /* loop over all data */
  for(n=0; n<N; n++) {
    /*
     * tdata = data{n}; tdata = tdata';
     * tT = size(tdata,2);
    */
    mxArray *tdata_mx = mxGetCell(prhs[0], n);
    tdata_t = parseMatrix(tdata_mx, 0, dim);  /* the original tdata [T x dim] (without transposing) */
    tT = mxGetM(tdata_mx);
            
#ifdef DEBUG
    printf("-- n=%d; tT=%d --\n", n, tT);
#endif
    
    /* [PRE-ALLOCATED]
     * delta = zeros(K,tT);
    */
    
    /* 
     * for k = 1:K
    */
    for(k=0; k<K; k++) {
      /*
       *   diff = bsxfun(@minus, tdata, m(:,k));
      */
#ifndef USEPTRS
      for(tj=0; tj<tT; tj++)
        for(ti=0; ti<dim; ti++)
          diff[IX(ti,tj,dim,maxT)] = tdata_t[IX(tj,ti,tT,dim)] - m[IX(ti,k,dim,K)];
#else     
      ptr1 = diff;
      for(tj=0; tj<tT; tj++) {        
        cptr5 = m+k*dim;
        cptr4 = tdata_t+tj;
        for(ti=0; ti<dim; ti++, ++ptr1, ++cptr5, cptr4+=tT) {
          *ptr1 = *cptr4 - *cptr5;
        }
      }
#endif
      
      
#ifdef DEBUG      
      printf("pHere 1");
#endif
      /*
       *   mterm = sum((W(:,:,k)*diff).*diff,1);
      */
      /* (W(:,:,k)*diff) */
#ifndef USEPTRS
      for(tj=0; tj<tT; tj++) {
        for(ti=0; ti<dim; ti++) {
          tmp = 0.0;
          for(int tk=0; tk<dim; tk++)
            tmp += W[IX3(ti,tk,k,dim,dim,K)]*diff[IX(tk,tj,dim,maxT)];
          diff2[IX(ti,tj,dim,maxT)] = tmp;
        }
      }
#else
      ptr1 = diff2;
      for(tj=0; tj<tT; tj++) {
        cptr4 = W+k*dim2; 
        for(ti=0; ti<dim; ti++, ++ptr1) {
          ptr2 = diff+dim*tj;
          tmp = 0.0;
          for(tk=0; tk<dim; tk++, ++cptr4, ++ptr2) {
            /* use the fact that W is symmetric when doing  */
            tmp += *cptr4 * *ptr2; 
          }
          *ptr1 = tmp;
        }
      }
#endif
      
#ifdef DEBUG            
      printf("pHere 2");
#endif
      
      /* .*diff */
#ifndef USEPTRS
      for(tj=0; tj<tT; tj++)
        for(ti=0; ti<dim; ti++)
          diff2[IX(ti,tj,dim,maxT)] *= diff[IX(ti,tj,dim,maxT)];
#else
      ptr1 = diff2;
      ptr2 = diff;
      for(tj=0; tj<tT; tj++)
        for(ti=0; ti<dim; ti++, ++ptr1, ++ptr2)
          *ptr1 *= *ptr2;      
#endif
      
#ifdef DEBUG
      printf("pHere 3");
#endif

      /*
       *   mterm = sum(...,1);
       *   delta(k,:) = dim/beta(k) + v(k) * mterm;
      */
#ifndef USEPTRS
      for(tj=0; tj<tT; tj++) {
        tmp = 0.0;
        for(ti=0; ti<dim; ti++)
          tmp += diff2[IX(ti,tj,dim,maxT)];
        delta[IX(k,tj,K,maxT)] = dim/beta[k] + v[k]*tmp;
      }
#else
      ptr1 = delta+k;
      ptr2 = diff2;
      for(tj=0; tj<tT; tj++, ptr1+=K) {
        tmp = 0.0;
        for(ti=0; ti<dim; ti++, ++ptr2)
          tmp += *ptr2;
        *ptr1 = dim/beta[k] + v[k]*tmp;
      }
#endif
      
#ifdef DEBUG
      printf("pHere 4");
#endif
    }
    
    /*
     * logrho = bsxfun(@minus, 0.5*logLambdaTilde(:), 0.5*delta) - const_denominator;
     * logrho_Saved(:,n,1:tT) = logrho;
    */
#ifndef USEPTRS
    for(tj=0; tj<tT; tj++) {
      for(ti=0; ti<K; ti++) {
        logrho_Saved[IX3(ti,n,tj,K,N,maxT)] = 0.5*(logLambdaTilde[ti] - delta[IX(ti,tj,K,maxT)]) - const_denominator;
      }
    }
#else
    ptr1 = logrho_Saved + n*K;
    ptr2 = delta;
    for(tj=0; tj<tT; tj++, ptr1+=NK1) {
      cptr4 = logLambdaTilde;
      for(ti=0; ti<K; ti++, ++ptr2, ++cptr4, ++ptr1) {
        *ptr1 = 0.5*(*cptr4 - *ptr2) - const_denominator;
      }
    }
#endif
    
    
#ifdef DEBUG          
    printf("pHere 5");
#endif
    
    /** forward_backward **********************************/
    
    /* [NOT NEEDED - use t_gamma and t_sumxi]
     * gamma = zeros(K,tT);
     * sumxi = zeros(K,K);  % [row]
     */
    
    /* [SEE BELOW]
     * fb_logrho = logrho';
     */
    
    /* [PRECOMPUTED]
     * if ~usegroups
     *   fb_logPiTilde = logPiTilde';
     *   fb_logATilde = logATilde;
     * else
     *   fb_logPiTilde = logPiTilde{group_map(n)}';
     *   fb_logATilde = logATilde{group_map(n)};
     * end
     */
    
    /* [PRE-ALLOCATED]
     * t_alpha = zeros(tT, K); % alpha hat
     * t_Delta = zeros(tT, K); % Delta (unnormalized version)
     * t_beta  = zeros(tT, K);
     * t_Eta   = zeros(tT, K); % Eta (unnormalized version)
     * t_c     = zeros(1,tT);
     */
    
    /* [PRE-COMPUTED]
     * t_pz1 = exp(fb_logPiTilde);    % priors p(z1)
     * t_tpztzt1 = exp(fb_logATilde); % transitions p(zt|z(t-1)) [row format]
     */
    
    /*
     * %t_pxtzt = exp(fb_logrho);      % emissions p(xt|zt)
     *
     * % emission p(xt|zt)/max(p(xt|zt)) (more stable)
     * max_fb_logrho = max(fb_logrho, [], 2);
     * t_pxtzt = exp(bsxfun(@minus, fb_logrho, max_fb_logrho));    
     */
#ifndef USEPTRS    
    for(tj=0; tj<tT; tj++) {
      /* find max */
      tmp = logrho_Saved[IX3(0,n,tj,K,N,maxT)];
      for(ti=1; ti<K; ti++) {
        if (logrho_Saved[IX3(ti,n,tj,K,N,maxT)] > tmp)
          tmp = logrho_Saved[IX3(ti,n,tj,K,N,maxT)];
      }
      max_fb_logrho[tj] = tmp;
      
      for(ti=0; ti<K; ti++)
        t_pxtzt[IX(tj,ti,maxT,K)] = exp(logrho_Saved[IX3(ti,n,tj,K,N,maxT)] - tmp);
    }
    
#else
    ptr2 = logrho_Saved + n*K;
    ptr4 = max_fb_logrho;
    for(tj=0; tj<tT; tj++, ptr2+=NK, ++ptr4) {
      ptr3 = ptr2;      
      tmp = *ptr3;
      ++ptr3;
      for(ti=1; ti<K; ti++, ++ptr3)
        if (*ptr3 > tmp)
          tmp = *ptr3;
      
      *ptr4 = tmp;
      
      ptr1 = t_pxtzt + tj;
      ptr3 = ptr2;
      for(ti=0; ti<K; ti++, ++ptr3, ptr1+=maxT)
        *ptr1 = exp(*ptr3 - tmp);
    }
#endif
    
#ifdef DEBUG      
    printf("pHere 6");
#endif
    
    /* [UNNEEDED - use tT]
     * t_x = tdata';
     * t_T = size(t_x,1);
     */
    
    /*
     * if t_T >= 1
     */
    if (tT >= 1) {
      
      /** FORWARD PASS **/
      
      /* [PRE-ALLOCATED]
       * t_gamma = zeros(K,t_T);
       * t_sumxi = zeros(K,K);
       */
      
      /*
       * %before normalizing
       * t_Delta(1,:) = t_pz1.*t_pxtzt(1,:);
       */      
#ifndef USEPTRS
      for(ti=0; ti<K; ti++)
        t_Delta[IX(0,ti,maxT,K)] = t_pz1[ti]*t_pxtzt[IX(0,ti,maxT,K)];
#else
      ptr1  = t_Delta;
      ptr2  = t_pxtzt;
      cptr4 = t_pz1;
      for(ti=0; ti<K; ti++, ++cptr4, ptr2+=maxT, ptr1+=maxT)
        *ptr1 = *cptr4 * *ptr2;
#endif
      
#ifdef DEBUG
      printf("pHere 7");
#endif

      /*
       * % 2016-04-29 ABC: rescale for numerical stability (otherwise values get too small)
       * t_c(1) = sum(t_Delta(1,:));
       */
#ifndef USEPTRS
      tmp = 0.0;
      for(ti=0; ti<K; ti++)
        tmp += t_Delta[IX(0,ti,maxT,K)];
      t_c[0] = tmp;
#else      
      ptr1 = t_Delta;
      tmp = 0.0;      
      for(ti=0; ti<K; ti++, ptr1+=maxT)
        tmp += *ptr1;
      t_c[0] = tmp;
#endif
      
#ifdef DEBUG
      printf("pHere 8");
#endif

      /*
       * % normalize
       * t_alpha(1,:) = t_Delta(1,:) / t_c(1);
       */
#ifndef USEPTRS
      for(ti=0; ti<K; ti++)
        t_alpha[IX(0,ti,maxT,K)] = t_Delta[IX(0,ti,maxT,K)] / t_c[0];
#else
      ptr1 = t_alpha;
      ptr2 = t_Delta;
      for(ti=0; ti<K; ti++, ptr1+=maxT, ptr2+=maxT)
        *ptr1 = *ptr2 / t_c[0];
#endif
      
#ifdef DEBUG
      printf("pHere 9");
#endif

      /* 
       * if t_T > 1
       *   for i=2:t_T
       */
      if (tT > 1) {
        for(t=1; t<tT; t++) { /* NOTE: "i" changed to "t" */
#ifdef DEBUG
          printf("front t=%d\n", t);
#endif
          
          /*
           * % before normalizing
           * t_Delta(i,:) = t_alpha(i-1,:)*t_tpztzt1.*t_pxtzt(i,:);
           */
#ifndef USEPTRS
          for(tj=0; tj<K; tj++) {
            tmp = 0.0;
            for(ti=0; ti<K; ti++) {
              tmp += t_alpha[IX(t-1,ti,maxT,K)]*t_tpztzt1[IX(ti,tj,K,K)];
            }
            t_Delta[IX(t,tj,maxT,K)] = tmp * t_pxtzt[IX(t,tj,maxT,K)];
          }
#else
          ptr1  = t_Delta+t;
          ptr3  = t_pxtzt+t;
          cptr4 = t_tpztzt1;
          for(tj=0; tj<K; tj++, ptr1+=maxT, ptr3+=maxT) {
            ptr2  = t_alpha+t-1;
            tmp = 0.0;
            for(ti=0; ti<K; ti++, ptr2+=maxT, ++cptr4) {
              tmp += *ptr2 * *cptr4;
            }
            *ptr1 = tmp * *ptr3;
          }
#endif
          
#ifdef DEBUG
          printf("pHere 10");
#endif
          
          /*
           * % 2016-04-29 ABC: rescale for numerical stability
           * t_c(i) = sum(t_Delta(i,:));
           */
#ifndef USEPTRS
          tmp = 0.0;
          for(ti=0; ti<K; ti++)
            tmp += t_Delta[IX(t,ti,maxT,K)];
          t_c[t] = tmp;
#else
          ptr1 = t_Delta+t;
          tmp = 0.0;
          for(ti=0; ti<K; ti++, ptr1+=maxT)
            tmp += *ptr1;
          t_c[t] = tmp;
#endif
          
#ifdef DEBUG
          printf("pHere 11");
#endif
          
          /* 
           * % normalize
           * t_alpha(i,:) = t_Delta(i,:) / t_c(i);
           */
#ifndef USEPTRS
          for(ti=0; ti<K; ti++)
            t_alpha[IX(t,ti,maxT,K)] = t_Delta[IX(t,ti,maxT,K)] / t_c[t];
#else
          ptr1 = t_alpha+t;
          ptr2 = t_Delta+t;
          for(ti=0; ti<K; ti++, ptr1+=maxT, ptr2+=maxT)
            *ptr1 = *ptr2 / t_c[t];
#endif
          
#ifdef DEBUG
          printf("pHere 12");
#endif

        }
      }

      
      /*
       * %% backward %%%%%%%%%%%%%%%%%%%%%%
       * % BUG FIX 2017-02-09: ones(1,K)./K;
       * % doesn't affect HMM results, since we had normalized anyways.
       * t_beta(t_T,:) = ones(1,K);
       */
#ifndef USEPTRS
      for(ti=0; ti<K; ti++)
        t_beta[IX(tT-1,ti,maxT,K)] = 1.0;
#else
      ptr1 = t_beta+tT-1;
      for(ti=0; ti<K; ti++, ptr1+=maxT)
        *ptr1 = 1.0;
#endif
      
#ifdef DEBUG
      printf("Here 1");
#endif
      
      /*
       * t_gamma(:,t_T) = (t_alpha(t_T,:).*t_beta(t_T,:))';
      */
#ifndef USEPTRS
      for(ti=0; ti<K; ti++)
        t_gamma[IX(ti,tT-1,K,maxT)] = t_alpha[IX(tT-1,ti,maxT,K)]*t_beta[IX(tT-1,ti,maxT,K)];
#else
      ptr1 = t_gamma+(tT-1)*K;
      ptr2 = t_alpha+(tT-1);
      ptr3 = t_beta+(tT-1);
      for(ti=0; ti<K; ti++, ++ptr1, ptr2+=maxT, ptr3+=maxT)
        *ptr1 = *ptr2 * *ptr3;
#endif
      
#ifdef DEBUG      
      printf("Here 2");
#endif
      
      /* initialize accumulator */
#ifndef USEPTRS
      for(tj=0; tj<K; tj++)
        for(ti=0; ti<K; ti++)
          t_sumxi[IX(ti,tj,K,K)] = 0.0;
#else
      ptr1 = t_sumxi;
      for(tj=0; tj<K; tj++)
        for(ti=0; ti<K; ti++, ++ptr1)
          *ptr1 = 0.0;
#endif
      
       /*
        * if t_T > 1
        *   for i=(t_T-1):-1:1
        */
      if (tT>1) {
#ifdef DEBUG        
        printf("Here 3");
#endif

        for(t=tT-2; t>=0; t--) {  /* note i changed to t */
#ifdef DEBUG
          printf("back t=%d\n", t);
#endif
          
          /*
           * % before normalizing
           * bpi = (t_beta(i+1,:).*t_pxtzt(i+1,:));
           */
#ifndef USEPTRS          
          for(ti=0; ti<K; ti++)
            bpi[ti] = t_beta[IX(t+1,ti,maxT,K)]*t_pxtzt[IX(t+1,ti,maxT,K)];
#else
          ptr1 = bpi;
          ptr2 = t_beta+t+1;
          ptr3 = t_pxtzt+t+1;
          for(ti=0; ti<K; ti++, ++ptr1, ptr2+=maxT, ptr3+=maxT)
            *ptr1 = *ptr2 * *ptr3;
#endif
          
#ifdef DEBUG          
          printf("Here 4");
#endif
          
          /*
           * t_Eta(i,:) = bpi*t_tpztzt1';
           */
#ifndef USEPTRS
          for(ti=0; ti<K; ti++) {
            tmp = 0.0;
            for (tj=0; tj<K; tj++)
              tmp += bpi[tj]*t_tpztzt1[IX(ti,tj,K,K)];
            t_Eta[IX(t,ti,maxT,K)] = tmp;
          } 
#else
          
          ptr1 = t_Eta+t;
          for(ti=0; ti<K; ti++, ptr1+=maxT) {
            tmp = 0.0;
            ptr2 = bpi;
            cptr4 = t_tpztzt1+ti;
            for (tj=0; tj<K; tj++, ++ptr2, cptr4+=K)
              tmp += *ptr2 * *cptr4;
            *ptr1 = tmp;
          } 
#endif
          
#ifdef DEBUG
          printf("Here 5");
#endif

          /*
           * % 2016-04-29 ABC: rescale
           * t_beta(i,:) = t_Eta(i,:)/t_c(i+1);
           */
#ifndef USEPTRS
          for(ti=0; ti<K; ti++)
            t_beta[IX(t,ti,maxT,K)] = t_Eta[IX(t,ti,maxT,K)] / t_c[t+1];
#else
          ptr1 = t_beta+t;
          ptr2 = t_Eta+t;
          for(ti=0; ti<K; ti++, ptr1+=maxT, ptr2+=maxT)
            *ptr1 = *ptr2 / t_c[t+1];
#endif
          
#ifdef DEBUG
          printf("Here 6");
#endif
          
          /*
           * % responsibility
           * t_gamma(:,i) = (t_alpha(i,:).*t_beta(i,:))';
           */
#ifndef USEPTRS
          for(ti=0; ti<K; ti++)
            t_gamma[IX(ti,t,K,maxT)] = t_alpha[IX(t,ti,maxT,K)]*t_beta[IX(t,ti,maxT,K)];
#else
          ptr1 = t_gamma+t*K;
          ptr2 = t_alpha+t;
          ptr3 = t_beta+t;
          for(ti=0; ti<K; ti++, ++ptr1, ptr2+=maxT, ptr3+=maxT)
            *ptr1 = *ptr2 * *ptr3;
#endif
          
#ifdef DEBUG          
          printf("Here 7");
#endif
          
          /*
           * %t_sumxi = t_sumxi + (t_logATilde.*(t_alpha(i,:)'*bpi));
           *
           * % 2016-04-29 ABC: rescale xi
           * tmp_xi = (t_tpztzt1.*(t_alpha(i,:)'*bpi)) / t_c(i+1);
           *
           * % 2016-04-29 ABC BUG FIX: normalize xi matrix to sum to 1
           * % (it's a joint probability matrix)
           * % 2017-02-09 - not necessary anymore after fixing beta initialization bug.
           * %tmp_xi = tmp_xi / sum(tmp_xi(:));
           *
           * % accumulate
           * t_sumxi = t_sumxi + tmp_xi;
           */
          
          /* computation of tmp_xi is merged with accumulator */
#ifndef USEPTRS
          for(tj=0; tj<K; tj++) 
            for(ti=0; ti<K; ti++)
              t_sumxi[IX(ti,tj,K,K)] += t_tpztzt1[IX(ti,tj,K,K)]*
                                        t_alpha[IX(t,ti,maxT,K)]*bpi[tj] / t_c[t+1];
#else
          ptr1 = t_sumxi;
          ptr2 = bpi;
          cptr4 = t_tpztzt1;
          for(tj=0; tj<K; tj++, ++ptr2) {
            ptr3 = t_alpha+t;
            for(ti=0; ti<K; ti++, ++ptr1, ++cptr4, ptr3+=maxT)
              *ptr1 += *cptr4 * *ptr3 * *ptr2/ t_c[t+1];
          }
#endif
          
#ifdef DEBUG          
          printf("Here 8");
#endif
          /* [UNSUPPORTED]
           * if savexi
           *   xi_Saved{n}(:,:,i) = tmp_xi;
           * end
           */
        }
      }
      
      /* [UNNECESSARY]
       * %for i = 1:size(t_gamma,2)
       * %  gamma(:,i) = t_gamma(:,i);
       * %end
       * gamma(:,1:tT) = t_gamma;
       * sumxi = t_sumxi;
       */
    
      /*
       * gamma_all(:,n,1:tT) = gamma;
       */
#ifndef USEPTRS
      for(tj=0; tj<tT; tj++)
        for(ti=0; ti<K; ti++)
          gamma_all[IX3(ti,n,tj,K,N,maxT)] = t_gamma[IX(ti,tj,K,maxT)];
#else
      ptr1 = gamma_all+K*n;
      ptr2 = t_gamma;
      for(tj=0; tj<tT; tj++, ptr1+=NK1)
        for(ti=0; ti<K; ti++, ++ptr1, ++ptr2)
          *ptr1 = *ptr2;
#endif
      
#ifdef DEBUG
      printf("Here 9");
#endif
      
      /*
       * xi_sum(:,:,n) = sumxi;
       */
#ifndef USEPTRS
      for(tj=0; tj<K; tj++)
        for(ti=0; ti<K; ti++)
          xi_sum[IX3(ti,tj,n,K,K,N)] = t_sumxi[IX(ti,tj,K,K)];
#else
      ptr1 = xi_sum + n*K2;
      ptr2 = t_sumxi;
      for(tj=0; tj<K; tj++)
        for(ti=0; ti<K; ti++, ++ptr1, ++ptr2)
          *ptr1 = *ptr2;
#endif
      
#ifdef DEBUG
      printf("Here 10");
#endif

      /*
       * % from scaling constants
       * %phi_norm(n) = sum(log(t_c));
       *
       * % from scaling constants (more stable)
       * phi_norm(n) = sum(log(t_c)) + sum(max_fb_logrho);
       */
#ifndef USEPTRS
      tmp = 0.0;
      for(ti=0; ti<tT; ti++)
        tmp += log(t_c[ti]) + max_fb_logrho[ti];
      phi_norm[n] = tmp;
#else
      ptr1 = t_c;
      ptr2 = max_fb_logrho;
      tmp = 0.0;
      for(ti=0; ti<tT; ti++, ++ptr1, ++ptr2)
        tmp += log(*ptr1) + *ptr2;
      phi_norm[n] = tmp;
#endif
      
#ifdef DEBUG
      printf("Here 11");
#endif
    }
    
    /* break early */
    /*printf("break\n");
    break; */
  }
#ifdef DEBUG
  printf("Here 12");
#endif

  /* clean up */
  mxFree(delta);
  mxFree(diff);
  mxFree(diff2);
  mxFree(t_alpha);
  mxFree(t_Delta);
  mxFree(t_beta);
#ifdef DEBUG
  printf("Here 13\n");
#endif

  mxFree(t_Eta);
  mxFree(t_c);
  mxFree(t_pxtzt);
  mxFree(t_gamma);
  mxFree(t_sumxi);
  mxFree(bpi);
  mxFree(max_fb_logrho);
}
