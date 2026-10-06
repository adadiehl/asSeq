/*
 *  ase.c
 *
 *  Created by Wei Sun on 5/25/2010.
 *  Modified by Vasyl Zhabotynsky on 04/05/2011 
 *
 */

#include "ase.h"

/**********************************************************************
 *
 * Closed forms of the sums in the beta-binomial likelihood
 *
 * The ASE likelihood with mean pi and over-dispersion theta has terms
 * prod_{k=0}^{m-1} (a + k theta) for a = pi, 1 - pi or 1, where m is a
 * read count. Summing log(a + k theta) and its derivatives term by term
 * costs O(m), which is prohibitive for highly expressed genes with
 * hundreds of thousands of allele-specific reads per sample. With
 * x = a/theta,
 *
 *   sum log(a + k theta) = m log(a) + sum log(1 + k/x)
 *                        = m log(a) + lgamma(x + m) - lgamma(x) - m log(x)
 *   sum 1/(a + k theta)  = (digamma(x + m) - digamma(x))/theta
 *   sum k/(a + k theta)  = sum k/(x + k) / theta
 *   sum 1/(a + k theta)^2 = (trigamma(x) - trigamma(x + m))/theta^2
 *
 * For x < 10 the gamma-function differences are evaluated directly. For
 * x >= 10 (small theta) they are written with Stirling-type expansions
 * so that the large terms cancel analytically, e.g.
 *
 *   sum log(1 + k/x) = x g(u) - log(1 + u)/2 + c(x + m) - c(x),
 *
 * with u = m/x, g(u) = (1 + u) log(1 + u) - u and c the Stirling
 * correction of lgamma; g and u - log(1 + u) are evaluated by their
 * Taylor series for small u. The expansions are accurate to about 1e-16
 * for x >= 10. If theta <= 0 (not used by asSeq) the sums are computed
 * term by term.
 *
 **********************************************************************/

#define BB_XMIN 10.0

/* lgamma(y) - [(y - 1/2) log(y) - y + log(2 pi)/2], for y >= 10 */
static double bb_stirling_c(double y)
{
  double r = 1.0/y, r2 = r*r;
  return r*(1.0/12 + r2*(-1.0/360 + r2*(1.0/1260 + r2*(-1.0/1680 + 
         r2*(1.0/1188 + r2*(-691.0/360360 + r2*(1.0/156)))))));
}

/* log(y) - 1/(2y) - digamma(y), for y >= 10 */
static double bb_digamma_e(double y)
{
  double r2 = 1.0/(y*y);
  return r2*(1.0/12 + r2*(-1.0/120 + r2*(1.0/252 + r2*(-1.0/240 + 
         r2*(1.0/132 + r2*(-691.0/32760 + r2*(1.0/12)))))));
}

/* trigamma(y) - 1/y - 1/(2 y^2), for y >= 10 */
static double bb_trigamma_r(double y)
{
  double r = 1.0/y, r2 = r*r;
  return r2*r*(1.0/6 + r2*(-1.0/30 + r2*(1.0/42 + r2*(-1.0/30 + 
         r2*(5.0/66 + r2*(-691.0/2730 + r2*(7.0/6)))))));
}

/* g(u) = (1 + u) log(1 + u) - u = sum_{j>=2} (-1)^j u^j/(j(j-1)) */
static double bb_g(double u)
{
  int j;
  double t, s;
  if (u >= 0.1) return (1.0 + u)*log1p(u) - u;
  s = 0.0;
  t = -u;
  for (j=2; j<60; j++) {
    t *= -u;
    s += t/(j*(j - 1.0));
    if (fabs(t) < 1e-18*fabs(s)) break;
  }
  return s;
}

/* h(u) = u - log(1 + u) = sum_{j>=2} (-1)^j u^j/j */
static double bb_h(double u)
{
  int j;
  double t, s;
  if (u >= 0.1) return u - log1p(u);
  s = 0.0;
  t = -u;
  for (j=2; j<60; j++) {
    t *= -u;
    s += t/j;
    if (fabs(t) < 1e-18*fabs(s)) break;
  }
  return s;
}

static int bb_closed_form_ok(double a, double theta, double m)
{
  return (theta > 0.0 && a > 0.0 && R_FINITE(theta) && R_FINITE(a) && 
          R_FINITE(m));
}

double bb_sum_log(double a, double theta, double m)
{
  int k;
  double x, u, s;
  if (m <= 0.0) return 0.0;
  if (!bb_closed_form_ok(a, theta, m)) {
    for (s=0.0, k=0; k<m; k++) s += log(a + k*theta);
    return s;
  }
  x = a/theta;
  if (x < BB_XMIN) {
    s = lgammafn(x + m) - lgammafn(x) - m*log(x);
  }else {
    u = m/x;
    s = x*bb_g(u) - 0.5*log1p(u) + bb_stirling_c(x + m) - bb_stirling_c(x);
  }
  return m*log(a) + s;
}

double bb_sum_inv(double a, double theta, double m)
{
  int k;
  double x, u, d;
  if (m <= 0.0) return 0.0;
  if (!bb_closed_form_ok(a, theta, m)) {
    for (d=0.0, k=0; k<m; k++) d += 1.0/(a + k*theta);
    return d;
  }
  x = a/theta;
  if (x < BB_XMIN) {
    d = digamma(x + m) - digamma(x);
  }else {
    u = m/x;
    d = log1p(u) + m/(2.0*x*(x + m)) - (bb_digamma_e(x + m) - bb_digamma_e(x));
  }
  return d/theta;
}

double bb_sum_kinv(double a, double theta, double m)
{
  int k;
  double x, u, q;
  /* the only term for m = 1 is k = 0 */
  if (m <= 1.0) return 0.0;
  if (!bb_closed_form_ok(a, theta, m)) {
    for (q=0.0, k=0; k<m; k++) q += k/(a + k*theta);
    return q;
  }
  x = a/theta;
  /* q = sum k/(x + k) = m - x (digamma(x + m) - digamma(x)) */
  if (x < BB_XMIN) {
    q = m - x*(digamma(x + m) - digamma(x));
  }else {
    u = m/x;
    q = x*bb_h(u) - m/(2.0*(x + m)) + x*(bb_digamma_e(x + m) - bb_digamma_e(x));
  }
  return q/theta;
}

double bb_sum_inv2(double a, double theta, double m)
{
  int k;
  double x, t;
  if (m <= 0.0) return 0.0;
  if (!bb_closed_form_ok(a, theta, m)) {
    for (t=0.0, k=0; k<m; k++) t += 1.0/((a + k*theta)*(a + k*theta));
    return t;
  }
  x = a/theta;
  if (x < BB_XMIN) {
    t = trigamma(x) - trigamma(x + m);
  }else {
    t = m/(x*(x + m)) + m*(2.0*x + m)/(2.0*x*x*(x + m)*(x + m)) 
        + bb_trigamma_r(x) - bb_trigamma_r(x + m);
  }
  return t/(theta*theta);
}

/**********************************************************************
 *
 * negative log likelihood and gradient function under H0
 *
 * H0: pi = 0.5
 *
 * The first parameter (n) is the number of paramters, n=1
 * Sample size N is included in the parameter "ex".
 *
 **********************************************************************/

double negLogH0 (int n, double* para, void* ex, SEXP x1){
  int i, N, h;
  double sumL, ni, ni0, pi0, piI, theta;
  double *exPara, *nA, *nTotal, *zeta;
  
  exPara = (double *) ex;
  N      = ceil(exPara[0]-0.5);
  h      = ceil(exPara[1]-0.5);
  pi0    = exPara[6];
  nA     = exPara + 7;
  nTotal = nA + N;
  zeta   = nTotal + N;
  
  theta  = para[0];
  
  sumL = 0.0;
  
  for(i=0; i<h; i++){
    ni   = nTotal[i];
    ni0  = nA[i];
    
    sumL += lchoose(ni, ni0);
    
    if(zeta[i] > 0){ piI = pi0; }else { piI = 0.5; }
    
    if(ni0 > 0){
      sumL += bb_sum_log(piI, theta, ni0);
    }
    
    if(ni0 < ni){
      sumL += bb_sum_log(1.0 - piI, theta, ni - ni0);
    }
    
    sumL -= bb_sum_log(1.0, theta, ni);
  }
  
  return(-sumL);
}


/**********************************************************************
 *
 * negGradLogH0A: gradient negLogH0
 *
 **********************************************************************/

void negGradLogH0(int n, double* para, double* gr, void* ex, SEXP x1)
{
  double grad, ni, ni0, pi0, piI, theta;
  int i, N, h;
  double *exPara, *nA, *nTotal, *zeta;
  
  exPara = (double *) ex;
  N = ceil(exPara[0]-0.5);
  h = ceil(exPara[1]-0.5);
  pi0    = exPara[6];
  nA     = exPara + 7;
  nTotal = nA + N;
  zeta   = nTotal + N;
  
  theta = para[0];

  grad = 0.0;
  
  for (i=0; i<h; i++) {
    ni  = nTotal[i];
    ni0 = nA[i];
    
    if(zeta[i] > 0){ piI = pi0; }else { piI = 0.5; }

    if(ni0 > 0){
      grad += bb_sum_kinv(piI, theta, ni0);
    }
    
    if(ni0 < ni){
      grad += bb_sum_kinv(1.0 - piI, theta, ni - ni0);
    }

    grad -= bb_sum_kinv(1.0, theta, ni);
  }
  
  gr[0] = -grad;
}

/**********************************************************************
 *
 * negative log likelihood and gradient function under H1
 *
 * H1: with alleleic imbalance, and pi_0 = 0.5
 *
 * Note the first parameter (n) is the number of paramters, n=2
 * Sample size N is included in the parameter "ex".
 *
 **********************************************************************/

double negLogH1 (int n, double* para, void* ex, SEXP x1){
  int i, N, h;
  double sumL, ni, ni0;
  double pi1, piI, theta;
  double *exPara, *nA, *nTotal, *zeta;
  
  exPara = (double *) ex;
  N      = ceil(exPara[0]-0.5);
  h      = ceil(exPara[1]-0.5);

  nA     = exPara + 7;
  nTotal = nA + N;
  zeta   = nTotal + N;
  
  theta = para[0];
  pi1   = para[1];

  sumL = 0.0;
  
  for(i=0; i<h; i++){
    ni   = nTotal[i];
    ni0  = nA[i];
    
    sumL += lchoose(ni, ni0);
    
    if(zeta[i] > 0){ piI = pi1; }else { piI = 0.5; }
  
    if(ni0 > 0){
      sumL += bb_sum_log(piI, theta, ni0);
    }
    
    if(ni0 < ni){
      sumL += bb_sum_log(1.0 - piI, theta, ni - ni0);
    }
    
    sumL -= bb_sum_log(1.0, theta, ni);
    
  }
  
  return(-sumL);
}

/**********************************************************************
 *
 * negGradLogH1: gradient negLogH1
 *
 **********************************************************************/

void negGradLogH1 (int n, double* para, double* gr, void* ex, SEXP x1)
{
  int i, N, h;
  double gradPi1, gradTh, pi1, piI, theta;
  double *exPara, *nA, *nTotal, *zeta, ni, ni0;
  
  exPara = (double *) ex;
  N      = ceil(exPara[0]-0.5);
  h      = ceil(exPara[1]-0.5);

  nA     = exPara + 7;
  nTotal = nA + N;
  zeta   = nTotal + N;

  theta = para[0];
  pi1   = para[1];
  
  gradPi1 = 0.0;
  gradTh  = 0.0;
  
  for(i=0; i<h; i++){
    ni   = nTotal[i];
    ni0  = nA[i];
    
    if(zeta[i] > 0){
      piI = pi1;
      
      if(ni0 > 0){
        gradPi1 += bb_sum_inv(piI, theta, ni0);
        gradTh  += bb_sum_kinv(piI, theta, ni0);
      }
      
      if(ni0 < ni){
        gradPi1 -= bb_sum_inv(1.0 - piI, theta, ni - ni0);
        gradTh  += bb_sum_kinv(1.0 - piI, theta, ni - ni0);
      }
      
    }else {
      piI = 0.5;
      
      if(ni0 > 0){
        gradTh += bb_sum_kinv(piI, theta, ni0);
      }
      
      if(ni0 < ni){
        gradTh += bb_sum_kinv(1.0 - piI, theta, ni - ni0);
      }
      
    }
    
    gradTh -= bb_sum_kinv(1.0, theta, ni);

  }
  
  gr[0] = -gradTh;
  gr[1] = -gradPi1;
}

/**********************************************************************
 *
 * ase
 *
 * allele specific expression. 
 *
  Input:
  
  Y1, Y2  Expression on two haplotypes
  Z       Covariates of interest
 **********************************************************************/


void ase (int* dims, double* Y1, double* Y2, double* Z, char** output, 
             double* RP_cut, int* cis_only, int* cis_distance, 
             int* eChr, int* ePos, int* mChr, int* mPos,  
             int* trace, int* succeed)
{
  int i, j, k, h0, h1;
  double chisq, pval, loglik0, loglik1, nT1, dfr;
  double pi0, th0, pi1, th1;
  double *exPara, *nA, *nTotal, *zeta;
  int nY = dims[0];
  int nZ = dims[1];
  int N  = dims[2];

  /* 
   * minimum of total number of reads to be considered 
   * If one sample has less than min_nT reads, it is disgarded
   */
  int min_nT = dims[3];
  
  /* 
   * minimum of sample size for testing
   */  
  int min_N  = dims[4];
  
  /* 
   * minimum of sample size for heterzygous genotypes
   */  
  int min_Nhet  = dims[5];
  
  /* 
   * p-value cutoff
   */  
  double P_cut = *RP_cut;
  
  /** 
   * we have to combine nA, nTotal, and zeta into one vector for it usage in 
   * sovling for MLE under H1
   */
  
  exPara = (double *) R_alloc(3*N+7, sizeof(double));
  nA     = exPara + 7;
  nTotal = nA + N;
  zeta   = nTotal + N;
  
  /* initial zeta */
  for (k=0; k<N; k++) {
    zeta[k] = -1.0;
  }
  
  /** 
   * exPara[0]  is sample size, fixed throught the computation 
   * exPara[1]  is sample size used in actual computation 
   *            which will be updated for each possible pair of gene and marker
   * exPara[6]  is pi, which will also be updated in each run
   */
  exPara[0] = (double) N; 
  exPara[1] = 0.0; 
  exPara[6] = 0.0; 

  /* pointers to Y and Z */
  double *pY1, *pY2, *pZ;
    
  /* position difference between gene and marker */
  int pos_diff;
  
  /* grid used to output frequency */
  double grid;
  
  /* 
   * p-value frequencies
   * freqs[100] = #{ pval < P_cut }
   * i = 0:99
   * freqs[i]   = #{ pval < [i/100, (i+1)/100) }
   */
  
  unsigned long freqs[101];
  for(i=0; i<=100; i++){ freqs[i] = 0; }
  
  /* 
   * variables for testing H0
   */
  double lower0, upper0, grH0;
  
  /* **********************************************************
   * parameters for function lbfgsb, which will be used to
     obtain MLE for H1: with allelic imbalance
   
   void lbfgsb(int n, int lmm, double *x, double *lower,
          double *upper, int *nbd, double *Fmin, optimfn fn,
          optimgr gr, int *fail, void *ex, double factr,
          double pgtol, int *fncount, int *grcount,
          int maxit, char *msg, int trace, int nREPORT);
   
   n:       the number of parameters
   
   lmm:     is an integer giving the number of BFGS updates 
            retained in the "L-BFGS-B" method, It defaults to 5.
   
   x:       starting parameters on entry and the final parameters on exit
   
   lower:   lower bounds
   
   upper:   upper bounds
   
   nbd:     specifies which bounds are to be used. 
            nbd(i)=0 if x(i) is unbounded,
            1 if x(i) has only a lower bound,
            2 if x(i) has both lower and upper bounds, and
            3 if x(i) has only an upper bound.
            On exit nbd is unchanged.
   
   Fmin:    final value of the function
   
   fn:      the function to be minimized
   
   gr:      the gradient function
   
   fail:    integer code, 0 for success, 51 for warning and 52 for error
   
   ex:      extra parameters for the function to be minimized
   
   factr:   controls the convergence of the "L-BFGS-B" method. 
            Convergence occurs when the reduction in the objective is 
            within this factor of the machine tolerance. Default is 1e7, 
            that is a tolerance of about 1e-8.
   
   pgtol:   helps control the convergence of the "L-BFGS-B" method. 
            It is a tolerance on the projected gradient in the current 
            search direction. This defaults to zero, when the check 
            is suppressed.
   
   fncount: the number of calls to fn 
   
   grcount: the number of calls to gr
   
   maxit:   maximum of iterations
   
   msg:     A character string giving any additional information 
            returned by the optimizer, or NULL
   
   trace:   Non-negative integer. If positive, tracing information 
            on the progress of the optimization is produced. 
            Higher values may produce more tracing information: 
            for method "L-BFGS-B" there are six levels of tracing. 
   
   nREPORT: The frequency of reports for the "BFGS", "L-BFGS-B" 
            and "SANN" methods if control$trace is positive. 
            Defaults to every 10 iterations for "BFGS" and "L-BFGS-B"
   
   * **********************************************************/
  
  int npara, lmm, fail, failA, failB, fncount, grcount, maxit, nREPORT;
  int nbd[2];

  npara   = 2;
  lmm     = 5;
  fail    = 0;
  failA   = 0;
  failB   = 0;
  fncount = 0;
  grcount = 0;
  maxit   = 100;
  nREPORT = 5;
  nbd[0]  = 1;
  nbd[1]  = 2;
  //technical parameters below:
  double *wa, *g1;
  int *iwa;
  SEXP x1;
  PROTECT(x1 = allocVector(REALSXP, npara));
  //consider replacing with simple Calloc
  wa  = (double *) S_alloc(2*lmm*npara+4*npara+11*lmm*lmm+8*lmm,sizeof(double));
  iwa = (int*) R_alloc(3*npara,sizeof(int));
  g1 = (double *)R_alloc(npara, sizeof(double));

  double gr[2];
  double initPara[2];
  double lower[2];
  double upper[2];
  double Fmin, factr, pgtol;
  
  /* initPara = c(theta, pi0, pi1) */
  
  initPara[0] = 0.1; 
  initPara[1] = 0.5;

  lower[0] = 0.0 + 1e-16;
  upper[0] = 1e16;
  
  lower[1] = 0.0 + 1e-16;
  upper[1] = 1.0 - 1e-16;

  factr = 1e7;
  pgtol = 0.0;
  
  char msg[1023];
  
  /**********************************************************/

  if(*trace){
    Rprintf("\n--------------------------------------------------\n");
    Rprintf("(nY, nZ, N) = (%d, %d, %d), trace=%d\n", nY, nZ, N, *trace);
    Rprintf("p.cut=%.4e, min.nTotal=%d, min.N=%d, min.Nhet=%d", 
            P_cut, min_nT, min_N, min_Nhet);
    Rprintf("\n--------------------------------------------------\n");
  }
  
  /* output file handles */
  FILE *fo, *ff;
  
  /* time records */
  time_t sec_s;
  time_t sec_e;
  
  /* starting time */
  sec_s = time(NULL);
  
  /* output file for the eQTL mapping results */
  fo = fopen (output[0], "w");
  
  /**
   * write out the header in the output file
   */
  fprintf(fo, "GeneRowID\tMarkerRowID\tTheta0\tlogLik0\t");
  fprintf(fo, "Theta1\tPi1\tlogLik1\tChisq\tPvalue\tdf\tN\tNhet\n");
  
  /***
   * identifify eQTL gene by gene
   */
  
  /* pY1/pY2 is the pointer to gene expression data */
  pY1 = Y1;
  pY2 = Y2;
  
  for(i=0; i<nY; i++,pY1+=N,pY2+=N){

    if(*trace > 1){
      Rprintf("\ni=%d\n", i);
    }
    
    /* *****************************************************
     * calculate nA and nTotal
     * *****************************************************/
    h0 = 0;
    
    for (k=0; k<N; k++) {
      nT1 = pY1[k] + pY2[k];
      
      if(nT1 < min_nT){ continue; }
      
      nTotal[h0] = nT1;
      nA[h0]     = pY2[k];
      
      h0++;
    }
    
    /* if sample size is not enough */
    if(h0 < min_N){ continue; }
        
    if(*trace > 1){
      Rprintf("\ni=%d, h0=%d\n", i, h0);
    }
    
    /* *****************************************************
     * obtain MLE for H0
     * situation A: assume pi_0 = 0.5
     * *****************************************************/
    
    exPara[1]   = (double) h0;
    exPara[6]   = 0.5; /* fixed pi0 */
    initPara[0] = 0.1; 
    initPara[1] = 0.5; /* value for pi0, not used for H0A */
    npara       = 1;
    
    //lbfgsb(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
    //       negLogH0, negGradLogH0, &failA, (void*)exPara, factr, pgtol,  
    //       &fncount, &grcount, maxit, msg, 0, nREPORT);
    lbfgsb1(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
           negLogH0, negGradLogH0, &failA, (void*)exPara, factr, pgtol,  
           &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1, x1);

    if (failA) {
      if (*trace)
        Rprintf("  i=%d, fail to fit baseline ASE model @ situation A\n", i);
    
      continue;
    }else{
      th0 = initPara[0];
      loglik0 = -Fmin;
      
      if(*trace > 1){
        Rprintf("\n  Obtained MLE for H0: twologlik0=%.4e", loglik0);
        Rprintf(" theta0A=%.4e", th0);
      }
    }
        
    /* *****************************************************
     * start eQTL mapping
     * *****************************************************/
    
    pZ = Z;
    
    for(j=0; j<nZ; j++,pZ+=N){
      
      if(*cis_only){
        if(eChr[i] != mChr[j]) continue;
        
        pos_diff = abs(ePos[i] - mPos[j]);
        
        if(pos_diff > *cis_distance) continue;
      }
            
      /* *****************************************************
       * calculate nA and nTotal
       * *****************************************************/
      h0 = 0;
      h1 = 0;
      
      for (k=0; k<N; k++) {
        nT1 = pY1[k] + pY2[k];
        
        if(nT1 < min_nT){ continue; }
        
        if(nTotal[h0] != nT1){
          error("mismatch ;( \n");
        }
        
        /* pZ[k] = 0 or 4 if homozygous */
        if (fabs(pZ[k] - 2.0) > 1.99) {
          zeta[h0] = -1.0;
          nA[h0]   = pY2[k];
        }else {
          zeta[h0] = 1.0;
          
          if(fabs(pZ[k] - 1.0) < 0.01){
            nA[h0] = pY2[k];
          }else if(fabs(pZ[k] - 3.0) < 0.01){
            nA[h0] = pY1[k];
          }else {
            error("invalid values for Z\n");
          }
          
          h1++;
        }
        h0 ++;
      }
            
      /* if sample size of heterzygous genotype is not enough */
      if(h1 < min_Nhet){ continue; }

      if(*trace > 1){
        Rprintf("\ni=%d, j=%d, h0=%d, h1=%d\n", i, j, h0, h1);
      }
      
      /* *****************************************************
       * obtain MLE for H1A: with allelic imbalance
       * *****************************************************/
      
      exPara[1]   = (double) h0;
      initPara[0] = th0;  /* theta */
      initPara[1] = 0.5;  /* pi1   */
      npara       = 2;
      
      //lbfgsb(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
      //       negLogH1, negGradLogH1, &fail, (void*)exPara, factr, pgtol,  
      //       &fncount, &grcount, maxit, msg, 0, nREPORT);
      lbfgsb1(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
             negLogH1, negGradLogH1, &fail, (void*)exPara, factr, pgtol,  
             &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1, x1);
      
      if (fail) {
        if(*trace){
          Rprintf("\n  i=%d, j=%d, h0=%d, h1=%d, fail for MLE of H1, fail=%d\n", 
                  i, j, h0, h1, fail);
        }
        continue;
      }        
      
      th1 = initPara[0];
      pi1 = initPara[1]; 
      loglik1 = -Fmin;
      
      if(*trace > 1){
        Rprintf("\n  Obtained MLE for H1: loglik1=%.4e fail=%d", loglik1, fail);
        Rprintf("\n  theta1=%.4e, pi1=%.4e", initPara[0], initPara[1]);
        negGradLogH1(2, initPara, gr, (void*)exPara, x1);
        Rprintf("\n  grad()=c(%.4e, %.4e)", gr[0], gr[1]);          
      }
            
      chisq = 2.0*(loglik1 - loglik0);

      if (chisq < -1e-5) {
        error("wrong loglik! i=%d, j=%d, loglik=(%.4e, %.4e)\n", 
              i, j, loglik0, loglik1);
      }
            
      dfr = 1.0;
      if (fabs(th0) < 1e-7)  dfr += 1.0;
      if (fabs(th1) < 1e-7)  dfr -= 1.0;

      if(fabs(dfr) < 0.01){
        if(*trace > 1){
          Rprintf("\n  i=%d, j=%d, h0=%d, h1=%d ", i, j, h0, h1);
          Rprintf("dfr=0, theta1=%.4e, pi1=%.4e\n", th1, pi1);
        }
        continue;
      }
      
      if (chisq < 1e-5) { 
        pval = 1.0; 
      }else{
        pval = pchisq(chisq, dfr, 0, 0);
      }
      
      k = (int) (pval / 0.01);
      freqs[k] += 1;
      
      if(pval < P_cut){
        freqs[100] += 1;
        
        /* gene ID and SNP ID */
        fprintf(fo, "%d\t%d\t%e\t%e\t", i+1, j+1, th0, loglik0);
        fprintf(fo, "%e\t%e\t%e\t", th1, pi1, loglik1);
        fprintf(fo, "%.3f\t%.2e\t%.0f\t%d\t%d\n", chisq, pval, dfr, h0, h1);
      }
      
    }
  }
  
  fclose(fo);
  
  // print out the frequencies
  ff   = fopen(output[1], "w");
  grid = 0.0;
  
  fprintf(ff, "<%.2e", P_cut);
  for(i=0; i<100; i++){
    fprintf(ff, "\t%.2f-%.2f", grid, grid+0.01);
    grid += 0.01;
  }
  fprintf(ff, "\n");
  
  fprintf(ff, "%lu", freqs[100]);
  for(i=0; i<100; i++){
    fprintf(ff, "\t%lu", freqs[i]);
  }
  fprintf(ff, "\n");  
  fclose(ff);
  UNPROTECT(1);
  /* end time */
  sec_e  = time(NULL);
  if(*trace){
    Rprintf("\ntotal time spent in ase is %ld secs\n", sec_e-sec_s);
    Rprintf("------------------------------------------------------\n");
  }
    
  *succeed = 1;
}
