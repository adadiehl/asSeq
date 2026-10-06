/*
 *  ase.h
 *
 *  Created by Wei Sun on 5/25/2010.
 *  Modified by Vasyl Zhabotynsky on 04/05/2011
 *
 */

#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include "utility.h"
#include "lbfgsb1.h"

/* closed forms of the beta-binomial sums over k = 0, ..., m-1 */
double bb_sum_log(double a, double theta, double m);   /* sum log(a + k theta) */
double bb_sum_inv(double a, double theta, double m);   /* sum 1/(a + k theta)  */
double bb_sum_kinv(double a, double theta, double m);  /* sum k/(a + k theta)  */
double bb_sum_inv2(double a, double theta, double m);  /* sum 1/(a + k theta)^2 */

double negLogH0 (int n, double* para, void* ex, SEXP x1);
double negLogH1 (int n, double* para, void* ex, SEXP x1);

void negGradLogH0 (int n, double* para, double* gr, void* ex, SEXP x1);
void negGradLogH1 (int n, double* para, double* gr, void* ex, SEXP x1);

void ase (int* dims, double* Y1, double* Y2, double* Z, char** output, 
          double* RP_cut, int* cis_only, int* cis_distance, 
          int* eChr, int* ePos, int* mChr, int* mPos,  
          int* trace, int* succeed);
