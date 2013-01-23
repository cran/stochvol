#ifndef _PROGUTILS_H_
#define _PROGUTILS_H_

/* Contains two code modules:

   a) Some helper functions such as progress bar tools and return value
      prettifier

   b) Some functions related to the Cholesky decomposition used for 
      sampling AWOL and efficiently solving the systems of linear
      equations

   c) function for inverse transform sampling
*/

#include <Rcpp.h>

// a)
// Sums up results and prepares return value
Rcpp::List cleanUp(Rcpp::NumericVector mu, Rcpp::NumericVector phi, Rcpp::NumericVector sigma, Rcpp::NumericMatrix hstore, Rcpp::NumericVector h0store);

// sets up the progress bar
int progressbar_init(int N);

// adds one '+' to progress bar
inline void progressbar_print() {
 Rprintf("+");
}

// finalizes progress bar
void progressbar_finish(int N);

inline void store_h(double * h, double * hstore, int timethin, int hstorelength, double h0, double * h0store) {
 for (int j = 0; j < hstorelength; j++) hstore[j] = h[timethin*j];
 *h0store = h0;
}


// b)
// Cholesky factor for a tridiagonal matrix with constant off-diagonal
void cholTridiag(Rcpp::NumericVector omega_diag, double omega_offdiag, double * chol_diag, double * chol_offdiag);

// Solves Chol*x = covector ("forward algorithm")
void forwardAlg(Rcpp::NumericVector chol_diag, Rcpp::NumericVector chol_offdiag, Rcpp::NumericVector covector, double * htmp);

// Solves (Chol')*x = htmp ("backward algorithm")
void backwardAlg(Rcpp::NumericVector chol_diag, Rcpp::NumericVector chol_offdiag, Rcpp::NumericVector htmp, double * h);

// c)
// draws length(r) RVs, expects the non-normalized CDF mixprob
void invTransformSampling(Rcpp::NumericMatrix mixprob, int * r);

#endif
