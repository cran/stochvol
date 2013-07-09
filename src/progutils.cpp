#include "progutils.h"

/* Contains the following code modules:

   a) Some helper functions such as progress bar tools and return value
      prettifier

   b) Some functions related to the Cholesky decomposition used for 
      sampling AWOL and efficiently solving the systems of linear
      equations

   c) function for inverse transform sampling
*/

// a)
// Sums up results and prepares return value
Rcpp::List cleanUp(Rcpp::NumericVector mu, Rcpp::NumericVector phi, Rcpp::NumericVector sigma, Rcpp::NumericMatrix hstore, Rcpp::NumericVector h0store) {
 Rcpp::NumericMatrix res(mu.length(), 3); 
 res(Rcpp::_,0) = mu;
 res(Rcpp::_,1) = phi;
 res(Rcpp::_,2) = sigma;

/* res.attr("dimnames") = Rcpp::List::create(
   R_NilValue, 
   Rcpp::CharacterVector::create("mu", "phi", "sigma")); */

 Rcpp::List val = Rcpp::List::create(
   Rcpp::_["para"] = res,
   Rcpp::_["latent"] = hstore,
   Rcpp::_["latent0"] = h0store);

 return val;
}

// sets up the progress bar
int progressbar_init(int N) {
 int show, tmp;
 REprintf("\n      ");
 if (N >= 2500) {
  for (int i = 0; i < 50+1; i++) REprintf(" ");
  show = N/50;
 }
 else {
  for (int i = 0; i < (N-1)/50+1; i++) REprintf(" ");
  show = 50;
 }
 REprintf("] 100%");
 REprintf("\r  0% [");
// R_FlushConsole();
 return show;
}

// finalizes progress bar
void progressbar_finish(int N) {
 if (!(N % 50) && N >= 2500) REprintf("+");
 REprintf("] 100%\n\n");
// R_FlushConsole();
}

// b)
// Cholesky factor for a tridiagonal matrix with constant off-diagonal
void cholTridiag(Rcpp::NumericVector omega_diag, double omega_offdiag, double * chol_diag, double * chol_offdiag)
{
 chol_diag[0] = sqrt(omega_diag[0]);  // maybe speed up via iterators?
 for (int j = 1; j < omega_diag.length(); j++) {
  chol_offdiag[j-1] = omega_offdiag/chol_diag[j-1];
  chol_diag[j] = sqrt(omega_diag[j]-chol_offdiag[j-1]*chol_offdiag[j-1]);
 }
}

// Solves Chol*x = covector ("forward algorithm")
void forwardAlg(Rcpp::NumericVector chol_diag, Rcpp::NumericVector chol_offdiag, Rcpp::NumericVector covector, double * htmp)
{
 htmp[0] = covector[0]/chol_diag[0];
 for (int j = 1; j < chol_diag.length(); j++) {
  htmp[j] = (covector[j] - chol_offdiag[j-1]*htmp[j-1])/chol_diag[j];
 }
}

// Solves (Chol')*x = htmp ("backward algorithm")
void backwardAlg(Rcpp::NumericVector chol_diag, Rcpp::NumericVector chol_offdiag, Rcpp::NumericVector htmp, double * h)
{
 int T = chol_diag.length();
 h[T-1] = htmp[T-1]/chol_diag[T-1];
 for (int j = T-2; j >= 0; j--) {
  h[j] = (htmp[j] - chol_offdiag[j]*h[j+1])/chol_diag[j];
 }
}

// c)
// draws length(r) RVs, expects the non-normalized CDF mixprob
void invTransformSampling(Rcpp::NumericMatrix mixprob, int * r) {
 int T = mixprob.ncol(), rows = mixprob.nrow(), index;
 Rcpp::NumericVector innov = Rcpp::runif(T); 
 double temp;
 bool larger, smaller;
 for (int j = 0; j < T; j++) {
  index = (rows-1)/2;  // start searching in the middle
  temp = innov[j]*mixprob(rows-1, j);  // current (non-normalized) value
  larger = false;  // indicates that we already went up
  smaller = false; // indicates that we already went down
  while(true) {
   if (temp > mixprob(index, j)) {
    if (smaller == true) {
     index++;
     break;
    }
    else {
    index++;
    larger = true;
    }
   }
   else {
    if (larger == true) {
     break;
    }
    else {
     if (index == 0) {
      break;
     }
     else {
      index--;
      smaller = true;
     }
    } 
   }
  }
 r[j] = index;
 }
}
