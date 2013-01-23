#ifndef _SAMPLER_H_
#define _SAMPLER_H_

// Main sampling steps and helper functions

#include <Rcpp.h>
#include "auxmix.h"
#include "progutils.h"
#include "densities.h"

// Main sampler (as called from R):
RcppExport SEXP sampler(const SEXP, const SEXP, const SEXP,
  const SEXP, const SEXP, const SEXP, const SEXP, const SEXP,
  const SEXP, const SEXP, const SEXP, const SEXP, const SEXP,
  const SEXP, const SEXP, const SEXP, const SEXP, const SEXP,
  const SEXP, const SEXP);

// Step (b): sample mu, phi, sigma - __CENTERED__ version:
Rcpp::NumericVector regressionCentered(
       double h0, Rcpp::NumericVector h,
       double mu, double phi, double sigma,
       double C0, double cT, double Bsigma,
       double a0, double b0,
       double bmu, double Bmu,
       double B011inv, double B022inv,
       bool gammaprior, bool truncnormal, double MHcontrol, int MHsteps);

// Step (b): sample mu, phi, sigma - __NONCENTERED__ version:
Rcpp::NumericVector regressionNoncentered(
       Rcpp::NumericVector data,
       double h0, Rcpp::NumericVector h, Rcpp::IntegerVector r,
       double mu, double phi, double sigma,
       double Bsigma, double a0, double b0,
       double bmu, double Bmu,
       bool truncnormal, int MHsteps);

#endif
