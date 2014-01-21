#include <R_ext/Rdynload.h>

#if defined(COMPILING_STOCHVOL)

// just declare, will be defined in sampler.cpp
// a single MCMC update:
void update(const Rcpp::NumericVector &, Rcpp::NumericVector &,
            Rcpp::NumericVector &, double &, Rcpp::NumericMatrix &,
	    Rcpp::IntegerVector &, const bool, const double,
	    const double, const double, const double, const double,
	    const double, const double, const double, const double,
	    const bool, const bool, const double, const int, const int);

#else

inline void update(const Rcpp::NumericVector &data, Rcpp::NumericVector &curpara,
            Rcpp::NumericVector &h, double &h0, Rcpp::NumericMatrix &mixprob,
	    Rcpp::IntegerVector &r, const bool centered_baseline, const double C0,
	    const double cT, const double Bsigma, const double a0, const double b0,
	    const double bmu, const double Bmu, const double B011inv, const double B022inv,
	    const bool Gammaprior, const bool truncnormal, const double MHcontrol,
	    const int MHsteps, const int parameterization) {
 
 typedef void(*Update)(const Rcpp::NumericVector &, Rcpp::NumericVector &,
            Rcpp::NumericVector &, double &, Rcpp::NumericMatrix &,
	    Rcpp::IntegerVector &, const bool, const double,
	    const double, const double, const double, const double,
	    const double, const double, const double, const double,
	    const bool, const bool, const double, const int, const int);
 
 static Update update = (Update)R_GetCCallable("stochvol", "update");

 update(data, curpara, h, h0, mixprob, r, centered_baseline, C0,
             cT, Bsigma, a0, b0, bmu, Bmu, B011inv, B022inv, Gammaprior,
             truncnormal, MHcontrol, MHsteps, parameterization);
}

#endif
