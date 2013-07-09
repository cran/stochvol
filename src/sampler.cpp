// Main sampling steps and helper functions

#include "sampler.h"

using namespace Rcpp; // avoid to type Rcpp:: every time
RNGScope scope;       // just in case no seed has been set at R level

// RcppExport is an alias for 'extern "C"'
RcppExport SEXP sampler(const SEXP y_in, const SEXP draws_in,
  const SEXP burnin_in, const SEXP bmu_in, const SEXP Bmu_in,
  const SEXP a0_in, const SEXP b0_in, const SEXP Bsigma_in,
  const SEXP thin_in, const SEXP timethin_in, const SEXP startpara_in,
  const SEXP startvol_in, const SEXP quiet_in, const SEXP para_in,
  const SEXP MHsteps_in, const SEXP B011_in, const SEXP B022_in,
  const SEXP mhcontrol_in, const SEXP gammaprior_in,
  const SEXP truncnormal_in) {

 // convert SEXP into Rcpp-structures (AFAIK no copy at this point)
 NumericVector y(y_in), startvol(startvol_in);
 List startpara(startpara_in);

 // length of time series
 int T = y.length();
 
 // number of MCMC draws
 int burnin = as<int>(burnin_in);
 int draws  = as<int>(draws_in);
 int N 	    = burnin + draws;
 
 // prior parameters
 double bmu    = as<double>(bmu_in);
 double Bmu    = as<double>(Bmu_in);
 double a0     = as<double>(a0_in);
 double b0     = as<double>(b0_in);
 double Bsigma = as<double>(Bsigma_in);
 
 // thinning parameters
 int timethin = as<int>(timethin_in);
 int thin     = as<int>(thin_in);

 // verbosity control
 bool verbose = !as<bool>(quiet_in);

 // "expert" settings:
 double B011inv         = 1/as<double>(B011_in);
 double B022inv         = 1/as<double>(B022_in);
 bool Gammaprior        = as<bool>(gammaprior_in);
 bool truncnormal       = as<bool>(truncnormal_in);
 double MHcontrol       = as<double>(mhcontrol_in);
 int MHsteps            = as<int>(MHsteps_in);
 int parameterization   = as<int>(para_in);
 bool centered_baseline = as<int>(para_in) % 2; // 0 for C, 1 for NC baseline

 // moment-matched IG-prior
 double c0 = 2.5;
 double C0 = 1.5*Bsigma;
 
 // pre-calculation of a posterior parameter
 double cT; 
 if (Gammaprior) {  
  if (MHsteps == 2 || MHsteps == 3) cT = T/2.0; // we want IG(-.5,0) as proposal
  else if (MHsteps == 1) cT = (T-1)/2.0; // we want IG(-.5,0) as proposal
 }
 else {
  if (MHsteps == 2) cT = c0 + (T+1)/2.0;  // pre-calculation outside the loop
  else return LogicalVector::create(NA_LOGICAL);  // not implemented!
 }
 
 // initialize the variables:
 NumericVector sigma2inv(N+1, pow(as<double>(startpara["sigma"]), -2));
 NumericVector phi(N+1, as<double>(startpara["phi"]));
 NumericVector mu(N+1, as<double>(startpara["mu"]));
 
 NumericVector h(T);  // contains h1 to hT, but not h0!
 for (int i = 0; i < T; i++) h[i] = startvol[i];  // need to make a manual copy (!)
 if (!centered_baseline) h = (h-mu[0])*sqrt(sigma2inv[0]);
 
 NumericVector h_alter(T);  // contains h1 to hT in alternative para
 double h0, h0_alter;
 
 IntegerVector r(T);  // mixture indicators
 
 int hstorelength = (T-1)/timethin+1;
 NumericMatrix hstore(hstorelength, draws/thin);
 NumericVector h0store(draws/thin);

 NumericMatrix mixprob(10, T);  // mixture probabilities
 NumericVector data = log(y*y);  // commonly used transformation
 NumericVector omega_diag(T);  // contains diagonal elements of precision matrix
 double omega_offdiag;  // contains off-diag element of precision matrix (const)
 NumericVector chol_offdiag(T-1), chol_diag(T);  // Cholesky-factor of Omega
 NumericVector covector(T);  // holds covector (see McCausland et al. 2011)
 NumericVector htmp(T);  // intermediate vector for sampling h
 NumericVector newvals(3);  // holds mu, phi, sigma in every iteration
 newvals[0] = mu[0];
 newvals[1] = phi[0];
 newvals[2] = 1/sqrt(sigma2inv[0]);

 double sigmainvtmp, phi2tmp;

 // initializes the progress bar
 // "show" holds the number of iterations per progress sign
 int show;
 if (verbose) show = progressbar_init(N);
 
 for (int i = 0; i < N; i++) {  // main MCMC loop

  // print a progress sign every "show" iterations
  if (verbose) if (!(i % show)) progressbar_print();

  /*
   * Step (c): sample indicators
   */

  // calculate non-normalized CDF of posterior indicator probs
  if (centered_baseline) findMixCDF(&mixprob(0,0), data-h);
  else findMixCDF(&mixprob(0,0), data-newvals[0]-newvals[2]*h); 

  // find correct indicators (currently by inversion method)
  invTransformSampling(mixprob, &r(0));

  /*
   * Step (a): sample the latent volatilities h:
   */
  
  if (centered_baseline) { // fill precision matrix omega and covector c for CENTERED para:
   omega_diag[0] = mix_varinv[r[0]] + sigma2inv[i];
   covector[0] = (data[0] - mix_mean[r[0]])*mix_varinv[r[0]]
               + mu[i]*(1-phi[i])*sigma2inv[i];
   for (int j = 1; j < (T-1); j++) {
    omega_diag[j] = mix_varinv[r[j]] + (1+phi[i]*phi[i])*sigma2inv[i]; 
    covector[j] = (data[j] - mix_mean[r[j]])*mix_varinv[r[j]]
                + mu[i]*pow((1-phi[i]),2)*sigma2inv[i];
   }
   omega_diag[T-1] = mix_varinv[r[T-1]] + sigma2inv[i];
   covector[T-1] = (data[T-1] - mix_mean[r[T-1]])*mix_varinv[r[T-1]]
                 + mu[i]*(1-phi[i])*sigma2inv[i];
   omega_offdiag = -phi[i]*sigma2inv[i];  // omega_offdiag is constant
  
  } else { // fill precision matrix omega and covector c for NONCENTERED para:

   sigmainvtmp = sqrt(sigma2inv[i]);
   phi2tmp = phi[i]*phi[i]; 
   omega_diag[0] = mix_varinv[r[0]]/sigma2inv[i] + 1;
   covector[0] = mix_varinv[r[0]]/sigmainvtmp*(data[0] - mix_mean[r[0]] - mu[i]);
   for (int j = 1; j < (T-1); j++) {
    omega_diag[j] = mix_varinv[r[j]]/sigma2inv[i] + 1 + phi2tmp; 
    covector[j] = mix_varinv[r[j]]/sigmainvtmp*(data[j] - mix_mean[r[j]] - mu[i]);
   }
   omega_diag[T-1] = mix_varinv[r[T-1]]/sigma2inv[i] + 1;
   covector[T-1] = mix_varinv[r[T-1]]/sigmainvtmp*(data[T-1] - mix_mean[r[T-1]] - mu[i]);
   omega_offdiag = -phi[i];  // omega_offdiag is constant
  } 

  // Cholesky decomposition
  cholTridiag(omega_diag, omega_offdiag, &chol_diag(0), &chol_offdiag(0));
 
  // Solution of Chol*x = covector ("forward algorithm")
  forwardAlg(chol_diag, chol_offdiag, covector, &htmp(0));

  htmp = htmp + rnorm(T);

  // Solution of (Chol')*x = htmp ("backward algorithm")
  backwardAlg(chol_diag, chol_offdiag, htmp, &h(0));

  // sample h0 | h1, phi, mu, sigma
  if (centered_baseline) h0 = as<double>(rnorm(1, mu[i] + phi[i]*(h[0]-mu[i]),
                                         1/sqrt(sigma2inv[i])));
  else h0 = as<double>(rnorm(1, phi[i]*h[0], 1)); 
  /*
   * Step (b): sample mu, phi, sigma
   */
  
  if (centered_baseline) {  // this means we have C as base
   newvals = regressionCentered(h0, h, mu[i], phi[i],
     sqrt(1/sigma2inv[i]), C0, cT, Bsigma, a0, b0, bmu, Bmu, B011inv,
     B022inv, Gammaprior, truncnormal, MHcontrol, MHsteps);
   
   if (parameterization == 3) {  // this means we should interweave
    h_alter = (h-newvals[0])/newvals[2];
    h0_alter = (h0-newvals[0])/newvals[2];
    newvals = regressionNoncentered(data, h0_alter, h_alter, r,
      newvals[0], newvals[1], newvals[2], Bsigma, a0, b0, bmu, Bmu,
      truncnormal, MHsteps);
    h = newvals[0] + newvals[2]*h_alter;
    h0 = newvals[0] + newvals[2]*h0_alter;
   }
   
   if (!((i+1) % thin)) if (i >= burnin) {  // this means we should store h
    store_h(&h(0), &hstore(0, (i-burnin)/thin), timethin, hstorelength,
            h0, &h0store((i-burnin)/thin));
   }
  
  } else {  // NC as base
  
   newvals = regressionNoncentered(data, h0, h, r, mu[i], phi[i], sqrt(1/sigma2inv[i]),
                                   Bsigma, a0, b0, bmu, Bmu, truncnormal, MHsteps);
   if (parameterization == 4) {  // this means we should interweave
    h_alter = newvals[0] + newvals[2]*h;
    h0_alter = newvals[0] + newvals[2]*h0;
    newvals = regressionCentered(h0_alter, h_alter, newvals[0], newvals[1], newvals[2],
                                 C0, cT, Bsigma, a0, b0, bmu, Bmu, B011inv, B022inv,
				 Gammaprior, truncnormal, MHcontrol, MHsteps);
    h = (h_alter-newvals[0])/newvals[2];
    h0 = (h0_alter-newvals[0])/newvals[2];
   }
  
   if (!((i+1) % thin)) if (i >= burnin) {  // this means we should store h
    h_alter = newvals[0] + newvals[2]*h;
    h0_alter = newvals[0] + newvals[2]*h0;
    store_h(&h_alter(0), &hstore(0, (i-burnin)/thin), timethin, hstorelength,
            h0_alter, &h0store((i-burnin)/thin));
   }
  }
  
  mu[i+1] = newvals[0];
  phi[i+1] = newvals[1];
  sigma2inv[i+1] = 1/(newvals[2]*newvals[2]);
 }

 if (verbose) progressbar_finish(N);  // finalize progress bar

 // Prepare return value and return
 return cleanUp(mu, phi, sqrt(1/sigma2inv), hstore, h0store);
}


// Step (b): sample mu, phi, sigma - __CENTERED__ version:
Rcpp::NumericVector regressionCentered(
       double h0, Rcpp::NumericVector h,
       double mu, double phi, double sigma,
       double C0, double cT, double Bsigma,
       double a0, double b0,
       double bmu, double Bmu,
       double B011inv, double B022inv,
       bool gammaprior, bool truncnormal, double MHcontrol, int MHsteps) {
 
 int T = h.length();
 double R, z, CT, sigma2_prop, sigma_prop, sum1, sum2, sum3, sum4, tmp1,
       	tmp2, BT11, BT12, BT22, bT1, bT2, chol11, chol12, chol22, phi_prop,
       	gamma_prop, tmpR, tmpR2, logR;
 Rcpp::NumericVector innov(2);
 Rcpp::NumericVector quant(2);

 // first calculate bT and BT:
 sum1 = h[0];
 sum3 = h0*h[0];
 sum4 = h[0]*h[0];
 for (int j = 1; j < T-1; j++) {
  sum1 += h[j];
  sum3 += h[j-1]*h[j];
  sum4 += h[j]*h[j];
 }
 sum2 = sum1 + h[T-1];  // h_1 + h_2 + ... + h_T
 sum1 += h0;            // h_0 + h_1 + ... + h_{T-1}
 sum3 += h[T-2]*h[T-1]; // h_0*h_1 + h_1*h_2 + ... + h_{T-1}*h_T
 sum4 += h0*h0;         // h_0^2 + h_1^2 + ... + h_{T-1}^2

 tmp1 = 1/(((sum4 + B011inv)*(T+B022inv)-sum1*sum1));
 BT11 = (T + B022inv)*tmp1;
 BT12 = -sum1*tmp1;
 BT22 = (sum4+B011inv)*tmp1;
 
 bT1 = BT11*sum3 + BT12*sum2;
 bT2 = BT12*sum3 + BT22*sum2;

 // draw sigma^2 
 if (MHsteps == 2 || MHsteps == 3) { // draw sigma^2 from full conditional
  z = pow(((h[0]-mu)-phi*(h0-mu)),2);  // TODO: more efficiently via sum1, sum2, etc.
  for (int j = 0; j < (T-1); j++) {
   z += pow((h[j+1]-mu)-phi*(h[j]-mu),2);
  }
  z = (z+(h0-mu)*(h0-mu)*(1-phi*phi));
  if (MHcontrol > 0) {  // let's do a log normal random walk
   sigma2_prop = exp(Rcpp::rnorm(1,log(sigma*sigma), MHcontrol)[0]);
   logR = logacceptrateRW(sigma2_prop, sigma*sigma, Bsigma, T, z);

//  REprintf("\nold: %f, new: %f, R: %f", sigma, sqrt(sigma2_prop), exp(logR));
   if (log(Rcpp::runif(1)[0]) < logR) sigma = sqrt(sigma2_prop);
  }
  else {  // either IG(-.5,0)-proposal or IG(1.5,1.5*Bsigma)-prior
   if (gammaprior) {
    CT = .5*z;
    sigma2_prop = 1/Rcpp::as<double>(Rcpp::rgamma(1, cT, 1/CT));
    if (log(Rcpp::as<double>(Rcpp::runif(1))) <
        logacceptrateGamma(sigma2_prop, sigma*sigma, Bsigma)) {
     sigma = sqrt(sigma2_prop);
    }
   } else {
     CT = C0+.5*z;
     sigma = sqrt(1/Rcpp::as<double>(Rcpp::rgamma(1, cT, 1/CT)));
   }
  }
 } else if (MHsteps == 1) {  // draw sigma^2 marginalized over gamma, phi
  if (gammaprior) {
   CT = .5*((sum4 - h0*h0 + h[T-1]*h[T-1]) - bT1*sum3 - bT2*sum2);
   sigma2_prop = 1/Rcpp::as<double>(Rcpp::rgamma(1, cT, 1/CT));
  }
 }
 
 // sampling of "betas" (i.e. phi and gamma)
 if (MHsteps == 3) {
  
  // sampling of phi from full conditional:
  double gamma = (1-phi)*mu;
  double BTsqrt = sigma/sqrt(sum4+B011inv);
  double bT = (sum3-gamma*sum1)/(sum4+B011inv);
  phi_prop = as<double>(Rcpp::rnorm(1, bT, BTsqrt));
  
  R = logdnorm(h0, mu, sigma/sqrt(1-phi_prop*phi_prop));
  R -= logdnorm(h0, mu, sigma/sqrt(1-phi*phi));
  R += logdbeta((phi_prop+1)/2, a0, b0);
  R -= logdbeta((phi+1)/2, a0, b0);
  R += logdnorm(phi, 0, sigma/sqrt(B011inv));
  R -= logdnorm(phi_prop, 0, sigma/sqrt(B011inv));

  if (log(Rcpp::as<double>(Rcpp::runif(1))) < R) {
   phi = phi_prop;
  }
 
  // sampling of gamma from full conditional:
  gamma = (1-phi)*mu;
  BTsqrt = sigma/sqrt(T+B022inv);
  bT = (sum2-phi*sum1)/(T+B022inv);
  gamma_prop = as<double>(Rcpp::rnorm(1, bT, BTsqrt));
  
  R = logdnorm(h0, gamma_prop/(1-phi), sigma/sqrt(1-phi*phi));
  R -= logdnorm(h0, gamma/(1-phi), sigma/sqrt(1-phi*phi));
  R += logdnorm(gamma_prop, bmu*(1-phi), sqrt(Bmu)*(1-phi));
  R -= logdnorm(gamma, bmu*(1-phi), sqrt(Bmu)*(1-phi));
  R += logdnorm(gamma, 0, sigma/sqrt(B022inv));
  R -= logdnorm(gamma_prop, 0, sigma/sqrt(B022inv));
  
  if (log(Rcpp::as<double>(Rcpp::runif(1))) < R) {
   mu = gamma_prop/(1-phi);
  }
 }
 else { 
  // Some more temps needed for sampling the betas
  chol11 = sqrt(BT11);
  chol12 = (BT12/chol11);
  chol22 = sqrt(BT22-chol12*chol12);
  chol11 *= sigma;
  chol12 *= sigma;
  chol22 *= sigma;

  if (truncnormal) { // draw from truncated normal via inversion method
   quant = Rcpp::pnorm(Rcpp::NumericVector::create(-1,1), bT1, chol11);
   phi_prop = (Rcpp::qnorm((quant[0] + Rcpp::runif(1)*(quant[1]-quant[0])),
               bT1, chol11))[0];
   gamma_prop = (Rcpp::rnorm(1, bT2 + chol12*((phi_prop-bT1)/chol11),
                 chol22))[0];
  }
  else { // draw from normal and reject (return) if outside
   innov = Rcpp::rnorm(1);
   phi_prop = bT1 + chol11*innov[0];
   if ((phi_prop >= 1) || (phi_prop <= -1)) { // outside the unit sphere
    Rcpp::NumericVector ret = Rcpp::NumericVector::create(mu, phi, sigma);
    return ret;
   }
   else gamma_prop = bT2 + chol12*innov[0] + chol22*Rcpp::rnorm(1)[0];
  }

  // acceptance probability exp(R) calculated on a log scale
  tmpR = 1-phi_prop;  // some temps used for acceptance probability
  tmpR2 = 1-phi;
 
  if (MHsteps == 2) {
   sigma_prop = sigma;  // sigma was accepted/rejected independently
   R = 0;  // initialize R
  } else if (MHsteps == 1) {
   sigma_prop = sqrt(sigma2_prop);  // accept sigma jointly with "betas"
   R = logacceptrateGamma(sigma2_prop, sigma*sigma, Bsigma);  // initialize R
  }
 
  R += logdnorm(h0, gamma_prop/tmpR, sigma_prop/sqrt(1-phi_prop*phi_prop));
  R -= logdnorm(h0, mu, sigma/sqrt(1-phi*phi));
  R += logdnorm(gamma_prop, bmu*tmpR, sqrt(Bmu)*tmpR);
  R -= logdnorm(mu*tmpR2, bmu*tmpR2, sqrt(Bmu)*tmpR2);
  R += logdbeta((phi_prop+1)/2, a0, b0);
  R -= logdbeta((phi+1)/2, a0, b0);
  R += logdnorm(phi, 0, sigma/sqrt(B011inv));
  R -= logdnorm(phi_prop, 0, sigma_prop/sqrt(B011inv));
  R += logdnorm(mu*tmpR2, 0, sigma/sqrt(B011inv));
  R -= logdnorm(gamma_prop, 0, sigma_prop/sqrt(B011inv));

  // accept/reject
  if (log(Rcpp::as<double>(Rcpp::runif(1))) < R) {
   mu = gamma_prop/(1-phi_prop);
   phi = phi_prop;
   if (MHsteps == 1) sigma = sigma_prop;
  }
 }

 Rcpp::NumericVector ret = Rcpp::NumericVector::create(mu, phi, sigma);
 return ret;
}

// Step (b): sample mu, phi, sigma - __NONCENTERED__ version:
Rcpp::NumericVector regressionNoncentered(
       Rcpp::NumericVector data,
       double h0, Rcpp::NumericVector h, Rcpp::IntegerVector r,
       double mu, double phi, double sigma,
       double Bsigma, double a0, double b0,
       double bmu, double Bmu,
       bool truncnormal, int MHsteps) {
 
 int T = h.length();
 double sumtmp1, sumtmp2, expR, phi_prop, BT11, BT12, BT22, bT1, bT2, tmp1,
       	tmp2, tmp3, tmp, chol11, chol12, chol22, tmpmean, tmpsd;
 Rcpp::NumericVector innov(2);
 Rcpp::NumericVector quant(2);

 if (MHsteps == 3) {  // Gibbs-sample mu|sigma,... and sigma|mu,...
  
  // first, draw sigma from the full conditional posterior:
  tmp1 = 0; 
  tmp2 = 0; 
  for (int j = 0; j < T; j++) {
   tmp1 += h[j]*h[j]*mix_varinv[r[j]];
   tmp2 += h[j]*(data[j]-mix_mean[r[j]]-mu)*mix_varinv[r[j]];
  }
  BT11 = 1/(tmp1+1/Bsigma);
  bT1 = BT11*tmp2; 
//  REprintf("old: %f, new: mean %f and sd %f\n", sigma, bT1, sqrt(BT11));
  sigma = as<double>(Rcpp::rnorm(1, bT1, sqrt(BT11)));

  // second, draw mu from the full conditional posterior:
  tmp1 = 0; 
  tmp2 = 0; 
  for (int j = 0; j < T; j++) {
   tmp1 += mix_varinv[r[j]];
   tmp2 += (data[j]-mix_mean[r[j]]-sigma*h[j])*mix_varinv[r[j]];
  }
  BT22 = 1/(tmp1+1/Bmu);
  bT2 = BT22*(tmp2 + bmu/Bmu);
//  REprintf("old: %f, new: mean %f and sd %f\n\n", mu, bT2, sqrt(BT22));
  mu = as<double>(Rcpp::rnorm(1, bT2, sqrt(BT22)));

 } else {  // Gibbs-sample mu and sigma jointly (regression) 
  BT11 = 1/Bmu;
  BT12 = 0;
  BT22 = 1/Bsigma;
  bT1 = 0;
  bT2 = bmu/Bmu;

  for (int j = 0; j < T; j++) {
   tmp1 = mix_varinv[r[j]];
   tmp2 = (data[j]-mix_mean[r[j]])*tmp1;
   tmp3 = h[j]*tmp1;
   BT11 += tmp1;
   BT12 -= tmp3;
   BT22 += h[j]*tmp3;
   bT1 += h[j]*tmp2;
   bT2 += tmp2;
  }

  tmp = BT11*BT22-BT12*BT12;
  BT11 /= tmp;
  BT12 /= tmp;
  BT22 /= tmp;
 
  tmp = bT1;
  bT1 = BT11*tmp + BT12*bT2;
  bT2 = BT12*tmp + BT22*bT2;
 
  chol11 = sqrt(BT11);
  chol12 = (BT12/chol11);
  chol22 = sqrt(BT22-chol12*chol12);
 
  innov = Rcpp::rnorm(2);
  sigma = bT1 + chol11*innov[0];
  mu = bT2 + chol12*innov[0] + chol22*innov[1];
 }


 // Sampling phi: find posterior mean muT and posterior variance sigma2T
 
 sumtmp1 = h0*h[0];
 sumtmp2 = h0*h0;
 for (int j = 0; j < T-1; j++) {
  sumtmp1 += h[j]*h[j+1];
  sumtmp2 += h[j]*h[j];
 }
 tmpmean = sumtmp1/sumtmp2;
 tmpsd = 1/sqrt(sumtmp2);

 // actual sampling
 if (truncnormal) {  // draw from truncated normal via inversion method
  quant = Rcpp::pnorm(Rcpp::NumericVector::create(-1,1), tmpmean, tmpsd);
  phi_prop = (Rcpp::qnorm((quant[0] + Rcpp::runif(1)*(quant[1]-quant[0])),
                          tmpmean, tmpsd))[0];
 }
 else {  // draw from normal and reject (return) if outside
  phi_prop = (Rcpp::rnorm(1, tmpmean, tmpsd))[0]; 
  if ((phi_prop >= 1) || (phi_prop <= -1)) { // outside the unit sphere
   Rcpp::NumericVector ret = Rcpp::NumericVector::create(mu, phi, sigma);
   return ret;
  }
 }
 
 // now for the MH step, acceptance prob expR
 expR  = exp(logdnorm(h0, 0, 1/sqrt(1-phi_prop*phi_prop))
         - logdnorm(h0, 0, 1/sqrt(1-phi*phi)));
 expR *= propBeta((phi_prop+1)/2, (phi+1)/2, a0, b0);
 // ^^note that factor 1/2 from transformation of densities cancels

 // accept/reject
 if (Rcpp::as<double>(Rcpp::runif(1)) < expR) phi = phi_prop;

 Rcpp::NumericVector ret = Rcpp::NumericVector::create(mu, phi, fabs(sigma));
 return ret;
}
