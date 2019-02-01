#include <RcppArmadillo.h>
#include "auxmix.h"
#include "mixture-state-sampler.h"
#include "parameterization.h"
#define _USE_MATH_DEFINES
#include <cmath>

arma::mat mixture_state_post_dist(
    const arma::vec& eps_star,
    const arma::vec& eta,
    const arma::ivec& d,
    const double mu,
    const double sigma2,
    const double rho,
    const Parameterization centering) {
  
  const int n = eps_star.size();
  const int mix_count = sizeof(mix_prob)/sizeof(mix_prob[0]);
  const double sigma2_used = centering == Parameterization::CENTERED ? sigma2 : 1.0;
  arma::mat result(mix_count, n);
  
  for (int r = 0; r < n; r++) {
    for (int c = 0; c < mix_count; c++) {
      const double a = mix_a[c];
      const double b = mix_b[c];
      const double m = mix_mean[c];
      const double v2 = mix_var[c];
      const double log_prior = log(mix_prob[c]);
      
      double log_eps_star_lik = -0.5 * pow((eps_star[r] - m), 2)/v2 - log(sqrt(2*M_PI*v2));
      double log_eta_lik;
      if (r < n-1) {
        log_eta_lik = -0.5/(sigma2_used*(1-rho*rho)) * pow(eta[r] - rho*sqrt(sigma2_used)*d[r]*exp(m/2)*(a+b*(eps_star[r]-m)), 2) - 0.5*log(2*M_PI*(1-rho*rho)*sigma2_used);
      } else {
        log_eta_lik = 0.0;
      }
      /*log_*/result.at(c, r) = log_prior + log_eps_star_lik + log_eta_lik;
    }
    const double max_log_result_col = arma::max(/*log_*/result.col(r));
    /*log_*/result.col(r) = /*log_*/result.col(r) - (max_log_result_col + log(arma::sum(arma::exp(/*log_*/result.col(r)-max_log_result_col))));
    result.col(r) = arma::exp(/*log_*/result.col(r));
    result.col(r) = arma::cumsum(result.col(r));
    result.col(r) = result.col(r) / result.at(mix_count-1, r);
  }
  
  return result;
}

arma::vec draw_s_auxiliary(
    const arma::vec& y_star,
    const arma::ivec& d,
    const arma::vec& h,
    const arma::vec& ht,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const Parameterization centering) {

  const int n = y_star.size();
  arma::vec eps_star;
  arma::vec eta;
  arma::mat post_dist;
  arma::vec unif_vec;
  arma::vec new_states(n);
  
  switch (centering) {
    case Parameterization::CENTERED:
    eps_star = y_star - h;
    eta = (h.tail(n-1) - mu) - phi*(h.head(n-1) - mu);
    break;
    case Parameterization::NONCENTERED:
    eps_star = y_star - mu - sqrt(sigma2)*ht;
    eta = ht.tail(n-1) - phi*ht.head(n-1);
    break;
  }
  post_dist = mixture_state_post_dist(eps_star, eta, d, mu, sigma2, rho, centering);

  for (int r = 0; r < n; r++) {
    const arma::vec post_dist_col_r = post_dist.unsafe_col(r);
    auto binary_search_result = std::lower_bound(post_dist_col_r.cbegin(), post_dist_col_r.cend(), R::runif(0, 1));
    new_states[r] = std::distance(post_dist_col_r.cbegin(), binary_search_result);
  }
  
  return new_states;
}
