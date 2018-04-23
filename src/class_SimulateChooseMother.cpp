#include "mitolina_types.h"

#include <RcppArmadillo.h> // FIXME: Avoid Rcpp here? Only in api_* files?


/*****************************************
WFRandomMother
******************************************/
WFRandomMother::WFRandomMother(size_t population_size) {
  m_population_size = (double)population_size;
}

void WFRandomMother::update_state_new_generation() {  
  //Rcpp::Rcout << "WFRandomMother: update_state_new_generation NOOP" << std::endl;
}

int WFRandomMother::get_mother_i() {
  //Rcpp::Rcout << "WFRandomMohter: get_mother_i" << std::endl;
  return R::runif(0, 1)*m_population_size;
}


/*****************************************
GammaVarianceRandomMother
******************************************/
GammaVarianceRandomMother::GammaVarianceRandomMother(size_t population_size, double gamma_parameter_shape, double gamma_parameter_scale) {
  m_population_size = population_size;
  m_gamma_parameter_shape = gamma_parameter_shape;
  m_gamma_parameter_scale = gamma_parameter_scale;
}

// modified from 
// https://github.com/RcppCore/RcppArmadillo/blob/master/inst/include/RcppArmadilloExtensions/sample.h
// ProbSampleReplace for size = 1
void GammaVarianceRandomMother::update_state_new_generation() {    
  //Rcpp::Rcout << "GammaVarianceRandomMother: update_state_new_generation" << std::endl;
  
  Rcpp::NumericVector mothers_prob_tmpl = Rcpp::rgamma(m_population_size, m_gamma_parameter_shape, m_gamma_parameter_scale);
  //Rcpp::Rcout << "mean[mothers_prob_tmpl] = " << Rcpp::mean(mothers_prob_tmpl) << ", var[mothers_prob_tmpl] = " << Rcpp::var(mothers_prob_tmpl) << std::endl;
  mothers_prob_tmpl = mothers_prob_tmpl / sum(mothers_prob_tmpl);    

  arma::vec mothers_prob(mothers_prob_tmpl.begin(), mothers_prob_tmpl.size(), false); // false means no copy
  arma::uvec mothers_prob_perm = arma::sort_index(mothers_prob, "descend"); //descending sort of index
  mothers_prob = arma::sort(mothers_prob, "descend");  // descending sort of prob
  mothers_prob = arma::cumsum(mothers_prob);
  
  m_mothers_prob_cum = mothers_prob;
  m_mothers_prob_perm = mothers_prob_perm;
}

int GammaVarianceRandomMother::get_mother_i() {
  //Rcpp::Rcout << "GammaVarianceRandomMother: get_mother_i" << std::endl;
  
  double rU = R::unif_rand(); // R's internal random number generation
  //double rU = R::runif(0, 1);

  int jj;
  size_t population_size_1 = m_population_size - 1;
  
  for (jj = 0; jj < population_size_1; ++jj) {
    if (rU <= m_mothers_prob_cum[jj]) {
      break;
    }
  }
  
  return m_mothers_prob_perm[jj];
}


