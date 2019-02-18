#include <RcppArmadillo.h> // FIXME: Avoid Rcpp here? Only in api_* files?

class SimulateChooseMother {
  public:
    virtual void update_state_new_generation() = 0;
    virtual int get_mother_i() = 0;
};

class WFRandomMother: public SimulateChooseMother {
  private:
    double m_population_size;
    
  public:
    WFRandomMother(size_t population_size);
    void update_state_new_generation();
    int get_mother_i();
};


class GammaVarianceRandomMother: public SimulateChooseMother {
  private:
    // fixed for entire simulate
    size_t m_population_size;
    double m_gamma_parameter_shape;
    double m_gamma_parameter_scale;
    
    // new for each generation
    arma::vec m_mothers_prob_cum;
    arma::uvec m_mothers_prob_perm;
    
  public:
    GammaVarianceRandomMother(size_t population_size, double gamma_parameter_shape, double gamma_parameter_scale);
    void update_state_new_generation();
    int get_mother_i();
 };
