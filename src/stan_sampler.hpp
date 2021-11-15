#ifndef STAN_SAMPLER_HPP
#define STAN_SAMPLER_HPP

#include <ext/Rinternals.h> // SEXP

#include <fstream> // fstream, ostream
#include <memory>  // unique_ptr
#include <string>  // string, c_str
#include <sstream> // stringstream
#include <vector>

// Uses pragmas to silence warnings generated by included headers. Has a special
// case for gcc where it will be treated as a system header, which is why it is
// spun off into its own file.
#include "stan_sampler_includes.hpp"

#include "double_writer.hpp"
#include "interruptable_sampler.hpp"

// #include "stan_files/continuous.hpp"

namespace continuous_model_namespace {
  class continuous_model;
}

namespace stan4bart {
  
  struct StanControl {
    unsigned int random_seed;
    double init_radius;
    int skip;
    double adapt_gamma;
    double adapt_delta;
    double adapt_kappa;
    unsigned int adapt_init_buffer;
    unsigned int adapt_term_buffer;
    unsigned int adapt_window;
    double adapt_t0;
    double stepsize;
    double stepsize_jitter;
    int max_treedepth;
  };
  
  struct R_CheckUserInterrupt_Functor : public stan::callbacks::interrupt {
    void operator()() {
      R_CheckUserInterrupt();
    }
  };
  
  typedef continuous_model_namespace::continuous_model StanModel;
  
  struct StanSampler {
    std::ostream& c_out;
    std::ostream& c_err;
    stan::callbacks::stream_logger logger;
    R_CheckUserInterrupt_Functor interrupt;
    
    std::fstream sample_stream;
    std::fstream diagnostic_stream;
    std::stringstream comment_stream;
    
    stan::callbacks::stream_writer diagnostic_writer;
    
    std::unique_ptr<stan::io::var_context> init_context_ptr;
    
    std::vector<std::string> constrained_param_names;
    
    
    std::vector<std::string> sample_names;
    std::vector<std::string> sampler_names;
    
    size_t sample_writer_offset; 
    
    stan4bart::double_writer init_writer;
    stan4bart::double_writer sample_writer;
  
    stan::io::dump dmp;
    stan::io::var_context& unit_e_metric;
    
    int num_pars;
    
    stan4bart::interruptable_sampler<StanModel>* sampler;
    
    StanSampler(StanModel& stanModel, const StanControl& stanControl, int chain_id, int num_warmup, int verbose);
    void run(bool isWarmup);
    
    void getParametricMean(const StanModel& model, double* result) const;
    void getParametricMean(const StanModel& model, double* result,
                           bool includeFixed, bool includeRandom) const;
    double getSigma(const StanModel& model) const;
    
    // bool isDivergentTransition() const;
    // int getTreeDepth() const;
  };
  
  struct double_writer;
   
  StanModel* createStanModelFromExpression(SEXP dataExpr);
  void deleteStanModel(StanModel* stanModel);
  void initializeStanControlFromExpression(StanControl& control, SEXP controlExpr);
  interruptable_sampler<StanModel>* createStanSampler();
  
  void setStanOffset(StanModel& model, const double* offset);
  void setResponse(StanModel& model, const double* y);
  
  
  bool isDivergentTransition(const StanSampler& sampler);
  
  SEXP createStanResultsExpr(const double_writer& sample_writer);
  void printStanControl(const StanControl& control);
}

#endif // STAN_SAMPLER_HPP
