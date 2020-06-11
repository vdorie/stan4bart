#ifndef STAN_SAMPLER_HPP
#define STAN_SAMPLER_HPP

#include <ext/Rinternals.h> // SEXP

#include <fstream> // fstream
#include <memory>  // unique_ptr
#include <string>  // string, c_str
#include <sstream> // stringstream
#include <vector>

#include <stan/callbacks/interrupt.hpp>
#include <stan/callbacks/stream_logger.hpp>
#include <stan/callbacks/stream_writer.hpp>

#include <stan/io/dump.hpp>
#include <stan/io/empty_var_context.hpp>
#include <stan/io/var_context.hpp>

#include "double_writer.hpp"
#include "interruptable_sampler.hpp"

// #include "stan_files/continuous.hpp"

namespace model_continuous_namespace {
  struct model_continuous;
}

namespace stan4bart {
  
  struct StanArgs {
    unsigned int random_seed;
    double init_radius;
    int thin;
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
  
  typedef model_continuous_namespace::model_continuous StanModel;
  
  struct StanSampler {
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
    
    StanSampler(StanModel& stanModel, const StanArgs& stanArgs, int chain_id, int num_warmup);
    void run(bool isWarmup);
  };
  
  struct double_writer;
   
  StanModel* createStanModelFromExpression(SEXP dataExpr);
  void deleteStanModel(StanModel* stanModel);
  void initializeStanArgsFromExpression(StanArgs& args, SEXP argsExpr);
  interruptable_sampler<StanModel>* createStanSampler();
  
  void setStanOffset(StanModel& model, const double* offset);
  void getZb(const StanSampler& sampler, const StanModel& model, double* zb);
  double getSigma(const StanSampler& sampler, const StanModel& model);
  
  SEXP createStanResultsExpr(const double_writer& sample_writer);
  void printStanArgs(const StanArgs& args);
}

#endif // STAN_SAMPLER_HPP
