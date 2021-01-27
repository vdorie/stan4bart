#ifndef STAN_SAMPLER_HPP
#define STAN_SAMPLER_HPP

#include <ext/Rinternals.h> // SEXP

#include <fstream> // fstream
#include <memory>  // unique_ptr
#include <string>  // string, c_str
#include <sstream> // stringstream
#include <vector>

#if defined(__GNUC__) && (\
  (!defined(__clang__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6))) || \
  ( defined(__clang__) && (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 7))))
#  define SUPPRESS_DIAGNOSTIC 1
#endif

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS 1
#ifdef SUPPRESS_DIAGNOSTIC
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunknown-pragmas"
#  pragma GCC diagnostic ignored "-Wunused-variable"
#  pragma GCC diagnostic ignored "-Wunused-parameter"
#  pragma GCC diagnostic ignored "-Wunused-local-typedef"
#  pragma GCC diagnostic ignored "-Wunneeded-internal-declaration"
#  pragma GCC diagnostic ignored "-Wunused-function"
#  pragma GCC diagnostic ignored "-Wsign-compare"
#  pragma GCC diagnostic ignored "-Wlanguage-extension-token"
#  pragma GCC diagnostic ignored "-Winfinite-recursion"
#  pragma GCC diagnostic ignored "-Wignored-qualifiers"
#endif
#include <stan/callbacks/interrupt.hpp>
#include <stan/callbacks/stream_logger.hpp>
#include <stan/callbacks/stream_writer.hpp>

#include <stan/io/dump.hpp>
#include <stan/io/empty_var_context.hpp>
#include <stan/io/var_context.hpp>
#ifdef SUPPRESS_DIAGNOSTIC
#  pragma GCC diagnostic pop
#endif

#include "double_writer.hpp"
#include "interruptable_sampler.hpp"

// #include "stan_files/continuous.hpp"

namespace model_continuous_namespace {
  class model_continuous;
}

namespace stan4bart {
  
  struct StanControl {
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
    
    StanSampler(StanModel& stanModel, const StanControl& stanControl, int chain_id, int num_warmup);
    void run(bool isWarmup);
  };
  
  struct double_writer;
   
  StanModel* createStanModelFromExpression(SEXP dataExpr);
  void deleteStanModel(StanModel* stanModel);
  void initializeStanControlFromExpression(StanControl& control, SEXP controlExpr);
  interruptable_sampler<StanModel>* createStanSampler();
  
  void setStanOffset(StanModel& model, const double* offset);
  void setResponse(StanModel& model, const double* y);
  void getParametricMean(const StanSampler& sampler, const StanModel& model, double* result);
  void getParametricMean(const StanSampler& sampler, const StanModel& model, double* result,
                         bool includeFixed, bool includeRandom);
  double getSigma(const StanSampler& sampler, const StanModel& model);
  
  SEXP createStanResultsExpr(const double_writer& sample_writer);
  void printStanControl(const StanControl& control);
}

#endif // STAN_SAMPLER_HPP
