#ifndef STAN_SAMPLER_HPP
#define STAN_SAMPLER_HPP

#include <ext/Rinternals.h> // SEXP

#include <fstream> // fstream, ostream
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
#  ifdef __clang__
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Wunknown-pragmas"
#    pragma clang diagnostic ignored "-Wunused-variable"
#    pragma clang diagnostic ignored "-Wunused-parameter"
#    pragma clang diagnostic ignored "-Wunused-local-typedef"
#    pragma clang diagnostic ignored "-Wunused-function"
#    pragma clang diagnostic ignored "-Wsign-compare"
#    pragma clang diagnostic ignored "-Wlanguage-extension-token"
#    pragma clang diagnostic ignored "-Winfinite-recursion"
#    pragma clang diagnostic ignored "-Wignored-qualifiers"
#    pragma clang diagnostic ignored "-Wunneeded-internal-declaration"
#    pragma clang diagnostic ignored "-Wdeprecated-copy"
#  else
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wunknown-pragmas"
#    pragma GCC diagnostic ignored "-Wunused-variable"
#    pragma GCC diagnostic ignored "-Wunused-parameter"
#    pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#    pragma GCC diagnostic ignored "-Wunused-function"
#    pragma GCC diagnostic ignored "-Wsign-compare"
#    pragma GCC diagnostic ignored "-Wignored-qualifiers"
#    pragma GCC diagnostic ignored "-Wignored-attributes"
#    pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#    pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
// This is for gcc under -Wpedantic, since some some warnings can't
// be silenced. It is highly undesirably, as it can suppress other,
// useful warnings with later code.
#    include <boost/math/tools/config.hpp>
#    ifdef BOOST_MATH_USE_FLOAT128
#      pragma GCC system_header
#    endif
#  endif
#endif

#include <stan/callbacks/interrupt.hpp>
#include <stan/callbacks/stream_logger.hpp>
#include <stan/callbacks/stream_writer.hpp>

#if defined(_WIN32)
#  define BOOST_MATH_DISABLE_DEPRECATED_03_WARNING 1
#endif
#include <stan/io/dump.hpp>
#include <stan/io/empty_var_context.hpp>
#include <stan/io/var_context.hpp>

#ifdef SUPPRESS_DIAGNOSTIC
#  ifdef __clang__
#    pragma clang diagnostic pop
#  else
#    pragma GCC diagnostic pop
#  endif
#endif

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
