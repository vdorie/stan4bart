#ifndef INTERRUPTABLE_SAMPLER_HPP
#define INTERRUPTABLE_SAMPLER_HPP

#include <cmath>  // log
#include <ctime> // clock_t, CLOCKS_PER_SEC
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
#include <stan/callbacks/logger.hpp>
#include <stan/callbacks/stream_writer.hpp>
#include <stan/callbacks/writer.hpp>

#include <stan/mcmc/hmc/nuts/adapt_diag_e_nuts.hpp>
#include <stan/mcmc/sample.hpp>

#include <stan/io/var_context.hpp>

#include <stan/services/util/create_rng.hpp>
#include <stan/services/util/initialize.hpp>
#include <stan/services/util/generate_transitions.hpp>
#include <stan/services/util/mcmc_writer.hpp>
#include <stan/services/util/read_diag_inv_metric.hpp>
#include <stan/services/util/validate_diag_inv_metric.hpp>
#ifdef SUPPRESS_DIAGNOSTIC
#  pragma GCC diagnostic pop
#endif

#include <ext/io.h>

namespace stan4bart {



template <class Model>
struct interruptable_sampler {
  Model& model;
  int num_thin;
  stan::callbacks::interrupt& interrupt;
  stan::callbacks::logger& logger;
  stan::callbacks::writer& sample_writer;
  double warm_delta_t;
  double sample_delta_t;
  
  boost::ecuyer1988 rng;
  std::vector<int> disc_vector;
  std::vector<double> cont_vector;
  Eigen::Map<Eigen::VectorXd> cont_params;
  
  Eigen::VectorXd inv_metric;
  stan::mcmc::adapt_diag_e_nuts<Model, boost::ecuyer1988> sampler;
  stan::services::util::mcmc_writer writer;
  stan::mcmc::sample s;
  
  // throws at least std::domain_error& e, best to catch std::exception
  // TODO: lookup exceptions thrown by:
  //   sampler.z().q = cont_params;
  //   sampler.init_stepsize(logger);
  interruptable_sampler(Model& model,
                        stan::io::var_context& init,
                        stan::io::var_context& init_inv_metric,
                        unsigned int random_seed,
                        unsigned int chain,
                        double init_radius,
                        int num_warmup,
                        int num_thin,
                        double stepsize,
                        double stepsize_jitter,
                        int max_depth,
                        double delta,
                        double gamma,
                        double kappa,
                        double t0, 
                        unsigned int init_buffer,
                        unsigned int term_buffer,
                        unsigned int window,
                        stan::callbacks::interrupt& interrupt,
                        stan::callbacks::logger& logger,
                        stan::callbacks::writer& init_writer,
                        stan::callbacks::writer& sample_writer,
                        stan::callbacks::writer& diagnostic_writer) :
    model(model),
    num_thin(num_thin),
    interrupt(interrupt),
    logger(logger),
    sample_writer(sample_writer),
    warm_delta_t(0),
    sample_delta_t(0),
    rng(stan::services::util::create_rng(random_seed, chain)),
    disc_vector(),
    cont_vector(stan::services::util::initialize(model, init, rng, init_radius, true, logger, init_writer)),
    cont_params(cont_vector.data(), cont_vector.size()),
    inv_metric(stan::services::util::read_diag_inv_metric(init_inv_metric, model.num_params_r(), logger)),
    sampler(model, rng),
    writer(sample_writer, diagnostic_writer, logger),
    s(cont_params, 0, 0)
    
  {
    stan::services::util::validate_diag_inv_metric(inv_metric, logger);
    
    sampler.set_metric(inv_metric);
    sampler.set_nominal_stepsize(stepsize);
    sampler.set_stepsize_jitter(stepsize_jitter);
    sampler.set_max_depth(max_depth);
    
    sampler.get_stepsize_adaptation().set_mu(std::log(10 * stepsize));
    sampler.get_stepsize_adaptation().set_delta(delta);
    sampler.get_stepsize_adaptation().set_gamma(gamma);
    sampler.get_stepsize_adaptation().set_kappa(kappa);
    sampler.get_stepsize_adaptation().set_t0(t0);

    sampler.set_window_params(num_warmup * num_thin, init_buffer, term_buffer, window, logger);
   
    sampler.engage_adaptation();
    
    sampler.z().q = cont_params;
    sampler.init_stepsize(logger);
    
    writer.write_sample_names(s, sampler, model);
    writer.write_diagnostic_names(s, sampler, model);
  }
  
  // implement a derived adaptation to access protected member if you need to get
  // the current adaptation window
  // reinterpret_cast<var_adaptation_derived&>(sampler.get_var_adaptation()).adapt_next_window_
  void run(bool warmup) {
    // handle thin internally - generate_transitions seems to always keep
    // first iteration, and instead we want to keep the last
    clock_t start = clock();
    if (num_thin > 1) {
      stan::services::util::generate_transitions(sampler, num_thin - 1, 0, 0,
                                                 1, 0, false, warmup, writer, s,
                                                 model, rng, interrupt, logger);
    }
    
    stan::services::util::generate_transitions(sampler, 1, 0, 0,
                                               1, 0, true, warmup, writer, s,
                                               model, rng, interrupt, logger);
    clock_t end = clock();
    if (warmup)
      warm_delta_t   += static_cast<double>(end - start) / CLOCKS_PER_SEC;
    else
      sample_delta_t += static_cast<double>(end - start) / CLOCKS_PER_SEC;
  }
  
  
  void disengage_adaptation() {
    sampler.disengage_adaptation();
    writer.write_adapt_finish(sampler);
    sampler.write_sampler_state(sample_writer);
  }
};

}

#endif // INTERRUPTABLE_SAMPLER_HPP

