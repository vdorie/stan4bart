#include <ext/Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include <exception>
#include <memory> // unique_ptr
#include <set> // external pointers set

#include <rc/util.h>
#include <rc/bounds.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/control.hpp>
#include <dbarts/data.hpp>
#include <dbarts/model.hpp>
#include <dbarts/results.hpp>

#include <rstan/io/r_ostream.hpp>

#include "bart_util.hpp"
#include "stan_sampler.hpp"

namespace {
  // stuff to handle external pointers
  typedef bool(*ExternalPointerComparator)(const SEXP&lhs, const SEXP& rhs);
  typedef std::set<SEXP, ExternalPointerComparator> PointerSet;
  
  PointerSet* activeSamplers;
  bool compareExternalPointers(const SEXP& lhs, const SEXP& rhs) {
    return R_ExternalPtrAddr(const_cast<SEXP>(lhs)) < R_ExternalPtrAddr(const_cast<SEXP>(rhs));
  }
  
  struct BARTFunctionTable {
    void (*initializeFit)(dbarts::BARTFit* fit, dbarts::Control* control, dbarts::Model* model, dbarts::Data* data);
    void (*invalidateFit)(dbarts::BARTFit* fit);
    void (*initializeControl)(dbarts::Control* control, SEXP controlExpr);
    void (*initializeData)(dbarts::Data* data, SEXP dataExpr);
    void (*invalidateData)(dbarts::Data* data);
    void (*initializeModel)(dbarts::Model* model, SEXP modelExpr, const dbarts::Control* control);
    void (*invalidateModel)(dbarts::Model* model);

    void (*runSamplerWithResults)(dbarts::BARTFit* fit, std::size_t numBurnIn, dbarts::Results* results);
    void (*setResponse)(dbarts::BARTFit* fit, const double* response);
    void (*setOffset)(dbarts::BARTFit* fit, const double* offset, bool updateState);
    void (*setSigma)(dbarts::BARTFit* fit, const double* sigma);
    void (*sampleTreesFromPrior)(dbarts::BARTFit* fit);
    void (*printInitialSummary)(const dbarts::BARTFit* fit);
  };
  BARTFunctionTable bartFunctions;
  
  double getRange(const double* x, size_t len) {
    double min, max;
    min = max = x[0];
    for (size_t i = 1; i < len; ++i) {
      if (x[i] < min) min = x[i];
      else if (x[i] > max) max = x[i];
    }
    return max - min;
  }
  
  enum UserOffsetType {
    OFFSET_DEFAULT = 0,
    OFFSET_FIXEF,
    OFFSET_RANEF
  };
  
  enum ResultsType {
    RESULTS_BOTH = 0,
    RESULTS_BART,
    RESULTS_STAN
  };
  
  struct Sampler {
    int numWarmup;
    int verbose;
    int refresh;
    const double* userOffset;
    UserOffsetType offsetType;
    
    model_continuous_namespace::model_continuous* stanModel;
    stan4bart::StanArgs stanArgs;
    stan4bart::StanSampler* stanSampler;
    
    dbarts::Control bartControl;
    dbarts::Data bartData;
    dbarts::Model bartModel;
    dbarts::BARTFit* bartSampler;
    
    double* bartOffset;
    double* stanOffset;
        
    Sampler() :
      stanModel(NULL), stanSampler(NULL), bartModel(false), bartSampler(NULL), bartOffset(NULL), stanOffset(NULL)
    {
    }
    ~Sampler() {
      delete [] stanOffset;
      delete [] bartOffset;
      
      // delete bart model
      if (bartSampler != NULL) {
        bartFunctions.invalidateFit(bartSampler);
        ::operator delete(bartSampler);
        bartSampler = NULL;
      }
      bartFunctions.invalidateModel(&bartModel);
      bartFunctions.invalidateData(&bartData);
      
      delete stanSampler;
      stan4bart::deleteStanModel(stanModel);
      stanModel = NULL;
    }
  };
  
  void initializeSamplerFromExpression(Sampler& sampler, SEXP commonControlExpr);
}

extern "C" {
  static void samplerFinalizer(SEXP samplerExpr);
  
  static SEXP createSampler(SEXP bartControlExpr, SEXP bartDataExpr, SEXP bartModelExpr,
                            SEXP stanDataExpr, SEXP stanArgsExpr,
                            SEXP commonControlExpr)
  {
    std::unique_ptr<Sampler> samplerPtr(new Sampler);
    Sampler& sampler(*samplerPtr);
    
    initializeSamplerFromExpression(sampler, commonControlExpr);
    
    const double* ranef_init = REAL(rc_getListElement(commonControlExpr, "ranef_init"));
    double sigma_init = rc_getDouble(rc_getListElement(commonControlExpr, "sigma_init"), "sigma_init",
      RC_VALUE | RC_GT, 0.0, RC_VALUE | RC_DEFAULT, 1.0,
      RC_END);
    
    
    sampler.stanModel = stan4bart::createStanModelFromExpression(stanDataExpr);
    stan4bart::initializeStanArgsFromExpression(sampler.stanArgs, stanArgsExpr);
    if (sampler.stanArgs.thin == R_NaInt) {
      sampler.stanArgs.thin = (2000 - sampler.numWarmup) / 1000;
      if (sampler.stanArgs.thin < 1) sampler.stanArgs.thin = 1;
    }
    
    int chain_id = 1;
    sampler.stanSampler = new stan4bart::StanSampler(*sampler.stanModel, sampler.stanArgs, chain_id, sampler.numWarmup);
    
    bartFunctions.initializeControl(&sampler.bartControl, bartControlExpr);
    bartFunctions.initializeData(&sampler.bartData, bartDataExpr);
    bartFunctions.initializeModel(&sampler.bartModel, bartModelExpr, &sampler.bartControl);
    
    sampler.bartSampler = static_cast<dbarts::BARTFit*>(::operator new (sizeof(dbarts::BARTFit)));
    bartFunctions.initializeFit(sampler.bartSampler, &sampler.bartControl, &sampler.bartModel, &sampler.bartData);
    
    size_t n = sampler.bartData.numObservations;
    sampler.bartOffset = new double[n];
    sampler.stanOffset = new double[n];
    
    if (ranef_init != NULL) {
      std::memcpy(sampler.bartOffset, ranef_init, n * sizeof(double));
    } else {
      for (size_t i = 0; i < n; ++i) sampler.bartOffset[i] = 0.0;
    }
    if (sampler.userOffset != NULL) {
      if (sampler.offsetType == OFFSET_FIXEF)
        for (size_t i = 0; i < n; ++i) sampler.bartOffset[i] = 0.0;
      else
        for (size_t i = 0; i < n; ++i) sampler.bartOffset[i] += sampler.userOffset[i];
    }
    
    bartFunctions.setOffset(sampler.bartSampler, sampler.bartOffset, true);
    bartFunctions.setSigma(sampler.bartSampler, &sigma_init);

    bartFunctions.sampleTreesFromPrior(sampler.bartSampler);
    
    for (size_t i = 0; i < sampler.bartData.numObservations; ++i) sampler.stanOffset[i] = 0.0;
    // sampler.stanModel->set_offset(sampler.stanOffset);
    stan4bart::setStanOffset(*sampler.stanModel, sampler.stanOffset);
    
    
    
    SEXP result = PROTECT(R_MakeExternalPtr(samplerPtr.get(), R_NilValue, R_NilValue));
    samplerPtr.release();
    R_RegisterCFinalizerEx(result, samplerFinalizer, static_cast<Rboolean>(FALSE));
    

    activeSamplers->insert(result);

    UNPROTECT(1);
    
    return result;
  }
  
  static SEXP run(SEXP samplerExpr, SEXP numIterExpr, SEXP isWarmupExpr, SEXP resultsTypeExpr)
  {
    Sampler* samplerPtr = static_cast<Sampler*>(R_ExternalPtrAddr(samplerExpr));
    if (samplerPtr == NULL) Rf_error("run called on NULL external pointer");
    Sampler& sampler(*samplerPtr);
    
    int numIter   = rc_getInt(numIterExpr, "num_iter", RC_VALUE | RC_GEQ, 1, RC_END);
    bool isWarmup = rc_getBool(isWarmupExpr, "warmup", RC_NA | RC_NO, RC_END);
    
    ResultsType resultsType = static_cast<ResultsType>(rc_getInt(resultsTypeExpr, "results_type",
     RC_VALUE | RC_DEFAULT, static_cast<int>(RESULTS_BOTH),
     RC_END));
    
    stan4bart::IterableBartResults* bartSamples = NULL;
    // allocate storage for results
    if (resultsType == RESULTS_BOTH || resultsType == RESULTS_BART)
      bartSamples = new stan4bart::IterableBartResults(
        sampler.bartData.numObservations, sampler.bartData.numPredictors, 0, numIter, 1 /* num chains */, false /* TODO: binary */);
    if (resultsType == RESULTS_BOTH || resultsType == RESULTS_STAN)
      sampler.stanSampler->sample_writer.resize(sampler.stanSampler->num_pars, numIter);
    
    size_t n = sampler.bartData.numObservations;
    
    
    if (sampler.verbose > 0)
      Rprintf("starting %s, %d draws, %s:\n", isWarmup ? "warmup" : "sampling", numIter,
              resultsType == RESULTS_BOTH ? "both BART and Stan" : (resultsType == RESULTS_BART ? "BART only" : "Stan only"));
    
    // TODO: add in capacity for fitting against true ranef and true fixef
    for (size_t iter = 0; iter < numIter; ++iter) {
      if (sampler.refresh > 0 && sampler.verbose > 1 && (iter + 1) % sampler.refresh == 0)
          Rprintf("  iter %.3d / %.3d\n", iter + 1, numIter);
      
      if (resultsType == RESULTS_BOTH || resultsType == RESULTS_BART) {
        
        bartFunctions.runSamplerWithResults(sampler.bartSampler, 0, bartSamples);
        
        // bart with an offset will produce predictions that have the offset added;
        // in order to just get the tree predictions, subtract out that offset
        for (size_t j = 0; j < sampler.bartData.numObservations; ++j) bartSamples->trainingSamples[j] -= sampler.bartOffset[j];
        
        std::memcpy(sampler.stanOffset, const_cast<const double*>(bartSamples->trainingSamples), n * sizeof(double));
        
        if (sampler.userOffset != NULL) {
          if (sampler.offsetType == OFFSET_DEFAULT)
            for (size_t j = 0; j < n; ++j) sampler.stanOffset[j] += sampler.userOffset[j];
          else if (sampler.offsetType == OFFSET_FIXEF)
            std::memcpy(sampler.stanOffset, sampler.userOffset, n * sizeof(double));
        }
        stan4bart::setStanOffset(*sampler.stanModel, sampler.stanOffset);
        
        bartSamples->incrementPointers();
      }
            
      if (resultsType == RESULTS_BOTH || resultsType == RESULTS_STAN) {
        sampler.stanSampler->run(isWarmup);
        
        stan4bart::getZb(*sampler.stanSampler, *sampler.stanModel, sampler.bartOffset);
        if (sampler.userOffset != NULL) {
          if (sampler.offsetType == OFFSET_DEFAULT)
            for (size_t j = 0; j < n; ++j) sampler.bartOffset[j] += sampler.userOffset[j];
          else if (sampler.offsetType == OFFSET_RANEF)
            std::memcpy(sampler.bartOffset, sampler.userOffset, n * sizeof(double));
        }
        double sigma = getSigma(*sampler.stanSampler, *sampler.stanModel);
        sampler.stanSampler->sample_writer.increment();
        
        bartFunctions.setOffset(sampler.bartSampler, sampler.bartOffset, isWarmup);
        bartFunctions.setSigma(sampler.bartSampler, &sigma);
      }
    }
    if (bartSamples != NULL) bartSamples->resetPointers();
        
    SEXP resultExpr = PROTECT(rc_newList(resultsType == RESULTS_BOTH ? 2 : 1));
    int pos = 0;
    if (resultsType == RESULTS_BOTH || resultsType == RESULTS_STAN)
      SET_VECTOR_ELT(resultExpr, pos++, createStanResultsExpr(sampler.stanSampler->sample_writer));
    if (resultsType == RESULTS_BOTH || resultsType == RESULTS_BART)
      SET_VECTOR_ELT(resultExpr, pos, stan4bart::createBartResultsExpr(*sampler.bartSampler, *bartSamples));
    
    SEXP namesExpr = PROTECT(rc_newCharacter(XLENGTH(resultExpr)));
    pos = 0;
    if (resultsType == RESULTS_BOTH || resultsType == RESULTS_STAN)
      SET_STRING_ELT(namesExpr, pos++, Rf_mkChar("stan"));
    if (resultsType == RESULTS_BOTH || resultsType == RESULTS_BART)
      SET_STRING_ELT(namesExpr, pos, Rf_mkChar("bart"));
    
    rc_setNames(resultExpr, namesExpr);
    UNPROTECT(1);
    
    delete bartSamples;
    
    UNPROTECT(1);
    
    return(resultExpr);
  }
    
        /* bartFunctions.runSamplerWithResults(fixefModel, 0, fixefWarmup);
        
        // bart with an offset will produce predictions that have the offset added;
        // in order to just get the tree predictions, subtract out that offset
        for (size_t j = 0; j < bartData.numObservations; ++j) fixefWarmup->trainingSamples[j] -= fixef_offset[j];
        
        std::memcpy(ranef_offset, const_cast<const double*>(fixefWarmup->trainingSamples), bartData.numObservations * sizeof(double));
        
        if (userOffset != NULL) {
          if (offset_type == OFFSET_DEFAULT)
            for (size_t j = 0; j < n; ++j) sampler.stanOffset[j] += sampler.userOffset[j];
          else if (offset_type == OFFSET_FIXEF)
            std::memcpy(ranef_offset, userOffset, bartData.numObservations * sizeof(double));
        }
        
        ranefModel->set_offset(sampler.stanOffset);
        
        if (save_warmup) {
          fixefWarmup->incrementPointers();
          if ((fixefWarmup->trainingSamples - fixefWarmup->base_trainingSamples) / bartData.numObservations != (iter + 1))
            Rf_error("samples beyond extent");
        }
      
        ranefSampler->run(true);
        
        ranefModel->get_Zb(sample_writer.x_curr + sample_writer_offset, fixef_offset);
        if (userOffset != NULL) {
          if (offset_type == OFFSET_DEFAULT)
            for (size_t j = 0; j < bartData.numObservations; ++j) fixef_offset[j] += userOffset[j];
          else if (offset_type == OFFSET_RANEF)
            std::memcpy(fixef_offset, userOffset, bartData.numObservations * sizeof(double));
        }
        sigma = ranefModel->get_aux(sample_writer.x_curr + sample_writer_offset);
        sample_writer.increment();
        
        std::memcpy(fixef_response, y_orig, bartData.numObservations * sizeof(double));
        for (size_t j = 0; j < bartData.numObservations; ++j) {
          fixef_response[j] -= fixef_offset[j];
          fixef_offset[j] = 0.0;
        }
        
        bartFunctions.setOffset(fixefModel, fixef_offset);
        bartFunctions.setResponse(fixefModel, fixef_response);
        y_scale = getRange(fixef_response, bartData.numObservations);
        bartModel.sigmaSqPrior->setScale((sigma / y_scale) * (sigma / y_scale));
        fixefModel->state[0].sigma = sigma / y_scale;
      }
      fixefWarmup->resetPointers();
    }
    
    ranefSampler->disengage_adaptation();
    

    delete bartSamples;
  } */
  
  static SEXP printInitialSummary(SEXP samplerExpr) {
    Sampler* samplerPtr = static_cast<Sampler*>(R_ExternalPtrAddr(samplerExpr));
    if (samplerPtr == NULL) Rf_error("printInitialSummary called on NULL external pointer");
    Sampler& sampler(*samplerPtr);
    
    Rprintf("stan args:\n");
    printStanArgs(sampler.stanArgs);
    Rprintf("bart init:\n");
    bartFunctions.printInitialSummary(sampler.bartSampler);
    
    if (sampler.userOffset != NULL) {
      Rprintf("\nuser offset: %f", sampler.userOffset[0]);
      for (size_t i = 1; i < (sampler.bartData.numObservations < 5 ? sampler.bartData.numObservations : 5); ++i)
        Rprintf(", %f", sampler.userOffset[i]);
      if (sampler.bartData.numObservations > 5) Rprintf("...");
      Rprintf("\n");
      if (sampler.offsetType != OFFSET_DEFAULT) Rprintf("  type: %s\n", sampler.offsetType == OFFSET_RANEF ? "ranef" : "fixef");
    }
    
    return R_NilValue;
  }
  
  static SEXP disengageAdaptation(SEXP samplerExpr)
  {
    Sampler* samplerPtr = static_cast<Sampler*>(R_ExternalPtrAddr(samplerExpr));
    if (samplerPtr == NULL) Rf_error("disengageAdaptation called on NULL external pointer");
    Sampler& sampler(*samplerPtr);
    
    sampler.stanSampler->sampler->disengage_adaptation();
    
    return R_NilValue;
  }

} // extern "C"

namespace {

void initializeSamplerFromExpression(Sampler& sampler, SEXP commonControlExpr)
{
  sampler.numWarmup = rc_getInt(rc_getListElement(commonControlExpr, "warmup"), "warmup",
    RC_VALUE | RC_DEFAULT, 1000,
    RC_END);
  sampler.verbose = rc_getInt(rc_getListElement(commonControlExpr, "verbose"), "verbose",
    RC_VALUE | RC_DEFAULT, 0,
    RC_END);
  sampler.refresh = rc_getInt(rc_getListElement(commonControlExpr, "refresh"), "refresh",
    RC_VALUE | RC_GEQ, 0,
    RC_NA | RC_YES, RC_END);
  
  SEXP offsetExpr = rc_getListElement(commonControlExpr, "offset");
  sampler.userOffset = offsetExpr == R_NilValue || rc_getLength(offsetExpr) == 0 || !Rf_isReal(offsetExpr) ? NULL : REAL(offsetExpr);
  
  sampler.offsetType = static_cast<UserOffsetType>(
    rc_getInt(rc_getListElement(commonControlExpr, "offset_type"), "offset_type",
     RC_VALUE | RC_DEFAULT, static_cast<int>(OFFSET_DEFAULT),
     RC_END));
  
  if (sampler.refresh == R_NaInt)
    sampler.refresh = 200;
}

}

extern "C" {
  
  
  // TODO: test observations
  /* 
  static SEXP mstan4bart(SEXP bartControlExpr, SEXP bartDataExpr, SEXP bartModelExpr,
                         SEXP stanDataExpr, SEXP stanArgsExpr,
                         SEXP commonControlExpr)
  {
    
    
    model_continuous_namespace::model_continuous* ranefModel = createStanModelFromExpression(stanDataExpr);
    
    StanArgs stanArgs;
    initializeStanArgsFromExpression(stanArgs, stanArgsExpr);
    
        
        
    
    std::ostream nullout(nullptr);
    // std::ostream& c_out = refresh > 0 ? rstan::io::rcout : nullout;
    // std::ostream& c_err = refresh > 0 ? rstan::io::rcerr : nullout;
    std::ostream& c_out = nullout;
    std::ostream& c_err = nullout;
    
    stan::callbacks::stream_logger 
      logger(c_out, c_out, c_out, c_err, c_err);
  
    R_CheckUserInterrupt_Functor interrupt;
  
    std::fstream sample_stream;
    std::fstream diagnostic_stream;
    std::stringstream comment_stream;
    bool append_samples(false);
    
    
    stan::callbacks::stream_writer diagnostic_writer(diagnostic_stream, "# ");
    
    std::unique_ptr<stan::io::var_context> init_context_ptr;
    init_context_ptr.reset(new stan::io::empty_var_context());
    
    std::vector<std::string> constrained_param_names;
    ranefModel->constrained_param_names(constrained_param_names);
    // init_writer allocation moved below
    
    int return_code = stan::services::error_codes::CONFIG;
    
    // start sampling block
    std::vector<std::string> sample_names;
    stan::mcmc::sample::get_sample_param_names(sample_names);
    std::vector<std::string> sampler_names;
    
    // start NUTS block
    sampler_names.resize(5);
    sampler_names[0] = "stepsize__";
    sampler_names[1] = "treedepth__";
    sampler_names[2] = "n_leapfrog__";
    sampler_names[3] = "divergent__";
    sampler_names[4] = "energy__";
    size_t sample_writer_offset = sample_names.size() + sampler_names.size();
   
    stan4bart::double_writer init_writer("init", ranefModel->num_params_r(), 1);
    stan4bart::double_writer sample_writer("sample", sample_names.size() + sampler_names.size() + constrained_param_names.size(), numWarmup + num_iter);
    
    stan::io::dump dmp
      = stan::services::util::create_unit_e_diag_inv_metric(ranefModel->num_params_r());
    stan::io::var_context& unit_e_metric(dmp);
    
    
    stan4bart::interruptable_sampler<model_continuous_namespace::model_continuous>* ranefSampler = NULL;
    
    try {
      ranefSampler = new stan4bart::interruptable_sampler<model_continuous_namespace::model_continuous>(
        *ranefModel,
        *init_context_ptr,
        unit_e_metric,
        stan_args.random_seed,
        chain_id,
        stan_args.init_radius,
        numWarmup,
        stan_args.thin,
        0, // stan_args.refresh,
        stan_args.stepsize,
        stan_args.stepsize_jitter,
        stan_args.max_treedepth,
        stan_args.adapt_delta,
        stan_args.adapt_gamma,
        stan_args.adapt_kappa,
        stan_args.adapt_t0,
        stan_args.adapt_init_buffer,
        stan_args.adapt_term_buffer,
        stan_args.adapt_window,
        interrupt, 
        logger,
        init_writer,
        sample_writer,
        diagnostic_writer);
    } catch (std::exception& e) {
      // TODO: delete stuff, maybe goto error
      std::stringstream errorMessage;
      errorMessage << "Error allocating interruptable sample: " << e.what();
      Rf_error(errorMessage.str().c_str());
    } */
    
    
    // create bart model
    /* dbarts::Control bartControl;
    dbarts::Data bartData;
    dbarts::Model bartModel(false);
    
    bartFunctions.initializeControl(&bartControl, bartControlExpr);
    bartFunctions.initializeData(&bartData, bartDataExpr);
    bartFunctions.initializeModel(&bartModel, bartModelExpr, &bartControl);
    
    dbarts::BARTFit* fixefModel = static_cast<dbarts::BARTFit*>(::operator new (sizeof(dbarts::BARTFit)));
    bartFunctions.initializeFit(fixefModel, &bartControl, &bartModel, &bartData);
    
    stan4bart::IterableBartResults* fixefSamples = new stan4bart::IterableBartResults( */
     // bartData.numObservations, bartData.numPredictors, 0, num_iter, 1 /* num chains */, false /* TODO: binary */);
    /* stan4bart::IterableBartResults* fixefWarmup = new stan4bart::IterableBartResults(
      bartData.numObservations, bartData.numPredictors, 0, save_warmup ? numWarmup : 1, 1, false);
    //Rprintf("allocating bart results with %lu obs, %lu preds, %lu fixef samples, %lu fixef warmup\n",
    //        bartData.numObservations, bartData.numPredictors, num_iter, save_warmup ? numWarmup : 1);
    
    bartFunctions.sampleTreesFromPrior(fixefModel);
    
    // y.scale an unfortunate side effect of how sigma priors are done in dbarts
    double y_scale, y_min, y_max;
    y_min = y_max = bartData.y[0];
    for (size_t i = 1; i < bartData.numObservations; ++i) {
      if (bartData.y[i] > y_max) y_max = bartData.y[i];
      if (bartData.y[i] < y_min) y_min = bartData.y[i];
    }
    y_scale = y_max - y_min; */
    /* 
    if (verbose > 0) {
      Rprintf("stan args:\n");
      printStanArgs(stan_args);
      Rprintf("bart init:\n");
      bartFunctions.printInitialSummary(fixefModel);
      Rprintf("  thin: %u\n", bartControl.treeThinningRate);
      
      if (userOffset != NULL) {
        Rprintf("\nuser offset: %f", userOffset[0]);
        for (size_t i = 1; i < (bartData.numObservations < 5 ? bartData.numObservations : 5); ++i)
          Rprintf(", %f", userOffset[i]);
        Rprintf("\n");
        if (offset_type != OFFSET_DEFAULT) Rprintf("  type: %s\n", offset_type == OFFSET_RANEF ? "ranef" : "fixef");
      }
    }
      
    
    double* fixef_offset = new double[bartData.numObservations];
    double* fixef_response = new double[bartData.numObservations];
    const double* y_orig = bartData.y;
    double* ranef_offset = new double[bartData.numObservations];
    
    std::memcpy(fixef_offset, ranef_init, bartData.numObservations * sizeof(double));
    if (ranef_true != NULL)
      std::memcpy(fixef_offset, ranef_true, bartData.numObservations * sizeof(double));
    
    if (userOffset != NULL) {
      if (offset_type == OFFSET_FIXEF)
        for (size_t i = 0; i < bartData.numObservations; ++i) fixef_offset[i] = 0.0;
      else
        for (size_t i = 0; i < bartData.numObservations; ++i) fixef_offset[i] += userOffset[i];
    }
    
    bartFunctions.setOffset(fixefModel, fixef_offset);
    bartModel.sigmaSqPrior->setScale((sigma_init / y_scale) * (sigma_init / y_scale));
    
    if (numWarmup > 0) {
      if (verbose > 0) Rprintf("starting warmup, %d iters\n", numWarmup);
      
      unsigned int iter = 0;
      unsigned int warmup_end_1 = numWarmup > 1 ? numWarmup / 2 - 1 : 0;
      
      // warmup just bart first
      for ( ; iter < warmup_end_1; ++iter) {
        bartFunctions.runSamplerWithResults(fixefModel, 0, fixefWarmup);
        
        // bart with an offset will produce predictions that have the offset added;
        // in order to just get the tree predictions, subtract out that offset
        for (size_t j = 0; j < bartData.numObservations; ++j) fixefWarmup->trainingSamples[j] -= fixef_offset[j];
        
        if (save_warmup) fixefWarmup->incrementPointers();
      }
      
      bartFunctions.runSamplerWithResults(fixefModel, 0, fixefWarmup);
      for (size_t j = 0; j < bartData.numObservations; ++j) fixefWarmup->trainingSamples[j] -= fixef_offset[j];
      
      std::memcpy(ranef_offset, const_cast<const double*>(fixefWarmup->trainingSamples), bartData.numObservations * sizeof(double));
      if (fixef_true != NULL)
        std::memcpy(ranef_offset, fixef_true, bartData.numObservations * sizeof(double));
      
      if (userOffset != NULL) {
        if (offset_type == OFFSET_DEFAULT)
          for (size_t j = 0; j < bartData.numObservations; ++j) ranef_offset[j] += userOffset[j];
        else if (offset_type == OFFSET_FIXEF)
          std::memcpy(ranef_offset, userOffset, bartData.numObservations * sizeof(double));
      }
      
      ranefModel->set_offset(ranef_offset);
      
      if (save_warmup) fixefWarmup->incrementPointers();
      
      // now just stan
      for (iter = 0 ; iter < warmup_end_1; ++iter) {
        if (refresh > 0 && verbose > 1 && (iter + 1) % refresh == 0)
          Rprintf("  iter %.3d / %.3d\n", iter + 1, numWarmup);
        
        ranefSampler->run(true);
        
        sample_writer.increment();
      }
      if (refresh > 0 && verbose > 1 && (iter + 1) % refresh == 0)
          Rprintf("  iter %.3d / %.3d\n", iter + 1, numWarmup);
      
      ranefSampler->run(true);
      
      ranefModel->get_Zb(sample_writer.x_curr + sample_writer_offset, fixef_offset);
      if (userOffset != NULL) {
        if (offset_type == OFFSET_DEFAULT)
          for (size_t j = 0; j < bartData.numObservations; ++j) fixef_offset[j] += userOffset[j];
        else if (offset_type == OFFSET_RANEF)
          std::memcpy(fixef_offset, userOffset, bartData.numObservations * sizeof(double));
      }
      double sigma = ranefModel->get_aux(sample_writer.x_curr + sample_writer_offset);
      sample_writer.increment();
      
      std::memcpy(fixef_response, y_orig, bartData.numObservations * sizeof(double));
      for (size_t j = 0; j < bartData.numObservations; ++j) {
        fixef_response[j] -= fixef_offset[j];
        fixef_offset[j] = 0.0;
      }
      bartFunctions.setOffset(fixefModel, fixef_offset);
      bartFunctions.setResponse(fixefModel, fixef_response);
      y_scale = getRange(fixef_response, bartData.numObservations);
      bartModel.sigmaSqPrior->setScale((sigma / y_scale) * (sigma / y_scale));
      fixefModel->state[0].sigma = sigma / y_scale;
      ++iter;
      
      // run for rest of warmup swapping between two
      for (  ; iter < numWarmup; ++iter) {
        if (refresh > 0 && verbose > 1 && (iter + 1) % refresh == 0)
          Rprintf("  iter %.3d / %.3d\n", iter + 1, numWarmup);
        
        bartFunctions.runSamplerWithResults(fixefModel, 0, fixefWarmup);
        
        // bart with an offset will produce predictions that have the offset added;
        // in order to just get the tree predictions, subtract out that offset
        for (size_t j = 0; j < bartData.numObservations; ++j) fixefWarmup->trainingSamples[j] -= fixef_offset[j];
        
        std::memcpy(ranef_offset, const_cast<const double*>(fixefWarmup->trainingSamples), bartData.numObservations * sizeof(double));
        
        if (userOffset != NULL) {
          if (offset_type == OFFSET_DEFAULT)
            for (size_t j = 0; j < bartData.numObservations; ++j) ranef_offset[j] += userOffset[j];
          else if (offset_type == OFFSET_FIXEF)
            std::memcpy(ranef_offset, userOffset, bartData.numObservations * sizeof(double));
        }
        
        ranefModel->set_offset(ranef_offset);
        
        if (save_warmup) {
          fixefWarmup->incrementPointers();
          if ((fixefWarmup->trainingSamples - fixefWarmup->base_trainingSamples) / bartData.numObservations != (iter + 1))
            Rf_error("samples beyond extent");
        }
      
        ranefSampler->run(true);
        
        ranefModel->get_Zb(sample_writer.x_curr + sample_writer_offset, fixef_offset);
        if (userOffset != NULL) {
          if (offset_type == OFFSET_DEFAULT)
            for (size_t j = 0; j < bartData.numObservations; ++j) fixef_offset[j] += userOffset[j];
          else if (offset_type == OFFSET_RANEF)
            std::memcpy(fixef_offset, userOffset, bartData.numObservations * sizeof(double));
        }
        sigma = ranefModel->get_aux(sample_writer.x_curr + sample_writer_offset);
        sample_writer.increment();
        
        std::memcpy(fixef_response, y_orig, bartData.numObservations * sizeof(double));
        for (size_t j = 0; j < bartData.numObservations; ++j) {
          fixef_response[j] -= fixef_offset[j];
          fixef_offset[j] = 0.0;
        }
        
        bartFunctions.setOffset(fixefModel, fixef_offset);
        bartFunctions.setResponse(fixefModel, fixef_response);
        y_scale = getRange(fixef_response, bartData.numObservations);
        bartModel.sigmaSqPrior->setScale((sigma / y_scale) * (sigma / y_scale));
        fixefModel->state[0].sigma = sigma / y_scale;
      }
      fixefWarmup->resetPointers();
    }
    
    ranefSampler->disengage_adaptation();
    
    if (num_iter > 0) {
      if (verbose > 0) Rprintf("starting sample, %d iters\n", num_iter);
      
      for (unsigned int i = 0; i < num_iter; ++i) {
        if (refresh > 0 && verbose > 1 && (i + 1) % refresh == 0)
          Rprintf("  iter %.3d / %.3d\n", i + 1, num_iter);
        
        bartFunctions.runSamplerWithResults(fixefModel, 0, fixefSamples);
        
        // bart with an offset will produce predictions that have the offset added;
        // in order to just get the tree predictions, subtract out that offset
        for (size_t j = 0; j < bartData.numObservations; ++j) fixefSamples->trainingSamples[j] -= fixef_offset[j];
        
        std::memcpy(ranef_offset, const_cast<const double*>(fixefSamples->trainingSamples), bartData.numObservations * sizeof(double));
        if (userOffset != NULL) {
          if (offset_type == OFFSET_DEFAULT)
            for (size_t j = 0; j < bartData.numObservations; ++j) ranef_offset[j] += userOffset[j];
          else if (offset_type == OFFSET_FIXEF)
            std::memcpy(ranef_offset, userOffset, bartData.numObservations * sizeof(double));
        }
        
        ranefModel->set_offset(ranef_offset);
        
        fixefSamples->incrementPointers();
      
        ranefSampler->run(false);
        
        ranefModel->get_Zb(sample_writer.x_curr + sample_writer_offset, fixef_offset);
        
        if (userOffset != NULL) {
          if (offset_type == OFFSET_DEFAULT)
            for (size_t j = 0; j < bartData.numObservations; ++j) fixef_offset[j] += userOffset[j];
          else if (offset_type == OFFSET_RANEF)
            std::memcpy(fixef_offset, userOffset, bartData.numObservations * sizeof(double));
        }
        double sigma = ranefModel->get_aux(sample_writer.x_curr + sample_writer_offset);
        sample_writer.increment();
        
        std::memcpy(fixef_response, y_orig, bartData.numObservations * sizeof(double));
        for (size_t j = 0; j < bartData.numObservations; ++j) {
          fixef_response[j] -= fixef_offset[j];
          fixef_offset[j] = 0.0;
        }
        bartFunctions.setOffset(fixefModel, fixef_offset);
        bartFunctions.setResponse(fixefModel, fixef_response);
        y_scale = getRange(fixef_response, bartData.numObservations);
        bartModel.sigmaSqPrior->setScale((sigma / y_scale) * (sigma / y_scale));
        fixefModel->state[0].sigma = sigma / y_scale;
      }
      fixefSamples->resetPointers();
    }
    
    delete [] fixef_offset;
    delete [] fixef_response;
    delete [] ranef_offset;
    */ 
    /* SEXP resultExpr = PROTECT(rc_newList(save_warmup ? 3 : 2));
    int pos = 0;
    if (save_warmup)
      SET_VECTOR_ELT(resultExpr, pos++, createBartResultsExpr(*fixefModel, *fixefWarmup));
    SET_VECTOR_ELT(resultExpr, pos++, createBartResultsExpr(*fixefModel, *fixefSamples));
    SET_VECTOR_ELT(resultExpr, pos, createStanResultsExpr(sample_writer));
    
    SEXP namesExpr = PROTECT(rc_newCharacter(XLENGTH(resultExpr)));
    pos = 0;
    if (save_warmup)
      SET_STRING_ELT(namesExpr, pos++, Rf_mkChar("bart_warmup"));
    SET_STRING_ELT(namesExpr, pos++, Rf_mkChar("bart_sample"));
    SET_STRING_ELT(namesExpr, pos, Rf_mkChar("stan_sample"));
    
    rc_setNames(resultExpr, namesExpr);
    UNPROTECT(1);
    
    delete fixefWarmup;
    delete fixefSamples;
    
    // delete bart model
    bartFunctions.invalidateFit(fixefModel);
    ::operator delete(fixefModel);
    bartFunctions.invalidateModel(&bartModel);
    bartFunctions.invalidateData(&bartData);
    
    // delete stan model
    delete ranefSampler;
    delete ranefModel;
    
    UNPROTECT(1);
    
    return resultExpr;
  } */
}

// this unusual set of declarations solves a rather obscure warning on Solaris
typedef void* (*C_voidPtrFunction)(void);
extern "C" typedef C_voidPtrFunction (*C_voidPtrFunctionLookup)(const char* _namespace, const char* name);

namespace {
  
  void lookupBARTFunctions()
  {
    bartFunctions.initializeFit         = reinterpret_cast<void (*)(dbarts::BARTFit*, dbarts::Control*, dbarts::Model*, dbarts::Data*)>(R_GetCCallable("dbarts", "initializeFit"));
    bartFunctions.invalidateFit         = reinterpret_cast<void (*)(dbarts::BARTFit*)>(R_GetCCallable("dbarts", "invalidateFit"));
    bartFunctions.initializeControl     = reinterpret_cast<void (*)(dbarts::Control*, SEXP)>(R_GetCCallable("dbarts", "initializeControl"));
    bartFunctions.initializeData        = reinterpret_cast<void (*)(dbarts::Data*, SEXP)>(R_GetCCallable("dbarts", "initializeData"));
    bartFunctions.invalidateData        = reinterpret_cast<void (*)(dbarts::Data*)>(R_GetCCallable("dbarts", "invalidateData"));
    bartFunctions.initializeModel       = reinterpret_cast<void (*)(dbarts::Model*, SEXP, const dbarts::Control*)>(R_GetCCallable("dbarts", "initializeModel"));
    bartFunctions.invalidateModel       = reinterpret_cast<void (*)(dbarts::Model*)>(R_GetCCallable("dbarts", "invalidateModel"));
    bartFunctions.setResponse           = reinterpret_cast<void (*)(dbarts::BARTFit*, const double*)>(R_GetCCallable("dbarts", "setResponse"));
    bartFunctions.setOffset             = reinterpret_cast<void (*)(dbarts::BARTFit*, const double*, bool)>(R_GetCCallable("dbarts", "setOffset"));
    bartFunctions.setSigma              = reinterpret_cast<void (*)(dbarts::BARTFit*, const double*)>(R_GetCCallable("dbarts", "setSigma"));
    
    bartFunctions.runSamplerWithResults = reinterpret_cast<void (*)(dbarts::BARTFit*, std::size_t, dbarts::Results*)>(R_GetCCallable("dbarts", "runSamplerWithResults"));
    bartFunctions.sampleTreesFromPrior  = reinterpret_cast<void (*)(dbarts::BARTFit*)>(R_GetCCallable("dbarts", "sampleTreesFromPrior"));
    bartFunctions.printInitialSummary   = reinterpret_cast<void (*)(const dbarts::BARTFit*)>(R_GetCCallable("dbarts", "printInitialSummary"));
  }
}

extern "C" {

static void samplerFinalizer(SEXP samplerExpr)
{
  Sampler* sampler = static_cast<Sampler*>(R_ExternalPtrAddr(samplerExpr));
  
  if (sampler == NULL) return;
  
  if (activeSamplers->find(samplerExpr) == activeSamplers->end()) return;
  
  activeSamplers->erase(samplerExpr);
  
  delete sampler;
  
  R_ClearExternalPtr(samplerExpr);
}

static SEXP isValidPointer(SEXP samplerExpr)
{
  Sampler* sampler = static_cast<Sampler*>(R_ExternalPtrAddr(samplerExpr));
  if (sampler == NULL) return Rf_ScalarLogical(FALSE);
  
  if (activeSamplers->find(samplerExpr) != activeSamplers->end())
    return Rf_ScalarLogical(TRUE);
  
  return Rf_ScalarLogical(FALSE);
}

static SEXP finalize(void)
{
  for (PointerSet::iterator it = activeSamplers->begin(); it != activeSamplers->end(); ) {
    SEXP samplerExpr = *it;
    Sampler* sampler = static_cast<Sampler*>(R_ExternalPtrAddr(samplerExpr));
    
    delete sampler;
    PointerSet::iterator prev = it;
    ++it;
    activeSamplers->erase(prev);
    R_ClearExternalPtr(samplerExpr);
  }
      
  delete activeSamplers;
  
  return R_NilValue;
}

#define DEF_FUNC(_N_, _F_, _A_) { _N_, reinterpret_cast<DL_FUNC>(&_F_), _A_ }

static R_CallMethodDef R_callMethods[] = {
  DEF_FUNC("stan4bart_create", createSampler, 6),
  DEF_FUNC("stan4bart_run", run, 4),
  DEF_FUNC("stan4bart_printInitialSummary", printInitialSummary, 1),
  DEF_FUNC("stan4bart_disengageAdaptation", disengageAdaptation, 1),
  {NULL, NULL, 0}
};

#undef DEF_FUNC

void attribute_visible R_init_stan4bart(DllInfo *info) {
  R_registerRoutines(info, NULL, R_callMethods, NULL, NULL);
  R_useDynamicSymbols(info, static_cast<Rboolean>(FALSE));
  
  lookupBARTFunctions();
  
  activeSamplers = new PointerSet(&compareExternalPointers);
}

} // extern "C"
