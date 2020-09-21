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
    
    const double* bart_offset_init = REAL(rc_getListElement(commonControlExpr, "bart_offset_init"));
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
    
    if (bart_offset_init != NULL) {
      std::memcpy(sampler.bartOffset, bart_offset_init, n * sizeof(double));
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
  
  static SEXP get_parametric_mean(SEXP samplerExpr)
  {
    Sampler* samplerPtr = static_cast<Sampler*>(R_ExternalPtrAddr(samplerExpr));
    if (samplerPtr == NULL) Rf_error("run called on NULL external pointer");
    Sampler& sampler(*samplerPtr);
    
    sampler.stanSampler->sample_writer.decrement();
    SEXP result = PROTECT(rc_newReal(sampler.bartData.numObservations));
    
    stan4bart::getParametricMean(*sampler.stanSampler, *sampler.stanModel, REAL(result));
    sampler.stanSampler->sample_writer.increment();
    
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
        sampler.bartData.numObservations, sampler.bartData.numPredictors, sampler.bartData.numTestObservations, numIter, 1 /* num chains */, false /* TODO: binary */);
    if (resultsType == RESULTS_BOTH || resultsType == RESULTS_STAN)
      sampler.stanSampler->sample_writer.resize(sampler.stanSampler->num_pars, numIter);
    
    size_t n = sampler.bartData.numObservations;
    
    if (sampler.verbose > 0)
      Rprintf("starting %s, %d draws, %s:\n", isWarmup ? "warmup" : "sampling", numIter,
              resultsType == RESULTS_BOTH ? "both BART and Stan" : (resultsType == RESULTS_BART ? "BART only" : "Stan only"));
    
    // TODO: add in capacity for fitting against true ranef and true fixef
    for (int iter = 0; iter < numIter; ++iter) {
      if (sampler.refresh > 0 && sampler.verbose > 1 && (iter + 1) % sampler.refresh == 0)
          Rprintf("  iter %.3d / %.3d\n", iter + 1, numIter);
      
      if (resultsType == RESULTS_BOTH || resultsType == RESULTS_BART) {
        
        bartFunctions.runSamplerWithResults(sampler.bartSampler, 0, bartSamples);
        
        // bart with an offset will produce predictions that have the offset added;
        // in order to just get the tree predictions, subtract out that offset
        for (size_t j = 0; j < n; ++j)
          bartSamples->trainingSamples[j] -= sampler.bartOffset[j];
        
        std::memcpy(sampler.stanOffset, const_cast<const double*>(bartSamples->trainingSamples), n * sizeof(double));
        
        if (sampler.userOffset != NULL) {
          if (sampler.offsetType == OFFSET_DEFAULT)
            for (size_t j = 0; j < n; ++j)
              sampler.stanOffset[j] += sampler.userOffset[j];
          else if (sampler.offsetType == OFFSET_FIXEF)
            std::memcpy(sampler.stanOffset, sampler.userOffset, n * sizeof(double));
        }
        stan4bart::setStanOffset(*sampler.stanModel, sampler.stanOffset);
        
        bartSamples->incrementPointers();
      }
      
      if (resultsType == RESULTS_BOTH || resultsType == RESULTS_STAN) {
        sampler.stanSampler->run(isWarmup);
        
        stan4bart::getParametricMean(*sampler.stanSampler, *sampler.stanModel, sampler.bartOffset);
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
  DEF_FUNC("stan4bart_finalize", finalize, 0),
  DEF_FUNC("stan4bart_get_parametric_mean", get_parametric_mean, 1),
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
