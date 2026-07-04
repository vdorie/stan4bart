#include <ext/Rinternals.h>
#include <R_ext/Arith.h> // R_NaReal
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#ifdef HAVE_STD_SNPRINTF
// snprintf in c++11, before that have to use C version
#  include <cstdio>
using std::snprintf;
#else
#  include <stdio.h>
#endif

#include <cstdint>
#include <exception>
#include <memory> // unique_ptr
#include <set> // external pointers set
#include <vector>

#include <rc/util.h>
#include <rc/bounds.h>

#include <dbarts/dbarts.h>

#include "rstan/io/r_ostream.hpp"

#include "bart_util.hpp"
#include "stan_sampler.hpp"

#if __cplusplus < 201112L
#  if defined(_WIN64) || SIZEOF_SIZE_T == 8
#    define SIZE_T_SPECIFIER "%lu"
#  else
#    define SIZE_T_SPECIFIER "%u"
#  endif
#else
#  define SIZE_T_SPECIFIER "%zu"
#endif

namespace {
  // stuff to handle external pointers
  typedef bool(*ExternalPointerComparator)(const SEXP&lhs, const SEXP& rhs);
  typedef std::set<SEXP, ExternalPointerComparator> PointerSet;
  
  PointerSet* activeSamplers;
  PointerSet* activeStoredBARTSamplers;
  
  bool compareExternalPointers(const SEXP& lhs, const SEXP& rhs) {
    return R_ExternalPtrAddr(const_cast<SEXP>(lhs)) < R_ExternalPtrAddr(const_cast<SEXP>(rhs));
  }
  
  // the flat C API of dbarts.h, resolved through R_GetCCallable at load
  struct BARTFunctionTable {
    int (*apiVersion)(void);
    dbarts_sampler* (*create)(SEXP control, SEXP model, SEXP data,
                              const char* family);
    void (*destroy)(dbarts_sampler* sampler);
    void (*run)(dbarts_sampler* sampler, std::size_t numBurnIn,
                std::size_t numSamples, dbarts_results* results);
    void (*sampleTreesFromPrior)(dbarts_sampler* sampler);
    void (*setOffset)(dbarts_sampler* sampler, const double* offset,
                      int updateScale);
    void (*setSigma)(dbarts_sampler* sampler, double sigma);
    int (*getLatents)(const dbarts_sampler* sampler, double* out);
    void (*predict)(dbarts_sampler* sampler, const double* x_test,
                    std::size_t numTestObservations,
                    const double* offset_test, double* out);
    void (*setTreeStorage)(dbarts_sampler* sampler, int keepTrees,
                           std::size_t numSamplesToStore);
    SEXP (*getTrees)(dbarts_sampler* sampler, const std::size_t*,
                     std::size_t, const std::size_t*, std::size_t,
                     const std::size_t*, std::size_t, int useLiveTrees);
    void (*printTrees)(dbarts_sampler* sampler, const std::size_t*,
                       std::size_t, const std::size_t*, std::size_t,
                       const std::size_t*, std::size_t);
    SEXP (*storeState)(dbarts_sampler* sampler);
    void (*setState)(dbarts_sampler* sampler, SEXP state);
    void (*setVerbose)(dbarts_sampler* sampler, int verbose,
                       std::uint32_t printEvery);
    std::size_t (*numObservations)(const dbarts_sampler* sampler);
    std::size_t (*numPredictors)(const dbarts_sampler* sampler);
    std::size_t (*numTestObservations)(const dbarts_sampler* sampler);
    std::size_t (*numChains)(const dbarts_sampler* sampler);
    std::size_t (*numTrees)(const dbarts_sampler* sampler);
    std::size_t (*numSavedSamples)(const dbarts_sampler* sampler);
    int (*kIsSampled)(const dbarts_sampler* sampler);
  };
  BARTFunctionTable bartFunctions;
  
  enum UserOffsetType {
    OFFSET_DEFAULT = 0,
    OFFSET_FIXEF,
    OFFSET_RANEF,
    OFFSET_BART,
    OFFSET_PARAMETRIC
  };
  
  const char* const userOffsetTypeNames[] = {
    "default",
    "fixef",
    "ranef",
    "bart",
    "parametric"
  };
  
  enum ResultsType {
    RESULTS_BOTH = 0,
    RESULTS_BART,
    RESULTS_STAN
  };
  
  // used for predict
  struct StoredBARTSampler {
    dbarts_sampler* fit;
    
    StoredBARTSampler() : fit(NULL) { }
    ~StoredBARTSampler() {
      if (fit != NULL) {
        bartFunctions.destroy(fit);
        fit = NULL;
      }
    }
  };
  
  struct Sampler {
    int defaultWarmup;
    int defaultIter;
    int verbose;
    int refresh;
    bool responseIsBinary;
    const double* userOffset;
    UserOffsetType offsetType;
    
    continuous_model_namespace::continuous_model* stanModel;
    stan4bart::StanControl stanControl;
    stan4bart::StanSampler* stanSampler;
    
    dbarts_sampler* bartSampler;
    std::size_t numObservations;
    std::size_t numTestObservations;
    bool kIsSampled;
    bool keepTrees;
    
    double* bartOffset;
    double* stanOffset;
    double* bartLatents;
    
    bool keepFits;
    SEXP callback;
    SEXP callbackEnv;
    
    Sampler() :
      stanModel(NULL), stanSampler(NULL), bartSampler(NULL),
      numObservations(0), numTestObservations(0), kIsSampled(false),
      keepTrees(false), bartOffset(NULL), stanOffset(NULL), bartLatents(NULL)
    {
    }
    ~Sampler() {
      delete [] bartLatents;
      delete [] stanOffset;
      delete [] bartOffset;
      
      if (bartSampler != NULL) {
        bartFunctions.destroy(bartSampler);
        bartSampler = NULL;
      }
      
      delete stanSampler;
      stan4bart::deleteStanModel(stanModel);
      stanModel = NULL;
    }
  };
  
  void initializeSamplerFromExpression(Sampler& sampler, SEXP commonControlExpr);
}

extern "C" {
  static void samplerFinalizer(SEXP samplerExpr);
  static void storedBARTSamplerFinalizer(SEXP samplerExpr);
  
#ifdef __clang__
#  if __has_warning("-Wenum-enum-conversion")
#    define SUPPRESS_ENUM_CONVERSION_WARNING 1
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Wenum-enum-conversion"
#  endif
#endif
  
  static SEXP createSampler(SEXP bartControlExpr, SEXP bartDataExpr, SEXP bartModelExpr,
                            SEXP stanDataExpr, SEXP stanControlExpr,
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
    stan4bart::initializeStanControlFromExpression(sampler.stanControl, stanControlExpr);
    if (sampler.stanControl.skip == R_NaInt) {
      sampler.stanControl.skip = (2000 - sampler.defaultWarmup) / 1000;
      if (sampler.stanControl.skip < 1) sampler.stanControl.skip = 1;
    }
    
    int chain_id = 1;
    sampler.stanSampler = new stan4bart::StanSampler(*sampler.stanModel, sampler.stanControl, chain_id, sampler.defaultWarmup, -1);
    sampler.stanSampler->setVerbose(sampler.verbose);
    
    // a verbose control prints the initial summary during creation; tree
    // storage is turned on only for the recorded sampling run
    sampler.keepTrees = rc_getBool(
      Rf_getAttrib(bartControlExpr, Rf_install("keepTrees")), "keepTrees",
      RC_NA | RC_NO, RC_END);
    sampler.bartSampler = bartFunctions.create(bartControlExpr,
                                               bartModelExpr, bartDataExpr,
                                               "");
    bartFunctions.setVerbose(sampler.bartSampler, 0, 100);
    bartFunctions.setTreeStorage(sampler.bartSampler, 0, 0);
    sampler.kIsSampled = bartFunctions.kIsSampled(sampler.bartSampler) != 0;
    sampler.numObservations =
      bartFunctions.numObservations(sampler.bartSampler);
    sampler.numTestObservations =
      bartFunctions.numTestObservations(sampler.bartSampler);
    
    size_t n = sampler.numObservations;
    sampler.bartOffset = new double[n];
    sampler.stanOffset = new double[n];
    if (sampler.responseIsBinary)
      sampler.bartLatents = new double[n];
    
    if (sampler.userOffset != NULL) {
      if (sampler.offsetType != OFFSET_BART) {
        std::memcpy(sampler.bartOffset, sampler.userOffset, n * sizeof(double));
        
        if (bart_offset_init != NULL && sampler.offsetType == OFFSET_DEFAULT)
          for (size_t i = 0; i < n; ++i) sampler.bartOffset[i] += bart_offset_init[i];
      } else if (bart_offset_init != NULL) {
        std::memcpy(sampler.bartOffset, bart_offset_init, n * sizeof(double));
      } else {
        for (size_t i = 0; i < n; ++i) sampler.bartOffset[i] = 0.0;
      }
    } else {
      if (bart_offset_init != NULL) {
        std::memcpy(sampler.bartOffset, bart_offset_init, n * sizeof(double));
      } else {
        for (size_t i = 0; i < n; ++i) sampler.bartOffset[i] = 0.0;
      }
    }
    
    bartFunctions.setOffset(sampler.bartSampler, sampler.bartOffset, true);
    if (!sampler.responseIsBinary)
      bartFunctions.setSigma(sampler.bartSampler, sigma_init);

    // the entry points that draw bracket R's RNG state internally
    bartFunctions.sampleTreesFromPrior(sampler.bartSampler);
    
    // draw once before running; only the training fits are needed
    std::vector<double> firstDraw(n);
    dbarts_results firstResults = {};
    firstResults.train = firstDraw.data();
    bartFunctions.run(sampler.bartSampler, 0, 1, &firstResults);
    
    for (size_t j = 0; j < n; ++j) firstDraw[j] -= sampler.bartOffset[j];
    
    if (sampler.userOffset != NULL && sampler.offsetType == OFFSET_BART) {
      // Override with user supplied bart offset
      std::memcpy(sampler.stanOffset, sampler.userOffset, n * sizeof(double));
    } else {
      std::memcpy(sampler.stanOffset, firstDraw.data(), n * sizeof(double));
      if (sampler.userOffset != NULL && sampler.offsetType == OFFSET_DEFAULT)
        for (size_t j = 0; j < n; ++j)
          sampler.stanOffset[j] += sampler.userOffset[j];
    }
    
    stan4bart::setStanOffset(*sampler.stanModel, sampler.stanOffset);
    if (sampler.responseIsBinary) {
      bartFunctions.getLatents(sampler.bartSampler, sampler.bartLatents);
      stan4bart::setResponse(*sampler.stanModel, sampler.bartLatents);
    }
    
    SEXP result = PROTECT(R_MakeExternalPtr(samplerPtr.get(), R_NilValue, R_NilValue));
    samplerPtr.release();
    R_RegisterCFinalizerEx(result, samplerFinalizer, static_cast<Rboolean>(FALSE));
    

    activeSamplers->insert(result);

    UNPROTECT(1);
    
    return result;
  }
  
#ifdef SUPPRESS_ENUM_CONVERSION_WARNING
#  pragma clang diagnostic pop
#endif
  
  static SEXP getParametricMean(SEXP samplerExpr)
  {
    Sampler* samplerPtr = static_cast<Sampler*>(R_ExternalPtrAddr(samplerExpr));
    if (samplerPtr == NULL) Rf_error("getParametricMean called on NULL external pointer");
    Sampler& sampler(*samplerPtr);
    
    sampler.stanSampler->sample_writer.decrement();
    SEXP result = PROTECT(rc_newReal(sampler.numObservations));
    
    sampler.stanSampler->getParametricMean(*sampler.stanModel, REAL(result));
    sampler.stanSampler->sample_writer.increment();
    
    UNPROTECT(1);
    
    return result;
  }
  
#ifdef SUPPRESS_ENUM_CONVERSION_WARNING
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wenum-enum-conversion"
#endif
  
  static SEXP predictBART(SEXP storedBARTSamplerExpr, SEXP x_testExpr, SEXP offset_testExpr)
  {
    StoredBARTSampler* samplerPtr = static_cast<StoredBARTSampler*>(R_ExternalPtrAddr(storedBARTSamplerExpr));
    if (samplerPtr == NULL) Rf_error("predictBART called on NULL external pointer");
    StoredBARTSampler& sampler(*samplerPtr);
    
    dbarts_sampler* fit(sampler.fit);
        
    if (Rf_isNull(x_testExpr)) return R_NilValue;
    
    if (!Rf_isReal(x_testExpr)) Rf_error("x.test must be of type real");
    
    rc_assertDimConstraints(x_testExpr, "dimensions of x_test", RC_LENGTH | RC_EQ, rc_asRLength(2),
                            RC_NA,
                            RC_VALUE | RC_EQ, static_cast<int>(bartFunctions.numPredictors(fit)),
                            RC_END);
    int* dims = INTEGER(Rf_getAttrib(x_testExpr, R_DimSymbol));
    
    size_t numChains = bartFunctions.numChains(fit);
    size_t numSavedSamples = bartFunctions.numSavedSamples(fit);
    bool usesSavedTrees = numSavedSamples > 0;
    size_t numSamples = usesSavedTrees ? numSavedSamples : 1;
    size_t numTestObservations = static_cast<size_t>(dims[0]);
    
    double* testOffset = NULL;
    if (!Rf_isNull(offset_testExpr)) {
      if (!Rf_isReal(offset_testExpr)) Rf_error("offset.test must be of type real");
      if (rc_getLength(offset_testExpr) != 1 || !ISNA(REAL(offset_testExpr)[0])) {
        if (rc_getLength(offset_testExpr) != numTestObservations) Rf_error("length of offset.test must equal number of rows in x.test");
        testOffset = REAL(offset_testExpr);
      }
    }
    
    SEXP result = PROTECT(Rf_allocVector(REALSXP, numTestObservations * numSamples * numChains));
    if (usesSavedTrees) {
      if (numChains <= 1)
        rc_setDims(result, static_cast<int>(numTestObservations), static_cast<int>(numSamples), -1);
      else
        rc_setDims(result, static_cast<int>(numTestObservations), static_cast<int>(numSamples), static_cast<int>(numChains), -1);
    } else {
      if (numChains > 1)
        rc_setDims(result, static_cast<int>(numTestObservations), static_cast<int>(numChains), -1);
    }
    
    // predictions arrive on the original response scale: the restored state
    // carries the fit's transform, so no R-side rescaling remains
    bartFunctions.predict(fit, REAL(x_testExpr), numTestObservations, testOffset, REAL(result));
    
    UNPROTECT(1);
    
    return result;
  }
  
#ifdef SUPPRESS_ENUM_CONVERSION_WARNING
#  pragma clang diagnostic pop
#endif
  
  static SEXP exportBARTState(SEXP samplerExpr)
  {
    Sampler* samplerPtr = static_cast<Sampler*>(R_ExternalPtrAddr(samplerExpr));
    if (samplerPtr == NULL) Rf_error("exportBARTState called on NULL external pointer");
    Sampler& sampler(*samplerPtr);
    
    return bartFunctions.storeState(sampler.bartSampler);
  }
  
  static SEXP createStoredBARTSampler(SEXP controlExpr, SEXP dataExpr, SEXP modelExpr, SEXP stateExpr)
  {
    std::unique_ptr<StoredBARTSampler> samplerPtr(new StoredBARTSampler);
    StoredBARTSampler& sampler(*samplerPtr);
    
    // the R side sizes the control for restoration (n.chains matching the
    // state, keepTrees with n.samples at the saved capacity); the state
    // carries the fit's response transform, so no scale pokes remain
    sampler.fit = bartFunctions.create(controlExpr, modelExpr, dataExpr, "");
    bartFunctions.setVerbose(sampler.fit, 0, 100);
    bartFunctions.setState(sampler.fit, stateExpr);
    
    SEXP result = PROTECT(R_MakeExternalPtr(samplerPtr.get(), R_NilValue, R_NilValue));
    samplerPtr.release();
    R_RegisterCFinalizerEx(result, storedBARTSamplerFinalizer, static_cast<Rboolean>(FALSE));

    activeStoredBARTSamplers->insert(result);

    UNPROTECT(1);
    
    return result;
  }

  static SEXP printTrees(SEXP storedBARTSamplerExpr, SEXP chainIndicesExpr, SEXP sampleIndicesExpr, SEXP treeIndicesExpr)
  {
    StoredBARTSampler* samplerPtr = static_cast<StoredBARTSampler*>(R_ExternalPtrAddr(storedBARTSamplerExpr));
    if (samplerPtr == NULL) Rf_error("printTrees called on NULL external pointer");
    StoredBARTSampler& sampler(*samplerPtr);
    
    dbarts_sampler* fit(sampler.fit);
    
    size_t numChains  = bartFunctions.numChains(fit);
    size_t numSamples = bartFunctions.numSavedSamples(fit);
    size_t numTrees   = bartFunctions.numTrees(fit);

    size_t numChainIndices  = Rf_isNull(chainIndicesExpr)  ? numChains  : rc_getLength(chainIndicesExpr);
    size_t numSampleIndices = Rf_isNull(sampleIndicesExpr) ? numSamples : rc_getLength(sampleIndicesExpr);
    size_t numTreeIndices   = Rf_isNull(treeIndicesExpr)   ? numTrees   : rc_getLength(treeIndicesExpr);

    if (numChainIndices > numChains)
      Rf_error(SIZE_T_SPECIFIER " chains specified but only " SIZE_T_SPECIFIER " in sampler", numChainIndices, numChains);
    if (numSampleIndices > numSamples)
      Rf_error(SIZE_T_SPECIFIER " samples specified but only " SIZE_T_SPECIFIER " in sampler", numSampleIndices, numSamples);
    if (numTreeIndices > numTrees)
      Rf_error(SIZE_T_SPECIFIER " trees specified but only " SIZE_T_SPECIFIER " in sampler", numTreeIndices, numTrees);
    
    size_t* chainIndices  = new size_t[numChainIndices];
    size_t* sampleIndices = new size_t[numSamples];
    size_t* treeIndices   = new size_t[numTreeIndices];
    
    if (Rf_isNull(chainIndicesExpr)) {
      for (size_t i = 0; i < numChains; ++i) chainIndices[i] = i;
    } else {
      int* i_chainIndices = INTEGER(chainIndicesExpr);
      for (size_t i = 0; i < numChainIndices; ++i) chainIndices[i] = static_cast<size_t>(i_chainIndices[i] - 1);
    }
    
    if (Rf_isNull(sampleIndicesExpr)) {
      for (size_t i = 0; i < numSamples; ++i) sampleIndices[i] = i;
    } else {
      int* i_sampleIndices = INTEGER(sampleIndicesExpr);
      for (size_t i = 0; i < numSampleIndices; ++i) sampleIndices[i] = static_cast<size_t>(i_sampleIndices[i] - 1);
    }
    
    if (Rf_isNull(treeIndicesExpr)) {
      for (size_t i = 0; i < numTrees; ++i) treeIndices[i] = i;
    } else {
      int* i_treeIndices = INTEGER(treeIndicesExpr);
      for (size_t i = 0; i < numTreeIndices; ++i) treeIndices[i] = static_cast<size_t>(i_treeIndices[i] - 1);
    }
    
    bartFunctions.printTrees(fit, chainIndices, numChainIndices, sampleIndices, numSampleIndices, treeIndices, numTreeIndices);

    delete [] treeIndices;
    delete [] sampleIndices;
    delete [] chainIndices;
    
    return R_NilValue;
  }
  
  static SEXP getTrees(SEXP storedBARTSamplerExpr, SEXP chainIndicesExpr, SEXP sampleIndicesExpr, SEXP treeIndicesExpr, SEXP currentExpr)
  {
    StoredBARTSampler* samplerPtr = static_cast<StoredBARTSampler*>(R_ExternalPtrAddr(storedBARTSamplerExpr));
    if (samplerPtr == NULL) Rf_error("getTrees called on NULL external pointer");
    StoredBARTSampler& sampler(*samplerPtr);
    
    dbarts_sampler* fit(sampler.fit);

    // when currentExpr is true, return the live working trees even for a
    // keepTrees sampler; there is then no sample dimension
    bool useLiveTrees = Rf_asLogical(currentExpr) == TRUE;
    bool treatAsSaved = bartFunctions.numSavedSamples(fit) > 0 && !useLiveTrees;
     
    size_t numChains  = bartFunctions.numChains(fit);
    size_t numSamples = treatAsSaved ? bartFunctions.numSavedSamples(fit) : 0;
    size_t numTrees   = bartFunctions.numTrees(fit);

    size_t numChainIndices  = Rf_isNull(chainIndicesExpr)  ? numChains  : rc_getLength(chainIndicesExpr);
    size_t numSampleIndices = Rf_isNull(sampleIndicesExpr) ? numSamples : rc_getLength(sampleIndicesExpr);
    size_t numTreeIndices   = Rf_isNull(treeIndicesExpr)   ? numTrees   : rc_getLength(treeIndicesExpr);
    
    if (numChainIndices > numChains)
      Rf_error(SIZE_T_SPECIFIER " chains specified but only " SIZE_T_SPECIFIER " in sampler", numChainIndices, numChains);
    if (numSampleIndices > numSamples)
      Rf_error(SIZE_T_SPECIFIER " samples specified but only " SIZE_T_SPECIFIER " in sampler", numSampleIndices, numSamples);
    if (numTreeIndices > numTrees)
      Rf_error(SIZE_T_SPECIFIER " trees specified but only " SIZE_T_SPECIFIER " in sampler", numTreeIndices, numTrees);
    
    size_t* chainIndices  = new size_t[numChainIndices];
    size_t* sampleIndices = treatAsSaved ? new size_t[numSamples] : NULL;
    size_t* treeIndices   = new size_t[numTreeIndices];
    
    if (Rf_isNull(chainIndicesExpr)) {
      for (size_t i = 0; i < numChains; ++i) chainIndices[i] = i;
    } else {
      int* i_chainIndices = INTEGER(chainIndicesExpr);
      for (size_t i = 0; i < numChainIndices; ++i) chainIndices[i] = static_cast<size_t>(i_chainIndices[i] - 1);
    }
    
    if (Rf_isNull(sampleIndicesExpr)) {
      for (size_t i = 0; i < numSamples; ++i) sampleIndices[i] = i;
    } else {
      int* i_sampleIndices = INTEGER(sampleIndicesExpr);
      for (size_t i = 0; i < numSampleIndices; ++i) sampleIndices[i] = static_cast<size_t>(i_sampleIndices[i] - 1);
    }
    
    if (Rf_isNull(treeIndicesExpr)) {
      for (size_t i = 0; i < numTrees; ++i) treeIndices[i] = i;
    } else {
      int* i_treeIndices = INTEGER(treeIndicesExpr);
      for (size_t i = 0; i < numTreeIndices; ++i) treeIndices[i] = static_cast<size_t>(i_treeIndices[i] - 1);
    }
        
    SEXP resultExpr = PROTECT(bartFunctions.getTrees(
      fit, chainIndices, numChainIndices, sampleIndices, numSampleIndices,
      treeIndices, numTreeIndices, useLiveTrees ? 1 : 0));
    
    delete [] treeIndices;
    delete [] sampleIndices;
    delete [] chainIndices;
    
    UNPROTECT(1);
    
    return resultExpr;
  }
  
#ifdef SUPPRESS_ENUM_CONVERSION_WARNING
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wenum-enum-conversion"
#endif
  
  static SEXP run(SEXP samplerExpr, SEXP numIterExpr, SEXP isWarmupExpr, SEXP resultsTypeExpr)
  {
    Sampler* samplerPtr = static_cast<Sampler*>(R_ExternalPtrAddr(samplerExpr));
    if (samplerPtr == NULL) Rf_error("run called on NULL external pointer");
    Sampler& sampler(*samplerPtr);
    
    int numIter    = rc_getInt(numIterExpr, "num_iter", RC_VALUE | RC_GEQ, 1, RC_END);
    bool isWarmup  = rc_getBool(isWarmupExpr, "warmup", RC_NA | RC_NO, RC_END);
    
    ResultsType resultsType = static_cast<ResultsType>(rc_getInt(resultsTypeExpr, "results_type",
     RC_VALUE | RC_DEFAULT, static_cast<int>(RESULTS_BOTH),
     RC_END));

    SEXP callbackClosure = R_NilValue;
    SEXP yhat_train = R_NilValue;
    SEXP yhat_test = R_NilValue;
    SEXP stan_pars = R_NilValue;
    SEXP callbackResults = R_NilValue;
    size_t callbackResultLength = 0;
    unsigned int protectCount = 0;
    if (sampler.callback != R_NilValue) {
      yhat_train = PROTECT(rc_newReal(sampler.numObservations));
      ++protectCount;
      if (sampler.numTestObservations > 0) {
        yhat_test = PROTECT(rc_newReal(sampler.numTestObservations));
        ++protectCount;
      }
      stan_pars = PROTECT(rc_newReal(sampler.stanSampler->num_pars));
      ++protectCount;
      SEXP stan_par_names = PROTECT(rc_newCharacter(sampler.stanSampler->num_pars));
      ++protectCount;
      size_t pos = 0;
      for (size_t i = 0; i < sampler.stanSampler->sample_names.size(); ++i)
        SET_STRING_ELT(stan_par_names, pos++, Rf_mkChar(sampler.stanSampler->sample_names[i].c_str()));
      for (size_t i = 0; i < sampler.stanSampler->sampler_names.size(); ++i)
        SET_STRING_ELT(stan_par_names, pos++, Rf_mkChar(sampler.stanSampler->sampler_names[i].c_str()));
      for (size_t i = 0; i < sampler.stanSampler->constrained_param_names.size(); ++i)
        SET_STRING_ELT(stan_par_names, pos++, Rf_mkChar(sampler.stanSampler->constrained_param_names[i].c_str()));
      rc_setNames(stan_pars, stan_par_names);
        
      callbackClosure = PROTECT(Rf_lang4(sampler.callback, yhat_train, yhat_test, stan_pars));
      ++protectCount;
    }
    
    stan4bart::IterableBartResults* bartSamples = NULL;
    // allocate storage for results
    
    size_t numStorageSamples = sampler.keepFits ? numIter : 1;

    if (resultsType == RESULTS_BOTH || resultsType == RESULTS_BART)
      bartSamples = new stan4bart::IterableBartResults(
        sampler.numObservations,
        bartFunctions.numPredictors(sampler.bartSampler),
        sampler.numTestObservations, numStorageSamples, sampler.kIsSampled);
    if (resultsType == RESULTS_BOTH || resultsType == RESULTS_STAN)
      sampler.stanSampler->sample_writer.resize(sampler.stanSampler->num_pars, numStorageSamples);
    
    size_t n = sampler.numObservations;
    
    if (sampler.keepTrees)
      bartFunctions.setTreeStorage(sampler.bartSampler, isWarmup ? 0 : 1,
                                   isWarmup ? 0 : static_cast<size_t>(numIter));
    
    if (sampler.verbose > 0)
      Rprintf("starting %s, %d draws, %s\n", isWarmup ? "warmup" : "sampling", numIter,
              resultsType == RESULTS_BOTH ? "both BART and Stan" : (resultsType == RESULTS_BART ? "BART only" : "Stan only"));
    
    for (int iter = 0; iter < numIter; ++iter) {
      if (sampler.refresh > 0 && sampler.verbose > 1 && (iter + 1) % sampler.refresh == 0)
          Rprintf("  iter %.3d / %.3d\n", iter + 1, numIter);
      
      // order of update matters - need to store a parametric components that go with a bart prediction
      // or else when they're added together they won't be consistent with `predict`
      if (resultsType == RESULTS_BOTH || resultsType == RESULTS_STAN) {
        // Rprintf("running stan\n");
        sampler.stanSampler->run(isWarmup);
        
        if (sampler.userOffset == NULL) {
          // Rprintf("getting stan para mean\n");
          sampler.stanSampler->getParametricMean(*sampler.stanModel, sampler.bartOffset);
        } else {
          // The user offset can be used to replace parts of the model for debugging purposes.
          // We only add in the remaining parts.
          switch (sampler.offsetType) {
            case OFFSET_DEFAULT:
            sampler.stanSampler->getParametricMean(*sampler.stanModel, sampler.bartOffset);
            for (size_t j = 0; j < n; ++j) sampler.bartOffset[j] += sampler.userOffset[j];
            break;
            
            case OFFSET_BART:
            sampler.stanSampler->getParametricMean(*sampler.stanModel, sampler.bartOffset);
            break;
            
            case OFFSET_RANEF:
            sampler.stanSampler->getParametricMean(*sampler.stanModel, sampler.bartOffset,
                                                   true, false);
            for (size_t j = 0; j < n; ++j) sampler.bartOffset[j] += sampler.userOffset[j];
            break;
            
            case OFFSET_FIXEF:
            sampler.stanSampler->getParametricMean(*sampler.stanModel, sampler.bartOffset,
                                                   false, true);
            for (size_t j = 0; j < n; ++j) sampler.bartOffset[j] += sampler.userOffset[j];
            break;
            
            case OFFSET_PARAMETRIC:
            // Replaces the whole Stan part.
            std::memcpy(sampler.bartOffset, sampler.userOffset, n * sizeof(double));
            break;
          }
        }
        if (!sampler.responseIsBinary) {
          // Rprintf("getting sigma\n");
          double sigma = sampler.stanSampler->getSigma(*sampler.stanModel);
          bartFunctions.setSigma(sampler.bartSampler, sigma);
        }
        
        /* if (!isWarmup) {
          if (stanSampler->isDivergentTransition()) {
            R_ShowMessage("bad things happened!\n")
          }
        } */
        
        // Rprintf("incrementing sampling\n");
        if (sampler.keepFits)
          sampler.stanSampler->sample_writer.increment();
        
        // Rprintf("setting bart offset\n");
        // this will adjusting the scale every iteration during the first 1/8 of warmup, 
        // ever other iteration during the second 1/8, every fourth iteration during the
        // third, and so forth
        int update_scale_mod = 1 << (8 * iter / numIter);
        bartFunctions.setOffset(sampler.bartSampler, sampler.bartOffset, isWarmup && iter % update_scale_mod == 0);
        // bartFunctions.setOffset(sampler.bartSampler, sampler.bartOffset, isWarmup);
      }
      
      if (resultsType == RESULTS_BOTH || resultsType == RESULTS_BART) {
        
        // Rprintf("sampling bart\n");
        bartFunctions.run(sampler.bartSampler, 0, 1, &bartSamples->current);
        double* trainingSample = const_cast<double*>(bartSamples->current.train);
        
        // bart with an offset will produce predictions that have the offset added;
        // in order to just get the tree predictions, subtract out that offset
        for (size_t j = 0; j < n; ++j)
          trainingSample[j] -= sampler.bartOffset[j];
        
        if (sampler.userOffset != NULL && sampler.offsetType == OFFSET_BART) {
          // Override with user supplied bart offset
          std::memcpy(sampler.stanOffset, sampler.userOffset, n * sizeof(double));
        } else {
          std::memcpy(sampler.stanOffset, trainingSample, n * sizeof(double));
          if (sampler.userOffset != NULL && sampler.offsetType == OFFSET_DEFAULT)
            for (size_t j = 0; j < n; ++j)
              sampler.stanOffset[j] += sampler.userOffset[j];
        }
        
        // Rprintf("setting stan offset\n");
        stan4bart::setStanOffset(*sampler.stanModel, sampler.stanOffset);
        if (sampler.responseIsBinary) {
          // Rprintf("getting latents\n");
          bartFunctions.getLatents(sampler.bartSampler, sampler.bartLatents);
          stan4bart::setResponse(*sampler.stanModel, sampler.bartLatents);
        }
        
        if (sampler.callback != R_NilValue) {
          std::memcpy(REAL(yhat_train), trainingSample, n * sizeof(double));
          if (yhat_test != R_NilValue)
            std::memcpy(REAL(yhat_test), bartSamples->current.test, sampler.numTestObservations * sizeof(double));
          sampler.stanSampler->copyOutParameters(REAL(stan_pars), sampler.keepFits ? -1 : 0);

          SEXP callbackIterResult = PROTECT(Rf_eval(callbackClosure, sampler.callbackEnv));
          bool resultAllocated = false;
          if (callbackIterResult != R_NilValue && rc_getLength(callbackIterResult) > 0) {
            callbackResultLength = rc_getLength(callbackIterResult);
            if (callbackResults == R_NilValue) {
              callbackResults = PROTECT(rc_newReal(callbackResultLength * numIter));
              resultAllocated = true;
              ++protectCount;
              SEXP existingDims = rc_getDims(callbackIterResult);
              if (existingDims != R_NilValue) {
                SEXP dims = PROTECT(rc_newInteger(rc_getLength(existingDims) + 1));
                for (size_t i = 0 ; i < rc_getLength(existingDims); ++i)
                  INTEGER(dims)[i] = INTEGER(existingDims)[i];
                INTEGER(dims)[rc_getLength(existingDims)] = numIter;
                Rf_setAttrib(callbackResults, R_DimSymbol, dims);
                UNPROTECT(1);
              } else {
                rc_setDims(callbackResults, static_cast<int>(callbackResultLength), numIter, -1);
              }
              if (rc_getNames(callbackIterResult) != R_NilValue) {
                SEXP dimNames = PROTECT(rc_newList(rc_getLength(rc_getDims(callbackResults))));
                rc_setDimNames(callbackResults, dimNames);
                UNPROTECT(1);
                SET_VECTOR_ELT(dimNames, 0, rc_getNames(callbackIterResult));
                SET_VECTOR_ELT(dimNames, 1, R_NilValue);
              } else if (rc_getDimNames(callbackIterResult) != R_NilValue) {
                SEXP dimNames = PROTECT(rc_newList(rc_getLength(rc_getDims(callbackResults))));
                rc_setDimNames(callbackResults, dimNames);
                UNPROTECT(1);
                SEXP existingDimNames = rc_getDimNames(callbackIterResult);
                for (size_t i = 0 ; i < rc_getLength(existingDimNames); ++i)
                  SET_VECTOR_ELT(dimNames, i, VECTOR_ELT(existingDimNames, i));
                SET_VECTOR_ELT(dimNames, rc_getLength(existingDimNames), R_NilValue);
                
                SEXP existingDimNamesNames = rc_getNames(existingDimNames);
                if (existingDimNamesNames != R_NilValue) {
                  SEXP dimNamesNames = PROTECT(rc_newCharacter(rc_getLength(rc_getDims(callbackResults))));
                  for (size_t i = 0; i < rc_getLength(existingDimNamesNames); ++i) {
                    SET_STRING_ELT(dimNamesNames, i, STRING_ELT(existingDimNamesNames, i));
                  }
                  SET_STRING_ELT(dimNamesNames, rc_getLength(existingDimNamesNames), PROTECT(Rf_mkChar("iterations")));
                  rc_setNames(dimNames, dimNamesNames);
                  UNPROTECT(2);
                }
              }
            }
            std::memcpy(REAL(callbackResults) + iter * callbackResultLength,
                        const_cast<const double*>(REAL(callbackIterResult)),
                        callbackResultLength * sizeof(double));
          }
          // We aren't balanced if it is our first result, so we have to unprotect
          // something that isn't the most recent on the stack.
          if (resultAllocated)
            UNPROTECT_PTR(callbackIterResult);
          else
            UNPROTECT(1);
        }

        // Rprintf("increment bart pointers\n");
        if (sampler.keepFits)
          bartSamples->incrementPointers();
      }
    }
    
    SEXP resultExpr = R_NilValue;
    // Rprintf("writing results\n");
    if (sampler.keepFits) {
      if (bartSamples != NULL) bartSamples->resetPointers();
      
      size_t resultsLength = (resultsType == RESULTS_BOTH ? 2 : 1) + 
                             (sampler.callback != R_NilValue ? 1 : 0);
      resultExpr = PROTECT(rc_newList(resultsLength));
      ++protectCount;
      int pos = 0;
      if (resultsType == RESULTS_BOTH || resultsType == RESULTS_STAN)
        SET_VECTOR_ELT(resultExpr, pos++, createStanResultsExpr(sampler.stanSampler->sample_writer));
      if (resultsType == RESULTS_BOTH || resultsType == RESULTS_BART)
        SET_VECTOR_ELT(resultExpr, pos++, stan4bart::createBartResultsExpr(*bartSamples));
      if (sampler.callback != R_NilValue)
        SET_VECTOR_ELT(resultExpr, pos, callbackResults);
      
      SEXP namesExpr = PROTECT(rc_newCharacter(rc_getLength(resultExpr)));
      pos = 0;
      if (resultsType == RESULTS_BOTH || resultsType == RESULTS_STAN)
        SET_STRING_ELT(namesExpr, pos++, Rf_mkChar("stan"));
      if (resultsType == RESULTS_BOTH || resultsType == RESULTS_BART)
        SET_STRING_ELT(namesExpr, pos++, Rf_mkChar("bart"));
      if (sampler.callback != R_NilValue)
        SET_STRING_ELT(namesExpr, pos, Rf_mkChar("callback"));
      
      rc_setNames(resultExpr, namesExpr);
      UNPROTECT(1);
    } else {
      resultExpr = PROTECT(rc_newList(1));
      ++protectCount;
      SET_VECTOR_ELT(resultExpr, 0, callbackResults);
      SEXP namesExpr = PROTECT(rc_newCharacter(rc_getLength(resultExpr)));
      SET_STRING_ELT(namesExpr, 0, Rf_mkChar("callback"));

      rc_setNames(resultExpr, namesExpr);
      UNPROTECT(1);
    }
    
    delete bartSamples;

    UNPROTECT(protectCount);
    
    return(resultExpr);
  }
  
#ifdef SUPPRESS_ENUM_CONVERSION_WARNING
#  pragma clang diagnostic pop
#endif
  
  static SEXP printInitialSummary(SEXP samplerExpr) {
    Sampler* samplerPtr = static_cast<Sampler*>(R_ExternalPtrAddr(samplerExpr));
    if (samplerPtr == NULL) Rf_error("printInitialSummary called on NULL external pointer");
    Sampler& sampler(*samplerPtr);
    
    // the bart initial summary prints during creation under a verbose
    // control; only the stan side and the user offset remain here
    Rprintf("stan control:\n");
    printStanControl(sampler.stanControl);
    Rprintf("stan model:\n");
    stan4bart::printStanModel(sampler.stanModel);
    
    if (sampler.userOffset != NULL) {
      Rprintf("\nuser offset: %f", sampler.userOffset[0]);
      for (size_t i = 1; i < (sampler.numObservations < 5 ? sampler.numObservations : 5); ++i)
        Rprintf(", %f", sampler.userOffset[i]);
      if (sampler.numObservations > 5) Rprintf("...");
      Rprintf("\n");
      if (sampler.offsetType != OFFSET_DEFAULT) Rprintf("  type: %s\n", userOffsetTypeNames[sampler.offsetType]);
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

#ifdef SUPPRESS_ENUM_CONVERSION_WARNING
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wenum-enum-conversion"
#endif

void initializeSamplerFromExpression(Sampler& sampler, SEXP commonControlExpr)
{
  sampler.defaultWarmup = rc_getInt(rc_getListElement(commonControlExpr, "warmup"), "warmup",
    RC_VALUE | RC_DEFAULT, 1000,
    RC_END);
  sampler.defaultIter = rc_getInt(rc_getListElement(commonControlExpr, "iter"), "iter",
    RC_VALUE | RC_DEFAULT, 2000,
    RC_END);
  sampler.verbose = rc_getInt(rc_getListElement(commonControlExpr, "verbose"), "verbose",
    RC_VALUE | RC_DEFAULT, 0,
    RC_END);
  sampler.refresh = rc_getInt(rc_getListElement(commonControlExpr, "refresh"), "refresh",
    RC_VALUE | RC_GEQ, 0,
    RC_NA | RC_YES, RC_END);
  sampler.responseIsBinary = rc_getBool(rc_getListElement(commonControlExpr, "is_binary"), "responseIsBinary",
    RC_NA | RC_NO, RC_END);
  
  SEXP offsetExpr = rc_getListElement(commonControlExpr, "offset");
  sampler.userOffset = offsetExpr == R_NilValue || rc_getLength(offsetExpr) == 0 || !Rf_isReal(offsetExpr) ? NULL : REAL(offsetExpr);
  
  sampler.offsetType = static_cast<UserOffsetType>(
    rc_getInt(rc_getListElement(commonControlExpr, "offset_type"), "offset_type",
     RC_VALUE | RC_DEFAULT, static_cast<int>(OFFSET_DEFAULT),
     RC_END));
  
  if (sampler.refresh == R_NaInt)
    sampler.refresh = 200;

  sampler.keepFits = rc_getBool(rc_getListElement(commonControlExpr, "keep_fits"), "keepFits",
    RC_NA | RC_NO, RC_END);
  sampler.callback = rc_getListElement(commonControlExpr, "callback");
  if (sampler.callback != R_NilValue && !Rf_isFunction(sampler.callback))
    Rf_error("callback must be a function or NULL");
  sampler.callbackEnv = rc_getListElement(commonControlExpr, "callbackEnv");
  if (sampler.callbackEnv != R_NilValue && !Rf_isEnvironment(sampler.callbackEnv))
    Rf_error("callbackEnv must be an environment or NULL");
}

#ifdef SUPPRESS_ENUM_CONVERSION_WARNING
#  pragma clang diagnostic pop
#  undef SUPPRESS_ENUM_CONVERSION_WARNING
#endif

}

#if __cplusplus >= 202002L
#  include <bit>
#else
#  include <cstring>

namespace std {

#  if __cplusplus >= 201103L
#    include <type_traits>

template <class To, class From>
typename std::enable_if<
  sizeof(To) == sizeof(From) &&
  std::is_trivially_copyable<From>::value &&
  std::is_trivially_copyable<To>::value,
  To>::type
// constexpr support needs compiler magic
bit_cast(const From& src) noexcept
{
  static_assert(std::is_trivially_constructible<To>::value,
    "This implementation additionally requires destination type to be trivially constructible");

  To dst;
  memcpy(&dst, &src, sizeof(To));
  return dst;
}

#  else

// We are only using this to cast function pointers, which are trivially copiable.
// is_trivially_copyable is compiler specific and isn't worth trying to drop in
// an implementation.
template <class To, class From>
To
bit_cast(const From& src)
{
  To dst;
  memcpy(&dst, &src, sizeof(To));
  return dst;
}

#  endif

}

#endif

// this unusual set of declarations solves a rather obscure warning on Solaris
typedef void* (*C_voidPtrFunction)(void);
extern "C" typedef C_voidPtrFunction (*C_voidPtrFunctionLookup)(const char* _namespace, const char* name);

namespace {
  
  void lookupBARTFunctions()
  {
    bartFunctions.apiVersion = std::bit_cast<int (*)(void)>(R_GetCCallable("dbarts", "dbarts_apiVersion"));
    if (bartFunctions.apiVersion() != DBARTS_C_API_VERSION)
      Rf_error("stan4bart was built against dbarts C API version %d but the installed dbarts provides %d; reinstall stan4bart",
               DBARTS_C_API_VERSION, bartFunctions.apiVersion());
    
    bartFunctions.create               = std::bit_cast<dbarts_sampler* (*)(SEXP, SEXP, SEXP, const char*)>(R_GetCCallable("dbarts", "dbarts_sampler_create"));
    bartFunctions.destroy              = std::bit_cast<void (*)(dbarts_sampler*)>(R_GetCCallable("dbarts", "dbarts_sampler_destroy"));
    bartFunctions.run                  = std::bit_cast<void (*)(dbarts_sampler*, std::size_t, std::size_t, dbarts_results*)>(R_GetCCallable("dbarts", "dbarts_sampler_run"));
    bartFunctions.sampleTreesFromPrior = std::bit_cast<void (*)(dbarts_sampler*)>(R_GetCCallable("dbarts", "dbarts_sampler_sampleTreesFromPrior"));
    bartFunctions.setOffset            = std::bit_cast<void (*)(dbarts_sampler*, const double*, int)>(R_GetCCallable("dbarts", "dbarts_sampler_setOffset"));
    bartFunctions.setSigma             = std::bit_cast<void (*)(dbarts_sampler*, double)>(R_GetCCallable("dbarts", "dbarts_sampler_setSigma"));
    bartFunctions.getLatents           = std::bit_cast<int (*)(const dbarts_sampler*, double*)>(R_GetCCallable("dbarts", "dbarts_sampler_getLatents"));
    bartFunctions.predict              = std::bit_cast<void (*)(dbarts_sampler*, const double*, std::size_t, const double*, double*)>(R_GetCCallable("dbarts", "dbarts_sampler_predict"));
    bartFunctions.setTreeStorage       = std::bit_cast<void (*)(dbarts_sampler*, int, std::size_t)>(R_GetCCallable("dbarts", "dbarts_sampler_setTreeStorage"));
    bartFunctions.getTrees             = std::bit_cast<SEXP (*)(dbarts_sampler*, const std::size_t*, std::size_t, const std::size_t*, std::size_t, const std::size_t*, std::size_t, int)>(R_GetCCallable("dbarts", "dbarts_sampler_getTrees"));
    bartFunctions.printTrees           = std::bit_cast<void (*)(dbarts_sampler*, const std::size_t*, std::size_t, const std::size_t*, std::size_t, const std::size_t*, std::size_t)>(R_GetCCallable("dbarts", "dbarts_sampler_printTrees"));
    bartFunctions.storeState           = std::bit_cast<SEXP (*)(dbarts_sampler*)>(R_GetCCallable("dbarts", "dbarts_sampler_storeState"));
    bartFunctions.setState             = std::bit_cast<void (*)(dbarts_sampler*, SEXP)>(R_GetCCallable("dbarts", "dbarts_sampler_setState"));
    bartFunctions.setVerbose           = std::bit_cast<void (*)(dbarts_sampler*, int, std::uint32_t)>(R_GetCCallable("dbarts", "dbarts_sampler_setVerbose"));
    bartFunctions.numObservations      = std::bit_cast<std::size_t (*)(const dbarts_sampler*)>(R_GetCCallable("dbarts", "dbarts_sampler_numObservations"));
    bartFunctions.numPredictors        = std::bit_cast<std::size_t (*)(const dbarts_sampler*)>(R_GetCCallable("dbarts", "dbarts_sampler_numPredictors"));
    bartFunctions.numTestObservations  = std::bit_cast<std::size_t (*)(const dbarts_sampler*)>(R_GetCCallable("dbarts", "dbarts_sampler_numTestObservations"));
    bartFunctions.numChains            = std::bit_cast<std::size_t (*)(const dbarts_sampler*)>(R_GetCCallable("dbarts", "dbarts_sampler_numChains"));
    bartFunctions.numTrees             = std::bit_cast<std::size_t (*)(const dbarts_sampler*)>(R_GetCCallable("dbarts", "dbarts_sampler_numTrees"));
    bartFunctions.numSavedSamples      = std::bit_cast<std::size_t (*)(const dbarts_sampler*)>(R_GetCCallable("dbarts", "dbarts_sampler_numSavedSamples"));
    bartFunctions.kIsSampled           = std::bit_cast<int (*)(const dbarts_sampler*)>(R_GetCCallable("dbarts", "dbarts_sampler_kIsSampled"));
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

static void storedBARTSamplerFinalizer(SEXP samplerExpr)
{
  StoredBARTSampler* sampler = static_cast<StoredBARTSampler*>(R_ExternalPtrAddr(samplerExpr));
  
  if (sampler == NULL) return;
  
  if (activeStoredBARTSamplers->find(samplerExpr) == activeStoredBARTSamplers->end()) return;
  
  activeStoredBARTSamplers->erase(samplerExpr);
  
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
  
  for (PointerSet::iterator it = activeStoredBARTSamplers->begin(); it != activeStoredBARTSamplers->end(); ) {
    SEXP samplerExpr = *it;
    StoredBARTSampler* sampler = static_cast<StoredBARTSampler*>(R_ExternalPtrAddr(samplerExpr));
    
    delete sampler;
    PointerSet::iterator prev = it;
    ++it;
    activeStoredBARTSamplers->erase(prev);
    R_ClearExternalPtr(samplerExpr);
  }
      
  delete activeStoredBARTSamplers;
  
  return R_NilValue;
}

#define DEF_FUNC(_N_, _F_, _A_) { _N_, std::bit_cast<DL_FUNC>(&_F_), _A_ }

static R_CallMethodDef R_callMethods[] = {
  DEF_FUNC("stan4bart_create", createSampler, 6),
  DEF_FUNC("stan4bart_run", run, 4),
  DEF_FUNC("stan4bart_printInitialSummary", printInitialSummary, 1),
  DEF_FUNC("stan4bart_disengageAdaptation", disengageAdaptation, 1),
  DEF_FUNC("stan4bart_finalize", finalize, 0),
  DEF_FUNC("stan4bart_exportBARTState", exportBARTState, 1),
  DEF_FUNC("stan4bart_createStoredBARTSampler", createStoredBARTSampler, 4),
  DEF_FUNC("stan4bart_predictBART", predictBART, 3),
  DEF_FUNC("stan4bart_getParametricMean", getParametricMean, 1),
  DEF_FUNC("stan4bart_printTrees", printTrees, 4),
  DEF_FUNC("stan4bart_getTrees", getTrees, 5),
  {NULL, NULL, 0}
};

#undef DEF_FUNC

void attribute_visible R_init_stan4bart(DllInfo *info) {
  R_registerRoutines(info, NULL, R_callMethods, NULL, NULL);
  R_useDynamicSymbols(info, static_cast<Rboolean>(FALSE));
  
  lookupBARTFunctions();
  
  activeSamplers = new PointerSet(&compareExternalPointers);
  activeStoredBARTSamplers = new PointerSet(&compareExternalPointers);
}

} // extern "C"

