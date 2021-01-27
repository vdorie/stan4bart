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

#include "rstan/io/r_ostream.hpp"

#include "bart_util.hpp"
#include "stan_sampler.hpp"

namespace {
  // stuff to handle external pointers
  typedef bool(*ExternalPointerComparator)(const SEXP&lhs, const SEXP& rhs);
  typedef std::set<SEXP, ExternalPointerComparator> PointerSet;
  
  PointerSet* activeSamplers;
  PointerSet* activeStoredBARTSamplers;
  
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
    SEXP (*createStateExpression)(const dbarts::BARTFit* fit);
    void (*initializeState)(dbarts::BARTFit* fit, SEXP stateExpr);
    void (*setControl)(dbarts::BARTFit* fit, const dbarts::Control* control);

    void (*runSamplerWithResults)(dbarts::BARTFit* fit, std::size_t numBurnIn, dbarts::Results* results);
    
    void (*predict)(const dbarts::BARTFit* fit, const double* x_test, std::size_t numTestObservations, const double* testOffset, double* result);
    void (*setResponse)(dbarts::BARTFit* fit, const double* response);
    void (*setOffset)(dbarts::BARTFit* fit, const double* offset, bool updateState);
    void (*setSigma)(dbarts::BARTFit* fit, const double* sigma);
    void (*sampleTreesFromPrior)(dbarts::BARTFit* fit);
    void (*printInitialSummary)(const dbarts::BARTFit* fit);
    
    void (*getLatentVariables)(const dbarts::BARTFit*, double*);
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
    dbarts::Control control;
    dbarts::Data data;
    dbarts::Model model;
    dbarts::BARTFit* fit;
    
    StoredBARTSampler() : model(false), fit(NULL) { }
    ~StoredBARTSampler() {
      if (fit != NULL) {
        bartFunctions.invalidateFit(fit);
        ::operator delete(fit);
        fit = NULL;
      }
      bartFunctions.invalidateModel(&model);
      bartFunctions.invalidateData(&data);
    }
  };
  
  struct Sampler {
    int defaultWarmup;
    int defaultIter;
    int verbose;
    int refresh;
    bool is_binary;
    const double* userOffset;
    UserOffsetType offsetType;
    
    model_continuous_namespace::model_continuous* stanModel;
    stan4bart::StanControl stanControl;
    stan4bart::StanSampler* stanSampler;
    
    dbarts::Control bartControl;
    dbarts::Data bartData;
    dbarts::Model bartModel;
    dbarts::BARTFit* bartSampler;
    bool keepTrees;
    
    double* bartOffset;
    double* stanOffset;
    double* bartLatents;
        
    Sampler() :
      stanModel(NULL), stanSampler(NULL), bartModel(false), bartSampler(NULL), bartOffset(NULL), stanOffset(NULL), bartLatents(NULL)
    {
    }
    ~Sampler() {
      delete [] bartLatents;
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
  static void storedBARTSamplerFinalizer(SEXP samplerExpr);
  
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
    if (sampler.stanControl.thin == R_NaInt) {
      sampler.stanControl.thin = (2000 - sampler.defaultWarmup) / 1000;
      if (sampler.stanControl.thin < 1) sampler.stanControl.thin = 1;
    }
    
    int chain_id = 1;
    sampler.stanSampler = new stan4bart::StanSampler(*sampler.stanModel, sampler.stanControl, chain_id, sampler.defaultWarmup);
    
    bartFunctions.initializeControl(&sampler.bartControl, bartControlExpr);
    sampler.keepTrees = sampler.bartControl.keepTrees;
    sampler.bartControl.keepTrees = false;
    if (sampler.keepTrees) {
      sampler.bartControl.defaultNumSamples = sampler.defaultIter - sampler.defaultWarmup;
      sampler.bartControl.defaultNumBurnIn  = sampler.defaultWarmup;
    }
    sampler.bartControl.responseIsBinary = sampler.is_binary;
    
    bartFunctions.initializeData(&sampler.bartData, bartDataExpr);
    bartFunctions.initializeModel(&sampler.bartModel, bartModelExpr, &sampler.bartControl);
    
    sampler.bartSampler = static_cast<dbarts::BARTFit*>(::operator new (sizeof(dbarts::BARTFit)));
    bartFunctions.initializeFit(sampler.bartSampler, &sampler.bartControl, &sampler.bartModel, &sampler.bartData);
    
    size_t n = sampler.bartData.numObservations;
    sampler.bartOffset = new double[n];
    sampler.stanOffset = new double[n];
    if (sampler.is_binary)
      sampler.bartLatents = new double[n];
    
    if (bart_offset_init != NULL) {
      std::memcpy(sampler.bartOffset, bart_offset_init, n * sizeof(double));
    } else {
      for (size_t i = 0; i < n; ++i) sampler.bartOffset[i] = 0.0;
    }
    if (sampler.userOffset != NULL && sampler.offsetType != OFFSET_BART) {
      for (size_t i = 0; i < n; ++i) sampler.bartOffset[i] += sampler.userOffset[i];
    }
    
    bartFunctions.setOffset(sampler.bartSampler, sampler.bartOffset, true);
    bartFunctions.setSigma(sampler.bartSampler, &sigma_init);

    bartFunctions.sampleTreesFromPrior(sampler.bartSampler);
    
    // draw once before running
    dbarts::Control bartControl = sampler.bartSampler->control;
    bool oldVerbose = bartControl.verbose;
    bartControl.verbose = false;
    bartFunctions.setControl(sampler.bartSampler, &bartControl);
    dbarts::Results* first_draw = new dbarts::Results(n,
                                                      sampler.bartSampler->data.numPredictors,
                                                      sampler.bartSampler->data.numTestObservations,
                                                      1, bartControl.numChains,
                                                      sampler.bartSampler->model.kPrior != NULL);
    bartFunctions.runSamplerWithResults(sampler.bartSampler, 0, first_draw);
    
    for (size_t j = 0; j < n; ++j) first_draw->trainingSamples[j] -= sampler.bartOffset[j];
    std::memcpy(sampler.stanOffset, const_cast<const double*>(first_draw->trainingSamples), n * sizeof(double));
    if (sampler.userOffset != NULL) {
      if (sampler.offsetType != OFFSET_BART)
        for (size_t j = 0; j < n; ++j)
          sampler.stanOffset[j] += sampler.userOffset[j];
      else if (sampler.offsetType == OFFSET_BART)
        std::memcpy(sampler.stanOffset, sampler.userOffset, n * sizeof(double));
    }
    stan4bart::setStanOffset(*sampler.stanModel, sampler.stanOffset);
    if (sampler.is_binary) {
      bartFunctions.getLatentVariables(sampler.bartSampler, sampler.bartLatents);
      stan4bart::setResponse(*sampler.stanModel, sampler.bartLatents);
    }
    
    delete first_draw;
    bartControl.verbose = oldVerbose;
    bartFunctions.setControl(sampler.bartSampler, &bartControl);
    
    
    SEXP result = PROTECT(R_MakeExternalPtr(samplerPtr.get(), R_NilValue, R_NilValue));
    samplerPtr.release();
    R_RegisterCFinalizerEx(result, samplerFinalizer, static_cast<Rboolean>(FALSE));
    

    activeSamplers->insert(result);

    UNPROTECT(1);
    
    return result;
  }
  
  static SEXP getBARTDataRange(SEXP samplerExpr)
  {
    Sampler* samplerPtr = static_cast<Sampler*>(R_ExternalPtrAddr(samplerExpr));
    if (samplerPtr == NULL) Rf_error("getBARTDataRange called on NULL external pointer");
    Sampler& sampler(*samplerPtr);
    
    SEXP resultExpr = PROTECT(rc_newReal(2));
    double* result = REAL(resultExpr);
    result[0] = sampler.bartSampler->sharedScratch.dataScale.min;
    result[1] = sampler.bartSampler->sharedScratch.dataScale.max;
    
    UNPROTECT(1);
    
    return resultExpr;
  }
  
  static SEXP getParametricMean(SEXP samplerExpr)
  {
    Sampler* samplerPtr = static_cast<Sampler*>(R_ExternalPtrAddr(samplerExpr));
    if (samplerPtr == NULL) Rf_error("getParametricMean called on NULL external pointer");
    Sampler& sampler(*samplerPtr);
    
    sampler.stanSampler->sample_writer.decrement();
    SEXP result = PROTECT(rc_newReal(sampler.bartData.numObservations));
    
    stan4bart::getParametricMean(*sampler.stanSampler, *sampler.stanModel, REAL(result));
    sampler.stanSampler->sample_writer.increment();
    
    UNPROTECT(1);
    
    return result;
  }
  
  static SEXP predictBART(SEXP storedBARTSamplerExpr, SEXP x_testExpr, SEXP offset_testExpr)
  {
    StoredBARTSampler* samplerPtr = static_cast<StoredBARTSampler*>(R_ExternalPtrAddr(storedBARTSamplerExpr));
    if (samplerPtr == NULL) Rf_error("predictBART called on NULL external pointer");
    StoredBARTSampler& sampler(*samplerPtr);
    
    const dbarts::Control& control(sampler.control);
    // const dbarts::Data& data(sampler.data);
    // const dbarts::Model& model(sampler.model);
    const dbarts::BARTFit* fit(sampler.fit);
        
    if (Rf_isNull(x_testExpr)) return R_NilValue;
    
    if (!Rf_isReal(x_testExpr)) Rf_error("x.test must be of type real");
    
    rc_assertDimConstraints(x_testExpr, "dimensions of x_test", RC_LENGTH | RC_EQ, rc_asRLength(2),
                            RC_NA,
                            RC_VALUE | RC_EQ, static_cast<int>(fit->data.numPredictors),
                            RC_END);
    int* dims = INTEGER(Rf_getAttrib(x_testExpr, R_DimSymbol));
    
    size_t numSamples = control.keepTrees ? fit->currentNumSamples : 1;
    size_t numTestObservations = static_cast<size_t>(dims[0]);
    
    double* testOffset = NULL;
    if (!Rf_isNull(offset_testExpr)) {
      if (!Rf_isReal(offset_testExpr)) Rf_error("offset.test must be of type real");
      if (rc_getLength(offset_testExpr) != 1 || !ISNA(REAL(offset_testExpr)[0])) {
        if (rc_getLength(offset_testExpr) != numTestObservations) Rf_error("length of offset.test must equal number of rows in x.test");
        testOffset = REAL(offset_testExpr);
      }
    }
    
    SEXP result = PROTECT(Rf_allocVector(REALSXP, numTestObservations * numSamples * control.numChains));
    if (control.keepTrees) {
      if (fit->control.numChains <= 1)
        rc_setDims(result, static_cast<int>(numTestObservations), static_cast<int>(numSamples), -1);
      else
        rc_setDims(result, static_cast<int>(numTestObservations), static_cast<int>(numSamples), static_cast<int>(control.numChains), -1);
    } else {
      if (fit->control.numChains > 1)
        rc_setDims(result, static_cast<int>(numTestObservations), static_cast<int>(control.numChains), -1);
    }
    
    bartFunctions.predict(fit, REAL(x_testExpr), numTestObservations, testOffset, REAL(result));
    
    UNPROTECT(1);
    
    return result;
  }
  
  static SEXP exportBARTState(SEXP samplerExpr)
  {
    Sampler* samplerPtr = static_cast<Sampler*>(R_ExternalPtrAddr(samplerExpr));
    if (samplerPtr == NULL) Rf_error("exportBARTState called on NULL external pointer");
    Sampler& sampler(*samplerPtr);
    
    return bartFunctions.createStateExpression(sampler.bartSampler);
  }
  
  static SEXP createStoredBARTSampler(SEXP controlExpr, SEXP dataExpr, SEXP modelExpr, SEXP stateExpr)
  {
    std::unique_ptr<StoredBARTSampler> samplerPtr(new StoredBARTSampler);
    StoredBARTSampler& sampler(*samplerPtr);
    
    bartFunctions.initializeControl(&sampler.control, controlExpr);
    sampler.control.numChains = XLENGTH(stateExpr);
    sampler.control.keepTrees = true;
    bartFunctions.initializeData(&sampler.data, dataExpr);
    bartFunctions.initializeModel(&sampler.model, modelExpr, &sampler.control);
    
    sampler.fit = static_cast<dbarts::BARTFit*>(::operator new (sizeof(dbarts::BARTFit)));
    bartFunctions.initializeFit(sampler.fit, &sampler.control, &sampler.model, &sampler.data);
    
    bartFunctions.initializeState(sampler.fit, stateExpr);
    sampler.fit->sharedScratch.dataScale.min = -0.5;
    sampler.fit->sharedScratch.dataScale.max = 0.5;
    sampler.fit->sharedScratch.dataScale.range = 1.0;
    
    SEXP result = PROTECT(R_MakeExternalPtr(samplerPtr.get(), R_NilValue, R_NilValue));
    samplerPtr.release();
    R_RegisterCFinalizerEx(result, storedBARTSamplerFinalizer, static_cast<Rboolean>(FALSE));

    activeStoredBARTSamplers->insert(result);

    UNPROTECT(1);
    
    return result;
  }
  
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
    
    stan4bart::IterableBartResults* bartSamples = NULL;
    // allocate storage for results
    if (resultsType == RESULTS_BOTH || resultsType == RESULTS_BART)
      bartSamples = new stan4bart::IterableBartResults(
        sampler.bartData.numObservations, sampler.bartData.numPredictors, sampler.bartData.numTestObservations, numIter, 1 /* num chains */, false /* TODO: binary */);
    if (resultsType == RESULTS_BOTH || resultsType == RESULTS_STAN)
      sampler.stanSampler->sample_writer.resize(sampler.stanSampler->num_pars, numIter);
    
    size_t n = sampler.bartData.numObservations;
    
    if (sampler.keepTrees) {
      if (!isWarmup) {
        sampler.bartControl.keepTrees = true;
      } else {
        sampler.bartControl.keepTrees = false;
      }
      bartFunctions.setControl(sampler.bartSampler, &sampler.bartControl);
    }
    
    if (sampler.verbose > 0)
      Rprintf("starting %s, %d draws, %s:\n", isWarmup ? "warmup" : "sampling", numIter,
              resultsType == RESULTS_BOTH ? "both BART and Stan" : (resultsType == RESULTS_BART ? "BART only" : "Stan only"));
    
    
    // TODO: add in capacity for fitting against true ranef and true fixef
    for (int iter = 0; iter < numIter; ++iter) {
      if (sampler.refresh > 0 && sampler.verbose > 1 && (iter + 1) % sampler.refresh == 0)
          Rprintf("  iter %.3d / %.3d\n", iter + 1, numIter);
      
      // order of update matters - need to store a parametric components that go with a bart prediction
      // or else when they're added together they won't be consistent with `predict`
      if (resultsType == RESULTS_BOTH || resultsType == RESULTS_STAN) {
        sampler.stanSampler->run(isWarmup);
        
        if (sampler.userOffset == NULL) {
          stan4bart::getParametricMean(*sampler.stanSampler, *sampler.stanModel, sampler.bartOffset);
        } else {
          switch (sampler.offsetType) {
            case OFFSET_DEFAULT:
            stan4bart::getParametricMean(*sampler.stanSampler, *sampler.stanModel, sampler.bartOffset);
            for (size_t j = 0; j < n; ++j) sampler.bartOffset[j] += sampler.userOffset[j];
            break;
            
            case OFFSET_BART:
            stan4bart::getParametricMean(*sampler.stanSampler, *sampler.stanModel, sampler.bartOffset);
            break;
            
            case OFFSET_RANEF:
            stan4bart::getParametricMean(*sampler.stanSampler, *sampler.stanModel, sampler.bartOffset,
                                         true, false);
            for (size_t j = 0; j < n; ++j) sampler.bartOffset[j] += sampler.userOffset[j];
            break;
            
            case OFFSET_FIXEF:
            stan4bart::getParametricMean(*sampler.stanSampler, *sampler.stanModel, sampler.bartOffset,
                                         false, true);
            for (size_t j = 0; j < n; ++j) sampler.bartOffset[j] += sampler.userOffset[j];
            break;
            
            case OFFSET_PARAMETRIC:
            std::memcpy(sampler.bartOffset, sampler.userOffset, n * sizeof(double));
            break;
          }
        }
        if (sampler.is_binary) {
          double sigma = getSigma(*sampler.stanSampler, *sampler.stanModel);
          bartFunctions.setSigma(sampler.bartSampler, &sigma);
        }
        
        sampler.stanSampler->sample_writer.increment();
        
        bartFunctions.setOffset(sampler.bartSampler, sampler.bartOffset, isWarmup);
      }
      
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
          else if (sampler.offsetType == OFFSET_BART)
            std::memcpy(sampler.stanOffset, sampler.userOffset, n * sizeof(double));
        }
        stan4bart::setStanOffset(*sampler.stanModel, sampler.stanOffset);
        if (sampler.is_binary) {
          bartFunctions.getLatentVariables(sampler.bartSampler, sampler.bartLatents);
          stan4bart::setResponse(*sampler.stanModel, sampler.bartLatents);
        }
        
        bartSamples->incrementPointers();
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
    
    Rprintf("stan control:\n");
    printStanControl(sampler.stanControl);
    Rprintf("bart init:\n");
    bartFunctions.printInitialSummary(sampler.bartSampler);
    
    if (sampler.userOffset != NULL) {
      Rprintf("\nuser offset: %f", sampler.userOffset[0]);
      for (size_t i = 1; i < (sampler.bartData.numObservations < 5 ? sampler.bartData.numObservations : 5); ++i)
        Rprintf(", %f", sampler.userOffset[i]);
      if (sampler.bartData.numObservations > 5) Rprintf("...");
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
  sampler.is_binary = rc_getBool(rc_getListElement(commonControlExpr, "is_binary"), "is_binary",
    RC_NA | RC_NO, RC_END);
  
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
    bartFunctions.setControl            = reinterpret_cast<void (*)(dbarts::BARTFit*, const dbarts::Control*)>(R_GetCCallable("dbarts", "setControl"));
    
    bartFunctions.predict               = reinterpret_cast<void (*)(const dbarts::BARTFit*, const double*, std::size_t, const double*, double*)>(R_GetCCallable("dbarts", "predict"));
    bartFunctions.setResponse           = reinterpret_cast<void (*)(dbarts::BARTFit*, const double*)>(R_GetCCallable("dbarts", "setResponse"));
    bartFunctions.setOffset             = reinterpret_cast<void (*)(dbarts::BARTFit*, const double*, bool)>(R_GetCCallable("dbarts", "setOffset"));
    bartFunctions.setSigma              = reinterpret_cast<void (*)(dbarts::BARTFit*, const double*)>(R_GetCCallable("dbarts", "setSigma"));
    
    bartFunctions.createStateExpression = reinterpret_cast<SEXP (*)(const dbarts::BARTFit*)>(R_GetCCallable("dbarts", "createStateExpression"));
    bartFunctions.initializeState       = reinterpret_cast<void (*)(dbarts::BARTFit*, SEXP)>(R_GetCCallable("dbarts", "initializeState"));
    
    bartFunctions.runSamplerWithResults = reinterpret_cast<void (*)(dbarts::BARTFit*, std::size_t, dbarts::Results*)>(R_GetCCallable("dbarts", "runSamplerWithResults"));
    bartFunctions.sampleTreesFromPrior  = reinterpret_cast<void (*)(dbarts::BARTFit*)>(R_GetCCallable("dbarts", "sampleTreesFromPrior"));
    bartFunctions.printInitialSummary   = reinterpret_cast<void (*)(const dbarts::BARTFit*)>(R_GetCCallable("dbarts", "printInitialSummary"));
    
    bartFunctions.getLatentVariables    = reinterpret_cast<void (*)(const dbarts::BARTFit*, double*)>(R_GetCCallable("dbarts", "storeLatents"));
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

#define DEF_FUNC(_N_, _F_, _A_) { _N_, reinterpret_cast<DL_FUNC>(&_F_), _A_ }

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
  DEF_FUNC("stan4bart_getBARTDataRange", getBARTDataRange, 1),
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
