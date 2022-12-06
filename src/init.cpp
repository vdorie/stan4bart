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
    void (*initializeModel)(dbarts::Model* model, SEXP modelExpr, const dbarts::Control* control, const dbarts::Data* data);
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
    void (*printTrees)(const dbarts::BARTFit*, const std::size_t*, std::size_t, const std::size_t*, std::size_t,
                       const std::size_t*, std::size_t);
    dbarts::FlattenedTrees* (*getTrees)(const dbarts::BARTFit*, const std::size_t*, std::size_t, const std::size_t*,
                                        std::size_t, const std::size_t*, std::size_t);
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
    bool responseIsBinary;
    const double* userOffset;
    UserOffsetType offsetType;
    
    continuous_model_namespace::continuous_model* stanModel;
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
    
    bool keepFits;
    SEXP callback;
    SEXP callbackEnv;
    
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
    
    bartFunctions.initializeControl(&sampler.bartControl, bartControlExpr);
    sampler.keepTrees = sampler.bartControl.keepTrees;
    sampler.bartControl.keepTrees = false;
    if (sampler.keepTrees) {
      sampler.bartControl.defaultNumSamples = sampler.defaultIter - sampler.defaultWarmup;
      sampler.bartControl.defaultNumBurnIn  = sampler.defaultWarmup;
    }
    sampler.bartControl.responseIsBinary = sampler.responseIsBinary;
    
    bartFunctions.initializeData(&sampler.bartData, bartDataExpr);
    bartFunctions.initializeModel(&sampler.bartModel, bartModelExpr, &sampler.bartControl, &sampler.bartData);
    
    sampler.bartSampler = static_cast<dbarts::BARTFit*>(::operator new (sizeof(dbarts::BARTFit)));
    bartFunctions.initializeFit(sampler.bartSampler, &sampler.bartControl, &sampler.bartModel, &sampler.bartData);
    
    size_t n = sampler.bartData.numObservations;
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
      bartFunctions.setSigma(sampler.bartSampler, &sigma_init);

    GetRNGstate();
   
    bartFunctions.sampleTreesFromPrior(sampler.bartSampler);
    
    // draw once before running
    dbarts::Control bartControl = sampler.bartSampler->control;
    bool oldVerbose = bartControl.verbose;
    bartControl.verbose = false;
    bartFunctions.setControl(sampler.bartSampler, &bartControl);
    dbarts::Results* firstDraw = new dbarts::Results(n,
                                                     sampler.bartSampler->data.numPredictors,
                                                     sampler.bartSampler->data.numTestObservations,
                                                     1, bartControl.numChains,
                                                     !sampler.bartSampler->model.kPrior->isFixed);
    bartFunctions.runSamplerWithResults(sampler.bartSampler, 0, firstDraw);
    
    for (size_t j = 0; j < n; ++j) firstDraw->trainingSamples[j] -= sampler.bartOffset[j];
    
    if (sampler.userOffset != NULL && sampler.offsetType == OFFSET_BART) {
      // Override with user supplied bart offset
      std::memcpy(sampler.stanOffset, sampler.userOffset, n * sizeof(double));
    } else {
      std::memcpy(sampler.stanOffset, const_cast<const double*>(firstDraw->trainingSamples), n * sizeof(double));
      if (sampler.userOffset != NULL && sampler.offsetType == OFFSET_DEFAULT)
        for (size_t j = 0; j < n; ++j)
          sampler.stanOffset[j] += sampler.userOffset[j];
    }
    
    stan4bart::setStanOffset(*sampler.stanModel, sampler.stanOffset);
    if (sampler.responseIsBinary) {
      bartFunctions.getLatentVariables(sampler.bartSampler, sampler.bartLatents);
      stan4bart::setResponse(*sampler.stanModel, sampler.bartLatents);
    }
    
    delete firstDraw;
    
    bartControl.verbose = oldVerbose;
    bartFunctions.setControl(sampler.bartSampler, &bartControl);
    
    PutRNGstate();
    
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
  
#ifdef SUPPRESS_ENUM_CONVERSION_WARNING
#  pragma clang diagnostic pop
#endif
  
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
    bartFunctions.initializeModel(&sampler.model, modelExpr, &sampler.control, &sampler.data);
    
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

  static SEXP printTrees(SEXP storedBARTSamplerExpr, SEXP chainIndicesExpr, SEXP sampleIndicesExpr, SEXP treeIndicesExpr)
  //static SEXP printTrees(SEXP samplerExpr, SEXP chainIndicesExpr, SEXP sampleIndicesExpr, SEXP treeIndicesExpr)
  {
    /* Sampler* samplerPtr = static_cast<Sampler*>(R_ExternalPtrAddr(samplerExpr));
    if (samplerPtr == NULL) Rf_error("printTrees called on NULL external pointer");
    Sampler& sampler(*samplerPtr); 

    const dbarts::BARTFit* fit(sampler.bartSampler);
    const dbarts::Control& control(sampler.bartControl); */
    
    StoredBARTSampler* samplerPtr = static_cast<StoredBARTSampler*>(R_ExternalPtrAddr(storedBARTSamplerExpr));
    if (samplerPtr == NULL) Rf_error("getTrees called on NULL external pointer");
    StoredBARTSampler& sampler(*samplerPtr);
    
    const dbarts::BARTFit* fit(sampler.fit);
    const dbarts::Control& control(sampler.control);
    
    size_t numChains  = control.numChains;
    size_t numSamples = control.keepTrees ? fit->currentNumSamples : 0;
    size_t numTrees   = control.numTrees;

    size_t numChainIndices  = Rf_isNull(chainIndicesExpr)  ? numChains  : rc_getLength(chainIndicesExpr);
    size_t numSampleIndices = Rf_isNull(sampleIndicesExpr) ? numSamples : rc_getLength(sampleIndicesExpr);
    size_t numTreeIndices   = Rf_isNull(treeIndicesExpr)   ? numTrees   : rc_getLength(treeIndicesExpr);

    if (numChainIndices > numChains)
      Rf_error("%lu chains specified but only %lu in sampler", numChainIndices, numChains);
    if (numSampleIndices > numSamples)
      Rf_error("%lu samples specified but only %lu in sampler", numSampleIndices, numSamples);
    if (numTreeIndices > numTrees)
      Rf_error("%lu trees specified but only %lu in sampler", numTreeIndices, numTrees);
    
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
  
  static SEXP getTrees(SEXP storedBARTSamplerExpr, SEXP chainIndicesExpr, SEXP sampleIndicesExpr, SEXP treeIndicesExpr)
  //static SEXP getTrees(SEXP samplerExpr, SEXP chainIndicesExpr, SEXP sampleIndicesExpr, SEXP treeIndicesExpr)
  {
    /*Sampler* samplerPtr = static_cast<Sampler*>(R_ExternalPtrAddr(samplerExpr));
    if (samplerPtr == NULL) Rf_error("getTrees called on NULL external pointer");
    Sampler& sampler(*samplerPtr);

    const dbarts::BARTFit* fit(sampler.bartSampler);
    const dbarts::Control& control(sampler.bartControl);
    */

    StoredBARTSampler* samplerPtr = static_cast<StoredBARTSampler*>(R_ExternalPtrAddr(storedBARTSamplerExpr));
    if (samplerPtr == NULL) Rf_error("getTrees called on NULL external pointer");
    StoredBARTSampler& sampler(*samplerPtr);
    
    const dbarts::Control& control(sampler.control);
    const dbarts::BARTFit* fit(sampler.fit);
     
    size_t numChains  = control.numChains;
    size_t numSamples = control.keepTrees ? fit->currentNumSamples : 0;
    size_t numTrees   = control.numTrees;

    size_t numChainIndices  = Rf_isNull(chainIndicesExpr)  ? numChains  : rc_getLength(chainIndicesExpr);
    size_t numSampleIndices = Rf_isNull(sampleIndicesExpr) ? numSamples : rc_getLength(sampleIndicesExpr);
    size_t numTreeIndices   = Rf_isNull(treeIndicesExpr)   ? numTrees   : rc_getLength(treeIndicesExpr);
    
    if (numChainIndices > numChains)
      Rf_error("%lu chains specified but only %lu in sampler", numChainIndices, numChains);
    if (numSampleIndices > numSamples)
      Rf_error("%lu samples specified but only %lu in sampler", numSampleIndices, numSamples);
    if (numTreeIndices > numTrees)
      Rf_error("%lu trees specified but only %lu in sampler", numTreeIndices, numTrees);
    
    size_t* chainIndices  = new size_t[numChainIndices];
    size_t* sampleIndices = control.keepTrees ? new size_t[numSamples] : NULL;
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
        
    dbarts::FlattenedTrees* flattenedTreesPtr = bartFunctions.getTrees(fit, chainIndices, numChainIndices, sampleIndices, numSampleIndices, treeIndices, numTreeIndices);
    
    delete [] treeIndices;
    delete [] sampleIndices;
    delete [] chainIndices;

    dbarts::FlattenedTrees& flattenedTrees(*flattenedTreesPtr);
    
    R_xlen_t numCols = 4 + (numChains > 1 ? 1 : 0) + (control.keepTrees ? 1 : 0);
    SEXP resultExpr = PROTECT(rc_newList(numCols));
        
    SEXP classExpr = PROTECT(rc_newCharacter(1));
    SET_STRING_ELT(classExpr, 0, Rf_mkChar("data.frame"));
    Rf_setAttrib(resultExpr, R_ClassSymbol, classExpr);
    UNPROTECT(1);
    
    SEXP resultRowNamesExpr;
    rc_allocateInSlot2(resultRowNamesExpr, resultExpr, R_RowNamesSymbol, STRSXP, flattenedTrees.totalNumNodes);
    
    SEXP resultNamesExpr;
    rc_allocateInSlot2(resultNamesExpr, resultExpr, R_NamesSymbol, STRSXP, numCols);
    
    int* chainNumber = NULL;
    int* sampleNumber = NULL;
    int* treeNumber, *numObservations, *variable;
    double* value;
        
    R_xlen_t colNum = 0;
    if (numChains > 1) {
      SET_VECTOR_ELT(resultExpr, colNum, PROTECT(rc_newInteger(flattenedTrees.totalNumNodes)));
      SET_STRING_ELT(resultNamesExpr, colNum, PROTECT(Rf_mkChar("chain")));
      UNPROTECT(2);
      chainNumber = INTEGER(VECTOR_ELT(resultExpr, colNum));
      ++colNum;
    }
    if (control.keepTrees) {
      SET_VECTOR_ELT(resultExpr, colNum, PROTECT(rc_newInteger(flattenedTrees.totalNumNodes)));
      SET_STRING_ELT(resultNamesExpr, colNum, PROTECT(Rf_mkChar("sample")));
      UNPROTECT(2);
      sampleNumber = INTEGER(VECTOR_ELT(resultExpr, colNum));
      ++colNum;
    }
    SET_VECTOR_ELT(resultExpr, colNum, PROTECT(rc_newInteger(flattenedTrees.totalNumNodes)));
    SET_STRING_ELT(resultNamesExpr, colNum, PROTECT(Rf_mkChar("tree")));
    treeNumber = INTEGER(VECTOR_ELT(resultExpr, colNum));
    ++colNum;
    SET_VECTOR_ELT(resultExpr, colNum, PROTECT(rc_newInteger(flattenedTrees.totalNumNodes)));
    SET_STRING_ELT(resultNamesExpr, colNum, PROTECT(Rf_mkChar("n")));
    numObservations = INTEGER(VECTOR_ELT(resultExpr, colNum));
    ++colNum;
    SET_VECTOR_ELT(resultExpr, colNum, PROTECT(rc_newInteger(flattenedTrees.totalNumNodes)));
    SET_STRING_ELT(resultNamesExpr, colNum, PROTECT(Rf_mkChar("var")));
    variable = INTEGER(VECTOR_ELT(resultExpr, colNum));
    ++colNum;
    SET_VECTOR_ELT(resultExpr, colNum, PROTECT(rc_newReal(flattenedTrees.totalNumNodes)));
    SET_STRING_ELT(resultNamesExpr, colNum, PROTECT(Rf_mkChar("value")));
    value = REAL(VECTOR_ELT(resultExpr, colNum));
    UNPROTECT(8);
    
    size_t numDigits = 1;
    size_t temp = flattenedTrees.totalNumNodes;
    while (temp >= 10) {
      temp /= 10;
      ++numDigits;
    }
    char* buffer = new char[numDigits + 1];
    for (size_t i = 0; i < flattenedTrees.totalNumNodes; ++i) {
      if (chainNumber != NULL)
        chainNumber[i] = static_cast<int>(flattenedTrees.chainNumber[i] + 1);
      if (sampleNumber != NULL)
        sampleNumber[i] = static_cast<int>(flattenedTrees.sampleNumber[i] + 1);
      treeNumber[i] = static_cast<int>(flattenedTrees.treeNumber[i] + 1);
      numObservations[i] = static_cast<int>(flattenedTrees.numObservations[i]);
      int variable_i = static_cast<int>(flattenedTrees.variable[i]);
      variable[i] = variable_i >= 0 ? variable_i + 1 : variable_i;
      value[i] = flattenedTrees.value[i];
#if defined(__MINGW32__) && __cplusplus < 201112L
#  ifdef _WIN64
      std::sprintf(buffer, "%lu", static_cast<unsigned long>(i + 1));
#  else
      std::sprintf(buffer, "%u", i + 1);
#  endif
#else
      std::sprintf(buffer, "%zu", i + 1);
#endif
      SET_STRING_ELT(resultRowNamesExpr, i, PROTECT(Rf_mkChar(buffer)));
      UNPROTECT(1);
    }
    
    delete [] buffer;
    // delete flattenedTreesPtr;
    delete [] flattenedTrees.value;
    delete [] flattenedTrees.variable;
    delete [] flattenedTrees.numObservations;
    delete [] flattenedTrees.treeNumber;
    delete [] flattenedTrees.sampleNumber;
    delete [] flattenedTrees.chainNumber;
    ::operator delete(flattenedTreesPtr);
    
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

    SEXP callbackClosure;
    SEXP yhat_train;
    SEXP yhat_test = R_NilValue;
    SEXP stan_pars;
    SEXP callbackResults = R_NilValue;
    size_t callbackResultLength = 0;
    unsigned int protectCount = 0;
    if (sampler.callback != R_NilValue) {
      yhat_train = PROTECT(rc_newReal(sampler.bartData.numObservations));
      ++protectCount;
      if (sampler.bartData.numTestObservations > 0L) {
        yhat_test = PROTECT(rc_newReal(sampler.bartData.numTestObservations));
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
        sampler.bartData.numObservations, sampler.bartData.numPredictors,
        sampler.bartData.numTestObservations, numStorageSamples,
        1 /* num chains */, !sampler.bartModel.kPrior->isFixed);
    if (resultsType == RESULTS_BOTH || resultsType == RESULTS_STAN)
      sampler.stanSampler->sample_writer.resize(sampler.stanSampler->num_pars, numStorageSamples);
    
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
      Rprintf("starting %s, %d draws, %s\n", isWarmup ? "warmup" : "sampling", numIter,
              resultsType == RESULTS_BOTH ? "both BART and Stan" : (resultsType == RESULTS_BART ? "BART only" : "Stan only"));
    
    GetRNGstate();
    
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
          bartFunctions.setSigma(sampler.bartSampler, &sigma);
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
        bartFunctions.runSamplerWithResults(sampler.bartSampler, 0, bartSamples);
        
        // bart with an offset will produce predictions that have the offset added;
        // in order to just get the tree predictions, subtract out that offset
        for (size_t j = 0; j < n; ++j)
          bartSamples->trainingSamples[j] -= sampler.bartOffset[j];
        
        if (sampler.userOffset != NULL && sampler.offsetType == OFFSET_BART) {
          // Override with user supplied bart offset
          std::memcpy(sampler.stanOffset, sampler.userOffset, n * sizeof(double));
        } else {
          std::memcpy(sampler.stanOffset, const_cast<const double*>(bartSamples->trainingSamples), n * sizeof(double));
          if (sampler.userOffset != NULL && sampler.offsetType == OFFSET_DEFAULT)
            for (size_t j = 0; j < n; ++j)
              sampler.stanOffset[j] += sampler.userOffset[j];
        }
        
        // Rprintf("setting stan offset\n");
        stan4bart::setStanOffset(*sampler.stanModel, sampler.stanOffset);
        if (sampler.responseIsBinary) {
          // Rprintf("getting latents\n");
          bartFunctions.getLatentVariables(sampler.bartSampler, sampler.bartLatents);
          stan4bart::setResponse(*sampler.stanModel, sampler.bartLatents);
        }
        
        if (sampler.callback != R_NilValue) {
          std::memcpy(REAL(yhat_train), const_cast<const double*>(bartSamples->trainingSamples), n * sizeof(double));
          if (yhat_test != R_NilValue)
            std::memcpy(REAL(yhat_test), const_cast<const double*>(bartSamples->testSamples), sampler.bartData.numTestObservations * sizeof(double));
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
    
    PutRNGstate();
    
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
        SET_VECTOR_ELT(resultExpr, pos++, stan4bart::createBartResultsExpr(*sampler.bartSampler, *bartSamples));
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
    
    Rprintf("stan control:\n");
    printStanControl(sampler.stanControl);
    Rprintf("stan model:\n");
    stan4bart::printStanModel(sampler.stanModel);
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
  std::memcpy(&dst, &src, sizeof(To));
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
  std::memcpy(&dst, &src, sizeof(To));
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
    bartFunctions.initializeFit         = std::bit_cast<void (*)(dbarts::BARTFit*, dbarts::Control*, dbarts::Model*, dbarts::Data*)>(R_GetCCallable("dbarts", "initializeFit"));
    bartFunctions.invalidateFit         = std::bit_cast<void (*)(dbarts::BARTFit*)>(R_GetCCallable("dbarts", "invalidateFit"));
    bartFunctions.initializeControl     = std::bit_cast<void (*)(dbarts::Control*, SEXP)>(R_GetCCallable("dbarts", "initializeControl"));
    bartFunctions.initializeData        = std::bit_cast<void (*)(dbarts::Data*, SEXP)>(R_GetCCallable("dbarts", "initializeData"));
    bartFunctions.invalidateData        = std::bit_cast<void (*)(dbarts::Data*)>(R_GetCCallable("dbarts", "invalidateData"));
    bartFunctions.initializeModel       = std::bit_cast<void (*)(
      dbarts::Model*,
      SEXP,
      const dbarts::Control*,
      const dbarts::Data*)>(R_GetCCallable("dbarts", "initializeModel"));
    bartFunctions.invalidateModel       = std::bit_cast<void (*)(dbarts::Model*)>(R_GetCCallable("dbarts", "invalidateModel"));
    bartFunctions.setControl            = std::bit_cast<void (*)(dbarts::BARTFit*, const dbarts::Control*)>(R_GetCCallable("dbarts", "setControl"));
    
    bartFunctions.predict               = std::bit_cast<void (*)(const dbarts::BARTFit*, const double*, std::size_t, const double*, double*)>(R_GetCCallable("dbarts", "predict"));
    bartFunctions.setResponse           = std::bit_cast<void (*)(dbarts::BARTFit*, const double*)>(R_GetCCallable("dbarts", "setResponse"));
    bartFunctions.setOffset             = std::bit_cast<void (*)(dbarts::BARTFit*, const double*, bool)>(R_GetCCallable("dbarts", "setOffset"));
    bartFunctions.setSigma              = std::bit_cast<void (*)(dbarts::BARTFit*, const double*)>(R_GetCCallable("dbarts", "setSigma"));
    
    bartFunctions.createStateExpression = std::bit_cast<SEXP (*)(const dbarts::BARTFit*)>(R_GetCCallable("dbarts", "createStateExpression"));
    bartFunctions.initializeState       = std::bit_cast<void (*)(dbarts::BARTFit*, SEXP)>(R_GetCCallable("dbarts", "initializeState"));
    
    bartFunctions.runSamplerWithResults = std::bit_cast<void (*)(dbarts::BARTFit*, std::size_t, dbarts::Results*)>(R_GetCCallable("dbarts", "runSamplerWithResults"));
    bartFunctions.sampleTreesFromPrior  = std::bit_cast<void (*)(dbarts::BARTFit*)>(R_GetCCallable("dbarts", "sampleTreesFromPrior"));
    bartFunctions.printInitialSummary   = std::bit_cast<void (*)(const dbarts::BARTFit*)>(R_GetCCallable("dbarts", "printInitialSummary"));
    
    bartFunctions.getLatentVariables    = std::bit_cast<void (*)(const dbarts::BARTFit*, double*)>(R_GetCCallable("dbarts", "storeLatents"));
    bartFunctions.printTrees            = std::bit_cast<void (*)(const dbarts::BARTFit*, const std::size_t*, std::size_t,
                                                                 const std::size_t*, std::size_t,
                                                                 const std::size_t*, std::size_t)>(R_GetCCallable("dbarts", "printTrees"));
    bartFunctions.getTrees              = std::bit_cast<dbarts::FlattenedTrees* (*)(
      const dbarts::BARTFit*, const std::size_t*, std::size_t,
      const std::size_t*, std::size_t, const std::size_t*, std::size_t)>(R_GetCCallable("dbarts", "getTrees"));
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
  DEF_FUNC("stan4bart_getBARTDataRange", getBARTDataRange, 1),
  DEF_FUNC("stan4bart_printTrees", printTrees, 4),
  DEF_FUNC("stan4bart_getTrees", getTrees, 4),
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

