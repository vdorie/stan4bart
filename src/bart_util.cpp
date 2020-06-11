#include "bart_util.hpp"

#include <misc/stddef.h>
#include <cstring>

#include <rc/util.h>

#include <dbarts/control.hpp>
#include <dbarts/data.hpp>

namespace stan4bart {

SEXP createBartResultsExpr(const dbarts::BARTFit& fit, const dbarts::Results& results)
{
  int protectCount = 0;
  
  SEXP resultExpr = PROTECT(rc_newList(results.kSamples == NULL ? 4 : 5));
  ++protectCount;
  SET_VECTOR_ELT(resultExpr, 0, rc_newReal(rc_asRLength(results.getNumSigmaSamples())));
  SET_VECTOR_ELT(resultExpr, 1, rc_newReal(rc_asRLength(results.getNumTrainingSamples())));
  if (fit.data.numTestObservations > 0)
    SET_VECTOR_ELT(resultExpr, 2, rc_newReal(rc_asRLength(results.getNumTestSamples())));
  else
    SET_VECTOR_ELT(resultExpr, 2, R_NilValue);
  SET_VECTOR_ELT(resultExpr, 3, rc_newInteger(rc_asRLength(results.getNumVariableCountSamples())));
  if (results.kSamples != NULL)
    SET_VECTOR_ELT(resultExpr, 4, rc_newReal(rc_asRLength(results.getNumSigmaSamples())));
  
  SEXP sigmaSamples = VECTOR_ELT(resultExpr, 0);
  if (fit.control.numChains > 1)
    rc_setDims(sigmaSamples, static_cast<int>(results.numSamples), static_cast<int>(fit.control.numChains), -1);
  std::memcpy(REAL(sigmaSamples), const_cast<const double*>(results.sigmaSamples), results.getNumSigmaSamples() * sizeof(double));
  
  SEXP trainingSamples = VECTOR_ELT(resultExpr, 1);
  if (fit.control.numChains <= 1)
    rc_setDims(trainingSamples, static_cast<int>(results.numObservations), static_cast<int>(results.numSamples), -1);
  else
    rc_setDims(trainingSamples, static_cast<int>(results.numObservations), static_cast<int>(results.numSamples), static_cast<int>(fit.control.numChains), -1);
  std::memcpy(REAL(trainingSamples), const_cast<const double*>(results.trainingSamples), results.getNumTrainingSamples() * sizeof(double));
  
  if (fit.data.numTestObservations > 0) {
    SEXP testSamples = VECTOR_ELT(resultExpr, 2);
    if (fit.control.numChains <= 1)
      rc_setDims(testSamples, static_cast<int>(results.numTestObservations), static_cast<int>(results.numSamples), -1);
    else
      rc_setDims(testSamples, static_cast<int>(results.numTestObservations), static_cast<int>(results.numSamples), static_cast<int>(fit.control.numChains), -1);
    std::memcpy(REAL(testSamples), const_cast<const double*>(results.testSamples), results.getNumTestSamples() * sizeof(double));
  }
  
  SEXP variableCountSamples = VECTOR_ELT(resultExpr, 3);
  if (fit.control.numChains <= 1)
    rc_setDims(variableCountSamples, static_cast<int>(results.numPredictors), static_cast<int>(results.numSamples), -1);
  else
    rc_setDims(variableCountSamples, static_cast<int>(results.numPredictors), static_cast<int>(results.numSamples), static_cast<int>(fit.control.numChains), -1);
  int* variableCountStorage = INTEGER(variableCountSamples);
  size_t length = results.getNumVariableCountSamples();
  // these likely need to be down-sized from 64 to 32 bits
  for (size_t i = 0; i < length; ++i) variableCountStorage[i] = static_cast<int>(results.variableCountSamples[i]);
  
  if (results.kSamples != NULL) {
    SEXP kSamples = VECTOR_ELT(resultExpr, 4);
    if (fit.control.numChains > 1)
      rc_setDims(kSamples, static_cast<int>(results.numSamples), static_cast<int>(fit.control.numChains), -1);
    std::memcpy(REAL(kSamples), const_cast<const double*>(results.kSamples), results.getNumSigmaSamples() * sizeof(double));
  }
      
  // create result storage and make it user friendly
  SEXP namesExpr;
  
  rc_setNames(resultExpr, namesExpr = rc_newCharacter(results.kSamples == NULL ? 4 : 5));
  SET_STRING_ELT(namesExpr, 0, Rf_mkChar("sigma"));
  SET_STRING_ELT(namesExpr, 1, Rf_mkChar("train"));
  SET_STRING_ELT(namesExpr, 2, Rf_mkChar("test"));
  SET_STRING_ELT(namesExpr, 3, Rf_mkChar("varcount"));
  if (results.kSamples != NULL)
    SET_STRING_ELT(namesExpr, 4, Rf_mkChar("k"));
  
  UNPROTECT(protectCount);
    
  return resultExpr;
}

}

