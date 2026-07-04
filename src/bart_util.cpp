#include "bart_util.hpp"

#include <cstring>

#include <rc/util.h>

namespace stan4bart {

SEXP createBartResultsExpr(const IterableBartResults& results)
{
  std::size_t numSamples = results.numSamples;

  SEXP resultExpr = PROTECT(rc_newList(results.kIsSampled ? 5 : 4));
  SET_VECTOR_ELT(resultExpr, 0, rc_newReal(rc_asRLength(numSamples)));
  SET_VECTOR_ELT(resultExpr, 1, rc_newReal(rc_asRLength(
    results.numObservations * numSamples)));
  if (results.numTestObservations > 0)
    SET_VECTOR_ELT(resultExpr, 2, rc_newReal(rc_asRLength(
      results.numTestObservations * numSamples)));
  else
    SET_VECTOR_ELT(resultExpr, 2, R_NilValue);
  SET_VECTOR_ELT(resultExpr, 3, rc_newInteger(rc_asRLength(
    results.numPredictors * numSamples)));
  if (results.kIsSampled)
    SET_VECTOR_ELT(resultExpr, 4, rc_newReal(rc_asRLength(numSamples)));

  SEXP sigmaSamples = VECTOR_ELT(resultExpr, 0);
  std::memcpy(REAL(sigmaSamples), results.sigmaSamples.data(),
              numSamples * sizeof(double));

  SEXP trainingSamples = VECTOR_ELT(resultExpr, 1);
  rc_setDims(trainingSamples, static_cast<int>(results.numObservations),
             static_cast<int>(numSamples), -1);
  std::memcpy(REAL(trainingSamples), results.trainingSamples.data(),
              results.numObservations * numSamples * sizeof(double));

  if (results.numTestObservations > 0) {
    SEXP testSamples = VECTOR_ELT(resultExpr, 2);
    rc_setDims(testSamples, static_cast<int>(results.numTestObservations),
               static_cast<int>(numSamples), -1);
    std::memcpy(REAL(testSamples), results.testSamples.data(),
                results.numTestObservations * numSamples * sizeof(double));
  }

  SEXP variableCountSamples = VECTOR_ELT(resultExpr, 3);
  rc_setDims(variableCountSamples, static_cast<int>(results.numPredictors),
             static_cast<int>(numSamples), -1);
  int* variableCountStorage = INTEGER(variableCountSamples);
  std::size_t length = results.numPredictors * numSamples;
  // these need to be down-sized from unsigned to int
  for (std::size_t i = 0; i < length; ++i)
    variableCountStorage[i] =
      static_cast<int>(results.variableCountSamples[i]);

  if (results.kIsSampled) {
    SEXP kSamples = VECTOR_ELT(resultExpr, 4);
    std::memcpy(REAL(kSamples), results.kSamples.data(),
                numSamples * sizeof(double));
  }

  // create result storage and make it user friendly
  SEXP namesExpr;

  rc_setNames(resultExpr, namesExpr = rc_newCharacter(
    results.kIsSampled ? 5 : 4));
  SET_STRING_ELT(namesExpr, 0, Rf_mkChar("sigma"));
  SET_STRING_ELT(namesExpr, 1, Rf_mkChar("train"));
  SET_STRING_ELT(namesExpr, 2, Rf_mkChar("test"));
  SET_STRING_ELT(namesExpr, 3, Rf_mkChar("varcount"));
  if (results.kIsSampled)
    SET_STRING_ELT(namesExpr, 4, Rf_mkChar("k"));

  UNPROTECT(1);

  return resultExpr;
}

}
