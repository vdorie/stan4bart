#include "bart_util.hpp"

#include <rc/util.h>

namespace stan4bart {

SEXP allocateBartResultsExpr(std::size_t numObservations,
                             std::size_t numPredictors,
                             std::size_t numTestObservations,
                             std::size_t numSamples, bool kIsSampled)
{
  SEXP resultExpr = PROTECT(rc_newList(kIsSampled ? 5 : 4));
  SET_VECTOR_ELT(resultExpr, 0, rc_newReal(rc_asRLength(numSamples)));
  SET_VECTOR_ELT(resultExpr, 1, rc_newReal(rc_asRLength(
    numObservations * numSamples)));
  if (numTestObservations > 0)
    SET_VECTOR_ELT(resultExpr, 2, rc_newReal(rc_asRLength(
      numTestObservations * numSamples)));
  else
    SET_VECTOR_ELT(resultExpr, 2, R_NilValue);
  SET_VECTOR_ELT(resultExpr, 3, rc_newInteger(rc_asRLength(
    numPredictors * numSamples)));
  if (kIsSampled)
    SET_VECTOR_ELT(resultExpr, 4, rc_newReal(rc_asRLength(numSamples)));

  rc_setDims(VECTOR_ELT(resultExpr, 1), static_cast<int>(numObservations),
             static_cast<int>(numSamples), -1);
  if (numTestObservations > 0)
    rc_setDims(VECTOR_ELT(resultExpr, 2),
               static_cast<int>(numTestObservations),
               static_cast<int>(numSamples), -1);
  rc_setDims(VECTOR_ELT(resultExpr, 3), static_cast<int>(numPredictors),
             static_cast<int>(numSamples), -1);

  // create result storage and make it user friendly
  SEXP namesExpr;

  rc_setNames(resultExpr, namesExpr = rc_newCharacter(kIsSampled ? 5 : 4));
  SET_STRING_ELT(namesExpr, 0, Rf_mkChar("sigma"));
  SET_STRING_ELT(namesExpr, 1, Rf_mkChar("train"));
  SET_STRING_ELT(namesExpr, 2, Rf_mkChar("test"));
  SET_STRING_ELT(namesExpr, 3, Rf_mkChar("varcount"));
  if (kIsSampled)
    SET_STRING_ELT(namesExpr, 4, Rf_mkChar("k"));

  UNPROTECT(1);

  return resultExpr;
}

}
