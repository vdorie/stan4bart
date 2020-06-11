#include "stan_util.hpp"

#include <misc/stddef.h>
#include <cstring>

#include <string>

#include <rc/util.h>

#include "double_writer.hpp"

namespace stan4bart {

SEXP createStanResultsExpr(const double_writer& sample_writer)
{
  SEXP resultExpr = PROTECT(rc_newReal(sample_writer.num_pars * sample_writer.num_samples));
  rc_setDims(resultExpr, static_cast<int>(sample_writer.num_pars), static_cast<int>(sample_writer.num_samples), -1);
  std::memcpy(REAL(resultExpr), sample_writer.x_base, sample_writer.num_pars * sample_writer.num_samples * sizeof(double));
  
  SEXP dimnamesExpr = PROTECT(rc_newList(2));
  SET_VECTOR_ELT(dimnamesExpr, 0, rc_newCharacter(sample_writer.num_pars));
  SET_VECTOR_ELT(dimnamesExpr, 1, R_NilValue);
  
  rc_setDimNames(resultExpr, dimnamesExpr);
      
  SEXP namesExpr = VECTOR_ELT(dimnamesExpr, 0);
  for (size_t i = 0; i < sample_writer.num_pars; ++i)
    SET_STRING_ELT(namesExpr, i, Rf_mkChar(sample_writer.names[i].c_str()));
  
  UNPROTECT(2);
  
  return resultExpr;
 }

}

