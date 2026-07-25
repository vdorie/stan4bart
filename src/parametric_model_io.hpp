#ifndef PARAMETRIC_MODEL_IO_HPP
#define PARAMETRIC_MODEL_IO_HPP

/// \file parametric_model_io.hpp
/// \brief Build a ParametricModel from the marshaled Stan data list (the same
///        named list continuous_model consumes). Shared by the gradient gate
///        surface (logdensity_export.cpp) and the WALNUTS sampler wiring
///        (walnuts_sampler.cpp) so both derive len_z_T/len_rho/delta identically.

#include <ext/Rinternals.h>  // SEXP

#include "parametric_model.hpp"

namespace stan4bart {

ParametricModel buildParametricModel(SEXP dataExpr);

}  // namespace stan4bart

#endif  // PARAMETRIC_MODEL_IO_HPP
