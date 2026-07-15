#ifndef WALNUTS_SAMPLER_HPP
#define WALNUTS_SAMPLER_HPP

/// \file walnuts_sampler.hpp
/// \brief The WALNUTS-backed parametric sampler for both response families.
///
/// Presents exactly the surface init.cpp's Gibbs bridge consumes (via the
/// ParametricSampler base), driving a WALNUTS AdaptiveWalnuts/WalnutsSampler
/// over the hand-derived ParametricModel target instead of Stan's NUTS. The
/// walnuts + Eigen machinery is hidden behind a PIMPL so init.cpp never sees
/// them and stays a plain R-API translation unit.
///
/// Lifecycle mirrors bairrtt's IrtSampler: construct in the adapting phase,
/// drive one adapting transition per warmup sweep (run(true)), freeze() at
/// disengageAdaptation, then one fixed-tuning transition per sampling sweep
/// (run(false)). offset_/response refreshes update the referenced functor data
/// in place, so the next transition conditions on the moved BART fit.

#include <ext/Rinternals.h>  // SEXP

#include "parametric_sampler.hpp"

namespace stan4bart {

struct WalnutsSampler : public ParametricSampler {
  /// \param dataExpr the marshaled Stan data list (built by the R side).
  /// \param random_seed the per-chain seed (threaded exactly as Stan's was).
  /// \param init_radius unconstrained init half-width (Stan's init_r).
  /// \param num_warmup the outer loop's warmup sweep count.
  WalnutsSampler(SEXP dataExpr, unsigned int random_seed, double init_radius,
                 int num_warmup);
  ~WalnutsSampler() override;

  void run(bool isWarmup) override;
  void freeze() override;

  void getParametricMean(double* result) const override;
  void getParametricMean(double* result, bool includeFixed,
                         bool includeRandom) const override;
  double getSigma() const override;

  void setOffset(const double* offset) override;
  void setResponse(const double* y) override;
  void setVerbose(int level) override;

 private:
  struct Impl;
  Impl* impl_;
};

}  // namespace stan4bart

#endif  // WALNUTS_SAMPLER_HPP
