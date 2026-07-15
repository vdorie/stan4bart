#ifndef PARAMETRIC_SAMPLER_HPP
#define PARAMETRIC_SAMPLER_HPP

/// \file parametric_sampler.hpp
/// \brief The common surface init.cpp's Gibbs bridge consumes from the
///        parametric conditional's sampler. Both the Stan-backed StanSampler
///        (binary) and the WALNUTS-backed WalnutsSampler (continuous) derive
///        from ParametricSampler, so the family switch is a single runtime
///        pointer and every call site goes through one shape.
///
/// The data members (sample_writer, the three name vectors, num_pars,
/// sample_writer_offset) are shared storage; the methods are virtual. The
/// sample_writer's row layout is Stan-identical for both backends:
///   [sample diag (lp__, accept_stat__)] +
///   [sampler diag (stepsize__, treedepth__, n_leapfrog__, divergent__,
///    energy__)] + [constrained params, at sample_writer_offset].

#include <cstddef>  // size_t
#include <cstring>  // memcpy
#include <string>
#include <vector>

#include "double_writer.hpp"

namespace stan4bart {

struct ParametricSampler {
  // shared draw storage + naming (Stan-identical layout for both backends)
  double_writer sample_writer;
  std::vector<std::string> sample_names;
  std::vector<std::string> sampler_names;
  std::vector<std::string> constrained_param_names;
  size_t sample_writer_offset;
  int num_pars;

  ParametricSampler() : sample_writer("sample"), sample_writer_offset(0), num_pars(0) {}
  virtual ~ParametricSampler() {}

  /// \brief Take one parametric transition (adapting when isWarmup).
  virtual void run(bool isWarmup) = 0;

  /// \brief Freeze adaptation (Stan: disengage_adaptation; WALNUTS: hand off
  ///        AdaptiveWalnuts -> WalnutsSampler).
  virtual void freeze() = 0;

  /// \brief The parametric mean X*beta (+ Z*b) of the current draw, offset-free.
  virtual void getParametricMean(double* result) const = 0;
  virtual void getParametricMean(double* result, bool includeFixed,
                                 bool includeRandom) const = 0;
  /// \brief The residual sd of the current draw (continuous only).
  virtual double getSigma() const = 0;

  /// \brief Refresh the conditioning data between Gibbs transitions.
  virtual void setOffset(const double* offset) = 0;
  virtual void setResponse(const double* y) = 0;

  virtual void setVerbose(int level) = 0;

  /// \brief Copy the full num_pars row at column `offset` from x_curr out
  ///        (shared: both backends store into sample_writer identically).
  void copyOutParameters(double* result, int offset) const {
    std::memcpy(result, const_cast<const double*>(sample_writer.x_curr) + offset * num_pars,
                num_pars * sizeof(double));
  }
};

}  // namespace stan4bart

#endif  // PARAMETRIC_SAMPLER_HPP
