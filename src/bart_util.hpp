#ifndef BART_UTIL_HPP
#define BART_UTIL_HPP

#include <cstddef>
#include <cstdint>
#include <vector>

#include <ext/Rinternals.h>

#include <dbarts/dbarts.h>

namespace stan4bart {

  /// Per-iteration landing buffers for dbarts_sampler_run: sized for a full
  /// run's draws, with the current-iteration pointers advanced by the
  /// caller between single-draw runs. The flat-C-API replacement for the
  /// retired dbarts::Results subclass; single chain.
  struct IterableBartResults {
    std::size_t numObservations, numPredictors, numTestObservations;
    std::size_t numSamples;
    bool kIsSampled;

    std::vector<double> sigmaSamples, trainingSamples, testSamples, kSamples;
    std::vector<std::uint32_t> variableCountSamples;

    std::size_t position;
    dbarts_results current = {};

    IterableBartResults(std::size_t numObservations_,
                        std::size_t numPredictors_,
                        std::size_t numTestObservations_,
                        std::size_t numSamples_, bool kIsSampled_)
      : numObservations(numObservations_), numPredictors(numPredictors_),
        numTestObservations(numTestObservations_), numSamples(numSamples_),
        kIsSampled(kIsSampled_), sigmaSamples(numSamples_),
        trainingSamples(numObservations_ * numSamples_),
        testSamples(numTestObservations_ * numSamples_),
        kSamples(kIsSampled_ ? numSamples_ : 0),
        variableCountSamples(numPredictors_ * numSamples_), position(0)
    {
      // the versioned-struct contract: dbarts_sampler_run fills only the
      // fields structSize says the caller's dbarts_results carries; left
      // unset, every output field is skipped and the train buffers stay zero
      current.structSize = sizeof(dbarts_results);
      setCurrentPointers();
    }

    /// Aims the run outputs at this iteration's slice.
    void setCurrentPointers() {
      current.sigma = sigmaSamples.data() + position;
      current.train = trainingSamples.data() + position * numObservations;
      current.test = numTestObservations > 0
        ? testSamples.data() + position * numTestObservations : NULL;
      current.varcount = variableCountSamples.data() +
                         position * numPredictors;
      current.k = kIsSampled ? kSamples.data() + position : NULL;
      current.varprobs = NULL;
    }

    void incrementPointers() {
      ++position;
      setCurrentPointers();
    }

    void resetPointers() {
      position = 0;
      setCurrentPointers();
    }
  };

  /// A named list of (sigma, train, test, varcount[, k]) in single-chain
  /// layout, as the retired dbarts::Results-based version produced.
  SEXP createBartResultsExpr(const IterableBartResults& results);
}

#endif // BART_UTIL_HPP
