#ifndef BART_UTIL_HPP
#define BART_UTIL_HPP

#include <cstddef>
#include <cstdint>
#include <vector>

#include <ext/Rinternals.h>

#include <dbarts/dbarts.h>

namespace stan4bart {

  /// Per-iteration landing buffers for dbarts_sampler_run, sized for a full
  /// run's draws, with the current-iteration pointers advanced by the caller
  /// between single-draw runs; single chain. The buffers are the R result
  /// list's own vectors (allocateBartResultsExpr), which the caller holds
  /// PROTECTed, so this struct owns nothing: a dbarts entry that raises - an
  /// interrupt, an engine error, a warning turned into one - longjmps past it
  /// and R reclaims the storage.
  struct IterableBartResults {
    std::size_t numObservations, numPredictors, numTestObservations;
    std::size_t numSamples;
    bool kIsSampled;

    double* sigmaSamples;
    double* trainingSamples;
    double* testSamples;
    double* kSamples;
    std::uint32_t* variableCountSamples;

    std::size_t position;
    dbarts_results current;

    /// resultsExpr is allocateBartResultsExpr's list, rooted by the caller
    /// for as long as this struct is used.
    IterableBartResults(SEXP resultsExpr, std::size_t numObservations_,
                        std::size_t numPredictors_,
                        std::size_t numTestObservations_,
                        std::size_t numSamples_, bool kIsSampled_)
      : numObservations(numObservations_), numPredictors(numPredictors_),
        numTestObservations(numTestObservations_), numSamples(numSamples_),
        kIsSampled(kIsSampled_),
        sigmaSamples(REAL(VECTOR_ELT(resultsExpr, 0))),
        trainingSamples(REAL(VECTOR_ELT(resultsExpr, 1))),
        testSamples(numTestObservations_ > 0
                      ? REAL(VECTOR_ELT(resultsExpr, 2)) : NULL),
        kSamples(kIsSampled_ ? REAL(VECTOR_ELT(resultsExpr, 4)) : NULL),
        // an R integer vector is 32 bits wide; a split count never
        // approaches 2^31, so the unsigned write reads back unchanged
        variableCountSamples(reinterpret_cast<std::uint32_t*>(
          INTEGER(VECTOR_ELT(resultsExpr, 3)))),
        position(0), current()
    {
      // the versioned-struct contract: dbarts_sampler_run fills only the
      // fields structSize says the caller's dbarts_results carries; left
      // unset, every output field is skipped and the train buffers stay zero.
      // The trailing pointers stay null from the value-init, so a field this
      // consumer does not ask for is never written.
      current.structSize = sizeof(dbarts_results);
      setCurrentPointers();
    }

    /// Aims the run outputs at this iteration's slice.
    void setCurrentPointers() {
      current.sigma = sigmaSamples + position;
      current.train = trainingSamples + position * numObservations;
      current.test = numTestObservations > 0
        ? testSamples + position * numTestObservations : NULL;
      current.varcount = variableCountSamples + position * numPredictors;
      current.k = kIsSampled ? kSamples + position : NULL;
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

  /// The named list of (sigma, train, test, varcount[, k]) in single-chain
  /// layout that a run's draws land in directly, unprotected on return.
  SEXP allocateBartResultsExpr(std::size_t numObservations,
                               std::size_t numPredictors,
                               std::size_t numTestObservations,
                               std::size_t numSamples, bool kIsSampled);
}

#endif // BART_UTIL_HPP
