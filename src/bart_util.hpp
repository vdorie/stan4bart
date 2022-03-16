#ifndef BART_UTIL_HPP
#define BART_UTIL_HPP

#include <ext/Rinternals.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/results.hpp>

namespace stan4bart {

  SEXP createBartResultsExpr(const dbarts::BARTFit& fit, const dbarts::Results& results);

  
  struct IterableBartResults : dbarts::Results {
    double* base_sigmaSamples;
    double* base_trainingSamples;
    double* base_testSamples;
    std::uint32_t* base_variableCountSamples;
    double* base_kSamples;
    
    size_t base_numSamples;
    
    IterableBartResults(std::size_t numObservations, std::size_t numPredictors,
                        std::size_t numTestObservations, std::size_t numSamples, std::size_t numChains,
                        bool kIsModeled) :
     Results(numObservations, numPredictors, numTestObservations, numSamples, numChains, kIsModeled),
     base_sigmaSamples(sigmaSamples),
     base_trainingSamples(trainingSamples),
     base_testSamples(testSamples),
     base_variableCountSamples(variableCountSamples),
     base_kSamples(kSamples),
     base_numSamples(numSamples)
     {
       this->numSamples = 1;
     }
    
     void incrementPointers() {
       sigmaSamples += 1;
       trainingSamples += numObservations;
       if (testSamples != NULL)
         testSamples += numTestObservations;
       variableCountSamples += numPredictors;
       if (kSamples != NULL)
         kSamples += 1;
     }
     
     void resetPointers() {
       kSamples = base_kSamples;
       variableCountSamples = base_variableCountSamples;
       testSamples = base_testSamples;
       trainingSamples = base_trainingSamples;
       sigmaSamples = base_sigmaSamples;
       numSamples = base_numSamples;
     }
     
     ~IterableBartResults() {
       resetPointers();
       base_kSamples = NULL;
       base_variableCountSamples = NULL;
       base_testSamples = NULL;
       base_trainingSamples = NULL;
       base_sigmaSamples = NULL;
     }
  };
}

#endif // BART_UTIL_HPP
