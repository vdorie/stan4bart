#ifndef DOUBLE_WRITER_HPP
#define DOUBLE_WRITER_HPP

#include <cstddef> // size_t
#include <cstring> // memcpy
#include <sstream> // ostringstream
#include <vector>

#include <stan/callbacks/writer.hpp>

// TODO: remove context and names
namespace stan4bart {

struct double_writer : public stan::callbacks::writer {
  size_t num_pars;
  size_t num_samples;
  const char* context;
  double* x_base;
  double* x_curr;
  std::vector<std::string> names;
  
  double_writer(const char* context) :
    num_pars(0), num_samples(0),
    context(context),
    x_base(NULL),
    x_curr(NULL)
  {
  }
  
  double_writer(const char* context, size_t num_pars, size_t num_samples) :
    num_pars(num_pars), num_samples(num_samples),
    context(context),
    x_base(new double[num_pars * num_samples]),
    x_curr(x_base)
  {
  }
   
  ~double_writer() {
    delete [] x_base;
  }
  
  // copies in current sample, if applicable
  void resize(size_t num_pars, size_t num_samples) {
    const double* x_curr_old = x_curr;
    const double* x_base_old = x_base;
    
    this->num_pars = num_pars;
    this->num_samples = num_samples;
    x_base = new double[num_pars * num_samples];
    x_curr = x_base;
    
    if (x_base_old != NULL) {
      std::memcpy(x_curr, x_curr_old, num_pars * sizeof(double));
      delete [] x_base_old;
    }
  }
   
  void operator()(const std::vector<std::string>& names) {
    this->names = names;
  }
  
  void printValues() const {
    Rprintf("%s values:\n", context);
    for (size_t i = 0; i < num_pars; ++i)
      Rprintf("  %-14s: %.8f\n", names[i].c_str(), x_curr[i]);
    Rprintf("\n");
  }
  
  void operator()(const std::vector<double>& state) {
    if (num_pars != state.size()) {
      std::ostringstream errorMessage;
      errorMessage << "double writer size mismatch: " << num_pars << " allocated, " << state.size() << " requested";
      throw std::out_of_range(errorMessage.str());
    }
    std::memcpy(x_curr, state.data(), num_pars * sizeof(double));
  }
  void increment() {
    x_curr += num_pars;
  }
  void reset() {
    x_curr = x_base;
  }
};

}

#endif // DOUBLE_WRITER_HPP

