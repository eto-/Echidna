#ifndef _MINUIT_WRAPPER_H
#define _MINUIT_WRAPPER_H

#include <vector>
#include <string>

class TMinuit;

/**
 * @brief Maximum number of independent variables to permit in a function
 * to be minimized with Minuit
 */
#define MAX_MINUIT_PARS 15

/**
 * @brief Wrapper class for ROOT's Minuit function minimizer implementation.
 */
class MinuitWrapper {
  TMinuit * mn;			    ///< Pointer to the wrapped TMinuit object
  double limits[2][MAX_MINUIT_PARS];///< Cached array of parameter limits
  double values[MAX_MINUIT_PARS];   ///< Cached array of parameter values (either initial settings or fit values)
  double errors[MAX_MINUIT_PARS];   ///< Cached array of estimated errors in parameter values
  std::string * names[MAX_MINUIT_PARS]; ///< Cached array of parameter labels
  bool is_set[MAX_MINUIT_PARS];     ///< Array specifying whether the user has set any properties for each given parameter

  public:
  MinuitWrapper(void (* fcn)(int32_t &, double *, double &, double *, int32_t));
  ///< Create a MinuitWrapper by supplying to it a function with this signature:
  ///< (* func)(int32_t & num_parameters, double * parameter_gradient,
  ///<		double & return_value, double * parameter_values, int32_t flag)
  ~MinuitWrapper();
  ///< Destroy the MinuitWrapper and the wrapped TMinuit object.

  void SetParameter(size_t index, double value);
  ///< Set parameter @a index to @a value.
  void SetParameterName(size_t index, std::string name);
  ///< Set the name of parameter @a index to @a name.
  void SetParLimits(size_t index, double min, double max);
  ///< Restrict parameter @a index to be in the range [@a min, @a max].
  /**
   * @brief Minimize the function supplied in the c'tor according to the given constraints.
   * @param errflag (output) Flag that tells whether there was a minimization error.
   * @param convergence (output) Flag that tells whether the minimization converged.
   */
  void Execute(int32_t & errflag, int32_t & convergence);
  
  double GetParameter(size_t index);
  ///< Get the value of parameter @a index.
  double GetParError(size_t index);
  ///< Get the estimated error on parameter @a (only makes sense after a minimization).
};

#endif
