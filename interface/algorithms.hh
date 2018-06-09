#ifndef ALGORITHMS_HH
#define ALGORITHMS_HH

#include <math.h>

template <typename input_iterator, typename mean_t>
mean_t compute_mean (input_iterator first, input_iterator last, mean_t start_value) {

  mean_t count = 0;
  for (input_iterator i = first; i!=last; i++) {
    count += *i;
    start_value += *i * (i-first);
  }
  return start_value/count;
}

template <typename input_iterator, typename rms_t, typename mean_t>
rms_t compute_rms (input_iterator first, input_iterator last, rms_t start_value, mean_t mean) {

  mean_t count = 0;
  for (input_iterator i = first; i!=last; i++) {
    count += *i;
    start_value += *i * ::pow((i-first) - mean, 2);
  }

  return ::sqrt(start_value/count);
}

#endif
