/**
 * @file BxMath.h
 * @brief Miscellaneous mathematical functions that may be used by all modules
*/
#ifndef _BX_MATH_H
#define _BX_MATH_H

#include <stdint.h>
#include <sys/types.h>
#include <vector>
#include <iterator> // for iterator_traits<T>::value_type
#include <cmath>

namespace BxMath {
  double npe_from_nhits(size_t nhits, size_t npmts);

  // No-op pass-thru functor
  template <typename T> class no_op {
   public:
    typedef T value_type;
    const T & operator() (const T & arg) const { return arg; }
  };

  // It seems to be impossible to supply a default template type or
  // default function argument to a templated function; thus supply two
  // versions of each of these.

  // Templated function to get hit distribution times central moments
  template <typename TIterator, typename GetValueFunctor>
  std::vector<double> get_distribution_moments(TIterator begin, TIterator end,
		      unsigned highest_moment,
		      const GetValueFunctor & get_value);
  template <typename TIterator>
  std::vector<double> get_distribution_moments(TIterator begin, TIterator end,
		      unsigned highest_moment);

  // Templated function to get hit distribution time approximate peak
  template <typename TIterator, typename GetValueFunctor>
  typename GetValueFunctor::value_type
  get_peak(TIterator begin, TIterator end, const GetValueFunctor & get_value);
	
	template <typename TIterator, typename GetValueFunctor>
  typename GetValueFunctor::value_type
	get_peak2(TIterator begin, TIterator end, const GetValueFunctor & get_value);
  
	template <typename TIterator>
  typename std::iterator_traits<TIterator>::value_type
  get_peak(TIterator begin, TIterator end);
  
	template <typename TIterator>
  typename std::iterator_traits<TIterator>::value_type
  get_peak2(TIterator begin, TIterator end);
	
  void DecodeFlag(int flag, bool field[8]);

  // Swap endianness from little-endian (format of raw data files) to host.
  inline void LE_to_32_inplace(uint32_t & val, bool is_little_endian)
    {
      if (is_little_endian) return;
      unsigned char * valbytes = (unsigned char *)&val;
      unsigned char dummy;
      dummy = valbytes[0]; valbytes[0] = valbytes[3]; valbytes[3] = dummy;
      dummy = valbytes[1]; valbytes[1] = valbytes[2]; valbytes[2] = dummy;
    }

  inline void LE_to_16_inplace(uint16_t & val, bool is_little_endian)
    {
      if (is_little_endian) return;
      unsigned char * valbytes = (unsigned char *)&val;
      unsigned char dummy;
      dummy = valbytes[0]; valbytes[0] = valbytes[1]; valbytes[1] = dummy;
    }
      
  inline uint32_t LE_to_32(uint32_t val, bool is_little_endian)
    { LE_to_32_inplace(val, is_little_endian); return val; }
  inline uint16_t LE_to_16(uint16_t val, bool is_little_endian)
    { LE_to_16_inplace(val, is_little_endian); return val; }

  /**
   * @brief Determine the distance between two points in space.
   * @param x1 X position of first point, in meters.
   * @param y1 Y position of first point.
   * @param z1 Z position of first point.
   * @param x2 X position of second point.
   * @param y2 Y position of second point.
   * @param z2 Z position of second point.
   * @return Distance between points, in meters.
   */
  inline double distance(double x1, double y1, double z1,
			 double x2, double y2, double z2)
  {
    double dx = x1 - x2, dy = y1 - y2, dz = z1 - z2;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  }
  
  void get_distances(double x1, double y1, double z1,
		     double x2, double y2, double z2,
		     double & dist_to_pmt, double & dist_in_scintillator,
		     double & dist_in_buffer);
  
  // add more if you need them
};


// Template functions have to be defined in the header file, otherwise
// it causes a linker error unless the specializations are specifically
// declared in the same translation unit as the main definition.  So here we go:

#include "mach4/BxMath.icc"

#endif
