#ifndef Distribution1D_h__
#define Distribution1D_h__

#include <algorithm>
#include <vector>
#include <stdint.h>
#include <assert.h>
#include "devices/device_singleray/default.h"

using std::vector;
using std::max;

namespace embree{

  template<typename T, typename U>
  struct __UToT__
  {
    static inline T convert(const U &u) { return static_cast<T>(u); }
  };

  template<typename T>
  struct __UToT__<T, Color >
  {
    static inline T convert(const Color &u) { return reduce_add(u)/3.0f; }
  };

  template<typename T>
  struct __UToT__<T, Vector3f >
  {
    static inline T convert(const Vector3f &u) { return reduce_add(u)/3.0f; }
  };

  template<typename T, typename H>
  struct __UToT__<T, Vec2<H> >
  {
    static inline T convert(const Vec2<H> &u) { return reduce_add(u) / 2.0f; }
  };


  template<typename T, typename U = T>
  class lsDistribution1D {
  public:
    // Distribution1Df Public Methods
    lsDistribution1D() { _funcInt = 0; _count = 0; }
    lsDistribution1D(const U *f, uint32_t n) 
      : _func(n), _cdf(n+1)
    {
      _count = n;
      _cdf[0] = 0.0f;

      for (uint64_t i = 0; i < n; ++i)
        _func[i] = std::max<T>((T)0, __UToT__<T, U>::convert(f[i]));

      for (uint64_t i = 1; i < _cdf.size(); ++i)
        _cdf[i] = _cdf[i-1] + _func[i-1] / n;

      // Transform step function integral into CDF
      _funcInt = _cdf[n];
      if (_funcInt > 0.0f)
      {
        for (uint32_t i = 1; i < _cdf.size(); ++i)
          _cdf[i] /= _funcInt;
      }
    }
    ~lsDistribution1D() {
    }

    void Set(const U *f, uint32_t n)
    {
      _func.resize(n);
      _cdf.resize(n+1);

      _count = n;
      _cdf[0] = 0.0f;

      for (uint64_t i = 1; i < n; ++i)
        _func[i] = std::max<T>((T)0, __UToT__<T, U>::convert(f[i]));

      for (uint64_t i = 1; i < _cdf.size(); ++i)
        _cdf[i] = _cdf[i-1] + _func[i-1] / n;

      // Transform step function integral into CDF
      _funcInt = _cdf[n];
      if (_funcInt > 0.0f)
      {
        for (uint64_t i = 1; i < _cdf.size(); ++i)
          _cdf[i] /= _funcInt;
      }
    }

    T SampleContinuous(T u, T *pdf) const {
      if (_funcInt == 0.0f)
      {
        *pdf = 1.0f / _count;
        return min((uint32_t)(u * _count), (uint32_t)(_count - 1));
      }
      // Find surrounding CDF segments and _offset_
      typename vector<T>::const_iterator ptr = std::lower_bound(_cdf.begin(), _cdf.end(), u);
      intptr_t offset = max(static_cast<intptr_t>(0), ptr - _cdf.begin() - 1);

      // Compute offset along CDF segment
      T du = (u - _cdf[offset]) / (_cdf[offset+1] - _cdf[offset]);

      // Compute PDF for sampled offset
      if (pdf) *pdf = _func[offset] / _funcInt;

      // Return $x \in [0,1)$ corresponding to sample
      return (offset + du) / _count;
    }
    uint32_t SampleDiscrete(T u, T *pdf) const {
      if (_funcInt == 0.0f)
      {
        *pdf = 1.0f / _count;
        return min((uint32_t)(u * _count), (uint32_t)(_count - 1));
      }
      // Find surrounding CDF segments and _offset_
      typename vector<T>::const_iterator ptr = std::lower_bound(_cdf.begin(), _cdf.end(), u);
      uint32_t offset = max<uint32_t>(0, static_cast<uint32_t>(ptr - _cdf.begin() - 1));

      // Compute PDF for sampled offset
      if (pdf) *pdf = _func[offset] / (_funcInt * _count);
      return offset;
    }

    bool IsValid(){ return _funcInt > 0.0f; } 
  private:
    //friend struct Distribution2D;
    // Distribution1Df Private Data
    vector<T> _func;
    vector<T> _cdf;
    T _funcInt;
    uint32_t _count;
  };


  typedef lsDistribution1D<float> lsDistribution1Df;
  typedef lsDistribution1D<double> lsDistribution1Dd;
}
#endif // Distribution1D_h_

