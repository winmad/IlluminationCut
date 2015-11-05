// Copyright goes to the authors who published LightSlice

#ifndef _FAST_CONET_H_
#define _FAST_CONET_H_

// INCLUDES ====================================================
#include "devices/device_singleray/default.h"
#include "common/math/bbox.h"

namespace embree{
  
  typedef unsigned int uint32_t;

  template<class T>
  class FastCone {
    // data ---------------------------------------------
  protected:
    BBox<Vector3f > _dirBox;
  public:
    // Constructors -------------------------------------
    FastCone(const Vector3f& ax = Vector3f(0,0,1));
    FastCone(const FastCone<T>& v);

    // Access ops ---------------------------------------
    Vector3f		GetAxis() const;
    T		    GetAngleCos() const;
    T		    GetAngleSin() const;
    void		Set(const Vector3f *v, const uint32_t count);
    bool		IsValid() const;

    // Comparison ops -----------------------------------
    int operator == (const FastCone& v) const;
    int operator != (const FastCone& v) const;

    // Binary ops ---------------------------------------

    // Test ops -----------------------------------------
    bool Contain(const Vector3f& v) const;

    // Assignment ops -----------------------------------
    FastCone& operator=(const FastCone& v);

    // Vector ops ---------------------------------------
    void Grow(const Vector3f& v);
    void Grow(const FastCone& b);

    static bool Overlap(const FastCone& a, const FastCone& b);
    static FastCone Union(const FastCone& a, const FastCone& b);
  };

  // MACROS USEFUL FOR CODE GENERATION ===========================
  template <class T>
  inline FastCone<T>::FastCone(const Vector3f& ax) : _dirBox(ax, ax) {}

  template <class T>
  inline FastCone<T>::FastCone(const FastCone<T>& v) : _dirBox(v._dirBox) {} 

  template <class T>
  inline Vector3f FastCone<T>::GetAxis() const 
  {
    return normalize(center(_dirBox));
  }

  template <class T>
  inline T FastCone<T>::GetAngleSin() const 
  {
    T c = GetAngleCos();
    return sqrt(1 - c * c);
  }

  template <class T>
  inline T FastCone<T>::GetAngleCos() const 
  {
    Vector3f center = embree::center(_dirBox);
    T r2 = squarelength(_dirBox.upper - _dirBox.lower) * static_cast<T>(0.25);
    T d2 = squarelength(center);
    if (d2 == static_cast<T>(0))
    {
      return static_cast<T>(0);
    }
    float d = sqrt(d2);
    float a = (d2 - r2 + 1) / (2 * d);
    if (a < static_cast<T>(0))
    {
      a = static_cast<T>(0);
    }
    return clamp(a, 0.0f, 1.0f);
  }

  template <class T>
  inline void FastCone<T>::Set(const Vector3f *v, const uint32_t count) 
  {
    _dirBox = BBox<Vector3f>(EmptyTy());
    for (uint32_t i = 0; i < count; i++)
    {
      _dirBox.grow(v[i]);
    }
  }

  template <class T>
  inline bool FastCone<T>::IsValid() const 
  { 
    Vector3f center = embree::center(_dirBox);
    T r2 = squarelength(_dirBox.upper - _dirBox.lower) * static_cast<T>(0.25);
    T d2 = squarelength(center);
    if (d2 == static_cast<T>(0))
    {
      return false;
    }
    float d = sqrt(d2);
    float a = (d2 - r2 + 1) / (2 * d);
    return (a >= static_cast<T>(0) && a <= static_cast<T>(1));
  }

  template <class T>
  inline int FastCone<T>::operator == (const FastCone<T>& v) const {
    return _dirBox == v._dirBox;
  }
  template <class T>
  inline int FastCone<T>::operator != (const FastCone<T>& v) const {
    return !operator==(v); 
  }

  template <class T>
  inline bool FastCone<T>::Contain(const Vector3f& v) const {
    return conjoint(_dirBox,v);
  }

  template <class T>
  inline FastCone<T>& FastCone<T>::operator=(const FastCone<T>& v) {
    _dirBox = v._dirBox;
    return *this;
  }

  template <class T>
  inline void FastCone<T>::Grow(const Vector3f& v) {
    _dirBox.grow(v);
  }

  template <class T>
  inline void FastCone<T>::Grow(const FastCone<T>& b) {
    *this = Union(*this, b);
  }


  template <class T>
  inline FastCone<T> FastCone<T>::Union(const FastCone<T>& a, const FastCone<T>& b) {
    FastCone<T> result;
    result._dirBox = merge(a._dirBox, b._dirBox);
    return result;
  }

  template <class T>
  inline bool FastCone<T>::Overlap(const FastCone<T>& a, const FastCone<T>& b)
  {
    return conjoint(a._dirBox, b._dirBox);
  }
  
  typedef FastCone<float> FastConef;
  typedef FastCone<double> FastConed;
}
#endif
