// Copyright goes to the authors of Lightslice who made this code available online.

// 
// Class ResizableArray 
//

#ifndef LIGHTSLICE_ARRAYS_H_
#define LIGHTSLICE_ARRAYS_H_

// LOOK AT VALARRAY as a substitute for this

// Include ---------------------------------------------------------------

//#include "compilerOptions.h"
#include <vector>
using std::vector;
#include <algorithm>
#include <stdint.h>
using std::copy;



namespace embree{

  typedef unsigned int uint32_t;

// Simple wrapper for C arrays
// allows easy creation, deletion, reallocation and keeps the size stored in
// use vector for more intelligent semantic than this
// use this only for simple stuff
template<class T>
class carray {
public:
  carray() { d = 0; l = 0; }
  carray(uint32_t ll) { l = ll; d = new T[l]; }
  carray(uint32_t ll, T* dd) { l = ll; d = new T[l]; copy(dd,dd+l,d); }
  carray(uint32_t ll, const T& dd) { l = ll; d = new T[l]; set(dd); }
  carray(const vector<T>& v) { l = (uint32_t)v.size(); d = new T[l]; copy(v.begin(), v.end(), d); }
  carray(const carray<T>& v) { l = v.size(); d = new T[l]; copy(v.d, v.d+l, d); }

  ~carray() { clear(); }

  carray<T>& operator = (const carray<T>& v) 
  {
    if (v.d == NULL || v.l == 0)
    {
      if(d) 
        delete[] d;
      l = 0;
    }
    else
    {
      resize(v.l); 
      copy(v.d, v.d+l, d); 
    }
    return *this; 
  }

  void tovector(vector<T>& v) { v.clear(); v.resize(l); copy(d, d+l, v.begin()); }

  bool operator==(const carray<T> &v) const {
    if(l != v.l) return false;
    for(uint32_t i = 0; i < l; i ++) {
      if(! ((*this)[i] == v[i]) ) return false;
    }
    return true;
  }

  void resize(uint32_t ll) { 
    if(l == ll) return;
    if(d != 0) delete[] d;
    l = ll;
    d = new T[l];
  }

  void clear() {
    l = 0;
    if(d) delete[] d;
    d = 0;
  }

  uint32_t size() const { return l; }
  bool empty() const { return l == 0; }

  T& operator[](uint32_t i) { return d[i]; }
  const T& operator[](uint32_t i) const { return d[i]; }
  const T& at(uint32_t i) const { return d[i]; }
  T& at(uint32_t i) { return d[i]; }

  T* data() { return d; }

  void set(const T& v) { for(uint32_t i = 0; i < l; i ++) d[i] = v; }

  uint32_t getMemoryAllocated() { return sizeof(uint32_t) + sizeof(T*) + l * sizeof(T); }

private:
  uint32_t l;
  T* d;
};

// Simple wrapper for C arrays accessed as 2D matrices / images
// NOTE: addressing is done as an image (i.e. this is a transposed matrix if looking at i,j)
// allows easy creation, deletion, reallocation and keeps the size stored in
// use vector for more intelligent semantic than this
// use this only for simple stuff
template<class T>
class carray2 {
public:
  carray2() : w(0), h(0) { }
  carray2(uint32_t ww, uint32_t hh) : w(ww), h(hh), d(ww*hh) { }
  carray2(uint32_t ww, uint32_t hh, T* dd) : w(ww), h(hh), d(ww*hh,dd) { }
  carray2(const carray2<T>& v) : w(v.w), h(v.h), d(v.d) { }

  carray2<T>& operator = (const carray2<T>& v) { w = v.w; h = v.h; d = v.d; return *this; }

  void resize(uint32_t ww, uint32_t hh) { 
    if(w == ww && h == hh) return;
    w = ww; h = hh;
    d.resize(w*h);
  }

  void clear() { w = 0; h = 0; d.clear(); }

  uint32_t size() const { return w*h; }
  uint32_t width() const { return w; }
  uint32_t height() const { return h; }
  bool empty() const { return size() == 0; }

  T& operator[](uint32_t i) { return d[i]; }
  const T& operator[](uint32_t i) const { return d[i]; }
  T& at(uint32_t i, uint32_t j) { return d[i+j*w]; }
  const T& at(uint32_t i, uint32_t j) const { return d[i+j*w]; }

  const T* row(uint32_t j) const { return &d[j*w]; } 
  T* row(uint32_t j) { return &d[j*w]; } 
  const T* data() const { return d.data(); }
  T* data() { return d.data(); }

  void set(const T& v) { for(uint32_t i = 0; i < w*h; i ++) d[i] = v; }

  long getMemoryAllocated() { return 2*sizeof(uint32_t) + sizeof(T*) + w*h * sizeof(T); }

private:
  uint32_t	w, h;
  carray<T>	d;
};
}
#endif // _RESIZABLEARRAY_H_
