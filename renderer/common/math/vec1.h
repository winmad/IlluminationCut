// ======================================================================== //
// Copyright 2009-2013 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#ifndef __EMBREE_VEC1_H__
#define __EMBREE_VEC1_H__

#include "sys/platform.h"
#include "math/math.h"

namespace embree
{
  ////////////////////////////////////////////////////////////////////////////////
  /// Generic 2D vector Class
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> struct Vec1
  {
    T x;

    typedef T Scalar;
    enum { N = 1 };

    ////////////////////////////////////////////////////////////////////////////////
    /// Construction
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vec1     ( )                  { }
    __forceinline Vec1     ( const Vec1& other ) { x = other.x; }
    template<typename T1> __forceinline Vec1( const Vec1<T1>& a ) : x(T(a.x)) {}
    template<typename T1> __forceinline Vec1& operator =( const Vec1<T1>& other ) { x = other.x; return *this; }

    __forceinline explicit Vec1( const T& a             ) : x(a) {}
    __forceinline explicit Vec1( const T* const a, const ssize_t stride = 1 ) : x(a[0]) {}

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vec1( ZeroTy   ) : x(zero) {}
    __forceinline Vec1( OneTy    ) : x(one) {}
    __forceinline Vec1( PosInfTy ) : x(pos_inf) {}
    __forceinline Vec1( NegInfTy ) : x(neg_inf) {}

    __forceinline const T& operator []( const size_t axis ) const { assert(axis < 1); return (&x)[axis]; }
    __forceinline       T& operator []( const size_t axis )       { assert(axis < 1); return (&x)[axis]; }

    __forceinline operator T(){return x;}
  };

  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline Vec1<T> operator +( const Vec1<T>& a ) { return Vec1<T>(+a.x); }
  template<typename T> __forceinline Vec1<T> operator -( const Vec1<T>& a ) { return Vec1<T>(-a.x); }
  template<typename T> __forceinline Vec1<T> abs       ( const Vec1<T>& a ) { return Vec1<T>(abs  (a.x)); }
  template<typename T> __forceinline Vec1<T> rcp       ( const Vec1<T>& a ) { return Vec1<T>(rcp  (a.x)); }
  template<typename T> __forceinline Vec1<T> rsqrt     ( const Vec1<T>& a ) { return Vec1<T>(rsqrt(a.x)); }
  template<typename T> __forceinline Vec1<T> sqrt      ( const Vec1<T>& a ) { return Vec1<T>(sqrt (a.x)); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline Vec1<T> operator +( const Vec1<T>& a, const Vec1<T>& b ) { return Vec1<T>(a.x + b.x); }
  template<typename T> __forceinline Vec1<T> operator -( const Vec1<T>& a, const Vec1<T>& b ) { return Vec1<T>(a.x - b.x); }
  template<typename T> __forceinline Vec1<T> operator *( const Vec1<T>& a, const Vec1<T>& b ) { return Vec1<T>(a.x * b.x); }
  template<typename T> __forceinline Vec1<T> operator *( const       T& a, const Vec1<T>& b ) { return Vec1<T>(a   * b.x); }
  template<typename T> __forceinline Vec1<T> operator *( const Vec1<T>& a, const       T& b ) { return Vec1<T>(a.x * b  ); }
  template<typename T> __forceinline Vec1<T> operator /( const Vec1<T>& a, const Vec1<T>& b ) { return Vec1<T>(a.x / b.x); }
  template<typename T> __forceinline Vec1<T> operator /( const Vec1<T>& a, const       T& b ) { return Vec1<T>(a.x / b  ); }
  template<typename T> __forceinline Vec1<T> operator /( const       T& a, const Vec1<T>& b ) { return Vec1<T>(a   / b.x); }

  template<typename T> __forceinline Vec1<T> min(const Vec1<T>& a, const Vec1<T>& b) { return Vec1<T>(min(a.x, b.x)); }
  template<typename T> __forceinline Vec1<T> max(const Vec1<T>& a, const Vec1<T>& b) { return Vec1<T>(max(a.x, b.x)); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline Vec1<T>& operator +=( Vec1<T>& a, const Vec1<T>& b ) { a.x += b.x; return a; }
  template<typename T> __forceinline Vec1<T>& operator -=( Vec1<T>& a, const Vec1<T>& b ) { a.x -= b.x; return a; }
  template<typename T> __forceinline Vec1<T>& operator *=( Vec1<T>& a, const       T& b ) { a.x *= b  ; return a; }
  template<typename T> __forceinline Vec1<T>& operator /=( Vec1<T>& a, const       T& b ) { a.x /= b  ; return a; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Reduction Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline T reduce_add( const Vec1<T>& a ) { return a.x ; }
  template<typename T> __forceinline T reduce_mul( const Vec1<T>& a ) { return a.x ; }
  template<typename T> __forceinline T reduce_min( const Vec1<T>& a ) { return a.x; }
  template<typename T> __forceinline T reduce_max( const Vec1<T>& a ) { return a.x; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline bool operator ==( const Vec1<T>& a, const Vec1<T>& b ) { return a.x == b.x; }
  template<typename T> __forceinline bool operator !=( const Vec1<T>& a, const Vec1<T>& b ) { return a.x != b.x; }
  template<typename T> __forceinline bool operator < ( const Vec1<T>& a, const Vec1<T>& b ) {
    if (a.x != b.x) return a.x < b.x;
    return false;
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Euclidian Space Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline T       dot      ( const Vec1<T>& a, const Vec1<T>& b ) { return a.x*b.x ; }
  template<typename T> __forceinline T       length   ( const Vec1<T>& a )                   { return sqrt(dot(a,a)); }
  template<typename T> __forceinline Vec1<T> normalize( const Vec1<T>& a )                   { return a*rsqrt(dot(a,a)); }
  template<typename T> __forceinline T       distance ( const Vec1<T>& a, const Vec1<T>& b ) { return length(a-b); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Select
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline Vec1<T> select ( bool s, const Vec1<T>& t, const Vec1<T>& f ) {
    return Vec1<T>(select(s,t.x,f.x));
  }

  template<typename T> __forceinline Vec1<T> select ( const typename T::Mask& s, const Vec1<T>& t, const Vec1<T>& f ) {
    return Vec1<T>(select(s,t.x,f.x));
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Output Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> inline std::ostream& operator<<(std::ostream& cout, const Vec1<T>& a) {
    return cout << "(" << a.x << ")";
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Default template instantiations
  ////////////////////////////////////////////////////////////////////////////////

  typedef Vec1<bool >         Vec1b;
  typedef Vec1<unsigned int > Vec1u;
  typedef Vec1<int  >         Vec1i;
  typedef Vec1<float>         Vec1f;
}

#endif
