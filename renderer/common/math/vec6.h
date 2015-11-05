// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#ifndef __EMBREE_VEC6_H__
#define __EMBREE_VEC6_H__

#include "sys/platform.h"
#include "vec3.h"
#include "vector3f_sse.h"

namespace embree
{
  ////////////////////////////////////////////////////////////////////////////////
  /// Generic 6D vector Class
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> struct Vec6
  {
    T x, y, z, d1, d2, d3;

    typedef T Scalar;
    enum { N = 6 };

    ////////////////////////////////////////////////////////////////////////////////
    /// Construction
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vec6    ( ) : x(zero),y(zero),z(zero),d1(zero),d2(zero),d3(zero)  { }
    __forceinline Vec6    ( const Vec6& other ) { x  = other.x ; y  = other.y ; z  = other.z ;
                                                  d1 = other.d1; d2 = other.d2; d3 = other.d3;}
    template<typename T1> __forceinline Vec6( const Vec6<T1>& a ) : x(T(a.x)), y(T(a.y)), z(T(a.z)), d1(T(a.d1),d2(T(a.d2)),d3(T(a.d3))) {}
    template<typename T1> __forceinline Vec6& operator =(const Vec6<T1>& other) { x = other.x; y = other.y; z = other.z; 
                                                                                 d1 = other.d1; d2 = other.d2; d3=other.d3; return *this; }

    __forceinline explicit Vec6( const T& a                                     ) : x(a), y(a), z(a), d1(a),d2(a),d3(a) {}
    //__forceinline explicit Vec6( const Vec3<T>& a, const Vec3<T>&   b           ) : x(a.x), y(a.y), z(a.z), d1(b.x),d2(b.y),d3(b.z) {}
    __forceinline explicit Vec6( const Vector3f& a, const Vector3f&   b           ) : x(a.x), y(a.y), z(a.z), d1(b.x),d2(b.y),d3(b.z) {}
    __forceinline explicit Vec6( const T& x, const T& y, const T& z, const T& d1,const T& d2,const T&d3) : x(x), y(y), z(z), d1(d1),d2(d2),d3(d3) {}
    __forceinline explicit Vec6( const T* const a, const size_t stride = 1     ) : x(a[0]), y(a[stride]), z(a[2*stride]), d1(a[3*stride]),d2(a[4*stride]),d3(a[5*stride]) {}

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vec6( ZeroTy   ) : x(zero), y(zero), z(zero), d1(zero),d2(zero),d3(zero) {}
    __forceinline Vec6( OneTy    ) : x(one),  y(one),  z(one),  d1(one),d2(one),d3(one) {}
    __forceinline Vec6( PosInfTy ) : x(pos_inf), y(pos_inf), z(pos_inf), d1(pos_inf),d2(pos_inf),d3(pos_inf) {}
    __forceinline Vec6( NegInfTy ) : x(neg_inf), y(neg_inf), z(neg_inf), d1(neg_inf),d2(neg_inf),d3(neg_inf) {}

    __forceinline const T& operator []( const size_t axis ) const { assert(axis < 6); return (&x)[axis]; }
    __forceinline       T& operator []( const size_t axis )       { assert(axis < 6); return (&x)[axis]; }

    //__forceinline       Vec3<T> pos(  )       { return Vec3<T>(x,y,z); }
    //__forceinline       Vec3<T> dir(  )       { return Vec3<T>(d1,d2,d3); }
    __forceinline       Vector3f pos(  )  const     { return Vector3f(x,y,z); }
    __forceinline       Vector3f dir(  )  const     { return Vector3f(d1,d2,d3); }
  };

  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline Vec6<T> operator +( const Vec6<T>& a ) { return Vec6<T>(+a.x, +a.y, +a.z, +a.d1, +a.d2,+a.d3); }
  template<typename T> __forceinline Vec6<T> operator -( const Vec6<T>& a ) { return Vec6<T>(-a.x, -a.y, -a.z, -a.d1, -a.d2,-a.d3); }
  template<typename T> __forceinline Vec6<T> abs       ( const Vec6<T>& a ) { return Vec6<T>(abs  (a.x), abs  (a.y), abs  (a.z), abs  (a.d1), abs  (a.d2), abs  (a.d3)); }
  template<typename T> __forceinline Vec6<T> rcp       ( const Vec6<T>& a ) { return Vec6<T>(rcp  (a.x), rcp  (a.y), rcp  (a.z), rcp  (a.d1), rcp  (a.d2), rcp  (a.d3)); }
  template<typename T> __forceinline Vec6<T> rsqrt     ( const Vec6<T>& a ) { return Vec6<T>(rsqrt(a.x), rsqrt(a.y), rsqrt(a.z), rsqrt(a.d1), rsqrt(a.d2), rsqrt(a.d3)); }
  template<typename T> __forceinline Vec6<T> sqrt      ( const Vec6<T>& a ) { return Vec6<T>(sqrt (a.x), sqrt (a.y), sqrt (a.z), sqrt (a.d1), sqrt (a.d2), sqrt (a.d3)); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline Vec6<T> operator +( const Vec6<T>& a, const Vec6<T>& b ) { return Vec6<T>(a.x + b.x, a.y + b.y, a.z + b.z, a.d1 + b.d1, a.d2 + b.d2, a.d3 + b.d3); }
  template<typename T> __forceinline Vec6<T> operator -( const Vec6<T>& a, const Vec6<T>& b ) { return Vec6<T>(a.x - b.x, a.y - b.y, a.z - b.z, a.d1 - b.d1, a.d2 - b.d2, a.d3 - b.d3); }
  template<typename T> __forceinline Vec6<T> operator *( const Vec6<T>& a, const Vec6<T>& b ) { return Vec6<T>(a.x * b.x, a.y * b.y, a.z * b.z, a.d1 * b.d1, a.d2 * b.d2, a.d3 * b.d3); }
  template<typename T> __forceinline Vec6<T> operator *( const       T& a, const Vec6<T>& b ) { return Vec6<T>(a   * b.x, a   * b.y, a   * b.z, a   * b.d1, a    * b.d2, a    * b.d3); }
  template<typename T> __forceinline Vec6<T> operator *( const Vec6<T>& a, const       T& b ) { return Vec6<T>(a.x * b  , a.y * b  , a.z * b  , a.d1 * b  , a.d2 * b   , a.d3 * b   ); }
  template<typename T> __forceinline Vec6<T> operator /( const Vec6<T>& a, const Vec6<T>& b ) { return Vec6<T>(a.x / b.x, a.y / b.y, a.z / b.z, a.d1 / b.d1, a.d2 / b.d2, a.d3 / b.d3); }
  template<typename T> __forceinline Vec6<T> operator /( const       T& a, const Vec6<T>& b ) { return Vec6<T>(a   / b.x, a   / b.y, a   / b.z, a   / b.d1, a    / b.d2, a    / b.d3); }
  template<typename T> __forceinline Vec6<T> operator /( const Vec6<T>& a, const       T& b ) { return Vec6<T>(a.x / b  , a.y / b  , a.z / b  , a.d1 / b  , a.d2 / b   , a.d3 / b   ); }

  template<typename T> __forceinline Vec6<T> min(const Vec6<T>& a, const Vec6<T>& b) { return Vec6<T>(min(a.x, b.x), min(a.y, b.y), min(a.z, b.z), min(a.d1, b.d1), min(a.d2,b.d2), min(a.d3,b.d3)); }
  template<typename T> __forceinline Vec6<T> max(const Vec6<T>& a, const Vec6<T>& b) { return Vec6<T>(max(a.x, b.x), max(a.y, b.y), max(a.z, b.z), max(a.d1, b.d1), max(a.d2,b.d2), max(a.d3,b.d3)); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline Vec6<T>& operator +=( Vec6<T>& a, const Vec6<T>& b ) { a.x += b.x; a.y += b.y; a.z += b.z; a.d1 += b.d1; a.d2 += b.d2; a.d3 += b.d3; return a; }
  template<typename T> __forceinline Vec6<T>& operator -=( Vec6<T>& a, const Vec6<T>& b ) { a.x -= b.x; a.y -= b.y; a.z -= b.z; a.d1 -= b.d1; a.d2 -= b.d2; a.d3 -= b.d3; return a; }
  template<typename T> __forceinline Vec6<T>& operator *=( Vec6<T>& a, const       T& b ) { a.x *= b  ; a.y *= b  ; a.z *= b  ; a.d1 *= b ; a.d2 *= b; a.d3*=b ; return a; }
  template<typename T> __forceinline Vec6<T>& operator /=( Vec6<T>& a, const       T& b ) { a.x /= b  ; a.y /= b  ; a.z /= b  ; a.d1 /= b ; a.d2 /= b; a.d3/=b ; return a; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Reduction Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline T reduce_add( const Vec6<T>& a ) { return a.x + a.y + a.z + a.d1 + a.d2 + a.d3; }
  template<typename T> __forceinline T reduce_mul( const Vec6<T>& a ) { return a.x * a.y * a.z * a.d1 * a.d2 * a.d3; }
  template<typename T> __forceinline T reduce_min( const Vec6<T>& a ) { return min(a.x, a.y, a.z, a.d1, a.d2, a.d3); }
  template<typename T> __forceinline T reduce_max( const Vec6<T>& a ) { return max(a.x, a.y, a.z, a.d1, a.d2, a.d3); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline bool operator ==( const Vec6<T>& a, const Vec6<T>& b ) { return a.x == b.x && a.y == b.y && a.z == b.z && a.d1 == b.d1 && a.d2 == b.d2 && a.d3 == b.d3; }
  template<typename T> __forceinline bool operator !=( const Vec6<T>& a, const Vec6<T>& b ) { return a.x != b.x || a.y != b.y || a.z != b.z || a.d1 != b.d1 || a.d2 != b.d2 || a.d3 != b.d3; }
  template<typename T> __forceinline bool operator < ( const Vec6<T>& a, const Vec6<T>& b ) {
    if (a.x != b.x) return a.x < b.x;
    if (a.y != b.y) return a.y < b.y;
    if (a.z != b.z) return a.z < b.z;
    if (a.d1 != b.d1) return a.d1 < b.d1;
    if (a.d2 != b.d2) return a.d2 < b.d2;
    if (a.d3 != b.d3) return a.d3 < b.d3;
    return false;
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Euclidian Space Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline T       dot      ( const Vec6<T>& a, const Vec6<T>& b ) { return a.x*b.x + a.y*b.y + a.z*b.z + a.d1*b.d1 + a.d2*b.d2 + a.d3*b.d3; }
  template<typename T> __forceinline T       length   ( const Vec6<T>& a )                   { return sqrt(dot(a,a)); }
  template<typename T> __forceinline T    sqrlength   ( const Vec6<T>& a )                   { return dot(a,a); }
  template<typename T> __forceinline Vec6<T> normalize( const Vec6<T>& a )                   { return a*rsqrt(dot(a,a)); }
  template<typename T> __forceinline T       distance ( const Vec6<T>& a, const Vec6<T>& b ) { return length(a-b); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Select
  ////////////////////////////////////////////////////////////////////////////////

  //template<typename T> __forceinline Vec4<T> select ( bool s, const Vec4<T>& t, const Vec4<T>& f ) {
  //  return Vec4<T>(select(s,t.x,f.x),select(s,t.y,f.y),select(s,t.z,f.z),select(s,t.w,f.w));
  //}


  template<typename T> __forceinline int maxDim ( const Vec6<T>& a ) 
  { 
    if(a.x >= a.y && a.x >= a.z && a.x >= a.d1 && a.x >= a.d2 && a.x >= a.d3) return 0;
    if(a.y >= a.x && a.y >= a.z && a.y >= a.d1 && a.y >= a.d2 && a.y >= a.d3) return 1;
    if(a.z >= a.x && a.z >= a.y && a.z >= a.d1 && a.z >= a.d2 && a.z >= a.d3) return 2;
    if(a.d1 >= a.x && a.d1 >= a.y && a.d1 >= a.z && a.d1 >= a.d2 && a.d1 >= a.d3) return 3;
    if(a.d2 >= a.x && a.d2 >= a.y && a.d2 >= a.z && a.d2 >= a.d1 && a.d2 >= a.d3) return 4;
    if(a.d3 >= a.x && a.d3 >= a.y && a.d3 >= a.z && a.d3 >= a.d2 && a.d3 >= a.d2) return 5;
    return -1; // error
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Output Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> inline std::ostream& operator<<(std::ostream& cout, const Vec6<T>& a) {
    return cout << "(" << a.x << ", " << a.y << ", " << a.z << ", " << a.d1 <<  ", " << a.d2 <<  ", " << a.d3 << ")";
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Default template instantiations
  ////////////////////////////////////////////////////////////////////////////////

  typedef Vec6<bool > Vec6b;
  typedef Vec6<int  > Vec6i;
  typedef Vec6<float> Vec6f;
}

#endif
