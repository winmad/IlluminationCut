#ifndef __EMBREE_VPLLIGHT_H__
#define __EMBREE_VPLLIGHT_H__

#include "devices/device_singleray/clustering/icParameters.h"
#include "devices/device_singleray/lights/light.h"
#include <iostream>

namespace embree
{
  /*! Implements a spot light source for VPL.  Always has a 180 degree falloff region*/
  class vplVPL
  {
  public:
    
    vplVPL(){}

    /*! Construction from position, direction, intensity and opening angles. */
    vplVPL(const Vector3f& P, const Vector3f& _D, const Color& I)
      : P(P), _D(normalize(_D)), I(I) {}

    __forceinline Color sample(const DifferentialGeometry& dg, Sample3f& wi, float& tMax) const
    {
      return sample(dg.P,wi,tMax);
    }

    __forceinline Color sample(const Vector3f& p, Sample3f& wi, float& tMax) const
    {
      Vector3f distVec = P - p;
      float distSquare = dot(distVec,distVec);
      float dist = sqrt(distSquare);
      wi = Sample3f(distVec * rcp(dist), max( (float) (CLAMPING_DISTANCE * CLAMPING_DISTANCE) , distSquare ) );
      tMax = dist;
      return I*clamp(dot(wi.value,_D));
    }

    __forceinline Color dirInt(const Vector3f wo) const
    {
      return I*clamp(dot(-wo,_D));
    }

    float pdf (const DifferentialGeometry& dg, const Vector3f& wi) const {
      return zero;
    }

    void read(FILE * ptr)
    {
      fread( &P.x, 1, sizeof(float), ptr ) ;
      fread( &P.y, 1, sizeof(float), ptr ) ;
      fread( &P.z, 1, sizeof(float), ptr ) ;
      fread( &_D.x, 1, sizeof(float), ptr ) ;
      fread( &_D.y, 1, sizeof(float), ptr ) ;
      fread( &_D.z, 1, sizeof(float), ptr ) ;
      fread( &I.r, 1, sizeof(float), ptr ) ;
      fread( &I.g, 1, sizeof(float), ptr ) ;
      fread( &I.b, 1, sizeof(float), ptr ) ;
    }
    void write(FILE * ptr)
    {
      fwrite( &P.x, 1, sizeof(float), ptr ) ;
      fwrite( &P.y, 1, sizeof(float), ptr ) ;
      fwrite( &P.z, 1, sizeof(float), ptr ) ;
      fwrite( &_D.x, 1, sizeof(float), ptr ) ;
      fwrite( &_D.y, 1, sizeof(float), ptr ) ;
      fwrite( &_D.z, 1, sizeof(float), ptr ) ;
      fwrite( &I.r, 1, sizeof(float), ptr ) ;
      fwrite( &I.g, 1, sizeof(float), ptr ) ;
      fwrite( &I.b, 1, sizeof(float), ptr ) ;
    }

  public:
    Vector3f  P;                        //!< Position of the spot light
    Vector3f _D;                        //!< Negative light direction of the spot light
    Color     I;                        //!< Radiant intensity (W/sr)
  };
}

#endif
