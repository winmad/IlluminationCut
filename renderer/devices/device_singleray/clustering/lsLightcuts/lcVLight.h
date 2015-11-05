//////////////////////////////////////////////////////////////////////////
//!
//!	\file    virtual_light.h
//!	\date    24:2:2010   15:08
//!	\author  Jiawei Ou
//!	
//!	\brief   virtual_light declaration
//!
//////////////////////////////////////////////////////////////////////////
#ifndef virtual_light_h__
#define virtual_light_h__

#include "common/math/vec4.h"
#include "common/math/affinespace.h"

namespace embree {

class VirtualLightCache
{
public:
  VirtualLightCache() {};
  ~VirtualLightCache() {};

  virtual void Clear() = 0;

  virtual void AddOmniDirLight(const Vec3f &position, const Vec3f &intensity) = 0;
  virtual void AddDirLight(const Vec3f &normal, const Vec3f &le) = 0;
  virtual void AddSpotLight(const Vec3f& position, const Vec3f &normal, const Vec3f &le,float cosinner,float cosouter)=0;
  virtual void AddOrientedLight(const Vec3f& position, const Vec3f &normal, const Vec3f &le) = 0;
  virtual void AddIndirectLight(const Vec3f& position, const Vec3f &normal, const Vec3f &le) = 0;

  virtual uint32_t DirectLightNum() const = 0;
  virtual uint32_t IndirectLightNum() const = 0;

  uint32_t		 LightNum() { return DirectLightNum() + IndirectLightNum(); }
};


struct VLightSample
{
  virtual void Transform(const AffineSpace3f &m) = 0;
  virtual bool IsValid() const = 0;
  virtual void AddToCache(VirtualLightCache *cache, uint32_t samples) const = 0;
};

struct OmniVLightSample : public VLightSample
{
  OmniVLightSample(const Vec3f &pos, const Vec3f &intensity);
  virtual void Transform(const AffineSpace3f &m);
  virtual bool IsValid() const;
  virtual void AddToCache(VirtualLightCache *cache, uint32_t samples) const;
  Vec3f   _position;
  Vec3f   _intensity;
};



struct DirVLightSample : public VLightSample
{
  DirVLightSample(const Vec3f &norm, const Vec3f &le);
  virtual void Transform(const AffineSpace3f &m);
  virtual bool IsValid() const;
  virtual void AddToCache(VirtualLightCache *cache, uint32_t samples) const;
  Vec3f   _le;
  Vec3f   _normal;
};

struct OrientdVLightSample : public VLightSample
{
  OrientdVLightSample(const Vec3f &pos, const Vec3f &norm, 
    const Vec3f &le);
  virtual void Transform(const AffineSpace3f &m);
  virtual bool IsValid() const;
  virtual void AddToCache(VirtualLightCache *cache, uint32_t samples) const;
  Vec3f   _position;
  Vec3f   _normal;
  Vec3f   _le;
};


struct SpotVLightSample : public VLightSample
{
  SpotVLightSample(const Vec3f &pos, const Vec3f& normal, const Vec3f &intensity, float cosIn, 
    float cosOut);
  virtual void Transform(const AffineSpace3f &m);
  virtual bool IsValid() const;
  virtual void AddToCache(VirtualLightCache *cache, uint32_t samples) const;
  Vec3f   _position;
  Vec3f   _normal;
  Vec3f   _le;
  float   _outterCosAngle;
  float   _innerCosAngle;
};
}
#endif // virtual_light_h__

