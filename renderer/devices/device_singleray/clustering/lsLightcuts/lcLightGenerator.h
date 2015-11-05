//////////////////////////////////////////////////////////////////////////
//!
//!	\file    LightGenerator.h
//!	\date    24:3:2010   13:13
//!	\author  Jiawei Ou
//!	
//!	\brief   LightGenerator declaration
//!
//////////////////////////////////////////////////////////////////////////
#ifndef _VIRTUAL_LIGHT_GENERATOR_H_
#define _VIRTUAL_LIGHT_GENERATOR_H_

#include <vector>
#include "devices/device_singleray/default.h"
#include "devices/device_singleray/clustering/lsLightcuts/lcVLight.h"
#include "devices/device_singleray/clustering/lightslice/lsDistribution1D.h"
#include "devices/device_singleray/clustering/lightslice/lsPathSampler.h"

namespace embree {
  class Scene;
  class RayEngine;
  class ReportHandler;
  class Light;
  class BackendScene;

  //////////////////////////////////////////////////////////////////////////
  // Virtual Light Declaration
  class VirtualLightGenerator
  {
  public:
    virtual bool Generate(uint32_t indirect, VirtualLightCache* cache, float time = 0.0f, ReportHandler *report = 0) = 0;
  };

  //////////////////////////////////////////////////////////////////////////
  // Virtual Point Light Declaration
  class VirtualPointLightGenerator : public VirtualLightGenerator
  {
  public:
    VirtualPointLightGenerator(Ref<BackendScene>& scene);
    virtual bool Generate(uint32_t indirect, VirtualLightCache* cache, float time = 0.0f );
  protected:
    virtual void _GenerateDirect(VirtualLightCache* cache, float time);
    virtual void _GenerateIndirect(uint32_t indirect, VirtualLightCache* cache, float time);
    Vector3f GetSceneCenter() const { return _sceneCenter; }
    float GetSceneRadius() const { return _sceneRadius; }
  protected:
    static void ComputeLightSamplingCDF(lsDistribution1Df &dist, const vector<Ref<Light> > &lights, float sceneRadius);

    Ref<BackendScene>							_scene;
    StratifiedPathSamplerStd		_sampler;
    Vector3f                           _sceneCenter;
    float                           _sceneRadius;
    lsDistribution1Df					_dist;
  private:
    float							ClampDistSqr(float distSqr);
  };
}
#endif //_VIRTUAL_LIGHT_GENERATOR_H_