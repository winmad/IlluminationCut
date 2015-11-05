#ifndef _LIGHT_TREE_BUILDER_H_
#define _LIGHT_TREE_BUILDER_H_

#include "devices/device_singleray/default.h"
#include "devices/device_singleray/clustering/lsLightcuts/lcLightTree.h"
#include "devices/device_singleray/api/scene.h"
#include "devices/device_singleray/clustering/lsLightcuts/lcVLight.h"
#include "devices/device_singleray/clustering/lsLightcuts/lcLightGenerator.h"
#include "devices/device_singleray/clustering/lsLightcuts/lcLightTreeBuilderUtil.h"
#include "devices/device_singleray/clustering/lightslice/lsRandom.h"

namespace embree{

class LightTreeBuilder
{
public:
  LightTreeBuilder() {} 
  virtual ~LightTreeBuilder(void) {}
protected:
  template<typename T>
  void _SampleLights(T* lightTree, BackendScene * scene, 
    uint32_t indirect);    
  float                               _c;
};


template<typename T>
void LightTreeBuilder::_SampleLights(T* lightTree, BackendScene* scene, 
  uint32_t indirect) 
{
  //BuilderLightCache<T> cache(lightTree);

  _c = scene->sceneradius / 16.0f;

  //_generator->Generate(indirect, &cache, 0.0f);

}
}

#endif
