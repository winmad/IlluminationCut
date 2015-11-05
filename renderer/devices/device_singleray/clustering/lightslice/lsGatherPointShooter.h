#ifndef GatherPointShooter_h__
#define GatherPointShooter_h__

#include "devices/device_singleray/default.h"
#include "devices/device_singleray/clustering/lightslice/lsArrays.h"
#include "devices/device_singleray/api/scene.h"
#include "devices/device_singleray/cameras/camera.h"
#include "devices/device_singleray/api/swapchain.h"

namespace embree{
  typedef unsigned long uint64_t;

  struct GatherPoint 
  {
    Color			emission;
    Color			weight;
    Vec2i           pixel;
    DifferentialGeometry    isect;
    Vector3f           wo;
    float			strength;
    uint32_t        index;
  };


  struct BackgroundPixel 
  {
    Color			background;
    Vec2i           pixel;
    float			strength;
  };


  namespace GatherPointShooter
  {
    void InitRandomSeeds(carray2<uint64_t> *randSeeds);
    void Shoot(Ref<BackendScene>& scene,Ref<Camera> camera, uint32_t i, uint32_t j, uint32_t samples, carray2<uint64_t> *randSeeds, vector<GatherPoint> &gatherPoint, vector<BackgroundPixel> &bgPixels);
    void Shoot(Ref<BackendScene>& scene,Ref<Camera> camera, uint32_t width, uint32_t height, uint32_t samples, vector<GatherPoint> &gatherPoint, vector<BackgroundPixel> &bgPixels,Ref<SwapChain>& image);
  };


}
#endif // GatherPointShooter_h__