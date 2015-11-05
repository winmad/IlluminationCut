#ifndef _GATHER_TREE_BUILDER_H_
#define _GATHER_TREE_BUILDER_H_

#include "devices/device_singleray/clustering/mdLightcuts/mdGatherTree.h"
#include "devices/device_singleray/api/scene.h"
#include "devices/device_singleray/cameras/camera.h"
#include "devices/device_singleray/api/swapchain.h"
#include "devices/device_singleray/clustering/lightslice/lsPathSampler.h"

namespace embree{
struct GpKdItem
{
  uint32_t    idx;
  Vec6f       point;
};


class GatherTreeBuilder
{
public:
  GatherTreeBuilder(Ref<BackendScene> &scene,Ref<Camera>& camera, uint32_t width, uint32_t height, uint32_t samples, float c);
  ~GatherTreeBuilder(void);

  GatherNode*				Build(uint32_t i, uint32_t j, uint64_t seed, std::vector<mdGatherPoint> &points, Color *background); 
protected:
  GatherNode*             _BuildTree( std::vector<mdGatherPoint> &points, PathSampler* sampler);
  GatherNode*             _Build( std::vector<GpKdItem>::iterator start, std::vector<GpKdItem>::iterator end, std::vector<mdGatherPoint> &points, PathSampler* sampler);
  GatherNode*             _MakeLeaf( std::vector<GpKdItem>::iterator it, std::vector<mdGatherPoint> &points );
  float                   _BoundGlossyCos( const BBox3f &gbox );

private:
  Ref<BackendScene>      &_scene;
  Ref<Camera>            &_camera;
  uint32_t                _width;
  uint32_t                _height;
  uint32_t                _samples;
  Vec2f                   _pixelSize;
  float                   _normalScale;
  float                   _radius;

  Vec6f                   _ComputeKdPoint( mdGatherPoint &gp );
  void					_UpdateNode(GatherNode* node, PathSampler *sampler);

};
}
#endif
