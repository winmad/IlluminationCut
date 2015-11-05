#ifndef _LIGHT_CUT_INTEGRATOR_H_
#define _LIGHT_CUT_INTEGRATOR_H_
#include "devices/device_singleray/clustering/lsLightcuts/lcLightTree.h"
#include "devices/device_singleray/clustering/lightslice/lsPathSampler.h"
#include "devices/device_singleray/api/swapchain.h"
#include "devices/device_singleray/cameras/camera.h"
#include "devices/device_singleray/api/scene.h"

namespace embree {
  class Lightcutter
  {
  public:
    Lightcutter(lcLightTree *lightTree, Ref<BackendScene> & scene,Ref<Camera>& camera, float error = 0.01f, uint32_t maxCutSize = CT_MAXCUT);
    ~Lightcutter(void);

    virtual void	Lightcut(Ref<SwapChain>& image, uint32_t samples);
    virtual Color   EvaluateLightcut( const Ray &ray, uint32_t &cutSize );
  protected:
    Color           _EvalutateNode(const OrientedLightTreeNode *node, DifferentialGeometry& dp, const Vector3f &wo, float rayEpsilon, CompositedBRDF *ms);
    Color           _EvalutateNode(const SpotLightTreeNode *node, DifferentialGeometry& dp, const Vector3f &wo, float rayEpsilon, CompositedBRDF *ms);
    Color           _EvalutateNode(const DirectionalLightTreeNode *node, DifferentialGeometry& dp, const Vector3f &wo, float rayEpsilon, CompositedBRDF *ms);

    Color           _ComputeUpperBound(const OrientedLightTreeNode *node, DifferentialGeometry& dp, const Vector3f &wo, CompositedBRDF *ms);
    Color           _ComputeUpperBound(const SpotLightTreeNode *node, DifferentialGeometry& dp, const Vector3f &wo, CompositedBRDF *ms);
    Color           _ComputeUpperBound(const DirectionalLightTreeNode *node, DifferentialGeometry& dp, const Vector3f &wo, CompositedBRDF *ms);
    Color           _BoundMaterial(const BBox3f &bbox, const Vector3f &wo, const DifferentialGeometry & dp, CompositedBRDF *m, bool dirLight );
    Color           _EvaluateCut(DifferentialGeometry &dp, const Vector3f &wo, float rayEpsilon, CompositedBRDF *m, uint32_t &cs );
    Ref<BackendScene>                   _scene;
    Ref<Camera>                   _camera;
    uint32_t                            _maxCutSize;
    lcLightTree                          *_lightTree;
    float                               _clamp;
    float                               _error;
  };
}
#endif // _LIGHT_CUT_INTEGRATOR_H_
