#ifndef _MULTIDIMENSIONAL_LIGHTCUTTER_H_
#define _MULTIDIMENSIONAL_LIGHTCUTTER_H_
#include "devices/device_singleray/clustering/lsLightcuts/lcLightTree.h"
#include "devices/device_singleray/clustering/mdLightcuts/mdGatherTree.h"
#include "devices/device_singleray/api/swapchain.h"
#include "devices/device_singleray/clustering/lightslice/lsPathSampler.h"
#include "devices/device_singleray/api/scene.h"
#include "devices/device_singleray/cameras/camera.h"


namespace embree{

  class MdLightcutter
  {
    friend class MTMdLightcutThread;
  public:
    MdLightcutter(MdLightTree *lightTree, Ref<BackendScene> &scene, Ref<Camera> & camera, uint32_t maxCutSize = CT_MAXCUT);
    ~MdLightcutter(void);
    void Lightcut(Ref<SwapChain>& image, uint32_t samples );
  protected:
    Color _EvaluateLightcut(GatherNode *node, uint32_t &cutSize);
    Color _EvalutateNode(const lsOrientedLight*l, const GatherNode* gn);
    Color _EvalutateNode(const lsDirLight* l, const GatherNode* gn);
    Color _ComputeUpperBound( MdOrientedLightTreeNode* ltNode, GatherNode *gpNode, bool &refineLight);
    Color _ComputeUpperBound( MdDirectionalLightTreeNode* ltNode, GatherNode *gpNode, bool &refineLight);
    Color _BoundMaterialOrientedLight( MdOrientedLightTreeNode* ltNode, GatherNode *gpNode );
    Color _BoundMaterialDirLight( MdDirectionalLightTreeNode* ltNode, GatherNode *gpNode );
    template<typename LIGHT> uint32_t _SelectRandomLight(LIGHT* lights);
    Ref<BackendScene>                   _scene;
    Ref<Camera>                   _camera;
    uint32_t                 _maxCutSize;
    MdLightTree               *_lightTree;
    float                    _clamp;
    RandomPathSamplerStd	 _sampler;
  };

  template<typename LIGHT>
  uint32_t MdLightcutter::_SelectRandomLight(LIGHT* lights)
  {
    int matchingComps = 0;
    for (uint32_t i = 0; i < MD_REP_SLOTS; i++)
    {
      if (lights[i])
        matchingComps++;
    }

    int which = min<uint32_t>((uint32_t)(_sampler.Next1D() * matchingComps), matchingComps-1);
    for (uint32_t i = 0; i < MD_REP_SLOTS; ++i)
    {
      if (lights[i] && which-- == 0) 
        return i;
    }
    //assert(0);
    return 0;
  }

}
#endif // _MULTIDIMENSIONAL_LIGHTCUTTER_H_
