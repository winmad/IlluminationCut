#ifndef _MT_MD_LIGHT_CUT_INTEGRATOR_H_
#define _MT_MD_LIGHT_CUT_INTEGRATOR_H_
#include "devices/device_singleray/clustering/lsLightcuts/lcLightTree.h"
#include "devices/device_singleray/clustering/mdLightcuts/mdLightcutter.h"
#include "devices/device_singleray/cameras/camera.h"

namespace embree {
class MTMdLightcutter : public MdLightcutter
{
  friend class MTMdLightcutThread;
public:
  MTMdLightcutter(MdLightTree *lightTree, Ref<BackendScene> &scene,Ref<Camera>& camera,  uint32_t maxCutSize = 10000);
  virtual ~MTMdLightcutter(void);

  virtual void Lightcut(Ref<SwapChain> &image, uint32_t samples);
protected:
  uint32_t                _nCore;
  uint32_t                _curLine;
};

}
#endif // _MT_MD_LIGHT_CUT_INTEGRATOR_H_
