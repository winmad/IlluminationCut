#ifndef _MT_LIGHT_CUT_INTEGRATOR_H_
#define _MT_LIGHT_CUT_INTEGRATOR_H_
//#include "devices/device_singleray/clustering/lsLightcuts/lcLightTree.h"
#include "devices/device_singleray/clustering/lsLightcuts/lcLightcutter.h"
#include "devices/device_singleray/cameras/camera.h"

namespace embree{
  class MTLightcutter : public Lightcutter
  {
    friend class MTLightcutThread;
  public:
    MTLightcutter(lcLightTree *lightTree, Ref<BackendScene> &scene,Ref<Camera>& camera, float error = 0.01f, uint32_t maxCutSize = CT_MAXCUT);
    virtual ~MTLightcutter(void);

    virtual void Lightcut(Ref<SwapChain>& image, uint32_t samples);
  protected:
    uint32_t                _nCore;
    uint32_t                _curLine;
  };
}
#endif // _MT_LIGHT_CUT_INTEGRATOR_H_
