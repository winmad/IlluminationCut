//////////////////////////////////////////////////////////////////////////
//!
//!	\file    VirtualLightData.h
//!	\date    31:3:2010   14:49
//!	\author  Jiawei Ou
//!	
//!	\brief   VirtualLightData declaration
//!
//////////////////////////////////////////////////////////////////////////
#ifndef _VIRTUAL_LIGHT_DATA_H_
#define _VIRTUAL_LIGHT_DATA_H_

#include "devices/device_singleray/default.h"
#include <vector>

namespace embree{
  typedef unsigned int uint32_t;
  struct lsVLight {};

  struct lsOmniDirLight : public lsVLight
  {
    lsOmniDirLight() {};
    lsOmniDirLight(const Vector3f &pos, const Color &it) 
      : position(pos), intensity(it) {}
    Vector3f position;
    Color intensity;
  };

  struct lsDirLight : public lsVLight
  {
    lsDirLight() {};
    lsDirLight(const Vector3f &norm, const Color &l) 
      : normal(norm), le(l) {} 
    Vector3f normal;
    Color le;
  };

  struct lsOrientedLight : public lsVLight
  {
    lsOrientedLight() {};
    lsOrientedLight(const Vector3f &pos, const Vector3f &norm, const Color &l) 
      : position(pos), normal(norm), le(l) {}
    Vector3f position;
    Vector3f normal;
    Color le;
  };

  struct lsSpotLight : public lsVLight
  {
    lsSpotLight() {};
    lsSpotLight(const Vector3f &pos, const Vector3f &norm, const Color &l,const float cosinner,float cosouter) 
      : position(pos), normal(norm), le(l) ,cosinner(cosinner),cosouter(cosouter){}
    Vector3f position;
    Vector3f normal;
    Color le;
    float cosinner;
    float cosouter;
  };

  enum VL_TYPE
  {
    OMNIDIR_LIGHT,
    DIRECTIONAL_LIGHT,
    ORIENTED_LIGHT,
    SPOT_LIGHT,
    UNKNOWN_LIGHT
  };


  struct LightList 
  {
    VL_TYPE         GetLightType(uint32_t idx) const;
    const lsVLight*   GetLight(uint32_t idx) const;
    lsVLight*         GetLight(uint32_t idx);
    Color           GetPower( uint32_t idx ) const;
    inline uint32_t GetSize() const { return static_cast<uint32_t>(_ptLights.size() + _dirLights.size() + _otrLights.size()+_spotLights.size()); }
    inline void     Clear() { _ptLights.clear(); _dirLights.clear(); _otrLights.clear(); _spotLights.clear();}
    BBox<Vec1u>     GetRange(VL_TYPE type) const;

    std::vector<lsOrientedLight>	_otrLights;
    std::vector<lsDirLight>		_dirLights;
    std::vector<lsOmniDirLight>	_ptLights;
    std::vector<lsSpotLight>	_spotLights;
  };
}
#endif // _VIRTUAL_LIGHT_DATA_H_
