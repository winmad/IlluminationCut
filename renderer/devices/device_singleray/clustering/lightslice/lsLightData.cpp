#include "devices/device_singleray/clustering/lightslice/lsLightData.h"

namespace embree{
VL_TYPE LightList::GetLightType(uint32_t idx) const
{
  if (idx < _ptLights.size())
    return OMNIDIR_LIGHT;

  if (idx < _ptLights.size() + _dirLights.size())
    return DIRECTIONAL_LIGHT;

  if (idx < _ptLights.size() + _dirLights.size() + _otrLights.size())
    return ORIENTED_LIGHT;

  assert(false);
  return UNKNOWN_LIGHT;
}

const lsVLight* LightList::GetLight(uint32_t idx) const
{
  if (idx < _ptLights.size())
    return &_ptLights[idx];

  if (idx < _ptLights.size() + _dirLights.size())
    return &_dirLights[idx - _ptLights.size()];

  if (idx < _ptLights.size() + _dirLights.size() + _otrLights.size())
    return &_otrLights[idx - _ptLights.size() - _dirLights.size()];;

  assert(false);
  return NULL;
}

lsVLight* LightList::GetLight(uint32_t idx)
{
  if (idx < _ptLights.size())
    return &_ptLights[idx];

  if (idx < _ptLights.size() + _dirLights.size())
    return &_dirLights[idx - _ptLights.size()];

  if (idx < _ptLights.size() + _dirLights.size() + _otrLights.size())
    return &_otrLights[idx - _ptLights.size() - _dirLights.size()];;

  assert(false);
  return NULL;
}

Color LightList::GetPower( uint32_t idx ) const
{
  if (idx < _ptLights.size())
    return _ptLights[idx].intensity;

  if (idx < _ptLights.size() + _dirLights.size())
    return _dirLights[idx - _ptLights.size()].le;

  if (idx < _ptLights.size() + _dirLights.size() + _otrLights.size())
    return _otrLights[idx - _ptLights.size() - _dirLights.size()].le;

  assert(false);
  return Color(zero);
}

BBox<Vec1u> LightList::GetRange( VL_TYPE type ) const
{
  switch(type)
  {
  case OMNIDIR_LIGHT:
    return BBox<Vec1u>(Vec1u(0u), 
                       Vec1u((uint32_t)_ptLights.size()));
  case DIRECTIONAL_LIGHT:
    return BBox<Vec1u>(Vec1u((uint32_t) _ptLights.size()),
                       Vec1u((uint32_t)(_ptLights.size() + _dirLights.size())));
  case ORIENTED_LIGHT:
    return BBox<Vec1u>(Vec1u((uint32_t)(_ptLights.size() + _dirLights.size())),
                       Vec1u((uint32_t)(_ptLights.size() + _dirLights.size() + _otrLights.size())));
  default:
    break;
  }
  return BBox<Vec1u>(EmptyTy());
}
}
