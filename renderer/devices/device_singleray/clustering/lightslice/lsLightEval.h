#ifndef _LIGHT_EVALUATION_H_
#define _LIGHT_EVALUATION_H_

#include "devices/device_singleray/default.h"
#include "devices/device_singleray/clustering/lightslice/lsLightData.h"
#include "devices/device_singleray/clustering/icParameters.h"
#include "devices/device_singleray/shapes/differentialgeometry.h"
#include "devices/device_singleray/api/scene.h"

namespace embree
{
  namespace LightEvalUtil{
    class EvalFunction
    {
    public:
      EvalFunction(float minGeoTerm = CLAMPING_DISTANCE*CLAMPING_DISTANCE) : _minGeoTerm(minGeoTerm) {}
      float _minGeoTerm;
    };
    /*

    class EvalIrrad : public EvalFunction
    {
    public:
    EvalIrrad(float minGeoTerm = DEFAULT_MIN_GEO_TERM) : EvalFunction(minGeoTerm) {}
    Vector3f operator()(const lsOrientedLight& light, const DifferentialGeometry& dp, const Vector3f &wo, Material *ms, RayEngine *engine, float rayEpsilon) const;
    Vector3f operator()(const lsDirLight& light, const DifferentialGeometry& dp, const Vector3f &wo, Material *ms, RayEngine *engine, float rayEpsilon) const;
    };
*/

    class EvalL : public EvalFunction
    {
    public:
      EvalL(float minGeoTerm = CLAMPING_DISTANCE*CLAMPING_DISTANCE) : EvalFunction(minGeoTerm) {}
      Color operator()(const lsOrientedLight& light, const DifferentialGeometry& dp, const Vector3f &wo,  Ref<BackendScene>& scene, float rayEpsilon) const;
      Color operator()(const lsDirLight& light, const DifferentialGeometry& dp, const Vector3f &wo,  Ref<BackendScene>& scene, float rayEpsilon) const;
    };

    /*
    class EvalShading : public EvalFunction
    {
    public:
    EvalShading(float minGeoTerm = DEFAULT_MIN_GEO_TERM) : EvalFunction(minGeoTerm) {}
    Vector3f operator()(const lsOrientedLight& light, const DifferentialGeometry& dp, const Vector3f &wo, Material *ms, RayEngine *engine, float rayEpsilon) const;
    Vector3f operator()(const lsDirLight& light, const DifferentialGeometry& dp, const Vector3f &wo, Material *ms, RayEngine *engine, float rayEpsilon) const;
    };

    class EvalVis : public EvalFunction
    {
    public:
    EvalVis(float minGeoTerm = DEFAULT_MIN_GEO_TERM) : EvalFunction(minGeoTerm) {}
    Vector3f operator()(const lsOrientedLight& light, const DifferentialGeometry& dp, const Vector3f &wo, Material *ms, RayEngine *engine, float rayEpsilon) const;
    Vector3f operator()(const lsDirLight& light, const DifferentialGeometry& dp, const Vector3f &wo, Material *ms, RayEngine *engine, float rayEpsilon) const;
    };
    */
    class EvalLight : public EvalFunction
    {
    public:
      EvalLight(float minGeoTerm = CLAMPING_DISTANCE*CLAMPING_DISTANCE) : EvalFunction(minGeoTerm) {}
      Color operator()(const lsOrientedLight& light, const DifferentialGeometry& dp, const Vector3f &wo,  Ref<BackendScene>& scene, float rayEpsilon) const;
      //Color operator()(const lsDirLight& light, const DifferentialGeometry& dp, const Vector3f &wo,  Ref<BackendScene>& scene, float rayEpsilon) const;
    };
  }
}
#endif // _LIGHT_EVALUATION_H_