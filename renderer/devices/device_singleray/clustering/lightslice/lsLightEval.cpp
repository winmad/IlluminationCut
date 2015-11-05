#include "devices/device_singleray/clustering/lightslice/lsLightEval.h"
#include "devices/device_singleray/clustering/icStatistics.h"

namespace embree{
  namespace LightEvalUtil{

    Color EvalL::operator()(const lsOrientedLight& light, const DifferentialGeometry& dp, const Vector3f &wo, Ref<BackendScene>& scene, float rayEpsilon) const
    {
      // STATS(stats.atomic.nbClusters++;)
      Vector3f wi = normalize(light.position - dp.P);
      if(max(0.0f, dot(wi, dp.Ns) ) <= 0.0f)
        return Color(zero);
      float maxDist = length(light.position - dp.P);
      float cosAngle = max(0.0f, dot(light.normal , -wi));
      float lenSqrEst = max(_minGeoTerm, squarelength(light.position - dp.P));
      Color L = light.le * max(0.0f, dot(wi, dp.Ns) ) * cosAngle / lenSqrEst;
      if( zero != L) 
      {
        STATS(stats.atomic.rays++);
        Ray shadowRay(dp.P, wi, rayEpsilon, maxDist-rayEpsilon, 0.0f);
        rtcOccluded(scene->scene,(RTCRay&)shadowRay);
        if( !shadowRay )
          return L;
      }
      return Color(zero);
    }

    Color EvalL::operator()(const lsDirLight& light, const DifferentialGeometry& dp, const Vector3f &wo, Ref<BackendScene>& scene, float rayEpsilon) const
    {
      Vector3f wi = -light.normal;
      if(max(0.0f, dot(wi, dp.Ns) ) <= 0.0f)
        return Color(zero);
      Color L = light.le * max(0.0f, dot(wi, dp.Ns) );
      if(zero != L) 
      {
        Ray shadowRay(dp.P, wi, rayEpsilon, inf, 0.0f);
        rtcOccluded(scene->scene,(RTCRay&)shadowRay);
        if( !shadowRay )
        {
          return L;
        }
      }
      return Color(zero);
    }



    Color EvalLight::operator()(const lsOrientedLight& light, const DifferentialGeometry& dp, const Vector3f &wo, Ref<BackendScene>& scene, float rayEpsilon) const
    {
      STATS(stats.atomic.clusters++);
      Vector3f wi = normalize(light.position - dp.P);
      if(max(0.0f, dot(wi, dp.Ns) ) <= 0.0f)
        return Color(zero);
      float maxDist = length(light.position - dp.P);
      float cosAngle = max(0.0f, dot(light.normal , -wi));
      float lenSqrEst = max(_minGeoTerm, squarelength(light.position - dp.P));
      // shading was already done once, restore the Ns value, because it might be changed
      Vector3f ns = dp.Ns;
      CompositedBRDF brdfs;
      Ray eyeRay(dp.P+wo, -wo, 128.0f*std::numeric_limits<float>::epsilon(), 1.0f, 0.0f);
      if (dp.material) dp.material->shade(eyeRay, Medium::Vacuum(), dp, brdfs);
      dp.Ns = ns;

      Color brdf = brdfs.eval(wo, dp,wi,ALL);
      Color L = light.le * brdf * cosAngle / lenSqrEst;

      if( zero != L) 
      {
        STATS(stats.atomic.rays++);
        Ray shadowRay(dp.P, wi, rayEpsilon, maxDist-rayEpsilon, 0.0f);
        rtcOccluded(scene->scene,(RTCRay&)shadowRay);
        if( !shadowRay )
        {
          return L;
        }
      }
      return Color(zero);
    }

    /*Color EvalLight::operator()(const lsDirLight& light, const DifferentialGeometry& dp, const Vector3f &wo, Ref<BackendScene>& scene, float rayEpsilon) const
    {
      Vector3f wi = -light.normal;
      if(max(0.0f, dot(wi, dp.Ns) ) <= 0.0f)
        return Color(zero);
      CompositedBRDF brdfs;
      Ray eyeRay(dp.P+wo, -wo, 128.0f*std::numeric_limits<float>::epsilon(), 1.0f, 0.0f);
      if (dp.material) dp.material->shade(eyeRay, Medium::Vacuum(), dp, brdfs);
      Color L = light.le * brdfs.eval(wo,dp, wi, ALL);
      if(zero != L) 
      {
        Ray shadowRay(dp.P, wi, rayEpsilon, inf, 0.0f);
        rtcOccluded(scene->scene,(RTCRay&)shadowRay);
        if( !shadowRay )
          return L;
      }
      return Color(zero);
    }*/
  }
}
