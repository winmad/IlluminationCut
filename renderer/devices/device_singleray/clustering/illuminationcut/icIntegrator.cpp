// ======================================================================== //
// Copyright 2009-2013 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#include "devices/device_singleray/clustering/illuminationcut/icIntegrator.h"

namespace embree
{
  icIntegrator::icIntegrator(const Parms& parms)
    : lightSampleID(-1), firstScatterSampleID(-1), firstScatterTypeSampleID(-1)
  {
    maxDepth        = parms.getInt  ("maxDepth"       ,10    );
    minContribution = parms.getFloat("minContribution",0.01f );
    epsilon         = parms.getFloat("epsilon"        ,32.0f)*float(ulp);
    backplate       = parms.getImage("backplate");
  }
  
  void icIntegrator::requestSamples(Ref<SamplerFactory>& samplerFactory, const Ref<BackendScene>& scene)
  {
    precomputedLightSampleID.resize(scene->allLights.size());

    lightSampleID = samplerFactory->request2D();
    for (size_t i=0; i<scene->allLights.size(); i++) {
      precomputedLightSampleID[i] = -1;
      if (scene->allLights[i]->precompute())
        precomputedLightSampleID[i] = samplerFactory->requestLightSample(lightSampleID, scene->allLights[i]);
    }
    firstScatterSampleID = samplerFactory->request2D((int)maxDepth);
    firstScatterTypeSampleID = samplerFactory->request1D((int)maxDepth);
  }

  Color icIntegrator::getPixel(Ray& ray , const Ref<BackendScene>& scene,IntegratorState& state, DifferentialGeometry& dg,Vector3f& wo,Color& kd, Color& ks,float& exp,bool& hit)
  {
    /*! Traverse ray. */
    LightPath lightPath(ray);
    rtcIntersect(scene->scene,(RTCRay&)lightPath.lastRay);
    scene->postIntersect(lightPath.lastRay,dg);
    state.numRays++;

    Color L = zero;
    wo = -lightPath.lastRay.dir;

    /*! Environment shading when nothing hit. */
    if (!lightPath.lastRay)
    {
      hit = false;
      if (backplate && lightPath.unbend) {
        const int x = clamp(int(state.pixel.x * backplate->width ), 0, int(backplate->width )-1);
        const int y = clamp(int(state.pixel.y * backplate->height), 0, int(backplate->height)-1);
        L = backplate->get(x, y);
      }
      else {
        if (!lightPath.ignoreVisibleLights)
          for (size_t i=0; i<scene->envLights.size(); i++)
            L += scene->envLights[i]->Le(wo);
      }
      return L;
    }
    hit = true;

    /*! face forward normals */
    bool backfacing = false;
    if (dot(dg.Ng, lightPath.lastRay.dir) > 0) {
      backfacing = true; dg.Ng = -dg.Ng; dg.Ns = -dg.Ns;
    }

    /*! Shade surface. */
    CompositedBRDF brdfs;
    if (dg.material) dg.material->shade(lightPath.lastRay, lightPath.lastMedium, dg, brdfs);
    // get Phong approximation
    brdfs.getPhong(kd,ks,exp);
   
    /*! Add light emitted by hit area light source. */
    if (!lightPath.ignoreVisibleLights && dg.light && !backfacing)
      L += dg.light->Le(dg,wo);

    return L;
  }

  Color icIntegrator::Li(Ray& ray, const Ref<BackendScene>& scene, IntegratorState& state) {
    throw std::runtime_error("Not implemented");
  }
}

