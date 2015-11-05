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

#include "devices/device_singleray/clustering/vpl/vplIntegrator.h"

namespace embree
{
  vplIntegrator::vplIntegrator(const Parms& parms)
    : lightSampleID(-1), firstScatterSampleID(-1), firstScatterTypeSampleID(-1)
  {
    maxDepth        = parms.getInt  ("maxDepth"       ,10    );
    minContribution = parms.getFloat("minContribution",0.01f );
    epsilon         = parms.getFloat("epsilon"        ,32.0f)*float(ulp);
    backplate       = parms.getImage("backplate");
  }
  
  void vplIntegrator::requestSamples(Ref<SamplerFactory>& samplerFactory, const Ref<BackendScene>& scene)
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

  Color vplIntegrator::Li(LightPath& lightPath, const Ref<BackendScene>& scene, IntegratorState& state)
  {
    /*! Traverse ray. */
    DifferentialGeometry dg;
    //scene->intersector->intersect(lightPath.lastRay);
    rtcIntersect(scene->scene,(RTCRay&)lightPath.lastRay);
    scene->postIntersect(lightPath.lastRay,dg);
    state.numRays++;
    
    Color L = zero;
    const Vector3f wo = -lightPath.lastRay.dir;
    BRDFType directLightingBRDFTypes = (BRDFType)(ALL); 

    /*! Environment shading when nothing hit. */
    if (!lightPath.lastRay)
    {
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

    /*! face forward normals */
    bool backfacing = false;
    if (dot(dg.Ng, lightPath.lastRay.dir) > 0) {
      backfacing = true; dg.Ng = -dg.Ng; dg.Ns = -dg.Ns;
    }

    /*! Shade surface. */
    CompositedBRDF brdfs;
    if (dg.material) dg.material->shade(lightPath.lastRay, lightPath.lastMedium, dg, brdfs);

    /*! Add light emitted by hit area light source. */
    if (!lightPath.ignoreVisibleLights && dg.light && !backfacing)
      L += dg.light->Le(dg,wo);

    /*! Check if any BRDF component uses direct lighting. */
    bool useDirectLighting = false;
    for (size_t i=0; i<brdfs.size(); i++)
      useDirectLighting |= (brdfs[i]->type & directLightingBRDFTypes) != NONE;

    /*! Direct lighting. Shoot shadow rays to all light sources. */
    if (useDirectLighting)
    {
      for (size_t i=0; i<scene->vpls.size(); i++)
      {
        LightSample ls;
        ls.L = scene->vpls[i].sample(dg, ls.wi, ls.tMax);

        /*! Ignore zero radiance or illumination from the back. */
        //if (ls.L == Color(zero) || ls.wi.pdf == 0.0f || dot(dg.Ns,Vector3f(ls.wi)) <= 0.0f) continue; 
        if (ls.L == Color(zero) || ls.wi.pdf == 0.0f) continue;

        /*! Evaluate BRDF */
        Color brdf = brdfs.eval(wo, dg, ls.wi, ALL);
        if (brdf == Color(zero)) continue;

        /*! Test for shadows. */
        Ray shadowRay(dg.P, ls.wi, dg.error*epsilon, ls.tMax-dg.error*epsilon, lightPath.lastRay.time,dg.shadowMask);
        TIMESTART(stats.atomic.timeRT,rt0);
        rtcOccluded(scene->scene,(RTCRay&)shadowRay);
        TIMESTOP(stats.atomic.timeRT,rt0);
        state.numRays++;
        if (shadowRay) continue;

        /*! Evaluate BRDF. */
        L += ls.L * brdf * rcp(ls.wi.pdf);
      }
    }
    
    return L;
  }

  Color vplIntegrator::Li(Ray& ray, const Ref<BackendScene>& scene, IntegratorState& state) {
    LightPath path(ray); return Li(path,scene,state);
  }
}

