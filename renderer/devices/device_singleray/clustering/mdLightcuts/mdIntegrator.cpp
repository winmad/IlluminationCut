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

#include "devices/device_singleray/clustering/mdLightcuts/mdIntegrator.h"

namespace embree
{
  mdIntegrator::mdIntegrator(const Parms& parms)
    : lightSampleID(-1), firstScatterSampleID(-1), firstScatterTypeSampleID(-1)
  {
    maxDepth        = parms.getInt  ("maxDepth"       ,10    );
    minContribution = parms.getFloat("minContribution",0.01f );
    epsilon         = parms.getFloat("epsilon"        ,32.0f)*float(ulp);
    backplate       = parms.getImage("backplate");
  }
  
  void mdIntegrator::requestSamples(Ref<SamplerFactory>& samplerFactory, const Ref<BackendScene>& scene)
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

  Color mdIntegrator::Li(Ray& ray, const Ref<BackendScene>& scene, IntegratorState& state) {
    throw std::runtime_error("Not implemented");
  }
}

