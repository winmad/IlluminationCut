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

#ifndef __EMBREE_LC_INTEGRATOR_H__
#define __EMBREE_LC_INTEGRATOR_H__

#include "devices/device_singleray/integrators/integrator.h"
#include "devices/device_singleray/renderers//renderer.h"
#include "common/image/image.h"

namespace embree
{

  class lcIntegrator : public Integrator
  {

  public:

    /*! Construction of integrator from parameters. */
    lcIntegrator(const Parms& parms);

    /*! Registers samples we need tom the sampler. */
    void requestSamples(Ref<SamplerFactory>& samplerFactory, const Ref<BackendScene>& scene);

    /*! Computes the radiance arriving at the origin of the ray from the ray direction. */
    Color Li(Ray& ray, const Ref<BackendScene>& scene, IntegratorState& state);

    /* Configuration. */
  private:
    size_t maxDepth;               //!< Maximal recursion depth (1=primary ray only)
    float minContribution;         //!< Minimal contribution of a path to the pixel.
    float epsilon;                 //!< Epsilon to avoid self intersections.
    Ref<Image> backplate;          //!< High resolution background.

    /*! Random variables. */
  private:
    int lightSampleID;            //!< 2D random variable to sample the light source.
    int firstScatterSampleID;     //!< 2D random variable to sample the BRDF.
    int firstScatterTypeSampleID; //!< 1D random variable to sample the BRDF type to choose.
    std::vector<int> precomputedLightSampleID;  //!< ID of precomputed light samples for lights that need precomputations.
  };
}

#endif
