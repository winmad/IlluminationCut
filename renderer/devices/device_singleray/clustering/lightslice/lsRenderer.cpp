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

#include "devices/device_singleray/clustering/lightslice/lsRenderer.h"

#include "devices/device_singleray/clustering/lightslice/lsIntegrator.h"

#include "devices/device_singleray/clustering/icStatistics.h"

/* include all samplers */
#include "devices/device_singleray/samplers/sampler.h"

/* include all image filters */
#include "devices/device_singleray/filters/boxfilter.h"
#include "devices/device_singleray/filters/bsplinefilter.h"

namespace embree
{
  lsRenderer::lsRenderer(const Parms& parms)
    : iteration(0)
  {
    this->knnMat = new KnnMatrix;
    /*! create integrator to use */
    std::string _integrator = parms.getString("integrator","lsintegrator");
    if (_integrator == "lsintegrator") integrator = new lsIntegrator(parms);
    else throw std::runtime_error("unknown integrator type: "+_integrator);

    /*! create sampler to use */
    std::string _samplers = parms.getString("sampler","multijittered");
    if (_samplers == "multijittered"   ) samplers = new SamplerFactory(parms);
    else throw std::runtime_error("unknown sampler type: "+_samplers);

    /*! create pixel filter to use */
    std::string _filter = parms.getString("filter","bspline");
    if      (_filter == "none"   ) filter = NULL;
    else if (_filter == "box"    ) filter = new BoxFilter;
    else if (_filter == "bspline") filter = new BSplineFilter;
    else throw std::runtime_error("unknown filter type: "+_filter);

    /*! get framebuffer configuration */
    gamma = parms.getFloat("gamma",1.0f);

    /*! show progress to the user */
    showProgress = parms.getInt("showprogress",0);

    budget = parms.getInt("budget",400);
    seed = parms.getInt("seed",400);

  }

  void lsRenderer::renderFrame(const Ref<Camera>& camera, const Ref<BackendScene>& scene, const Ref<ToneMapper>& toneMapper, Ref<SwapChain > swapchain, int accumulate) 
  {
    if (accumulate == 0) iteration = 0;
    new RenderJob(this,camera,scene,toneMapper,swapchain,accumulate,iteration);
    iteration++;
  }

  lsRenderer::RenderJob::RenderJob (Ref<lsRenderer> renderer, const Ref<Camera>& camera, const Ref<BackendScene>& scene, 
                                            const Ref<ToneMapper>& toneMapper, Ref<SwapChain > swapchain, int accumulate, int iteration)
    : renderer(renderer), camera(camera), scene(scene), toneMapper(toneMapper), swapchain(swapchain), 
      accumulate(accumulate), iteration(iteration), tileID(0), atomicNumRays(0)
  {
    numTilesX = ((int)swapchain->getWidth() +TILE_SIZE-1)/TILE_SIZE;
    numTilesY = ((int)swapchain->getHeight()+TILE_SIZE-1)/TILE_SIZE;
    rcpWidth  = rcp(float(swapchain->getWidth()));
    rcpHeight = rcp(float(swapchain->getHeight()));
    this->framebuffer = swapchain->buffer();
    renderer->samplers->reset();
    renderer->integrator->requestSamples(renderer->samplers, scene);
    renderer->samplers->init(iteration,renderer->filter);

    double tt = getSeconds();

    //framebuffer->clear(Vec2i(zero),Vec2i((int)framebuffer->getWidth(),(int)framebuffer->getHeight())-Vec2i(1));
    renderer->knnMat->init(this->scene,this->camera);
    renderer->knnMat->Render(this->swapchain, renderer->samplers->samplesPerPixel, scene->vpls,renderer->seed,renderer->budget);
    for (size_t i = 0 ; i < swapchain->getWidth(); ++i)
    {
      for (size_t j = 0 ; j < swapchain->getHeight(); ++j)
      {
        // HACK to get the color of the pixel
        Color L = swapchain->update(i, j,Color(zero),0,true);
        const Color L1 = toneMapper->eval(L,i,j,swapchain);
        framebuffer->set(i, j, L1);
      }
    }

    double dt = getSeconds()-tt;

    /*! print fps, render time, and rays per second */
    std::ostringstream stream;
    stream << "render  ";
    stream.setf(std::ios::fixed, std::ios::floatfield);
    stream.precision(2);
    stream << 1.0f/dt << " fps, ";
    stream.precision(0);
    stream << dt*1000.0f << " ms, ";
    stream.precision(3);
    stream << atomicNumRays/dt*1E-6 << " mrps";
    std::cout << stream.str() << std::endl;
    // delete this
    delete this;
  }
}
