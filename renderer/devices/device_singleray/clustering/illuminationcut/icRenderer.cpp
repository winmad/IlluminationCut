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

#include "devices/device_singleray/clustering/illuminationcut/icRenderer.h"

#include "devices/device_singleray/clustering/illuminationcut/icIntegrator.h"

#include "devices/device_singleray/clustering/icStatistics.h"

/* include all samplers */
#include "devices/device_singleray/samplers/sampler.h"

/* include all image filters */
#include "devices/device_singleray/filters/boxfilter.h"
#include "devices/device_singleray/filters/bsplinefilter.h"

namespace embree
{
  icRenderer::icRenderer(const Parms& parms)
    : iteration(0)
  {
    /*! create integrator to use */
    std::string _integrator = parms.getString("integrator","icintegrator");
    if (_integrator == "icintegrator") integrator = new icIntegrator(parms);
    else throw std::runtime_error("unknown integrator type: "+_integrator);

    /*! create sampler to use */
    std::string _samplers = parms.getString("sampler","multijittered");
    if (_samplers == "multijittered"   ) samplers = new SamplerFactory(parms);
    else throw std::runtime_error("unknown sampler type: "+_samplers);

    /*! create pixel filter to use */
    std::string _filter = parms.getString("filter","none");
    if      (_filter == "none"   ) filter = NULL;
    else if (_filter == "box"    ) filter = new BoxFilter;
    else if (_filter == "bspline") filter = new BSplineFilter;
    else throw std::runtime_error("unknown filter type: "+_filter);

    /*! get framebuffer configuration */
    gamma = parms.getFloat("gamma",1.0f);

    /*! show progress to the user */
    showProgress = parms.getInt("showprogress",0);

    float error           = parms.getFloat("error",0.02f);
    bool sampling           = parms.getBool("sampling",true);
    this->epsilon           = parms.getFloat("epsilon",0.03);

    this->ic = new icIlluminationCut(error,sampling,samplers->samplesPerPixel);

  }

  void icRenderer::renderFrame(const Ref<Camera>& camera, const Ref<BackendScene>& scene, const Ref<ToneMapper>& toneMapper, Ref<SwapChain > swapchain, int accumulate) 
  {
    if (accumulate == 0) iteration = 0;
    new RenderJob(this,camera,scene,toneMapper,swapchain,accumulate,iteration);
    iteration++;
  }

  icRenderer::RenderJob::RenderJob (Ref<icRenderer> renderer, const Ref<Camera>& camera, const Ref<BackendScene>& scene, 
                                            const Ref<ToneMapper>& toneMapper, Ref<SwapChain > swapchain, int accumulate, int iteration)
    : renderer(renderer), camera(camera), scene(scene), toneMapper(toneMapper), swapchain(swapchain), 
      accumulate(accumulate), iteration(iteration), tileID(0), atomicNumRays(0)
  {
    numTilesX = ((int)swapchain->getWidth() +TILE_SIZE-1)/TILE_SIZE;
    numTilesY = ((int)swapchain->getHeight()+TILE_SIZE-1)/TILE_SIZE;
    rcpWidth  = rcp(float(swapchain->getWidth()));
    rcpHeight = rcp(float(swapchain->getHeight()));
    this->framebuffer = swapchain->buffer();
    this->framebuffer->startRendering(numTilesX*numTilesY);
    if (renderer->showProgress) new (&progress) Progress(numTilesX*numTilesY);
    if (renderer->showProgress) progress.start();
    renderer->samplers->reset();
    renderer->integrator->requestSamples(renderer->samplers, scene);
    renderer->samplers->init(iteration,renderer->filter);

    double tt = getSeconds();
    renderer->ic->init(scene,swapchain->getHeight() * swapchain->getWidth() * renderer->samplers->samplesPerPixel,renderer->epsilon);
    // get pixels and add them to the illuminationcut structure
    TaskScheduler::EventSync event;
    TaskScheduler::Task task(&event,_getPixels,this,TaskScheduler::getNumThreads(),_donothing,this,"render::getPixels");
    TaskScheduler::addTask(-1,TaskScheduler::GLOBAL_BACK,&task);
    event.sync();
    renderer->ic->pixeltree.shrink();
    // render using illuminationcut
    renderer->ic->illuminationcut();
    // create the image from the calculated values
    this->setPixels();

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

  void icRenderer::RenderJob::donothing(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event)
  {
    if (renderer->showProgress) progress.end();
  }


  void icRenderer::RenderJob::getPixels(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event)
  {
    /*! create a new sampler */
    IntegratorState state;
    if (taskIndex == taskCount-1) t0 = getSeconds();

    /*! tile pick loop */
    while (true)
    {
      /*! pick a new tile */
      size_t tile = tileID++;
      if (tile >= numTilesX*numTilesY) break;

      /*! process all tile samples */
      const int tile_x = (tile%numTilesX)*TILE_SIZE;
      const int tile_y = (tile/numTilesX)*TILE_SIZE;
      Random randomNumberGenerator(tile_x * 91711 + tile_y * 81551 + 3433*swapchain->firstActiveLine());

      for (size_t dy=0; dy<TILE_SIZE; dy++)
      {
        size_t y = tile_y+dy;
        if (y >= swapchain->getHeight()) continue;

        if (!swapchain->activeLine(y)) continue;
        size_t _y = swapchain->raster2buffer(y);

        for (size_t dx=0; dx<TILE_SIZE; dx++)
        {
          size_t x = tile_x+dx;
          if (x >= swapchain->getWidth()) continue;

          const int set = randomNumberGenerator.getInt(renderer->samplers->sampleSets);

          Color L = zero;
          size_t spp = renderer->samplers->samplesPerPixel;
          STATS(stats.atomic.pixels++);
          for (size_t s=0; s<spp; s++)
          {
            PrecomputedSample& sample = renderer->samplers->samples[set][s];
            const float fx = (float(x) + sample.pixel.x)*rcpWidth;
            const float fy = (float(y) + sample.pixel.y)*rcpHeight;

            Ray primary; camera->ray(Vec2f(fx,fy), sample.getLens(), primary);
            primary.time = sample.getTime();

            state.sample = &sample;
            state.pixel = Vec2f(fx,fy);
            bool hit = true;
            DifferentialGeometry dg;
            Vector3f wo;
            Color kd(zero),ks(zero);
            float exp(0.0f);
            L = renderer->integrator->getPixel(primary, scene, state,dg,wo,kd,ks,exp,hit);
            STATS(stats.atomic.gatherpoints++);
            // add the pixels to the ic structure or if not hit directly to the image
            if(hit)
            {
              renderer->ic->pixeltree.addPixel(dg,wo,kd,ks,exp,x,y,s);
#ifndef UB_IMAGE
              swapchain->update(x, _y, L, 0, true);
#endif
            }
            else
            {
#ifdef UB_IMAGE
              float r,g,b;
              ColormapJet::mapToJet(0,r,g,b);
              swapchain->update(x, _y, Color(r,g,b), 1, true);
#else
              swapchain->update(x, _y, L, 1, true);
#endif
            }
          }
        }
      }

      /*! print progress bar */
      if (renderer->showProgress) progress.next();

      /*! mark one more tile as finished */
      framebuffer->finishTile();
    }

    /*! we access the atomic ray counter only once per tile */
    atomicNumRays += state.numRays;
  }
  void icRenderer::RenderJob::setPixels()
  {
    // accumulate pixel values
    // should have a nice method for this
    std::vector<icPixel>& pxs = renderer->ic->pixeltree.pixels;
    for (size_t i = 0 ; i < pxs.size(); ++i)
    {
      size_t x  = pxs[i].x;
      size_t y  = pxs[i].y;
      size_t _y = swapchain->raster2buffer(y);
      Color L   = pxs[i].L;
      swapchain->update(x, _y, L, 1, true);
    }
    // tonemapping and setting the framebuffer
    for (size_t i = 0 ; i < swapchain->getWidth(); ++i)
    {
      for (size_t j = 0 ; j < swapchain->getHeight(); ++j)
      {
        // HACK to get the color of the pixel
        Color L = swapchain->update(i, j,Color(zero),0,true);
#if defined(UB_IMAGE)
        float ubNumber = L.r;
        float r,g,b;
        ColormapJet::mapToJet(ubNumber,r,g,b);
        const Color L1 = Color(r,g,b);
#else
        const Color L1 = toneMapper->eval(L,i,j,swapchain);
#endif
        framebuffer->set(i, j, L1);
      }
    }
  }
}
