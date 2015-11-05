#include "devices/device_singleray/clustering/lsLightcuts/lcMTLightcutter.h"
#include "devices/device_singleray/api/scene.h"
#include "devices/device_singleray/api/swapchain.h"
#include "devices/device_singleray/clustering/lightslice/lsPathSampler.h"
#include "common/sys/taskscheduler.h"
#include "devices/device_singleray/clustering/lightslice/lsArrays.h"
#include "common/sys/sync/atomic.h"
#include "devices/device_singleray/clustering/icStatistics.h"

namespace embree{

  typedef unsigned long int uint64_t;
#define MAX_RAYTRACE_DEPTH 5

  MTLightcutter::MTLightcutter(lcLightTree *lightTree, Ref<BackendScene> &scene, Ref<Camera>& camera,float error, uint32_t maxCutSize) 
    : Lightcutter(lightTree, scene,camera, error, maxCutSize) {}

  MTLightcutter::~MTLightcutter(void)
  {
  }

  class MTLightcutThread
  {
  public:
    typedef uint32_t argument_type;

    MTLightcutThread(MTLightcutter *lcutter, Ref<SwapChain>& img,Ref<Camera>& camera, carray2<uint64_t> *seeds, uint32_t s) 
      : lightcutter(lcutter), image(img), _camera(camera), randSeeds(seeds), samples(s)
    {
      scene = lightcutter->_scene;
      pixelSize = Vec2f(1.0f/image->getWidth(), 1.0f/image->getHeight());
      pixelID = 0;
    }
    TASK_RUN_FUNCTION(MTLightcutThread,compute);
    TASK_COMPLETE_FUNCTION(MTLightcutThread,donothing);
    void compute()
    {
      TaskScheduler::EventSync event;
      TaskScheduler::Task task(&event,_compute,this,TaskScheduler::getNumThreads(),_donothing,this,"render::getPixels");
      TaskScheduler::addTask(-1,TaskScheduler::GLOBAL_BACK,&task);
      event.sync();
    }
  private:
    MTLightcutter                   *lightcutter;
    Ref<SwapChain>                image;
    Ref<Camera>                     _camera;
    carray2<uint64_t>               *randSeeds;
    uint32_t                        samples;
    Ref<BackendScene>               scene;
    Vec2f                           pixelSize;
    Atomic                          pixelID;
  };

  void MTLightcutThread::compute(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event)  
  {
    while (true)
    {
      /*! pick a new tile */
      uint32_t j =(uint32_t) pixelID++;
      if (j >= image->getHeight()) break;

      uint64_t seed = randSeeds->operator[](j);
      StratifiedPathSamplerStd sampler(seed);

      for(uint32_t i = 0; i < image->getWidth(); i ++) 
      {
        sampler.BeginPixel(samples);
        uint32_t cutSize = 0;
        Color L(zero);
        for (uint32_t s = 0; s < samples; s++)
        {
          STATS(stats.atomic.gatherpoints++);
          uint32_t cs = 0;
          Vec2f pixel(Vec2i(i,j));
          Vec2f puv = (pixel + ((samples == 1) ? Vec2f(0.5f, 0.5f) : sampler.Pixel())) * pixelSize;
          Ray ray ; this->_camera->ray(puv, sampler.Lens(), ray);  /// camera lens hack as 0.5
          L += lightcutter->EvaluateLightcut(ray, cs) ;
          cutSize += cs;
          //std::cout<<cs<<std::endl;
          sampler.NextPixelSample();
        }
        sampler.EndPixel();
        image->accumulate(i,j,L,static_cast<float>(samples));
      }
    }
  }
  void MTLightcutThread::donothing(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event)
  {
  }

  void MTLightcutter::Lightcut(Ref<SwapChain>& image, uint32_t samples)
  {
    carray2<uint64_t> randSeeds(image->getWidth(), image->getHeight());
    Random rnd; rnd.setSeed((int)time(0));
    for (uint32_t i = 0; i < randSeeds.size(); i++)
    {
      //static std::minstd_rand0 seeder;
      //uint64_t seed = seeder();
      uint64_t seed = rnd.getInt();
      randSeeds[i] = seed; 
    }

    MTLightcutThread thread(this, image, _camera, &randSeeds, samples);
    thread.compute();
  }
}
