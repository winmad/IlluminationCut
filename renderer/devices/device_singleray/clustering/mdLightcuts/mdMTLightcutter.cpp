#include "devices/device_singleray/clustering/mdLightcuts/mdMTLightcutter.h"
#include "devices/device_singleray/api/scene.h"
#include "devices/device_singleray/api/swapchain.h"
#include "devices/device_singleray/clustering/lightslice/lsPathSampler.h"
#include "devices/device_singleray/clustering/mdLightcuts/mdGatherTreeBuilder.h"
#include "devices/device_singleray/clustering/icStatistics.h"
#include "common/sys/taskscheduler.h"

namespace embree{

  typedef unsigned long int uint64_t;
#define MAX_RAYTRACE_DEPTH 5

  MTMdLightcutter::MTMdLightcutter(MdLightTree *lightTree, Ref<BackendScene> &scene,Ref<Camera>& camera, uint32_t maxCutSize) 
    : MdLightcutter(lightTree, scene, camera, maxCutSize) {}

  MTMdLightcutter::~MTMdLightcutter(void)
  {
  }

  class MTMdLightcutThread
  {
  public:
    typedef uint32_t argument_type;

    MTMdLightcutThread(MTMdLightcutter *lcutter, GatherTreeBuilder *gtBuilder, carray2<uint64_t> *seeds, Ref<SwapChain>& img,Ref<Camera>& camera , uint32_t s) 
      : lightcutter(lcutter), gatherTreeBuilder(gtBuilder) , randSeeds(seeds), image(img), camera(camera), samples(s){
        scene = lightcutter->_scene;
        pixelSize = Vec2f(1.0f/image->getWidth(), 1.0f/image->getHeight());
        this->pixelID = 0;
    }
    TASK_RUN_FUNCTION(MTMdLightcutThread,compute);
    TASK_COMPLETE_FUNCTION(MTMdLightcutThread,donothing);
    void compute()
    {
      TaskScheduler::EventSync event;
      TaskScheduler::Task task(&event,_compute,this,TaskScheduler::getNumThreads(),_donothing,this,"render::getPixels");
      TaskScheduler::addTask(-1,TaskScheduler::GLOBAL_BACK,&task);
      event.sync();
    }
  private:
    MTMdLightcutter                 *lightcutter;
    GatherTreeBuilder		    *gatherTreeBuilder;
    carray2<uint64_t>		    *randSeeds;
    Ref<SwapChain>                 image;
    Ref<Camera>                      camera;
    uint32_t                         samples;
    Ref<BackendScene>	             scene;
    Vec2f                            pixelSize;
    Atomic                          pixelID;
  };

  void MTMdLightcutThread::compute(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event)
  {
    while (true)
    {
      /*! pick a new tile */ 
      uint32_t j =(uint32_t) pixelID++;
      if (j >= image->getHeight()) break;

      uint64_t seed = randSeeds->operator[](j);
      const float time = 0.0f;
      for(uint32_t i = 0; i < image->getWidth(); i ++) 
      {
        uint32_t cutSize = 0;
        std::vector<mdGatherPoint> points;
        Color background(zero);
        GatherNode *gpRoot = gatherTreeBuilder->Build(i, j, seed, points, &background);
        uint32_t cs = 0;
        Color L(zero);
        if (gpRoot) L = lightcutter->_EvaluateLightcut(gpRoot, cutSize);
        STATS(stats.atomic.clusters+=cutSize);
        image->accumulate(i, j , (background + L), 1.0f);
        if (gpRoot) delete gpRoot;
      }
    }
  }

  void MTMdLightcutThread::donothing(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event)
{
}

  void MTMdLightcutter::Lightcut(Ref<SwapChain> &image, uint32_t samples)
  {
    float normScale = _scene->sceneradius / 8.0f;

    std::shared_ptr<GatherTreeBuilder> gatherTreeBuilder = 
      std::shared_ptr<GatherTreeBuilder>(new GatherTreeBuilder(_scene,_camera,(uint32_t) image->getWidth(),(uint32_t) image->getHeight(), samples, normScale));


    carray2<uint64_t> randSeeds(image->getWidth(), image->getHeight());
    Random rnd; rnd.setSeed((int)time(0));
    for (uint32_t i = 0; i < randSeeds.size(); i++)
    {
      //static std::minstd_rand0 seeder;
      //uint64_t seed = seeder();
      uint64_t seed = rnd.getInt();
      randSeeds[i] = seed; 
    }

    MTMdLightcutThread thread(this, gatherTreeBuilder.get(), &randSeeds, image, _camera, samples);
    thread.compute();

  }
}
