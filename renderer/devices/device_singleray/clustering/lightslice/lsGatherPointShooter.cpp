#include "devices/device_singleray/clustering/lightslice/lsGatherPointShooter.h"
#include "devices/device_singleray/clustering/lightslice/lsPathSampler.h"
#include "common/sys/taskscheduler.h"
#include <vector>

namespace embree{
  
  class ShootGatherPointThread
  {
  public:
    typedef uint32_t argument_type;

    ShootGatherPointThread(Ref<BackendScene>& scene,Ref<Camera>& camera, const carray2<uint64_t> *randomSeeds, uint32_t width, uint32_t height, uint32_t samples,Ref<SwapChain>& image);
    TASK_RUN_FUNCTION(ShootGatherPointThread,compute);
    TASK_COMPLETE_FUNCTION(ShootGatherPointThread,donothing);
    void compute()
    {
      TaskScheduler::EventSync event;
      TaskScheduler::Task task(&event,_compute,this,TaskScheduler::getNumThreads(),_donothing,this,"render::getPixels");
      TaskScheduler::addTask(-1,TaskScheduler::GLOBAL_BACK,&task);
      event.sync();
    }
    void Init(void);
    void Clean(void);
    void shrink(void);
    static std::vector<GatherPoint>      _gatherPoints;  
    static std::vector<BackgroundPixel>  _bkPixels;     
    MutexActive bkMutex;
  private:
    void _Trace(const Vec2i &pixel, const Vec2f &offset, const Vec2f &auv, const Vec2f &pixelSize, uint32_t index, float time ) ;
    uint32_t                            _width;
    uint32_t                            _height;
    uint32_t                            _samples;
    Vec2f                               _pixelSize;
    uint32_t                            _nss;
    const carray2<uint64_t>               *_randSeeds;
    Ref<BackendScene>                     &_scene;
    Ref<Camera>                           &_camera;
    Ref<SwapChain>                   _image;
    Atomic                   pixelID;
    Atomic                   gpIdx;
    Atomic                   bkIdx;
  };


  ShootGatherPointThread::ShootGatherPointThread(Ref<BackendScene>& scene,Ref<Camera>& camera, const carray2<uint64_t> *randomSeeds, uint32_t width, uint32_t height, uint32_t samples,Ref<SwapChain>& image)
    : _scene(scene), _camera(camera), _width(width), _height(height), _samples(samples), _randSeeds(randomSeeds),_image(image)
  {
    _pixelSize = Vec2f(1.0f / _width, 1.0f / _height); 
    pixelID = 0;
    gpIdx = 0;
    bkIdx = 0;
  }

  void ShootGatherPointThread::Init()
  {
    _gatherPoints.clear();
    _bkPixels.clear();
    _gatherPoints.resize(_width * _height * _samples);
    _bkPixels.resize(_width * _height * _samples);
  }

  void ShootGatherPointThread::Clean()
  {
    _gatherPoints.clear();
    _gatherPoints.shrink_to_fit();
    _bkPixels.clear();
    _bkPixels.shrink_to_fit();
  }
  void ShootGatherPointThread::shrink()
  {
    _gatherPoints.resize(gpIdx);
    _bkPixels.resize(bkIdx);
  }

  void ShootGatherPointThread::compute(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event)
  {
    while (true)
    {
      /*! pick a new tile */
      uint32_t j =(uint32_t) pixelID++;
      if (j >= _randSeeds->height()) break;

    uint64_t seed = _randSeeds->operator[](j);
    StratifiedPathSamplerStd sampler(seed);
    for (uint32_t i = 0; i < _width; i++)
    {
      Vec2i pixel(i, j);
      if (_samples == 1)
        _Trace(pixel, Vec2f((float)0.5f, (float)0.5f), Vec2f((float)0.5f, (float)0.5f), _pixelSize, 0, 0.0f);
      else
      {
        sampler.BeginPixel(_samples);
        for (uint32_t s = 0; s < _samples; s++)
        {
          Vec2i pixel(i, j);
          _Trace(pixel, sampler.Pixel(), sampler.Lens(), _pixelSize, s, sampler.Time());
          sampler.NextPixelSample();
        }
        sampler.EndPixel();
      }
    }
    }
  }

  void ShootGatherPointThread::donothing(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event)
  {
  }

  //concurrent_vector<GatherPoint>          ShootGatherPointThread::_gatherPoints;
  //concurrent_vector<BackgroundPixel>      ShootGatherPointThread::_bkPixels;
  // ERROR resolve concurrency
  std::vector<GatherPoint>          ShootGatherPointThread::_gatherPoints;
  std::vector<BackgroundPixel>      ShootGatherPointThread::_bkPixels;

  void ShootGatherPointThread::_Trace( const Vec2i &pixel, const Vec2f &offset, const Vec2f &auv, const Vec2f &pixelSize, uint32_t index, float time ) 
  {
    STATS(stats.atomic.gatherpoints++);
    Vec2f puv = (Vec2f(pixel) + offset) * pixelSize;
    Ray ray ; _camera->ray(puv, auv, ray);
    GatherPoint gp;

    Color throughput = Color(one);

    uint32_t maxSpecDepth = 10;
    rtcIntersect(_scene->scene,(RTCRay&)ray);
    _scene->postIntersect(ray,gp.isect);
    if (ray)
    {
      // face forward normals
      bool backfacing = false;
      if (dot(gp.isect.Ng, ray.dir) > 0) {
        gp.isect.Ng = -gp.isect.Ng; gp.isect.Ns = -gp.isect.Ns; backfacing = true;
      }
      CompositedBRDF msu;
      if (gp.isect.material) gp.isect.material->shade(ray, Medium::Vacuum(), gp.isect, msu);
      if ( gp.isect.light && !backfacing)
        gp.emission = gp.isect.light->Le(gp.isect,ray.dir);
      else
        gp.emission = Color(zero);
      gp.pixel = pixel;
      gp.wo = -ray.dir;
      gp.weight = throughput;
      gp.strength = 1.0f / _samples;
      gp.index = index;
      _gatherPoints[gpIdx++]=gp;
    }
    else
    {
      Color L(zero);
      for (size_t i=0; i<_scene->envLights.size(); i++)
        L += _scene->envLights[i]->Le(-ray.dir);
      BackgroundPixel bkPixel;
      bkPixel.pixel = pixel;
      bkPixel.background = throughput * L;
      bkPixel.strength = 1.0f / _samples;
      _bkPixels[bkIdx++]=bkPixel;
    }
  }


  void GatherPointShooter::InitRandomSeeds(carray2<uint64_t> *randSeeds)
  {
    Random rnd; rnd.setSeed((int)time(0));
    for (uint32_t i = 0; i < randSeeds->size(); i++)
    {
      //static std::minstd_rand0 seeder;
      //uint64_t seed = seeder();
      uint64_t seed = rnd.getInt();
      randSeeds->operator[](i) = seed; 
    }
  }

  void GatherPointShooter::Shoot(Ref<BackendScene>& scene,Ref<Camera> camera, uint32_t width, uint32_t height, uint32_t samples, 
    vector<GatherPoint> &gatherPoints, vector<BackgroundPixel> &bgPixels,Ref<SwapChain>& image)
  {
    carray2<uint64_t> randSeeds(width, height);
    InitRandomSeeds(&randSeeds);

    ShootGatherPointThread thread(scene, camera, &randSeeds, width, height, samples,image);
    thread.Init();
    thread.compute();
    thread.shrink();

    //bgPixels.assign(thread._bkPixels.begin(), thread._bkPixels.end());
    //gatherPoints.assign(thread._gatherPoints.begin(), thread._gatherPoints.end());
    bgPixels.swap(thread._bkPixels);
    gatherPoints.swap(thread._gatherPoints);
  }

  /*
  void _Trace(Scene* scene, RayEngine *engine, const Vec2i &pixel, const Vec2f &offset, const Vec2f &auv, const Vec2f &pixelSize, uint32_t index, uint32_t samples, float time, vector<GatherPoint> &gatherPoints, vector<BackgroundPixel> &backPxs)
  {
    Vec2f puv = (Vec2f(pixel) + offset) * pixelSize;
    Ray ray = scene->MainCamera()->GenerateRay(puv, auv, time);
    GatherPoint gp;

    Vec3f throughput = Vec3f::One();
    bool hit = false;

    uint32_t maxSpecDepth = 10;
    while(engine->Intersect(ray, &gp.isect))
    {
      BxdfUnion msu;
      gp.isect.m->SampleReflectance(gp.isect.dp, msu);
      if(msu.HasDelta())
      {
        if (maxSpecDepth == 0)
        {
          hit = true;
          break;
        }
        maxSpecDepth--;
        BxdfSample sample = msu.SampleCos(DELTA_BXDF, -ray.D, gp.isect.dp, Vec2f::Zero(), 0.0f);
        throughput *= sample.brdfCos;
        ray = Ray(gp.isect.dp.P, sample.wi, gp.isect.rayEpsilon, RAY_INFINITY, ray.time);
        continue;
      }
      else if (msu.HasSmooth())
      {
        gp.emission = msu.Emission(gp.wo, gp.isect.dp);
        gp.pixel = pixel;
        gp.wo = -ray.D;
        gp.weight = throughput;
        gp.strength = 1.0f / samples;
        gp.index = index;
        gatherPoints.push_back(gp);
        hit = true;
      }
      break;
    }
    if (!hit)
    {
      if(scene->MainBackground()) 
      {	
        BackgroundPixel bkPixel;
        bkPixel.pixel = pixel;
        bkPixel.background = 
          throughput * scene->MainBackground()->SampleBackground(ray.D, ray.time);
        backPxs.push_back(bkPixel);
      }
    }
  }

  void Shoot(Scene* scene, uint32_t i, uint32_t j, uint32_t samples, carray2<uint64_t> *randSeeds, vector<GatherPoint> &gatherPoint, vector<BackgroundPixel> &bgPixels)
  {
    Vec2f pixelSize(1.0f / randSeeds->Width(), 1.0f / randSeeds->Height()); 
    Vec2i pixel(i, j);
    if (samples == 1)
      _Trace(scene, engine, pixel, Vec2f((float)0.5f, (float)0.5f), Vec2f((float)0.5f, (float)0.5f), pixelSize, 0, 1, 0.0f, gatherPoint, bgPixels);
    else
    {
      uint64_t seed = randSeeds->ElementAt(i, j);
      StratifiedPathSamplerStd::Engine e(seed);
      StratifiedPathSamplerStd sampler(e);
      sampler.BeginPixel(samples);
      for (uint32_t s = 0; s < samples; s++)
      {
        Vec2i pixel(i, j);
        _Trace(scene, engine, pixel, sampler.Pixel(), sampler.Lens(), pixelSize, s, samples, sampler.Time(), gatherPoint, bgPixels);
        sampler.NextPixelSample();
      }
      sampler.EndPixel();
    }
  }
  */
}
