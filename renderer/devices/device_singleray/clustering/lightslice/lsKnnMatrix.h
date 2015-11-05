#ifndef _LIGHTSLICE_KNN_MATRIX_H_
#define _LIGHTSLICE_KNN_MATRIX_H_


#include "devices/device_singleray/default.h"
#include "devices/device_singleray/api/swapchain.h"
#include "devices/device_singleray/clustering/lightslice/lsFastCone.h"
#include "devices/device_singleray/clustering/vpl/vplVPL.h"
#include "devices/device_singleray/brdfs/compositedbrdf.h"
#include "devices/device_singleray/clustering/lightcuts/ctKdTree.h"
#include "devices/device_singleray/brdfs/lambertian.h"
#include "devices/device_singleray/brdfs/specular.h"
#include "devices/device_singleray/clustering/lightslice/lsGatherPointShooter.h"
#include "devices/device_singleray/clustering/lightslice/lsLightData.h"
#include "devices/device_singleray/clustering/lightslice/lsPathSampler.h"
#include <vector>
#include <algorithm>

namespace embree{


  struct GatherKdItem
  {
    GatherKdItem() : idx(0) {}
    GatherKdItem(uint32_t i, const Vec6f& point) : idx(i), p(point) {} 
    uint32_t    idx;
    Vec6f       p;
  };

  struct ScaleLight
  {
    std::vector<uint32_t>    indices;
    std::vector<Color>       weights;
  };

  struct GatherGroup
  {
    GatherGroup() : bbox(BBox3f(EmptyTy())),normal(zero) {}
    std::vector<uint32_t>	indices;
    uint32_t            seed;
    std::vector<uint32_t>    neighbors;
    BBox3f				bbox;
    Vector3f               normal;
  };

  struct LightGroup
  {
    Vector3f               L;
    BBox3f				bbox;
    Vector3f               normal;
    FastConef           cone;
  };


  struct CloseSeed
  {
    void Set(GatherKdItem* p, float d2) { seedKdItem = p; distanceSquared = d2; }
    bool operator<(const CloseSeed &p2) const 
    {
      return distanceSquared == p2.distanceSquared ? 
        (seedKdItem < p2.seedKdItem) : distanceSquared < p2.distanceSquared;
    }

    GatherKdItem*         seedKdItem;
    float               distanceSquared;
  };

#define MAX_FOUND_SEED 64
  struct SeedProcess
  {
    SeedProcess(uint32_t n) : nLookup(n), foundSeeds(0) { }

    void operator()(GatherKdItem& hpItem, float dist2, float& maxDistSquared) const 
    {
      if (foundSeeds < nLookup) 
      {
        // Add photon to unordered array of photons
        closeSeeds[foundSeeds].Set(&hpItem, dist2);
        foundSeeds++;
        if (foundSeeds == nLookup) 
        {
          std::make_heap(&closeSeeds[0], &closeSeeds[nLookup]);
          maxDistSquared = closeSeeds[0].distanceSquared;
        }
      }
      else 
      {
        // Remove most distant photon from heap and add new photon
        std::pop_heap(&closeSeeds[0], &closeSeeds[nLookup]);
        closeSeeds[nLookup-1].Set(&hpItem, dist2);
        std::push_heap(&closeSeeds[0], &closeSeeds[nLookup]);
        maxDistSquared = closeSeeds[0].distanceSquared;
      }
    }

    mutable CloseSeed   closeSeeds[MAX_FOUND_SEED];
    uint32_t            nLookup;
    mutable uint32_t    foundSeeds;
  };


  class KnnMatrix : public RefCount
  {
    //friend class RenderKNNImageThread;
    //friend class RenderReducedMatrixThread;
    //friend struct FinalRenderThread;
    //friend struct FinalGroupRenderThread;
  public:
    KnnMatrix(){};
    void init(const Ref<BackendScene> &scene,Ref<Camera>& camera);
    ~KnnMatrix(void) {};
    //void                                  RenderGatherGroup(Image<Vector3f> *image);
    // seed:350,budget:2200 gives good results on office machine, exactly fits the memory etc..
    // for the blade set it to seed:
    void                                    Render( Ref<SwapChain>& image, uint32_t spp, const std::vector<vplVPL>& indirect, uint32_t seedNum = 400, uint32_t budget = 400);

    template<typename T> Color              RenderCell(const T &t, uint32_t col, uint32_t row );
    template<typename T> Color              RenderCell(const T &t, uint32_t col, const GatherPoint &gp );

  protected:
    void                                    _GenerateLights(const std::vector<vplVPL>& indirect);
    void                                    _ShootGatherPoints(uint32_t width, uint32_t height, uint32_t sample,Ref<SwapChain>& image);
    void                                    _GroupGatherPoints(uint32_t seedNum);
    void                                    _KdGatherGroup( vector<GatherKdItem>::iterator start, vector<GatherKdItem>::iterator end);
    void                                    _FindGatherGroupNeighbors();
    void                                    _RenderReducedMatrix(carray2<Color> &matrix);
    void                                    _InitialClusters(vector<vector<uint32_t> > &clusters, carray2<Color> &matrix, uint32_t budget, bool randProj);
    void                                    _RefineClusters(vector<vector<uint32_t> > &clusters, carray2<Color> &matrix, uint32_t budget, uint32_t samples, Ref<SwapChain>& image);
    void                                    _SetBackground(Ref<SwapChain> &image);

    /*
    protected:
    void                                    _KmeanGatherGroup( vector<GatherKdItem> items, uint32_t seedNum );
    void                                    _BuildGatherGroupKdTree();

    void									_NewInitialClusters(vector<vector<uint32_t> > &clusters, carray2<Vector3f> &matrix, uint32_t budget);

    VirtualLightGenerator					*_generator;




    */
  public:
    uint32_t                                _maxGatherGroupSize;
    float                                    _diagonal;
    Ref<BackendScene>                       _scene;
    Ref<Camera>                             _camera;
    float                                   _normScale;
    float                                   _clamp;
    StratifiedPathSamplerStd                _sampler;
    std::vector<GatherPoint>                     _gatherPoints;
    std::vector<BackgroundPixel>                 _bkPixels;
    std::vector<std::vector<ScaleLight> >             _scaledLights;
    std::vector<GatherGroup>                     _gpGroups;
    LightList                               _lightList;
  };

  template<typename T>
  Color KnnMatrix::RenderCell(const T &t, uint32_t col, uint32_t row )
  {
    GatherPoint &gp = _gatherPoints[row];
    return RenderCell(t, col, gp);
  }

  template<typename T>
  Color KnnMatrix::RenderCell(const T &t, uint32_t col, const GatherPoint &gp )
  {
    Color L(zero);
    switch (_lightList.GetLightType(col))
    {
    case OMNIDIR_LIGHT:
      assert(false); 
      L = Color(zero);
      break;
    case DIRECTIONAL_LIGHT:
      throw std::runtime_error("Embree::lightslice does not support directional lights");
        //L = t(*reinterpret_cast<DirLight*>(_lightList.GetLight(col)), 
        //    gp.isect.dp, gp.wo, gp.isect.m, _engine, gp.isect.rayEpsilon);
        break;
    case ORIENTED_LIGHT:
      L = t(*reinterpret_cast<lsOrientedLight*>(_lightList.GetLight(col)), 
        gp.isect, gp.wo, _scene, gp.isect.error * 128.0f * std::numeric_limits<float>::epsilon());
      break;
    default:
      break;
    }
    return L;
  }
}
#endif // _KNN_MATRIX_H_

