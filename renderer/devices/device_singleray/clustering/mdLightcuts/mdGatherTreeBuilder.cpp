#include "devices/device_singleray/default.h"
#include "devices/device_singleray/clustering/mdLightcuts/mdGatherTreeBuilder.h"
#include "devices/device_singleray/cameras/camera.h"
#include "devices/device_singleray/api/scene.h"
#include "devices/device_singleray/brdfs/lambertian.h"
#include "devices/device_singleray/brdfs/specular.h"
#include <vector>


namespace embree{
  using std::vector;

  GatherTreeBuilder::GatherTreeBuilder(Ref<BackendScene> &scene,Ref<Camera>& camera, uint32_t width, uint32_t height, uint32_t samples, float c) 
    : _scene(scene), _camera(camera), _width(width), _height(height), _samples(samples), _normalScale(c)
  {
    _pixelSize = Vec2f(1.0f/_width, 1.0f/_height);
    _radius = _scene->sceneradius / 2.0f;
  }

  GatherTreeBuilder::~GatherTreeBuilder(void)
  {
  }

  GatherNode* GatherTreeBuilder::Build( uint32_t i, uint32_t j, uint64_t seed, vector<mdGatherPoint> &points, Color *background)
  { 
    STATS(stats.atomic.gatherpoints+=_samples);
    const float invSamples = 1.0f / _samples;
    StratifiedPathSamplerStd sampler(seed);
    sampler.BeginPixel(_samples);
    for (uint32_t s = 0; s < _samples; s++)
    {
      Vec2f pixel((float)i,(float)j);
      Vec2f puv = (pixel + ((_samples == 1) ? Vec2f(0.5f, 0.5f) : sampler.Pixel())) * _pixelSize;
      Ray ray ; this->_camera->ray(puv, sampler.Lens(), ray);  

      mdGatherPoint gp;
      Color throughput = Color(one);
      rtcIntersect(_scene->scene,(RTCRay&)ray);
      _scene->postIntersect(ray,gp.isect);

      if(!ray) {
        Color L(zero);
        for (size_t i=0; i<_scene->envLights.size(); i++)
          L += _scene->envLights[i]->Le(-ray.dir);
        *background += throughput * L * invSamples;
        continue;
      }
      // face forward normals
      bool backfacing = false;
      if (dot(gp.isect.Ng, ray.dir) > 0) {
        gp.isect.Ng = -gp.isect.Ng; gp.isect.Ns = -gp.isect.Ns; backfacing = true;
      }
      if ( gp.isect.light && !backfacing)
        gp.emission = gp.isect.light->Le(gp.isect,ray.dir);
      else
        gp.emission = Color(zero);
      gp.pixel = Vec2i(i,j);
      gp.wo = -ray.dir;
      gp.weight = throughput;
      points.push_back(gp);
      sampler.NextPixelSample();
    }
    sampler.EndPixel();
    return _BuildTree(points, &sampler);
  }

  GatherNode* GatherTreeBuilder::_BuildTree( vector<mdGatherPoint> &points, PathSampler* sampler)
  {
    if (points.size() == 0)
      return NULL;

    vector<GpKdItem> inputData;
    for (uint32_t s = 0; s < points.size(); s++)
    {
      GpKdItem item;
      item.idx = s;
      item.point = _ComputeKdPoint(points[s]);
      inputData.push_back(item);
    }

    GatherNode *node = _Build(inputData.begin(), inputData.end(), points, sampler);
    return node;
  }

  Vec6f GatherTreeBuilder::_ComputeKdPoint( mdGatherPoint &gp )
  {
    CompositedBRDF brdfs;
    // isect.Ns is modified in the shade function although it has been already set
    Vector3f ns = gp.isect.Ns;
    if (gp.isect.material) gp.isect.material->shade(Ray(), Medium::Vacuum(), gp.isect, brdfs);
    gp.isect.Ns = ns;
    float n=0;
    Color kd(zero),ks(zero);
    brdfs.getPhong(kd,ks,n);

    Vector3f direction;
    if (ks!=zero && n!=0.0f)
    {
      float scale = 0.0f;
      if(n < 20.0f)
        scale = 1.0f;
      else if(n < 40.0f)
        scale = 2.0f;
      else if(n < 80.0f)
        scale = 3.0f;
      else
        scale = 4.0f;
      direction =  (- gp.wo + gp.isect.Ns* 2*dot(gp.wo, gp.isect.Ns) ) * scale * _normalScale;
    }
    else
      direction = gp.isect.Ns * _normalScale;

    return Vec6f(gp.isect.P, direction);
  }

  void GatherTreeBuilder::_UpdateNode(GatherNode* node, PathSampler *sampler)
  {
    GatherNode *left = node->left;
    GatherNode *right = node->right;
    node->bbox = merge(left->bbox, right->bbox);
    node->normalCone = FastConef::Union(left->normalCone, right->normalCone);
    node->hasGlossy = left->hasGlossy || right->hasGlossy;
    if (left->hasGlossy && right->hasGlossy)
      node->mirrorCone = FastConef::Union(left->mirrorCone, right->mirrorCone);
    else if (right->hasGlossy == node->hasGlossy)
      node->mirrorCone = right->mirrorCone;
    else if (left->hasGlossy == node->hasGlossy)
      node->mirrorCone = left->mirrorCone;

    node->strength = left->strength + right->strength;
    node->emission = left->emission + right->emission;

    node->kd = max(left->kd, right->kd);
    node->ks = max(left->ks, right->ks);
    node->n = max(left->n, right->n);


    if (left->gp && right->gp)
      node->gp = sampler->Next1D() > 0.5f ? left->gp : right->gp;
    else if (left->gp) node->gp = left->gp;
    else if (right->gp) node->gp = right->gp;
  }


  GatherNode* GatherTreeBuilder::_Build( vector<GpKdItem>::iterator start, vector<GpKdItem>::iterator end, vector<mdGatherPoint> &points, PathSampler* sampler)
  {
    assert(end > start);
    if (end - start == 1)
      return _MakeLeaf(start, points);

    BBox<Vec6f> bbox = BBox<Vec6f>(EmptyTy());
    for (vector<GpKdItem>::iterator it = start; it != end; it++)
      bbox.grow(it->point);

    uint32_t dim = maxDim(bbox.size());
    vector<GpKdItem>::iterator mid = start + (end - start) / 2;
    std::nth_element(start, mid, end, [dim](const GpKdItem &d1, const GpKdItem &d2) -> bool {
      return d1.point[dim] < d2.point[dim];
    });

    GatherNode* node = new GatherNode();
    node->left = _Build(start, mid, points, sampler);
    node->right = _Build(mid, end, points, sampler);
    _UpdateNode(node, sampler);
    return node;
  }

  GatherNode* GatherTreeBuilder::_MakeLeaf( vector<GpKdItem>::iterator it, vector<mdGatherPoint> &points)
  {
    GatherNode *node = new GatherNode();
    mdGatherPoint &gp = points[it->idx];
    // the shading normal has been already computed, so restore it
    Vector3f ns = gp.isect.Ns;
    CompositedBRDF brdfs;
    if (gp.isect.material) gp.isect.material->shade(Ray(), Medium::Vacuum(), gp.isect, brdfs);
    gp.isect.Ns = ns;
    float n=0;
    Color kd(zero),ks(zero);
    brdfs.getPhong(kd,ks,n);

    node->gp = &gp;
    node->strength = 1.0f / _samples;
    node->hasGlossy = false;
    node->bbox = BBox3f(gp.isect.P, gp.isect.P);
    node->emission = gp.emission;
    node->kd = kd;
    node->ks = ks;
    node->n = n;

    if(zero != kd)
      node->normalCone = FastConef(gp.isect.Ns);
    if (zero != ks && n > 0.0f)
    {
      node->hasGlossy = true;
      node->mirrorCone =  FastConef(- gp.wo + gp.isect.Ns* 2*dot(gp.wo, gp.isect.Ns) );
    }
    return node;
  }

  float GatherTreeBuilder::_BoundGlossyCos( const BBox3f &gbox )
  {
    const Vector3f &gm = gbox.lower;
    const Vector3f &gM = gbox.upper;
    float glossyCosBound = 1.0f;
    if (gM.z > 0)
    {
      float minx2, miny2;
      if (gm.x * gM.x <= 0)
        minx2 = 0.0f;
      else
        minx2 = min(gm.x * gm.x, gM.x * gM.x);

      if (gm.y * gM.y <= 0)
        miny2 = 0.0f;
      else
        miny2 = min(gm.y * gm.y, gM.y * gM.y);
      float maxz2 = gM.z * gM.z;
      glossyCosBound = gM.z / sqrt(minx2 + miny2 + maxz2);
    }
    else
      glossyCosBound = 0.0f;

    glossyCosBound = clamp(glossyCosBound, 0.0f, 1.0f);
    assert(glossyCosBound <= 1.0f);
    return glossyCosBound;
  }
}
