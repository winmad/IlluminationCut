#include "devices/device_singleray/clustering/lsLightcuts/lcLightcutter.h"
#include "devices/device_singleray/cameras/camera.h"
#include "devices/device_singleray/api/scene.h"
#include "devices/device_singleray/brdfs/lambertian.h"
#include "devices/device_singleray/brdfs/specular.h"
#include "devices/device_singleray/clustering/icStatistics.h"

namespace embree {
#define MAX_RAYTRACE_DEPTH 10

  struct LightcutHeapItem
  {
    LightcutHeapItem(LightTreeNode *node, Color ub, Color est = zero) : _node(node), _upperBound(ub), _estimate(est)
    {
      _e = reduce_avg(_upperBound);
    }
    bool operator< (const LightcutHeapItem &p2) const 
    {
      return _e == p2._e ? (_node < p2._node) : _e < p2._e;
    }
    LightTreeNode   *_node;
    Color           _upperBound;
    float           _e;
    Color           _estimate;
  };


  Lightcutter::Lightcutter(lcLightTree *lightTree, Ref<BackendScene> &scene, Ref<Camera>& camera,  float error, uint32_t maxCutSize) 
    : _scene(scene),_camera(camera), _maxCutSize(maxCutSize), _lightTree(lightTree), _error(error)
  {
    //float radius = (length(_scene->accel->bounds().size()) / 2.0f) * 0.05f;
    //_clamp = 4 * radius * radius;
    this->_clamp = CLAMPING_DISTANCE*CLAMPING_DISTANCE;
  }

  Lightcutter::~Lightcutter(void)
  {
  }

  void Lightcutter::Lightcut(Ref<SwapChain>& image, uint32_t samples)
  {
    Vec2f pixelSize(1.0f/image->getWidth(), 1.0f/image->getHeight());

    StratifiedPathSamplerStd sampler;
    for(uint32_t j = 0; j < image->getHeight(); j ++) 
    {
      for(uint32_t i = 0; i < image->getWidth(); i ++) 
      {
        Vec2f pixel(Vec2i(i,j));
        sampler.BeginPixel(samples);
        uint32_t cutSize = 0;
        Color L(zero);
        for (uint32_t s = 0; s < samples; s++)
        {
          uint32_t cs = 0;
          Vec2f puv = (pixel + ((samples == 1) ? Vec2f(0.5f, 0.5f) : sampler.Pixel())) * pixelSize;
          Ray ray ; this->_camera->ray(puv, sampler.Lens(), ray);  /// camera lens hack as 0.5
          L += EvaluateLightcut(ray, cs) / static_cast<float>(samples);
          cutSize += cs;
          sampler.NextPixelSample();
        }
        sampler.EndPixel();
        image->accumulate(pixel.x,pixel.y,L,1.0f);
      }
    }
  }

  Color Lightcutter::EvaluateLightcut( const Ray &r, uint32_t &cutSize )
  {
    Ray ray = r;
    Color throughput = one;
    DifferentialGeometry dg;
    rtcIntersect(_scene->scene,(RTCRay&)ray);
    _scene->postIntersect(ray,dg);
    if(!ray) {
      Color L(zero);  // 
      for (size_t i=0; i<_scene->envLights.size(); i++)
        L += _scene->envLights[i]->Le(-ray.dir);
      return throughput * L;
    }

    // face forward normals
    bool backfacing = false;
    if (dot(dg.Ng, ray.dir) > 0) {
      dg.Ng = -dg.Ng; dg.Ns = -dg.Ns; backfacing = true;
    }

    CompositedBRDF msu;
    if (dg.material) dg.material->shade(ray, Medium::Vacuum(), dg, msu);

    Vector3f wo = -ray.dir;

    uint32_t cs = 0;
    Color Ld = zero;
    Ld += _EvaluateCut(dg, wo, dg.error * 128.0f * std::numeric_limits<float>::epsilon() , &msu, cs);
    cutSize += cs;
    STATS(stats.atomic.clusters++);

    if ( dg.light && !backfacing)
      Ld += dg.light->Le(dg,ray.dir) ;
    return throughput * Ld;
  }

  Color Lightcutter::_EvalutateNode(const OrientedLightTreeNode *node, DifferentialGeometry& dp, const Vector3f &wo, float rayEpsilon,CompositedBRDF* brdfs)
  {
    lsOrientedLight& light = *(node->light);
    Vector3f wi = light.position - dp.P;
    float maxDist = length(wi);
    wi *= rcp(maxDist);
    float cosAngle = max(0.0f, dot(light.normal , -wi));
    float lenSqrEst = max(_clamp, maxDist*maxDist);
    //float lenSqrEst = (light.position - dp.P).GetLengthSqr();
    Color brdf = brdfs->eval(wo, dp,wi,ALL);
    Color estimate = brdf * cosAngle / lenSqrEst;

    if(zero != estimate) 
    {
      STATS(stats.atomic.rays++);
      Ray shadowRay(dp.P, wi, rayEpsilon, maxDist-rayEpsilon, 0.0f);
      rtcOccluded(_scene->scene,(RTCRay&)shadowRay);
      if( !shadowRay ){
        return estimate;
      }
    }
    return Color(zero);
  }


  Color Lightcutter::_EvalutateNode(const SpotLightTreeNode *node, DifferentialGeometry& dp, const Vector3f &wo, float rayEpsilon,CompositedBRDF* brdfs)
  {
    lsSpotLight& light = *(node->light);
    Vector3f wi = normalize(light.position - dp.P);
    float maxDist = length(light.position - dp.P);
    float cosAngle = max(0.0f, dot(light.normal , -wi));
    if(cosAngle <= light.cosouter) cosAngle =0 ;
    else if(cosAngle >= light.cosinner) cosAngle = 1;
    else {
      float delta = (cosAngle - light.cosouter)/(light.cosinner-light.cosouter);
      cosAngle = delta*delta*delta*delta;
    }
    float lenSqrEst = max(_clamp, squarelength(light.position - dp.P));
    //float lenSqrEst = (light.position - dp.P).GetLengthSqr();
    Color brdf = brdfs->eval(wo, dp,wi,ALL);
    Color estimate = brdf * cosAngle / lenSqrEst;

    if(zero != estimate) 
    {
      Ray shadowRay(dp.P, wi, rayEpsilon, maxDist-rayEpsilon, 0.0f);
      rtcOccluded(_scene->scene,(RTCRay&)shadowRay);
      if( !shadowRay )
        return estimate;
    }
    return Color(zero);
  }

  Color Lightcutter::_EvalutateNode(const DirectionalLightTreeNode *node, DifferentialGeometry& dp, const Vector3f &wo, float rayEpsilon,CompositedBRDF* brdfs)
  {
    lsDirLight* light = node->light;
    Vector3f wi = -light->normal;
    Color brdf = brdfs->eval(wo, dp,wi,ALL);
    Color estimate = brdf ;

    if(zero != estimate) 
    {
      Ray shadowRay(dp.P, wi, rayEpsilon, inf, 0.0f);
      rtcOccluded(_scene->scene,(RTCRay&)shadowRay);
      if( !shadowRay )
        return estimate;
    }
    return Color(zero);
  }

  Color Lightcutter::_ComputeUpperBound(const OrientedLightTreeNode *node, DifferentialGeometry & dp, const Vector3f &wo,CompositedBRDF* brdfs)
  {
    const Vector3f &m = node->bbox.lower;
    const Vector3f &M = node->bbox.upper;
    Vector3f bp = clamp(dp.P,m,M);
    float lenSqr = max((float) _clamp, squarelength(dp.P - bp));

    float cosBound = 1.0f;
    float cosHalfAngle = node->cone.GetAngleCos();
    if (cosHalfAngle <= 0.0f)
      cosBound = 1.0f;
    else
    {
      Vector3f vv = cross(Vector3f(0,0,1), node->cone.GetAxis());
      if(vv == zero) vv = Vector3f(0,1,0); // proper rotation axis if crossprod was 0
      vv = normalize(vv);
      LinearSpace3f mt = LinearSpace3f::rotate(vv, -acosf(clamp(dot(Vector3f(0,0,1) , node->cone.GetAxis()),-1.0f,1.0f)));

      Vector3f corners[8];
      BBox3f xbox = BBox3f(EmptyTy());
      node->bbox.GetCorners(corners);
      for (int i = 0; i < 8; i++)
        xbox.grow(xfmPoint( mt , dp.P - corners[i] ) );


      Vector3f &xm = xbox.lower;
      Vector3f &xM = xbox.upper;

      float cosTheta;
      if (xM.z > 0)
      {
        float minx2, miny2;
        if (xm.x * xM.x <= 0)
          minx2 = 0.0f;
        else
          minx2 = min(xm.x * xm.x, xM.x * xM.x);

        if (xm.y * xM.y <= 0)
          miny2 = 0.0f;
        else
          miny2 = min(xm.y * xm.y, xM.y * xM.y);
        float maxz2 = xM.z * xM.z;
        cosTheta = xM.z / sqrt(minx2 + miny2 + maxz2);
      }
      else
        cosTheta = 0;

      cosTheta = clamp(cosTheta, 0.0f, 1.0f);

      if (cosTheta > cosHalfAngle)
        cosBound = 1.0f;
      else
      {
        float sinHalfAngle = sqrt(1 - cosHalfAngle * cosHalfAngle);
        float sinTheta = sqrt(1 - cosTheta * cosTheta);
        cosBound = clamp(cosTheta * cosHalfAngle + sinTheta * sinHalfAngle, 0.0f, 1.0f);
        assert(cosBound >= 0.0f && cosBound <= 1.0f);
      }
    }
    return abs(node->L * _BoundMaterial(node->bbox, wo, dp, brdfs, false) * (cosBound / lenSqr));
  }

  Color Lightcutter::_ComputeUpperBound(const SpotLightTreeNode *node, DifferentialGeometry & dp, const Vector3f &wo,CompositedBRDF* brdfs)
  {
    const Vector3f &m = node->bbox.lower;
    const Vector3f &M = node->bbox.upper;
    Vector3f bp = clamp(dp.P,m, M);
    float lenSqr = max((float ) _clamp, squarelength(dp.P - bp));

    float cosBound = 1.0f;
    float cosHalfAngle = node->cone.GetAngleCos();
    if (cosHalfAngle <= 0.0f)
      cosBound = 1.0f;
    else
    {
      Vector3f vv = cross(Vector3f(0,0,1), node->cone.GetAxis());
      if(vv == zero) vv = Vector3f(0,1,0); // proper rotation axis if crossprod was 0
      vv = normalize(vv);
      LinearSpace3f mt = LinearSpace3f::rotate(vv, -acosf(clamp(dot(Vector3f(0,0,1) , node->cone.GetAxis()),-1.0f,1.0f)));

      Vector3f corners[8];
      BBox3f xbox = BBox3f(EmptyTy());
      node->bbox.GetCorners(corners);
      for (int i = 0; i < 8; i++)
        xbox.grow(xfmPoint( mt , dp.P - corners[i] ) );

      Vector3f &xm = xbox.lower;
      Vector3f &xM = xbox.upper;

      float cosTheta;
      if (xM.z > 0)
      {
        float minx2, miny2;
        if (xm.x * xM.x <= 0)
          minx2 = 0.0f;
        else
          minx2 = min(xm.x * xm.x, xM.x * xM.x);

        if (xm.y * xM.y <= 0)
          miny2 = 0.0f;
        else
          miny2 = min(xm.y * xm.y, xM.y * xM.y);
        float maxz2 = xM.z * xM.z;
        cosTheta = xM.z / sqrt(minx2 + miny2 + maxz2);
      }
      else
        cosTheta = 0;

      cosTheta = clamp(cosTheta, 0.0f, 1.0f);

      if (cosTheta > cosHalfAngle)
        cosBound = 1.0f;
      else
      {
        float sinHalfAngle = sqrt(1 - cosHalfAngle * cosHalfAngle);
        float sinTheta = sqrt(1 - cosTheta * cosTheta);
        cosBound = clamp(cosTheta * cosHalfAngle + sinTheta * sinHalfAngle, 0.0f, 1.0f);
        assert(cosBound >= 0.0f && cosBound <= 1.0f);
      }
    }
    return abs(node->L * _BoundMaterial(node->bbox, wo, dp, brdfs, false) * (cosBound / lenSqr));
  }

  Color Lightcutter::_ComputeUpperBound(const DirectionalLightTreeNode *node, DifferentialGeometry & dp, const Vector3f &wo,CompositedBRDF* brdfs)
  {
    return abs(node->L * _BoundMaterial(node->bbox, wo, dp, brdfs, true));
  }


  Color Lightcutter::_BoundMaterial( const BBox3f &bbox, const Vector3f &wo, const DifferentialGeometry & dp,CompositedBRDF* brdfs, bool dirLight )
  {
    //GLPhongApproximation matApprox = m->ApprtoximateAsGLPhong(dp);
    float n=0;
    Color kd(zero),ks(zero);
    brdfs->getPhong(kd,ks,n);

    Vector3f vv = cross(Vector3f(0,0,1) , dp.Ns);
    if (vv == zero) vv = Vector3f(0,1,0);
    vv=normalize(vv);
    LinearSpace3f mt = LinearSpace3f::rotate(vv, -acosf(clamp(dot(Vector3f(0,0,1) , dp.Ns),-1.0f,1.0f)));
    Vector3f corners[8];
    BBox3f xbox = BBox3f(EmptyTy());
    bbox.GetCorners(corners);
    for (int i = 0; i < 8; i++)
    {
      if (!dirLight)
        corners[i] -= dp.P;
      else
        corners[i] *= -1;
      xbox.grow(xfmPoint(mt,corners[i]));
    }

    Vector3f &xm = xbox.lower;
    Vector3f &xM = xbox.upper;

    float cosTheta;
    if (xM.z > 0)
    {
      float minx2, miny2;
      if (xm.x * xM.x <= 0)
        minx2 = 0.0f;
      else
        minx2 = min(xm.x * xm.x, xM.x * xM.x);

      if (xm.y * xM.y <= 0)
        miny2 = 0.0f;
      else
        miny2 = min(xm.y * xm.y, xM.y * xM.y);
      float maxz2 = xM.z * xM.z;
      cosTheta = xM.z / sqrt(minx2 + miny2 + maxz2);
    }
    else
      cosTheta = 0.0f;

    assert(cosTheta >= 0.0f && cosTheta <= 1.0f);
    cosTheta = clamp(cosTheta, 0.0f, 1.0f);

    if (zero!=ks && n != 0)
    {
      //Phong material
      float glossyCosBound = 1.0f;
      Vector3f R = - wo + dp.Ns*(2*(dot(wo,dp.Ns))); 
      if (R == wo || dp.Ns == wo)
      {
        R = wo;
      }

      Vector3f gv = cross(Vector3f(0,0,1) , R);
      // proper rotation axis in the case of 0 crossprod
      if (gv == zero) gv = Vector3f(0,1,0);
      gv=normalize(gv);
      LinearSpace3f gmt = LinearSpace3f::rotate(gv, -acosf(clamp(dot(Vector3f(0,0,1) , R),-1.0f,1.0f)));
      BBox3f gbox = BBox3f(EmptyTy());
      for (int i = 0; i < 8; i++)
        gbox.grow(xfmPoint(gmt,corners[i]));

      Vector3f &gm = gbox.lower;
      Vector3f &gM = gbox.upper;

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

      return (float(one_over_pi) * kd + (float(one_over_two_pi) * ks * (n+2) * pow(glossyCosBound, n))) * cosTheta;
    }

    return float(one_over_pi) * kd * cosTheta;
  }

  Color Lightcutter::_EvaluateCut( DifferentialGeometry &dp, const Vector3f &wo, float rayEpsilon, CompositedBRDF * brdfs, uint32_t &cs )
  {
    vector<LightcutHeapItem> lightcut;
    lightcut.reserve(_maxCutSize);

    Color totalEstL=zero, estimate=zero, bound=zero;
    if (_lightTree->GetOrientedRoot())
    {
      estimate = _EvalutateNode(_lightTree->GetOrientedRoot(), dp, wo, rayEpsilon, brdfs);
      bound = _ComputeUpperBound(_lightTree->GetOrientedRoot(), dp, wo, brdfs);
      totalEstL += _lightTree->GetOrientedRoot()->L * estimate;
      cs++;
      lightcut.push_back(LightcutHeapItem(_lightTree->GetOrientedRoot(), bound, estimate));
      push_heap(lightcut.begin(), lightcut.end());
    }
    if (false && _lightTree->GetSpotRoot())
    {
      throw std::runtime_error("Spotlights are not handled.");
      // HACK: Norbert -> Spotlights can only be original lights so we force them to descend until the root
      estimate = _EvalutateNode(_lightTree->GetSpotRoot(), dp, wo, rayEpsilon, brdfs);
      bound = _ComputeUpperBound(_lightTree->GetSpotRoot(), dp, wo, brdfs);
      bound = Color(FLT_MAX/4.0f);   // HACK: Norbert as indicated above
      totalEstL += _lightTree->GetSpotRoot()->L * estimate;
      cs++;
      lightcut.push_back(LightcutHeapItem(_lightTree->GetSpotRoot(), bound, estimate));
      push_heap(lightcut.begin(), lightcut.end());
    }

    if (_lightTree->GetDirectionalRoot())
    {
      estimate = _EvalutateNode(_lightTree->GetDirectionalRoot(), dp, wo, rayEpsilon, brdfs);
      bound = _ComputeUpperBound(_lightTree->GetDirectionalRoot(), dp, wo, brdfs);
      totalEstL += _lightTree->GetDirectionalRoot()->L * estimate;
      cs++;
      lightcut.push_back(LightcutHeapItem(_lightTree->GetDirectionalRoot(), bound, estimate));
      push_heap(lightcut.begin(), lightcut.end());
    }

    vector<Color> leafContrib;
    //std::cout<<"start while"<<std::endl;
    while(lightcut.size() > 0 && cs < _maxCutSize)
    {
      pop_heap(lightcut.begin(), lightcut.end());
      LightcutHeapItem cutItem = lightcut.back();
      //float error = (cutItem._upperBound / totalEstL).Average();
      Color avoidZeroDiv;
      avoidZeroDiv.r = (totalEstL.r == 0.0f) ? 1e-8f : totalEstL.r;
      avoidZeroDiv.g = (totalEstL.g == 0.0f) ? 1e-8f : totalEstL.g;
      avoidZeroDiv.b = (totalEstL.b == 0.0f) ? 1e-8f : totalEstL.b;

      float error = (cutItem._upperBound==zero) ? 0.0f : reduce_avg(cutItem._upperBound / avoidZeroDiv);
      //std::cout<<cutItem._upperBound<<" - "<<totalEstL<<": "<<error<<std::endl;
      if (error <= _error)
        break;
      //if (average(cutItem._upperBound) <= _error*average(totalEstL))
      //  break;
      lightcut.pop_back();

      if (cutItem._node == NULL || cutItem._node->IsLeaf())
        continue;

      if(cutItem._node->GetNodeType() == ORIENTED_NODE)
      {
        STATS(stats.atomic.clusters++);
        OrientedLightTreeNode *node = static_cast<OrientedLightTreeNode*>(cutItem._node);
        OrientedLightTreeNode *leftNode = node->left;
        Color leftEst = (node->light == leftNode->light) ? cutItem._estimate : _EvalutateNode(leftNode, dp, wo, rayEpsilon, brdfs);
        Color lL = leftNode->L * leftEst;

        OrientedLightTreeNode *rightNode = node->right;
        Color rightEst = (node->light == rightNode->light) ? cutItem._estimate : _EvalutateNode(rightNode, dp, wo, rayEpsilon, brdfs);
        Color rL = rightNode->L * rightEst;

        totalEstL -= node->L * cutItem._estimate;
        totalEstL += lL + rL; 
        cs++;

        if (leftNode->IsLeaf())
          leafContrib.push_back(lL);
        else
        {
          Color leftErr = _ComputeUpperBound(leftNode, dp, wo, brdfs);
          lightcut.push_back(LightcutHeapItem(leftNode, leftErr, leftEst));
          push_heap(lightcut.begin(), lightcut.end());
        }

        if (rightNode->IsLeaf())
          leafContrib.push_back(rL);
        else
        {
          Color rightErr = _ComputeUpperBound(rightNode, dp, wo, brdfs);
          lightcut.push_back(LightcutHeapItem(rightNode, rightErr, rightEst));
          push_heap(lightcut.begin(), lightcut.end());
        }
      }
      else if(cutItem._node->GetNodeType() == SPOT_NODE)
      {
        // HACK -> make it descend no matter what since we want every spotlight to be evaluated
        SpotLightTreeNode *node = static_cast<SpotLightTreeNode*>(cutItem._node);
        SpotLightTreeNode *leftNode = node->left;
        Color leftEst = (node->light == leftNode->light) ? cutItem._estimate : _EvalutateNode(leftNode, dp, wo, rayEpsilon, brdfs);
        Color lL = leftNode->L * leftEst;

        SpotLightTreeNode *rightNode = node->right;
        Color rightEst = (node->light == rightNode->light) ? cutItem._estimate : _EvalutateNode(rightNode, dp, wo, rayEpsilon, brdfs);
        Color rL = rightNode->L * rightEst;

        totalEstL -= node->L * cutItem._estimate;
        totalEstL += lL + rL; 
        cs++;

        if (leftNode->IsLeaf())
          leafContrib.push_back(lL);
        else
        {
          Color leftErr = Color(FLT_MAX/4.0f);
          lightcut.push_back(LightcutHeapItem(leftNode, leftErr, leftEst));
          push_heap(lightcut.begin(), lightcut.end());
        }

        if (rightNode->IsLeaf())
          leafContrib.push_back(rL);
        else
        {
          Color rightErr = Color(FLT_MAX/4.0f);
          lightcut.push_back(LightcutHeapItem(rightNode, rightErr, rightEst));
          push_heap(lightcut.begin(), lightcut.end());
        }
      }
      else if(cutItem._node->GetNodeType() == DIRECT_NODE)
      {
        DirectionalLightTreeNode *node = static_cast<DirectionalLightTreeNode*>(cutItem._node);
        DirectionalLightTreeNode *leftNode = node->left;
        Color leftEst = (node->light == leftNode->light) ? cutItem._estimate : _EvalutateNode(leftNode, dp, wo, rayEpsilon, brdfs);
        Color lL = leftNode->L * leftEst;

        DirectionalLightTreeNode *rightNode = node->right;
        Color rightEst = (node->light == rightNode->light) ? cutItem._estimate : _EvalutateNode(rightNode, dp, wo, rayEpsilon, brdfs);
        Color rL = rightNode->L * rightEst;

        totalEstL -= node->L * cutItem._estimate;
        totalEstL += lL + rL; 
        cs++;

        if (leftNode->IsLeaf())
          leafContrib.push_back(lL);
        else
        {
          Color leftErr = _ComputeUpperBound(leftNode, dp, wo, brdfs);
          lightcut.push_back(LightcutHeapItem(leftNode, leftErr, leftEst));
          push_heap(lightcut.begin(), lightcut.end());
        }

        if (rightNode->IsLeaf())
          leafContrib.push_back(rL);
        else
        {
          Color rightErr = _ComputeUpperBound(rightNode, dp, wo, brdfs);
          lightcut.push_back(LightcutHeapItem(rightNode, rightErr, rightEst));
          push_heap(lightcut.begin(), lightcut.end());
        }
      }
    }

    return totalEstL;
  }
}
