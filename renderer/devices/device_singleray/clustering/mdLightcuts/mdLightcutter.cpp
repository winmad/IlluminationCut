#include "devices/device_singleray/clustering/mdLightcuts/mdLightcutter.h"
#include "devices/device_singleray/cameras/camera.h"
#include "devices/device_singleray/api/scene.h"
#include "devices/device_singleray/clustering/mdLightcuts/mdGatherTreeBuilder.h"
#include "devices/device_singleray/clustering/icStatistics.h"

namespace embree{

MdLightcutter::MdLightcutter(MdLightTree *lightTree, Ref<BackendScene>& scene,Ref<Camera>& camera, uint32_t maxCutSize)
  : _scene(scene), _camera(camera),  _maxCutSize(maxCutSize), _lightTree(lightTree) 
{
  _clamp = CLAMPING_DISTANCE*CLAMPING_DISTANCE;
}

MdLightcutter::~MdLightcutter(void) {}

void MdLightcutter::Lightcut(Ref<SwapChain> &image, uint32_t samples)
{
  float normScale = _scene->sceneradius / 8.0f;

  std::shared_ptr<GatherTreeBuilder> gatherTreeBuilder = 
    std::shared_ptr<GatherTreeBuilder>(new GatherTreeBuilder(_scene, _camera, (uint32_t) image->getWidth(), (uint32_t)image->getHeight(), samples, normScale));

  carray2<uint64_t> randSeeds(image->getWidth(), image->getHeight());
  Random rnd; rnd.setSeed((int)time(0));
  for (uint32_t i = 0; i < randSeeds.size(); i++)
  {
    //static std::minstd_rand0 seeder;
    //uint64_t seed = seeder();
    uint64_t seed = rnd.getInt();
    randSeeds[i] = seed; 
  }

  for(uint32_t j = 0; j < image->getHeight(); j ++) 
  {
    for(uint32_t i = 0; i < image->getWidth(); i ++) 
    {
      vector<mdGatherPoint> points;
      Color background(zero);
      GatherNode *gpRoot = gatherTreeBuilder->Build(i, j, randSeeds.at(i, j), points, &background);
      uint32_t cutSize = 0;
      Color L = _EvaluateLightcut(gpRoot, cutSize);
      image->accumulate(i, j , background + L, 1.0f);
      if(gpRoot) delete gpRoot;
    }
  }
}

struct MdLightcutHeapItem
{
  MdLightcutHeapItem(LightTreeNode *node, GatherNode *gpNode, uint32_t ltIdx, const Color &ub, const Color &est, bool refineLight, uint32_t refineSeq) 
    : _ltNode(node), _gpNode(gpNode), _ltIdx(ltIdx), _upperBound(ub), _estimate(est), _refineLight(refineLight), _refineSeq(refineSeq) {
      _bound = reduce_avg(_upperBound);
  }

  bool operator< (const MdLightcutHeapItem &p2) const {
    return _bound == p2._bound ? (_ltNode < p2._ltNode) : _bound < p2._bound;
  }
  LightTreeNode   *_ltNode;
  GatherNode      *_gpNode;
  uint32_t		_ltIdx;
  Color           _upperBound;
  float           _bound;
  Color           _estimate;
  bool            _refineLight;
  uint32_t        _refineSeq;
};

Color MdLightcutter::_EvaluateLightcut(GatherNode *gpRoot, uint32_t &cutSize)
{
  if (gpRoot->gp == 0)
    return Color(zero);

  vector<MdLightcutHeapItem> lightcut;
  lightcut.reserve(_maxCutSize);
  cutSize = 0;

  Color totalEstL(zero), estimate(zero), bound(zero);
  if (_lightTree->GetOrientedRoot())
  {
    bool heuristic;
    MdOrientedLightTreeNode *root = _lightTree->GetOrientedRoot();
#if REP_LIGHT_NUM > 1 
    uint32_t index = _SelectRandomLight(root->lights);
    estimate = _EvalutateNode(root->lights[index], gpRoot);
#else
    uint32_t index = 0;
    estimate = _EvalutateNode(root->light, gpRoot);
#endif // MULTI_REP

    bound = _ComputeUpperBound(root, gpRoot, heuristic);
    totalEstL += root->L * estimate * gpRoot->strength;
    cutSize++;
    lightcut.push_back(MdLightcutHeapItem(root, gpRoot, index, bound, estimate, heuristic, 0));
    push_heap(lightcut.begin(), lightcut.end());
  }

  if (_lightTree->GetDirectionalRoot())
  {
    MdDirectionalLightTreeNode *root = _lightTree->GetDirectionalRoot();
    bool heuristic;
#if REP_LIGHT_NUM > 1 
    uint32_t index = _SelectRandomLight(root->lights);
    estimate = _EvalutateNode(root->lights[index], gpRoot);
#else
    uint32_t index = 0;
    estimate = _EvalutateNode(root->light, gpRoot);
#endif
    bound = _ComputeUpperBound(root, gpRoot, heuristic);
    totalEstL += root->L * estimate * gpRoot->strength;
    cutSize++;
    lightcut.push_back(MdLightcutHeapItem(root, gpRoot, index, bound, estimate, heuristic, 0));
    push_heap(lightcut.begin(), lightcut.end());
  }


  while(lightcut.size() > 0 && cutSize < _maxCutSize)
  {
    MdLightcutHeapItem cutItem = lightcut.front();
    if ( cutItem._upperBound.r <= 0.01f * totalEstL.r &&
         cutItem._upperBound.g <= 0.01f * totalEstL.g &&
         cutItem._upperBound.b <= 0.01f * totalEstL.b)
      break;

    STATS(stats.atomic.clusters++);
    pop_heap(lightcut.begin(), lightcut.end());
    lightcut.pop_back();

    LightTreeNode* node = cutItem._ltNode;
    GatherNode* gpNode = cutItem._gpNode;
    uint32_t index = cutItem._ltIdx;

    if (node == NULL || gpNode == NULL)
      std::cout << "null error" << endl;

    const mdGatherPoint *gp = gpRoot->gp;

    // Refine heuristic
    bool refineLightNode = cutItem._refineLight;
    if (cutItem._refineSeq > 4 && !gpNode->IsLeaf())
      refineLightNode = false;

    if (refineLightNode && !node->IsLeaf())
    {
      if (node->GetNodeType() == ORIENTED_NODE)
      {
        MdOrientedLightTreeNode *ltNode = static_cast<MdOrientedLightTreeNode*>(node);
        Color leftEst(zero), rightEst(zero);
        uint32_t lIdx = 0;
        uint32_t rIdx = 0;
        bool lHeur, rHeur;

        uint32_t seq = cutItem._refineSeq >= 0 ? cutItem._refineSeq + 1 : 0; 
        MdOrientedLightTreeNode *leftNode = ltNode->left;
        MdOrientedLightTreeNode *rightNode = ltNode->right;

#if REP_LIGHT_NUM > 1 
        if(ltNode->lights[index] == leftNode->lights[index]) {
          lIdx = index;
          leftEst = cutItem._estimate;
        } else {
          lIdx = _SelectRandomLight(leftNode->lights);
          leftEst = _EvalutateNode(leftNode->lights[lIdx], gpNode);
        }
#else
        leftEst = (ltNode->light == leftNode->light) ? cutItem._estimate : _EvalutateNode(leftNode->light, gpNode);
#endif
        Color leftBound = _ComputeUpperBound(leftNode, gpNode, lHeur) * gpNode->strength;
        lightcut.push_back(MdLightcutHeapItem(leftNode, gpNode, lIdx, leftBound, leftEst, lHeur, seq));
        push_heap(lightcut.begin(), lightcut.end());
        Color lL = leftNode->L * leftEst;

#if REP_LIGHT_NUM > 1 

        if(ltNode->lights[index] == rightNode->lights[index]) {
          rIdx = index;
          rightEst = cutItem._estimate;
        } else {
          rIdx = _SelectRandomLight(rightNode->lights);
          rightEst = _EvalutateNode(rightNode->lights[rIdx], gpNode);
        }
#else
        rightEst = (ltNode->light == rightNode->light) ? cutItem._estimate : _EvalutateNode(rightNode->light, gpNode);
#endif

        Color rightBound = _ComputeUpperBound(rightNode, gpNode, rHeur) * gpNode->strength;
        lightcut.push_back(MdLightcutHeapItem(rightNode, gpNode, rIdx, rightBound, rightEst, rHeur, seq));
        push_heap(lightcut.begin(), lightcut.end());
        Color rL = rightNode->L * rightEst;

        totalEstL -= (ltNode->L * cutItem._estimate) * gpNode->strength;
        totalEstL += (lL + rL) * gpNode->strength;
        cutSize++;
      }
      else if (node->GetNodeType() == DIRECT_NODE)
      {
        MdDirectionalLightTreeNode *ltNode = static_cast<MdDirectionalLightTreeNode*>(node);
        Color leftEst(zero), rightEst(zero);
        uint32_t lIdx = 0, rIdx = 0;
        bool lHeur, rHeur;
        uint32_t seq = cutItem._refineSeq >= 0 ? cutItem._refineSeq + 1 : 0; 

        MdDirectionalLightTreeNode *leftNode = ltNode->left;
        MdDirectionalLightTreeNode *rightNode = ltNode->right;
#if REP_LIGHT_NUM > 1 
        if(ltNode->lights[index] == leftNode->lights[index]) {
          lIdx = index;
          leftEst = cutItem._estimate;
        } else {
          lIdx = _SelectRandomLight(leftNode->lights);
          leftEst = _EvalutateNode(leftNode->lights[lIdx], gpNode);
        }
#else
        leftEst = (ltNode->light == leftNode->light) ? cutItem._estimate : _EvalutateNode(leftNode->light, gpNode);
#endif

        Color leftBound = _ComputeUpperBound(leftNode, gpNode, lHeur) * gpNode->strength;
        lightcut.push_back(MdLightcutHeapItem(leftNode, gpNode, lIdx, leftBound, leftEst, lHeur, seq));
        push_heap(lightcut.begin(), lightcut.end());
        Color lL = leftNode->L * leftEst;

#if REP_LIGHT_NUM > 1 
        if(ltNode->lights[index] == rightNode->lights[index]) {
          rIdx = index;
          rightEst = cutItem._estimate;
        } else {
          rIdx = _SelectRandomLight(rightNode->lights);
          rightEst = _EvalutateNode(rightNode->lights[rIdx], gpNode);
        }
#else
        rightEst = (ltNode->light == rightNode->light) ? cutItem._estimate : _EvalutateNode(rightNode->light, gpNode);
#endif

        Color rightBound = _ComputeUpperBound(rightNode, gpNode, rHeur) * gpNode->strength;
        lightcut.push_back(MdLightcutHeapItem(rightNode, gpNode, rIdx, rightBound, rightEst, rHeur, seq));
        push_heap(lightcut.begin(), lightcut.end());
        Color rL = rightNode->L * rightEst;

        totalEstL -= (ltNode->L * cutItem._estimate) * gpNode->strength;
        totalEstL += (lL + rL) * gpNode->strength;
        cutSize++;
      }
    }
    else if (!gpNode->IsLeaf())
    {
      if (node->GetNodeType() == ORIENTED_NODE)
      {
        MdOrientedLightTreeNode *ltNode = static_cast<MdOrientedLightTreeNode*>(node);
        Color leftEst(zero), rightEst(zero);
        bool lHeur, rHeur;
        uint32_t seq = cutItem._refineSeq <= 0 ? cutItem._refineSeq - 1 : 0;
#if REP_LIGHT_NUM > 1 
        lsOrientedLight* light = ltNode->lights[index];
#else
        lsOrientedLight* light = ltNode->light;
#endif

        GatherNode *leftNode = gpNode->left;
        leftEst = gpNode->gp == leftNode->gp ? cutItem._estimate : _EvalutateNode(light, leftNode);
        Color leftBound = _ComputeUpperBound(ltNode, leftNode, lHeur) * leftNode->strength;
        lightcut.push_back(MdLightcutHeapItem(ltNode, leftNode, index, leftBound, leftEst, lHeur, 0));
        push_heap(lightcut.begin(), lightcut.end());
        Color lL = ltNode->L * leftEst;

        GatherNode *rightNode = gpNode->right;
        rightEst = gpNode->gp == rightNode->gp ? cutItem._estimate : _EvalutateNode(light, rightNode);
        Color rightBound = _ComputeUpperBound(ltNode, rightNode, rHeur) * rightNode->strength;
        lightcut.push_back(MdLightcutHeapItem(ltNode, rightNode, index, rightBound, rightEst, rHeur, 0));
        push_heap(lightcut.begin(), lightcut.end());
        Color rL = ltNode->L * rightEst;

        totalEstL -= (ltNode->L * cutItem._estimate) * gpNode->strength;
        totalEstL += lL * leftNode->strength + rL * rightNode->strength;
        cutSize++;
      }
      else if (node->GetNodeType() == DIRECT_NODE)
      {
        MdDirectionalLightTreeNode *ltNode = static_cast<MdDirectionalLightTreeNode*>(node);
        Color leftEst(zero), rightEst(zero);
        bool lHeur, rHeur;
        uint32_t seq = cutItem._refineSeq <= 0 ? cutItem._refineSeq - 1 : 0; 
#if REP_LIGHT_NUM > 1 
        lsDirLight* light = ltNode->lights[index];
#else
        lsDirLight* light = ltNode->light;
#endif

        GatherNode *leftNode = gpNode->left;
        leftEst = gpNode->gp == leftNode->gp ? cutItem._estimate : _EvalutateNode(light, leftNode);
        Color leftBound = _ComputeUpperBound(ltNode, leftNode, lHeur) * leftNode->strength;
        lightcut.push_back(MdLightcutHeapItem(ltNode, leftNode, index, leftBound, leftEst, lHeur, 0));
        push_heap(lightcut.begin(), lightcut.end());
        Color lL = ltNode->L * leftEst;

        GatherNode *rightNode = gpNode->right;
        rightEst = gpNode->gp == rightNode->gp ? cutItem._estimate : _EvalutateNode(light, rightNode);
        Color rightBound = _ComputeUpperBound(ltNode, rightNode, rHeur) * rightNode->strength;
        lightcut.push_back(MdLightcutHeapItem(ltNode, rightNode, index, rightBound, rightEst, rHeur, 0));
        push_heap(lightcut.begin(), lightcut.end());
        Color rL = ltNode->L * rightEst;

        totalEstL -= (ltNode->L * cutItem._estimate) * gpNode->strength;
        totalEstL += lL * leftNode->strength + rL * rightNode->strength;
        cutSize++;
      }
    }
    else
      std::cout << "leaf error" << endl;
  }
  return gpRoot->emission + totalEstL;
}

Color MdLightcutter::_EvalutateNode(const lsOrientedLight* l, const GatherNode* gn )
{
  const mdGatherPoint* g = gn->gp;

  if (l == NULL)
    return Color(zero);
  if (g == 0)
    return Color(zero);

  const lsOrientedLight& light = *l;
  const mdGatherPoint& gp = *(g);

  const DifferentialGeometry &dp = gp.isect;
  Vector3f wi = normalize(light.position - dp.P);
  float maxDist = length(light.position - dp.P);
  float cosAngle = max(0.0f, dot(light.normal , -wi));
  float lenSqrEst = max(_clamp, squarelength(light.position - dp.P));
  // the Ns value was already set, so do not update it (restore it after updating)
  Vector3f ns = gp.isect.Ns;
  CompositedBRDF msu;
  if (gp.isect.material) gp.isect.material->shade(Ray(), Medium::Vacuum(), gp.isect, msu);
  gp.isect.Ns = ns;
  Color brdf = msu.eval(gp.wo, gp.isect, wi,ALL);

  Color estimate = brdf * cosAngle / lenSqrEst * gp.weight;

  if(zero!=estimate) 
  {
    STATS(stats.atomic.rays++);
    float rayerror = gp.isect.error*128.0f*std::numeric_limits<float>::epsilon();
    Ray shadowRay(gp.isect.P, wi, rayerror, maxDist-rayerror, 0.0f);
    rtcOccluded(_scene->scene,(RTCRay&)shadowRay);
    if(!shadowRay)
      return estimate;
  }
  return Color(zero);
}

Color MdLightcutter::_EvalutateNode(const lsDirLight* l, const GatherNode* gn )
{
  const mdGatherPoint* g = gn->gp;

  if (l == NULL)
    return Color(zero);
  if (g == 0)
    return Color(zero);

  const lsDirLight& light = *l;
  const mdGatherPoint& gp = *(g);

  const DifferentialGeometry &dp = gp.isect;
  Vector3f wi = -light.normal;
  CompositedBRDF msu;
  if (gp.isect.material) gp.isect.material->shade(Ray(), Medium::Vacuum(), gp.isect, msu);
  Color brdf = msu.eval(gp.wo, gp.isect,wi,ALL);
  Color estimate = brdf * gp.weight;

  if(zero!=estimate) 
  {
    float rayerror = gp.isect.error*128.0f*std::numeric_limits<float>::epsilon();
    Ray shadowRay(dp.P, wi, rayerror, inf, 0.0f);
    rtcOccluded(_scene->scene,(RTCRay&)shadowRay);
    if(!shadowRay)
      return estimate;
  }
  return Color(zero);
}

Color MdLightcutter::_ComputeUpperBound( MdOrientedLightTreeNode* ltNode, GatherNode *gpNode, bool &refineLight )
{
  if (ltNode->IsLeaf() && gpNode->IsLeaf())
    return Color(zero);

  const mdGatherPoint& gp = *gpNode->gp;
  if (!gp.isect.material) {
    refineLight = true;
    return Color(zero);
  }

  const DifferentialGeometry &dp = gp.isect;

  const BBox3f &ltBBox = ltNode->bbox;
  const BBox3f &gpBBox = gpNode->bbox;

  refineLight = true;
  float lenSqr = distanceSqr(ltBBox, gpBBox);

  float cosBound = 1.0f;
  float cosHalfAngle = ltNode->cone.GetAngleCos();
  if (cosHalfAngle <= 0.0f)
    cosBound = 1.0f;
  else
  {
    Vector3f vv = cross(Vector3f(0,0,1),ltNode->cone.GetAxis());
    if(vv == zero) vv = Vector3f(0,1,0); // proper rotation axis if crossprod was 0
    vv = normalize(vv);
    LinearSpace3f mt = LinearSpace3f::rotate(vv, -acosf(clamp(dot(Vector3f(0,0,1) , ltNode->cone.GetAxis()),-1.0f,1.0f)));


    BBox3f xbox(gpBBox.lower - ltBBox.upper, gpBBox.upper - ltBBox.lower);
    BBox3f tbox = xfmBBox(mt,xbox);

    Vector3f &xm = tbox.lower;
    Vector3f &xM = tbox.upper;

    float cosTheta;
    if (xM.z > 0)
    {
      float minx2 = (xm.x * xM.x <= 0) ? 0.0f : min(xm.x * xm.x, xM.x * xM.x);
      float miny2 = (xm.y * xM.y <= 0) ? 0.0f : min(xm.y * xm.y, xM.y * xM.y);
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

  Color MatBound = _BoundMaterialOrientedLight(ltNode, gpNode);
  Color upperBound = abs(MatBound * (cosBound / max(_clamp, lenSqr)));

  // compute heuristic
  if (conjoint(ltBBox, gpBBox))
  {
    refineLight = length(ltBBox.size()) > length(gpBBox.size());
  }
  else
  {
    float ltHeur = 0.0f;
    float gpHeur = 0.0f;
    // light heuristic
    if (!ltNode->IsLeaf())
    {
      float lmin = distanceSqr(ltNode->left->bbox, gpBBox);
      float rmin = distanceSqr(ltNode->right->bbox, gpBBox);

      float redDist = min(lmin, rmin) / lenSqr;
      float redMat = (1.0f - ltNode->cone.GetAngleCos()) * length(ltNode->bbox.size()) / sqrt(lenSqr);
      ltHeur = redDist * redMat;
    }

    if (!gpNode->IsLeaf())
    {
      float lmin = distanceSqr(ltBBox, gpNode->left->bbox);
      float rmin = distanceSqr(ltBBox, gpNode->right->bbox);

      float redDist = min(lmin, rmin) / lenSqr;
      float redMat = (1.0f - gpNode->normalCone.GetAngleCos()) * length(gpNode->bbox.size()) / sqrt(lenSqr);

      gpHeur = redDist * redMat;
    }
    refineLight = gpHeur <= ltHeur;
  }

  return abs(ltNode->L * upperBound);
}

Color MdLightcutter::_ComputeUpperBound( MdDirectionalLightTreeNode* ltNode, GatherNode *gpNode, bool &refineLight )
{
  if (ltNode->IsLeaf() && gpNode->IsLeaf())
    return Color(zero);

  const mdGatherPoint& gp = *gpNode->gp;
  if (!gp.isect.material) {
    refineLight = true;
    return Color(zero);
  }

  const BBox3f &ltBBox = ltNode->bbox;
  const BBox3f &gpBBox = gpNode->bbox;

  refineLight = true;
  float lenSqr = distanceSqr(ltBBox, gpBBox);

  Color MatBound = _BoundMaterialDirLight(ltNode, gpNode);
  Color upperBound = abs(MatBound);

  // compute heuristic
  if (conjoint(ltBBox, gpBBox))
  {
    refineLight = length(ltBBox.size()) > length(gpBBox.size());
  }
  else
  {
    float ltHeur = 0.0f;
    float gpHeur = 0.0f;
    // light heuristic
    if (!ltNode->IsLeaf())
    {
      float redMat = length(ltNode->bbox.size()) / sqrt(lenSqr);
      ltHeur = redMat;
    }

    if (!gpNode->IsLeaf())
    {
      float redMat = length(gpNode->bbox.size()) / sqrt(lenSqr);
      gpHeur = redMat;
    }
    refineLight = gpHeur <= ltHeur;
  }

  return abs(ltNode->L * upperBound);
}

Color MdLightcutter::_BoundMaterialOrientedLight( MdOrientedLightTreeNode* ltNode, GatherNode *gpNode)
{
  const mdGatherPoint& gp = *(gpNode->gp);

  const BBox3f &ltBBox = ltNode->bbox;
  const BBox3f &gpBBox = gpNode->bbox;
  float n(gpNode->n);
  Color kd(gpNode->kd),ks(gpNode->ks);
  

  BBox3f xbox(ltBBox.lower - gpBBox.upper, ltBBox.upper - gpBBox.lower);
  float cosBound = 1.0f;
  float cosHalfAngle = gpNode->normalCone.GetAngleCos();
  if (cosHalfAngle <= 0.0f)
    cosBound = 1.0f;
  else
  {
    Vector3f vv = cross(Vector3f(0,0,1), gpNode->normalCone.GetAxis());
    if(vv == zero) vv = Vector3f(0,1,0); // proper rotation axis if crossprod was 0
    vv = normalize(vv);
    LinearSpace3f mt = LinearSpace3f::rotate(vv, -acosf(clamp(dot(Vector3f(0,0,1) , gpNode->normalCone.GetAxis()),-1.0f,1.0f)));


    BBox3f tbox = xfmBBox(mt,xbox);

    Vector3f &xm = tbox.lower;
    Vector3f &xM = tbox.upper;

    float cosTheta;
    if (xM.z > 0)
    {
      float minx2 = (xm.x * xM.x <= 0) ? 0.0f : min(xm.x * xm.x, xM.x * xM.x);
      float miny2 = (xm.y * xM.y <= 0) ? 0.0f : min(xm.y * xm.y, xM.y * xM.y);
      float maxz2 = xM.z * xM.z;
      cosTheta = xM.z / sqrt(minx2 + miny2 + maxz2);
    }
    else
      cosTheta = 0.0f;

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

  if (zero!=ks && n != 0.0f)
  {
    //Phong material
    float glossyCosBound;
    float glossyCosTheta = 1.0f;
    Vector3f R = gpNode->mirrorCone.GetAxis();

    float glossyCosHalfAngle = gpNode->mirrorCone.GetAngleCos();
    if (glossyCosHalfAngle <= 0.0f)
      glossyCosBound = 1.0f;
    else
    {
      Vector3f gv = cross(Vector3f(0,0,1), R);
      if(gv == zero) gv = Vector3f(0,1,0); // proper rotation axis if crossprod was 0
      gv = normalize(gv);
      LinearSpace3f gmt = LinearSpace3f::rotate(gv, -acosf(clamp(dot(Vector3f(0,0,1) , R),-1.0f,1.0f)));

      BBox3f gbox = xfmBBox(gmt,xbox);

      Vector3f &gm = gbox.lower;
      Vector3f &gM = gbox.upper;

      if (gM.z > 0)
      {
        float minx2 = (gm.x * gM.x <= 0) ? 0.0f : min(gm.x * gm.x, gM.x * gM.x);
        float miny2 = (gm.y * gM.y <= 0) ? 0.0f : min(gm.y * gm.y, gM.y * gM.y);
        float maxz2 = gM.z * gM.z;
        glossyCosTheta = gM.z / sqrt(minx2 + miny2 + maxz2);
      }
      else
        glossyCosTheta = 0.0f;

      glossyCosTheta = clamp(glossyCosTheta, 0.0f, 1.0f);
      assert(glossyCosTheta <= 1.0f);

      if (glossyCosTheta > glossyCosHalfAngle)
        glossyCosBound = 1.0f;
      else
      {
        float sinHalfAngle = sqrt(1 - glossyCosHalfAngle * glossyCosHalfAngle);
        float sinTheta = sqrt(1 - glossyCosTheta * glossyCosTheta);
        glossyCosBound = clamp(glossyCosTheta * glossyCosHalfAngle + sinTheta * sinHalfAngle, 0.0f, 1.0f);
        assert(glossyCosBound >= 0.0f && glossyCosBound <= 1.0f);
      }
    }
    return (float(one_over_pi) * kd + (float(one_over_two_pi) * ks * (n+2) * pow(glossyCosBound, n))) * cosBound ;
  }
  return float(one_over_pi) * kd * cosBound;
}

Color MdLightcutter::_BoundMaterialDirLight( MdDirectionalLightTreeNode* ltNode, GatherNode *gpNode)
{
  const mdGatherPoint& gp = *(gpNode->gp);

  const BBox3f &ltBBox = ltNode->bbox;

  Color kd(gpNode->kd);
  Color ks(gpNode->ks);
  float n(gpNode->n);

  BBox3f xbox = BBox3f(EmptyTy());
  xbox.grow(-ltBBox.lower);
  xbox.grow(-ltBBox.upper);
  float cosBound = 1.0f;
  float cosHalfAngle = gpNode->normalCone.GetAngleCos();
  if (cosHalfAngle <= 0.0f)
    cosBound = 1.0f;
  else
  {
    Vector3f vv = cross(Vector3f(0,0,1), gpNode->normalCone.GetAxis());
    if(vv == zero) vv = Vector3f(0,1,0); // proper rotation axis if crossprod was 0
    vv = normalize(vv);
    LinearSpace3f mt = LinearSpace3f::rotate(vv, -acosf(clamp(dot(Vector3f(0,0,1) , gpNode->normalCone.GetAxis()),-1.0f,1.0f)));

    BBox3f tbox = xfmBBox(mt,xbox);

    Vector3f &xm = tbox.lower;
    Vector3f &xM = tbox.upper;

    float cosTheta;
    if (xM.z > 0)
    {
      float minx2 = (xm.x * xM.x <= 0) ? 0.0f : min(xm.x * xm.x, xM.x * xM.x);
      float miny2 = (xm.y * xM.y <= 0) ? 0.0f : min(xm.y * xm.y, xM.y * xM.y);
      float maxz2 = xM.z * xM.z;
      cosTheta = xM.z / sqrt(minx2 + miny2 + maxz2);
    }
    else
      cosTheta = 0.0f;

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

  if (zero!=ks && n != 0.0f)
  {
    //Phong material
    float glossyCosBound;
    float glossyCosTheta = 1.0f;
    Vector3f R = gpNode->mirrorCone.GetAxis();

    float glossyCosHalfAngle = gpNode->mirrorCone.GetAngleCos();
    if (glossyCosHalfAngle <= 0.0f)
      glossyCosBound = 1.0f;
    else
    {
      Vector3f gv = cross(Vector3f(0,0,1), R);
      if(gv == zero) gv = Vector3f(0,1,0); // proper rotation axis if crossprod was 0
      gv = normalize(gv);
      LinearSpace3f gmt = LinearSpace3f::rotate(gv, -acosf(clamp(dot(Vector3f(0,0,1) , R),-1.0f,1.0f)));

      BBox3f gbox = xfmBBox(gmt,xbox);

      Vector3f &gm = gbox.lower;
      Vector3f &gM = gbox.upper;

      if (gM.z > 0)
      {
        float minx2 = (gm.x * gM.x <= 0) ? 0.0f : min(gm.x * gm.x, gM.x * gM.x);
        float miny2 = (gm.y * gM.y <= 0) ? 0.0f : min(gm.y * gm.y, gM.y * gM.y);
        float maxz2 = gM.z * gM.z;
        glossyCosTheta = gM.z / sqrt(minx2 + miny2 + maxz2);
      }
      else
        glossyCosTheta = 0.0f;

      glossyCosTheta = clamp(glossyCosTheta, 0.0f, 1.0f);
      assert(glossyCosTheta <= 1.0f);

      if (glossyCosTheta > glossyCosHalfAngle)
        glossyCosBound = 1.0f;
      else
      {
        float sinHalfAngle = sqrt(1 - glossyCosHalfAngle * glossyCosHalfAngle);
        float sinTheta = sqrt(1 - glossyCosTheta * glossyCosTheta);
        glossyCosBound = clamp(glossyCosTheta * glossyCosHalfAngle + sinTheta * sinHalfAngle, 0.0f, 1.0f);
        assert(glossyCosBound >= 0.0f && glossyCosBound <= 1.0f);
      }
    }

    return (kd + (ks * pow(glossyCosBound, n))) * cosBound;
  }
  return kd * cosBound;
}
}
