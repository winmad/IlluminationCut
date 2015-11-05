#ifndef MDLC_GATHER_TREE_H_
#define MDLC_GATHER_TREE_H_
#include "devices/device_singleray/default.h"
#include "devices/device_singleray/clustering/lightslice/lsFastCone.h"
#include "devices/device_singleray/clustering/lsLightcuts/lcLightTree.h"
#include "devices/device_singleray/clustering/lightslice/lsRandom.h"
#include "devices/device_singleray/clustering/icParameters.h"
#include "devices/device_singleray/shapes/differentialgeometry.h"
//#include "renderer/brdfs/compositedbrdf.h"
//#include "image/cubemap.h"
#include <vector>

namespace embree{
  struct mdGatherPoint
  {
    mdGatherPoint():emission(zero),weight(zero){}
    Color           emission;
    Color           weight;
    Vec2i           pixel;
    DifferentialGeometry    isect;
    Vector3f           wo;
  };

  struct GatherNode
  {
    GatherNode() : left(NULL), right(NULL), gp(NULL), strength(0.0f), emission(zero), kd(zero), ks(zero), n(0.0f){}
    ~GatherNode() {if(left) delete left; if(right) delete right;}
    inline bool             IsLeaf() const { return left == NULL && right == NULL; }
    GatherNode              *left;
    GatherNode              *right;
    mdGatherPoint             *gp;

    BBox3f                 bbox;
    FastConef               normalCone;
    FastConef               mirrorCone;
    bool                    hasGlossy;

    float                   strength;
    Color                   emission;

    Color    kd;
    Color    ks;
    float    n;
  };

}

#endif // _GATHER_TREE_H_
