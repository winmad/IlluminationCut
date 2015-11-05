#ifndef _LIGHT_TREE_BUILDER_UTIL_H_
#define _LIGHT_TREE_BUILDER_UTIL_H_

#include "devices/device_singleray/default.h"
#include "devices/device_singleray/clustering/lsLightcuts/lcLightTree.h"

using namespace embree;

namespace LightTreeBuilderUtil
{
  __forceinline Vec6f GetPoint(OrientedLightTreeNode *node, float c) {
    Vector3f center = embree::center(node->bbox);
    Vector3f direction = node->cone.GetAxis();
    return Vec6f(center, direction * c);
  }

  __forceinline Vec6f GetPoint(SpotLightTreeNode *node, float c) {
    Vector3f center = embree::center(node->bbox);
    Vector3f direction = node->cone.GetAxis();
    return Vec6f(center, direction * c);
  }

  __forceinline Vec6f GetPoint(DirectionalLightTreeNode *node, float c) {
    return Vec6f(center(node->bbox), Vector3f(zero));
  }

  __forceinline Vec6f GetPoint(OmniDirLightTreeNode *node, float c) {
    return Vec6f(center(node->bbox), Vector3f(zero));
  }

#if REP_LIGHT_NUM > 1 
  inline Vec6f GetPoint(MdOrientedLightTreeNode *node, float c) {
    Vector3f center = embree::center(node->bbox);
    Vector3f direction = node->cone.GetAxis();
    return Vec6f(center, direction * c);
  }

  inline Vec6f GetPoint(MdDirectionalLightTreeNode *node, float c) {
    return Vec6f(center(node->bbox), Vector3f(zero));
  }

  inline Vec6f GetPoint(MdOmniDirLightTreeNode *node, float c) {
    return Vec6f(center(node->bbox), Vector3f(zero));
  }
#endif // MULTI_REP

  template<typename NodeType>
  struct LightTreeBuildItem
  {
    LightTreeBuildItem() : node(NULL) {}
    LightTreeBuildItem(NodeType *n, float c)
    {
      assert(n);
      node = n;
      point = GetPoint(n, c);
    }
    NodeType                                *node;
    Vec6f									point;
  };


  inline float Distance(const OrientedLightTreeNode *o1, const OrientedLightTreeNode *o2, float c)
  {
    BBox3f bbox = merge(o1->bbox, o2->bbox);
    FastConef cone = FastConef::Union(o1->cone, o2->cone);
    float a2 = squarelength(bbox.size());
    float sinBeta = cone.GetAngleSin();
    float d = c * sinBeta;
    Color le = o1->L + o2->L;
    return (reduce_add(le)/3.0f * (a2 + d * d));
  }

  inline float Distance(const DirectionalLightTreeNode *o1, const DirectionalLightTreeNode *o2, float c)
  {
    BBox3f bbox = merge(o1->bbox, o2->bbox);
    float a2 = squarelength(bbox.size());
    Color le = o1->L + o2->L;
    return reduce_add(le)/3.0f * a2;
  }

  inline float Distance(const OmniDirLightTreeNode *o1, const OmniDirLightTreeNode *o2, float c)
  {
    BBox3f bbox = merge(o1->bbox, o2->bbox);
    float a2 = squarelength(bbox.size());
    Color le = o1->L + o2->L;
    return reduce_add(le)/3.0f * a2;
  }
};
#endif // _LIGHT_TREE_BUILDER_UTIL_H_
