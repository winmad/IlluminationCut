#ifndef _DIVISIVE_LIGHT_TREE_BUILDER_H_
#define _DIVISIVE_LIGHT_TREE_BUILDER_H_
#include "devices/device_singleray/clustering/lsLightcuts/lcLightTreeBuilder.h"
#include "devices/device_singleray/clustering/lightslice/lsPathSampler.h"

using namespace LightTreeBuilderUtil;

namespace embree {

class DivisiveLightTreeBuilder : public LightTreeBuilder 
{
public:
  DivisiveLightTreeBuilder() {};
  ~DivisiveLightTreeBuilder(void) {};
  lcLightTree* Build(BackendScene* scene, std::vector<vplVPL>& indirect);
private:
  void        _ConstructTree(lcLightTree* lightTree);
  template<typename LightType, typename NodeType>
  NodeType*   _ConstructSubTree(vector<LightType> &lights);
  template<typename NodeType>
  NodeType*   _BuildSubTree(typename vector<LightTreeBuildItem<NodeType> >::iterator start, typename vector<LightTreeBuildItem<NodeType> >::iterator end);
  template<typename NodeType>
  NodeType*   _MakeLeaf( typename vector<LightTreeBuildItem<NodeType> >::iterator start ) { return start->node; }
  template<typename NodeType>
  void		_UpdateNode(NodeType *node);
  RandomPathSamplerStd	sampler;
};


template<typename NodeType>
NodeType* DivisiveLightTreeBuilder::_BuildSubTree(typename vector<LightTreeBuildItem<NodeType> >::iterator start, typename vector<LightTreeBuildItem<NodeType> >::iterator end)
{
  typedef typename vector<LightTreeBuildItem<NodeType> >::iterator DataIter;
  assert(end > start);
  if (end - start == 1)
    return _MakeLeaf<NodeType>(start);

  BBox<Vec6f> bbox = BBox<Vec6f>(EmptyTy());
  for (DataIter it = start; it != end; it++)
    bbox.grow(it->point);

  uint32_t dim = maxDim(bbox.size());
  DataIter mid = start + (end - start) / 2;
  auto pred =  [dim](const LightTreeBuildItem<NodeType> &d1, const LightTreeBuildItem<NodeType> &d2)->bool {
    return d1.point[dim] == d2.point[dim] ? (d1.node < d2.node) : d1.point[dim] < d2.point[dim]; 
  };
  std::nth_element(start, mid, end, pred );

  NodeType *node = new NodeType();
  node->left = _BuildSubTree<NodeType>(start, mid);
  node->right = _BuildSubTree<NodeType>(mid, end);
  _UpdateNode(node);
  return node;
}

template<typename LightType, typename NodeType>
NodeType* DivisiveLightTreeBuilder::_ConstructSubTree( vector<LightType> &lights)
{
  if (lights.size() < 1)
    return NULL;
  vector<LightTreeBuildItem<NodeType> > inputData;
  for (uint32_t i = 0; i < lights.size(); i++)
  {
    NodeType *node = new NodeType(lights[i]);
    inputData.push_back(LightTreeBuildItem<NodeType> (node, _c));
  }

  NodeType *node = _BuildSubTree<NodeType>(inputData.begin(), inputData.end());
  return node;
}


class DivisiveMdLightTreeBuilder : public LightTreeBuilder 
{
public:
  DivisiveMdLightTreeBuilder() { };
  ~DivisiveMdLightTreeBuilder(void) {};
  MdLightTree* Build(BackendScene* scene, std::vector<vplVPL>& indirect);
private:
  void        _ConstructTree(MdLightTree* lightTree);
  template<typename LightType, typename NodeType>
  NodeType*   _ConstructSubTree(vector<LightType> &lights);
  template<typename NodeType>
  NodeType*   _BuildSubTree(typename vector<LightTreeBuildItem<NodeType> >::iterator start, typename vector<LightTreeBuildItem<NodeType> >::iterator end);
  template<typename NodeType>
  NodeType*   _MakeLeaf( typename vector<LightTreeBuildItem<NodeType> >::iterator start ) { return start->node; }
  template<typename NodeType>
  void		_UpdateNode(NodeType *node);
  RandomPathSamplerStd	sampler;
};


template<typename LightType, typename NodeType>
NodeType* DivisiveMdLightTreeBuilder::_ConstructSubTree( vector<LightType> &lights )
{
  if (lights.size() < 1)
    return NULL;

  vector<LightTreeBuildItem<NodeType> > inputData;
  for (uint32_t i = 0; i < lights.size(); i++)
  {
#if REP_LIGHT_NUM > 1 
    uint32_t index = static_cast<uint32_t>(MD_REP_SLOTS * sampler.Next1D());
    NodeType *node = new NodeType(lights[i], index);
#else
    NodeType *node = new NodeType(lights[i]);
#endif
    inputData.push_back(LightTreeBuildItem<NodeType> (node, _c));
  }

  NodeType *node = _BuildSubTree<NodeType>(inputData.begin(), inputData.end());
  return node;
}

template<typename NodeType>
NodeType* DivisiveMdLightTreeBuilder::_BuildSubTree(typename vector<LightTreeBuildItem<NodeType> >::iterator start, typename vector<LightTreeBuildItem<NodeType> >::iterator end)
{
  typedef typename vector<LightTreeBuildItem<NodeType> >::iterator DataIter;
  assert(end > start);
  if (end - start == 1)
    return _MakeLeaf<NodeType>(start);

  BBox<Vec6f> bbox = BBox<Vec6f>(EmptyTy());
  for (DataIter it = start; it != end; it++)
    bbox.grow(it->point);

  uint32_t dim = maxDim(bbox.size());
  DataIter mid = start + (end - start) / 2;
  std::nth_element(start, mid, end, 
    [dim](const LightTreeBuildItem<NodeType> &d1, const LightTreeBuildItem<NodeType> &d2)->bool {
      return d1.point[dim] == d2.point[dim] ? (d1.node < d2.node) : d1.point[dim] < d2.point[dim]; 
  });

  NodeType *node = new NodeType();
  node->left = _BuildSubTree<NodeType>(start, mid);
  node->right = _BuildSubTree<NodeType>(mid, end);
  _UpdateNode(node);
  return node;
}
}
#endif // _DIVISIVE_LIGHT_TREE_BUILDER_H_
