#include "devices/device_singleray//clustering/lsLightcuts/lcDivisiveLightTreeBuilder.h"

namespace embree{

lcLightTree* DivisiveLightTreeBuilder::Build(BackendScene *scene, std::vector<vplVPL>& indirect )
{
  lcLightTree *lightTree = new lcLightTree();
  // add lights to lighttree
  for (size_t i = 0 ; i < indirect.size() ; ++i)
  {
    vplVPL& vl = indirect[i];
    lightTree->_orientedLights.push_back(lsOrientedLight(vl.P, - vl._D, vl.I));
  }
  this->_SampleLights(lightTree,scene,indirect.size());
  lightTree->_orientedRoot = _ConstructSubTree<lsOrientedLight, OrientedLightTreeNode>(lightTree->_orientedLights);
  lightTree->_directionalRoot = _ConstructSubTree<lsDirLight, DirectionalLightTreeNode>(lightTree->_dirLights);
  lightTree->_spotRoot = _ConstructSubTree<lsSpotLight, SpotLightTreeNode>(lightTree->_spotLights);
  return lightTree;
}

template<>
void DivisiveLightTreeBuilder::_UpdateNode(OrientedLightTreeNode *node)
{
  OrientedLightTreeNode *l = node->left;
  OrientedLightTreeNode *r = node->right;
  node->bbox = merge(l->bbox, r->bbox);
  node->cone = FastConef::Union(l->cone, r->cone);
  node->L = l->L + r->L;
  node->light = (sampler.Next1D() < (reduce_avg(l->L) / reduce_avg(node->L))) ? l->light : r->light;
}

template<>
void DivisiveLightTreeBuilder::_UpdateNode(SpotLightTreeNode *node)
{
  SpotLightTreeNode *l = node->left;
  SpotLightTreeNode *r = node->right;
  node->bbox = merge(l->bbox, r->bbox);
  node->cone = FastConef::Union(l->cone, r->cone);
  node->L = l->L + r->L;
  node->light = (sampler.Next1D() < (reduce_avg(l->L) / reduce_avg(node->L))) ? l->light : r->light;
}

template<>
void DivisiveLightTreeBuilder::_UpdateNode(DirectionalLightTreeNode *node)
{
  DirectionalLightTreeNode *l = node->left;
  DirectionalLightTreeNode *r = node->right;
  node->bbox = merge(l->bbox, r->bbox);
  node->L = l->L + r->L;
  node->light = (sampler.Next1D() < (reduce_avg(l->L) / reduce_avg(node->L))) ? l->light : r->light;
}

template<>
void DivisiveLightTreeBuilder::_UpdateNode(OmniDirLightTreeNode *node)
{
  OmniDirLightTreeNode *l = node->left;
  OmniDirLightTreeNode *r = node->right;
  node->bbox = merge(l->bbox, r->bbox);
  node->L = l->L + r->L;
  node->light = (sampler.Next1D() < (reduce_avg(l->L) / reduce_avg(node->L))) ? l->light : r->light;
}


MdLightTree* DivisiveMdLightTreeBuilder::Build(BackendScene *scene, std::vector<vplVPL>& indirect )
{
  MdLightTree *lightTree = new MdLightTree();
  for (size_t i = 0 ; i < indirect.size() ; ++i)
  {
    vplVPL& vl = indirect[i];
    lightTree->_orientedLights.push_back(lsOrientedLight(vl.P, - vl._D, vl.I));
  }
  this->_SampleLights(lightTree,scene,indirect.size());
  lightTree->_orientedRoot = _ConstructSubTree<lsOrientedLight, MdOrientedLightTreeNode>(lightTree->_orientedLights);
  lightTree->_directionalRoot = _ConstructSubTree<lsDirLight, MdDirectionalLightTreeNode>(lightTree->_dirLights);
  return lightTree;
}


template<>
void DivisiveMdLightTreeBuilder::_UpdateNode(MdOrientedLightTreeNode *node)
{
  MdOrientedLightTreeNode *l = node->left;
  MdOrientedLightTreeNode *r = node->right;
  node->bbox = merge(l->bbox, r->bbox);
  node->cone = FastConef::Union(l->cone, r->cone);
  node->L = l->L + r->L;
  float prob = reduce_avg(l->L) / reduce_avg(node->L);
#if REP_LIGHT_NUM > 1 
  for (uint32_t i = 0; i < MD_REP_SLOTS; i++)
  {
    if(l->lights[i] && r->lights[i])
      node->lights[i] = (sampler.Next1D() < prob) ? l->lights[i] : r->lights[i];
    else
      node->lights[i] = l->lights[i] ? l->lights[i] : r->lights[i];
  }
#else
  node->light = (sampler.Next1D() < prob) ? l->light : r->light;
#endif // MULTI_REP

}

template<>
void DivisiveMdLightTreeBuilder::_UpdateNode(MdDirectionalLightTreeNode *node)
{
  MdDirectionalLightTreeNode *l = node->left;
  MdDirectionalLightTreeNode *r = node->right;
  node->bbox = merge(l->bbox, r->bbox);
  node->L = l->L + r->L;
  float prob = reduce_avg(l->L) / reduce_avg(node->L);
#if REP_LIGHT_NUM > 1 
  for (uint32_t i = 0; i < MD_REP_SLOTS; i++)
  {
    if(l->lights[i] && r->lights[i])
      node->lights[i] = (sampler.Next1D() < prob) ? l->lights[i] : r->lights[i];
    else
      node->lights[i] = l->lights[i] ? l->lights[i] : r->lights[i];
  }
#else
  node->light = (sampler.Next1D() < prob) ? l->light : r->light;
#endif
}

template<>
void DivisiveMdLightTreeBuilder::_UpdateNode(MdOmniDirLightTreeNode *node)
{
  MdOmniDirLightTreeNode *l = node->left;
  MdOmniDirLightTreeNode *r = node->right;
  node->bbox = merge(l->bbox, r->bbox);
  node->L = l->L + r->L;
  float prob = reduce_avg(l->L) / reduce_avg(node->L);
#if REP_LIGHT_NUM > 1 
  for (uint32_t i = 0; i < MD_REP_SLOTS; i++)
  {
    if(l->lights[i] && r->lights[i])
      node->lights[i] = (sampler.Next1D() < prob) ? l->lights[i] : r->lights[i];
    else
      node->lights[i] = l->lights[i] ? l->lights[i] : r->lights[i];
  }
#else
  node->light = (sampler.Next1D() < prob) ? l->light : r->light;
#endif
}

}
