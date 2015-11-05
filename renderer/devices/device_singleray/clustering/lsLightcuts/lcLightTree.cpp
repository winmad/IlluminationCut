#include "devices/device_singleray/clustering/lsLightcuts/lcLightTree.h"
#include "devices/device_singleray/clustering/lightslice/lsRandom.h"

namespace embree{
OrientedLightTreeNode::OrientedLightTreeNode(lsOrientedLight &l) : left(0), right(0)
{
  light = &l;
  L = l.le;
  bbox = BBox3f(l.position, l.position);
  cone = FastConef(l.normal);
}


SpotLightTreeNode::SpotLightTreeNode(lsSpotLight &l) : left(0), right(0)
{
  light = &l;
  L = l.le;
  bbox = BBox3f(l.position, l.position);
  cone = FastConef(l.normal);
}

DirectionalLightTreeNode::DirectionalLightTreeNode(lsDirLight &l) : left(0), right(0)
{
  light = &l;
  L = l.le;
  bbox = BBox3f(l.normal, l.normal);
}

OmniDirLightTreeNode::OmniDirLightTreeNode(lsOmniDirLight &l) : left(0), right(0)
{
  light = &l;
  L = l.intensity;
  bbox = BBox3f(l.position, l.position);
}

OrientedLightTreeNode::~OrientedLightTreeNode()
{
  if(left) delete left;
  if(right) delete right;
}


SpotLightTreeNode::~SpotLightTreeNode()
{
  if(left) delete left;
  if(right) delete right;
}

DirectionalLightTreeNode::~DirectionalLightTreeNode()
{
  if(left) delete left;
  if(right) delete right;
}

OmniDirLightTreeNode::~OmniDirLightTreeNode()
{
  if(left) delete left;
  if(right) delete right;
}

#if REP_LIGHT_NUM > 1 
MdOrientedLightTreeNode::MdOrientedLightTreeNode(lsOrientedLight &l, uint32_t idx) : left(0), right(0)
{
  memset(lights, 0, sizeof(lights));
  lights[idx] = &l;
  L = l.le;
  bbox = BBox3f(l.position, l.position);
  cone = FastConef(l.normal);
}

MdDirectionalLightTreeNode::MdDirectionalLightTreeNode(lsDirLight &l, uint32_t idx) : left(0), right(0)
{
  memset(lights, 0, sizeof(lights));
  lights[idx] = &l;
  L = l.le;
  bbox = BBox3f(l.normal, l.normal);
}

MdOmniDirLightTreeNode::MdOmniDirLightTreeNode(lsOmniDirLight &l, uint32_t idx) : left(0), right(0)
{
  memset(lights, 0, sizeof(lights));
  lights[idx] = &l;
  L = l.intensity;
  bbox = BBox3f(l.position, l.position);
}

MdOrientedLightTreeNode::~MdOrientedLightTreeNode()
{
  if(left) delete left;
  if(right) delete right;
}

MdDirectionalLightTreeNode::~MdDirectionalLightTreeNode()
{
  if(left) delete left;
  if(right) delete right;
}

MdOmniDirLightTreeNode::~MdOmniDirLightTreeNode()
{
  if(left) delete left;
  if(right) delete right;
}
#endif

}