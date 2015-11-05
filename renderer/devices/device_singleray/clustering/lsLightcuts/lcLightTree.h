#ifndef _LCLIGHT_TREE_H_
#define _LCLIGHT_TREE_H_

#include <vector>
#include "devices/device_singleray/default.h"
#include "devices/device_singleray/clustering/lightslice/lsFastCone.h"
#include "devices/device_singleray/clustering/lightslice/lsLightData.h"
#include "devices/device_singleray/clustering/lightslice/lsRandom.h"
#include "devices/device_singleray/clustering/icParameters.h"
//#define MULTI_REP

using std::vector;

namespace embree{

  enum LightTreeNodeType
  {
    ORIENTED_NODE,
    SPOT_NODE,
    DIRECT_NODE,
    POINT_NODE
  };

  struct LightTreeNode
  {
    LightTreeNode() {}
    virtual ~LightTreeNode() {}
    Color           L;
    virtual bool	IsLeaf() const = 0;
    virtual         LightTreeNodeType GetNodeType() const = 0;
  };

  struct OrientedLightTreeNode : public LightTreeNode
  {
    typedef lsOrientedLight lighttype;
    OrientedLightTreeNode() : LightTreeNode(), light(NULL), left(0), right(0) {}
    OrientedLightTreeNode(lsOrientedLight &light);
    virtual ~OrientedLightTreeNode();
    virtual LightTreeNodeType           GetNodeType() const { return ORIENTED_NODE; }
    bool								IsLeaf() const { return left == NULL && right == NULL; }
    FastConef		                    cone;
    BBox3f			                    bbox;
    lsOrientedLight                       *light;
    OrientedLightTreeNode				*left;
    OrientedLightTreeNode				*right;
  };

  struct SpotLightTreeNode : public LightTreeNode
  {
    typedef lsSpotLight lighttype;
    SpotLightTreeNode() : LightTreeNode(), light(NULL), left(0), right(0) {}
    SpotLightTreeNode(lsSpotLight &light);
    virtual ~SpotLightTreeNode();
    virtual LightTreeNodeType           GetNodeType() const { return SPOT_NODE; }
    bool								IsLeaf() const { return left == NULL && right == NULL; }
    FastConef		                    cone;
    BBox3f			                    bbox;
    lsSpotLight                       *light;
    SpotLightTreeNode				*left;
    SpotLightTreeNode				*right;
  };

  struct DirectionalLightTreeNode : public LightTreeNode
  {
    typedef lsDirLight lighttype;
    DirectionalLightTreeNode() : LightTreeNode(), light(NULL), left(0), right(0) {}
    DirectionalLightTreeNode(lsDirLight &light);
    virtual ~DirectionalLightTreeNode();
    virtual LightTreeNodeType           GetNodeType() const { return DIRECT_NODE; }
    bool								IsLeaf() const { return left == NULL && right == NULL; }
    BBox3f			                    bbox;
    lsDirLight                            *light;
    DirectionalLightTreeNode			*left;
    DirectionalLightTreeNode			*right;
  };

  struct OmniDirLightTreeNode : public LightTreeNode
  {
    typedef lsOmniDirLight lighttype;
    OmniDirLightTreeNode() : LightTreeNode(), left(0), right(0) {}
    OmniDirLightTreeNode(lsOmniDirLight &light);
    virtual ~OmniDirLightTreeNode();
    virtual LightTreeNodeType           GetNodeType() const { return POINT_NODE; }
    bool								IsLeaf() const { return left == NULL && right == NULL; }
    BBox3f			                    bbox;
    lsOmniDirLight                        *light;
    OmniDirLightTreeNode				*left;
    OmniDirLightTreeNode				*right;
  };

#if REP_LIGHT_NUM > 1 
  const unsigned int MD_REP_SLOTS = 9;
  struct MdOrientedLightTreeNode : public LightTreeNode
  {
    typedef lsOrientedLight lighttype;
    MdOrientedLightTreeNode() : LightTreeNode(), left(0), right(0) { memset(lights, 0, sizeof(lights));}
    MdOrientedLightTreeNode(lsOrientedLight &light, uint32_t idx);
    virtual ~MdOrientedLightTreeNode();
    virtual LightTreeNodeType           GetNodeType() const { return ORIENTED_NODE; }
    bool								IsLeaf() const { return left == NULL && right == NULL; }
    FastConef		                    cone;
    BBox3f			                    bbox;
    lsOrientedLight                       *lights[MD_REP_SLOTS];
    MdOrientedLightTreeNode				*left;
    MdOrientedLightTreeNode				*right;
  };

  struct MdSpotLightTreeNode : public LightTreeNode
  {
    typedef lsSpotLight lighttype;
    MdSpotLightTreeNode() : LightTreeNode(), left(0), right(0) {memset(lights,0,sizeof(lights));}
    MdSpotLightTreeNode(lsSpotLight &light,uint32_t idx);
    virtual ~MdSpotLightTreeNode();
    virtual LightTreeNodeType           GetNodeType() const { return SPOT_NODE; }
    bool								IsLeaf() const { return left == NULL && right == NULL; }
    FastConef		                    cone;
    BBox3f			                    bbox;
    lsSpotLight                       *lights[MD_REP_SLOTS];
    SpotLightTreeNode				*left;
    SpotLightTreeNode				*right;
  };

  struct MdDirectionalLightTreeNode : public LightTreeNode
  {
    typedef lsDirLight lighttype;
    MdDirectionalLightTreeNode() : LightTreeNode(), left(0), right(0) { memset(lights, 0, sizeof(lights)); }
    MdDirectionalLightTreeNode(lsDirLight &light, uint32_t idx);
    virtual ~MdDirectionalLightTreeNode();
    virtual LightTreeNodeType           GetNodeType() const { return DIRECT_NODE; }
    bool								IsLeaf() const { return left == NULL && right == NULL; }
    BBox3f			                    bbox;
    lsDirLight                            *lights[MD_REP_SLOTS];
    MdDirectionalLightTreeNode			*left;
    MdDirectionalLightTreeNode			*right;
  };

  struct MdOmniDirLightTreeNode : public LightTreeNode
  {
    typedef lsOmniDirLight lighttype;
    MdOmniDirLightTreeNode() : LightTreeNode(), left(0), right(0) { memset(lights, 0, sizeof(lights)); }
    MdOmniDirLightTreeNode(lsOmniDirLight &light, uint32_t idx);
    virtual ~MdOmniDirLightTreeNode();
    virtual LightTreeNodeType           GetNodeType() const { return POINT_NODE; }
    bool								IsLeaf() const { return left == NULL && right == NULL; }
    BBox3f			                    bbox;
    lsOmniDirLight                        *lights[MD_REP_SLOTS];
    MdOmniDirLightTreeNode				*left;
    MdOmniDirLightTreeNode				*right;
  };
#else
  const unsigned int MD_REP_SLOTS = 1;  // resolve compile error

  typedef OrientedLightTreeNode MdOrientedLightTreeNode;
  typedef DirectionalLightTreeNode MdDirectionalLightTreeNode;
  typedef OmniDirLightTreeNode MdOmniDirLightTreeNode;
#endif

  template<typename OT, typename DT, typename PT,typename ST>
  class GenericLightTree
  {
    //friend class DivisiveMdLightTreeBuilder;
    friend class DivisiveLightTreeBuilder;
  public:
    GenericLightTree(void):_orientedRoot(NULL),_pointRoot(NULL),_spotRoot(NULL) {}
    ~GenericLightTree(void) { Clear(); }
    inline OT				*GetOrientedRoot() const { return _orientedRoot; }
    inline DT				*GetDirectionalRoot() const { return _directionalRoot; }
    inline ST				*GetSpotRoot() const { return _spotRoot; }

    inline uint32_t         OrientedLightNum() const { return static_cast<uint32_t>(_orientedLights.size()); }
    inline uint32_t         DirectionalLightNum() const { return static_cast<uint32_t>(_dirLights.size()); }
    inline uint32_t         OmniDirLightNum() const { return static_cast<uint32_t>(_omniDirLights.size()); }
    inline uint32_t         SpotLightNum() const { return static_cast<uint32_t>(_spotLights.size()); }

    vector<lsOmniDirLight>&   OmniDirLights() { return _omniDirLights; }
    vector<lsDirLight>&       DirectionalLights() { return _dirLights; }
    vector<lsOrientedLight>&  OrientedLights() { return _orientedLights; }
    vector<lsSpotLight>&      SpotLights() { return _spotLights; }

    void                    Clear() { 
      if (_orientedRoot) delete _orientedRoot; 
      if (_directionalRoot) delete _directionalRoot;
      if (_spotRoot) delete _spotRoot;
    }

  public:
    vector<lsOmniDirLight>            _omniDirLights;
    vector<lsDirLight>                _dirLights;
    vector<lsOrientedLight>           _orientedLights;
    vector<lsSpotLight>               _spotLights;

    OT								*_orientedRoot;
    DT								*_directionalRoot;
    PT								*_pointRoot;
    ST								*_spotRoot;
  };


  typedef GenericLightTree<OrientedLightTreeNode, DirectionalLightTreeNode, OmniDirLightTreeNode,SpotLightTreeNode> lcLightTree;
#if REP_LIGHT_NUM > 1 
  typedef GenericLightTree<MdOrientedLightTreeNode, MdDirectionalLightTreeNode, MdOmniDirLightTreeNode,MdSpotLightTreeNode> MdLightTree;
#else
  typedef lcLightTree MdLightTree;
#endif // MULTI_REP

}

#endif // _LIGHT_TREE_H_
