/*
*  This program is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
*  Release date: 20/05/2015
*  Authors: Norbert Bus & Nabil Mustafa*  
*  If you find this code useful, please cite it by referencing the following paper:
*  @article{BMB15, title={Illumination{Cut}}, author={Bus, Norbert and Mustafa, Nabil H. and Biri, Venceslas},
*           journal = {Computer Graphics Forum (Proceedings of Eurographics 2015)},year={2015}}
*/
#ifndef icLighttree_h__
#define icLighttree_h__

#include "devices/device_singleray/default.h"
#include "devices/device_singleray/clustering/icParameters.h"
#include "devices/device_singleray/clustering/vpl/vplVPL.h"
#include "common/math/vec6.h"
#include "devices/device_singleray/clustering/lightcuts/ctKdTree.h"
#include <time.h>
#include <vector>
#include <stack>

namespace embree
{

  typedef BBox<Vec6f> BBox6f;
  class icLighttree;

  
  // Cone with axis and half apex angle + rotation to (0,0,1)
  class icCone
  {
  public:


    /*! \brief Default constructor.  */
    icCone(){}


    /*! \brief Constructor.  */
    icCone(Vector3f dir,float halfAngle=0);


    /*! \brief Merge two cones.  */
    __forceinline icCone mergeCones( icCone& other ) const
    {
      float angle = acos(dot(this->axis,other.axis));

      // check if one cone contains the other
      if (abs(this->halfAngle - other.halfAngle ) > angle ) 
        return this->halfAngle > other.halfAngle ? *this : other;
      else
      {
        // define new cone containing both cones
        float newHalfAngle = 0.5f * (this->halfAngle + other.halfAngle + angle );
        newHalfAngle = clamp(newHalfAngle,0.0f,float(pi));
        Vector3f newDir;
        // handle extreme cases when crossprod will return IND values
        float len = length(this->axis - other.axis) ;    
        if (len < 0.001  )
        {
          newDir = this->axis; 
        }
        else if (len > 1.999 )
        {
          newDir = this->axis;
          newHalfAngle = float(pi);
        }
        else
        {
          LinearSpace3f rot = LinearSpace3f::rotate(cross(this->axis,other.axis),
                              newHalfAngle - this->halfAngle); 
          newDir = rot * (this->axis);
        }

        return icCone(newDir,newHalfAngle);
      }
    }


  public:
    Vector3f axis;         // representing the axis of a cone
    float halfAngle;       // angle between an half line and the axis
    float cosHalfAngle;    // angle between an half line and the axis
    LinearSpace3f rotation;   // rotate the axis to the z axis
  };


  /*! Node in the lightree. */
  class icLightNode 
  {
  public:
    friend class icLighttree;
    enum { nChildren = 2 };
    typedef BBox6f BB;


    /*! Constructor. */
    icLightNode();


    /*! Constructor. */
    icLightNode(vplVPL light[REP_LIGHT_NUM],BBox6f bb,icCone cone, 
        bool right[REP_LIGHT_NUM], icLightNode* lc=NULL, icLightNode* rc=NULL);


    /*! Destructor. */
    ~icLightNode();


    /*! Subtree size. */
    size_t subTreeSize();


    /*! Subtree depth. */
    size_t subTreeDepth();


    /*! \brief Distance between two nodes for agglomerative clustering. Custom metric from Lightcuts.  */
    __forceinline float dist( icLightNode* node, float sceneradiusSqr)
    {
      if (this==node)
        return FLT_MAX;

      Color I = this->repLights[0].I;
      I += node->repLights[0].I;
      float diagSqr = squarelength(size(embree::merge(this->cell,node->cell)).pos()); 
      float coeffSqr = sceneradiusSqr;
      icCone mergedCones =  this->cone.mergeCones(node->cone);
      float directionTermSqr = (1-cos(mergedCones.halfAngle))*( 1-cos(mergedCones.halfAngle));
      return reduce_add(I) * (diagSqr + coeffSqr * directionTermSqr); 
    }


    /*! Merge two nodes */
    icLightNode* merge(icLightNode* node, Random& rnd);


    /*! Calculate the dimension and the position of the next subdividing plane */
    void getNextDividingPlanePosition(float& pos, size_t& dim);


    /*! Create bbox for the i-th child in the octree. */
    BBox6f createBboxForChild( size_t i);


    /*! Upper bound the cos for a bb from (0,0,1). */
    __forceinline float CosUpperBoundForBB(Vector3f& lower, Vector3f& upper) const {
      // calculating cos upper bound supposing we are looking from (0,0,1)
      float cosTheta;
      if (upper.z > 0)
      {
        float minx2, miny2;
        if (lower.x * upper.x <= 0)
          minx2 = 0.0f;
        else
          minx2 = min(lower.x * lower.x, upper.x * upper.x);

        if (lower.y * upper.y <= 0)
          miny2 = 0.0f;
        else
          miny2 = min(lower.y * lower.y, upper.y * upper.y);
        float maxz2 = upper.z * upper.z;
        cosTheta = upper.z / sqrt(minx2 + miny2 + maxz2);
      }
      else
        cosTheta = 0.0f;

      assert(cosTheta >= 0.0f && cosTheta <= 1.0f);

      cosTheta = clamp(cosTheta, 0.0f, 1.0f);
      return cosTheta;
    }

    __forceinline Color upperBoundIntensity( const BBox6f & bb ) const
    {
      // transform the bb such that the cone axis would be (0,0,1)
      BBox3f b = BBox3f(bb.lower.pos() - this->cell.upper.pos(),bb.upper.pos() - this->cell.lower.pos());
      Vector3f corners[8];
      b.GetCorners(corners);
      BBox3f newbox = BBox3f(EmptyTy());
      for (int i = 0; i < 8; i++)
      {
        newbox.grow(xfmPoint(this->cone.rotation, corners[i]));
      }

      // get cos upper bound in the transformed setup
      float cosTheta = this->CosUpperBoundForBB(newbox.lower , newbox.upper);
      if (this->cone.cosHalfAngle <= 0.0f)
        cosTheta = 1.0f;
      else if (cosTheta > this->cone.cosHalfAngle)
        cosTheta = 1.0f;
      else
      {
        float sinHalfAngle = sqrt(1 - this->cone.cosHalfAngle * this->cone.cosHalfAngle);
        float sinTheta = sqrt(1 - cosTheta * cosTheta);
        cosTheta = clamp(cosTheta * this->cone.cosHalfAngle + sinTheta * sinHalfAngle);
      }
      // FIXME: regarding position the replights differ... 
      if (conjoint(this->repLights[0].P,BBox3f(bb.lower.pos(),bb.upper.pos())))
        cosTheta = 1.0f;
      return this->repLights[0].I * cosTheta;
    } 

  public:

    vplVPL repLights[REP_LIGHT_NUM];
    bool repLightFromRightChild[REP_LIGHT_NUM];
    BBox6f cell;
    icCone cone;
    KdNode<icLightNode>* kdnode;
    icLightNode *children[nChildren], *parent;
    bool valid;
    float posRadius;
    int level;
  };



  class icLighttree: public RefCount 
  {
  public:

    enum{nChildren=icLightNode::nChildren};


    /*! Range in an array of vpls corresponding to a node in the tree. */
    class NodeRange
    {
    public:
      NodeRange(icLightNode* node,vplVPL* start, vplVPL* end):node(node),start(start),end(end){}
      size_t size(){return end-start;}
      icLightNode* node;
      vplVPL* start, *end;
    };


    /*! Constructor. */
    icLighttree(std::vector<vplVPL>& vpls,float sceneradius):root(NULL),vpls(vpls)
    {
      this->sceneradiusSqr = sceneradius * sceneradius;
    }


    /*! Copy constructor. NOT implemented.*/
    icLighttree(icLighttree& cp);


    /*! Destructor. */
    ~icLighttree()
    { 
      if(root)delete root; 
    }


    /*! Build lightree. */
    void buildTree();


    /*! Build lightcuts metric based lighttree. Does not delete existing children.*/
    void buildTreeAgglomeraticeLocal();


    /*! Build compressed octree. Does not delete existing children.*/
    void createOctreeNonRecursiveAndCompress();


    /*! Set bb, replights, cones etc.*/
    void setProperties(icLightNode* lnode);


  public:
    icLightNode* root;              //! root of the tree
    MutexActive mutex;               //! mutex for threadsafety
    std::vector<vplVPL>& vpls;       //! vpls
    Random rnd;                      //! rnd number generator
    float sceneradiusSqr;            //! square of the scene radius
  };
}

#endif // icLighttree_h__
