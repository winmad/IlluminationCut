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

#ifndef clusterLightCuts_h__
#define clusterLightCuts_h__

#include "devices/device_singleray/default.h"
#include "devices/device_singleray/clustering/vpl/vplVPL.h"
#include "devices/device_singleray/brdfs/compositedbrdf.h"
#include "devices/device_singleray/clustering/lightcuts/ctKdTree.h"
#include "devices/device_singleray/brdfs/lambertian.h"
#include "devices/device_singleray/brdfs/specular.h"
#include <algorithm>

#pragma warning ( disable : 4068 )   //  disable 'unknown pragma warning'

namespace embree
{
  
  // Class representing a circular cone with axis(dir) and halfangle
  class Cone
  {
  public:

    /*! \brief Constructor for Cone.  */
    Cone(Vector3f dir,float halfAngle=0);


    /*! \brief Merge two cones.  */
    __forceinline Cone merge( Cone& other ) const
    {
      // angle between the two axis
      float angle = acos(dot(this->dir,other.dir));
      // check if one cone contains the other
      if (abs(this->halfAngle - other.halfAngle ) > angle )   // one contains the other
      {
        return this->halfAngle > other.halfAngle ? *this : other;
      }else      // they have separate parts
      {
        // define new cone containing both cones
        float newHalfAngle = 0.5f * (this->halfAngle + other.halfAngle + angle );
        newHalfAngle = clamp(newHalfAngle,0.0f,float(pi));      // may be bigger than the full angle
        Vector3f newdir;
        float len = length(this->dir - other.dir) ;    // handle two extreme cases when cross will return IND values
        if (len < 0.001  ){
          newdir = this->dir; 
        }else if (len > 1.999 ){
          newdir = this->dir;  
          newHalfAngle = float(pi);
        }else{
          // get the new axis by rotation
          LinearSpace3f rot = LinearSpace3f::rotate(cross(this->dir,other.dir),newHalfAngle - this->halfAngle); 
          newdir = rot * (this->dir);
        }
        // return cone
        return Cone(newdir,newHalfAngle);
      }
    }


  public:
    Vector3f dir;        // representing the axis of a cone
    float halfAngle;     // angle between an half line and the dir
    float cosHalfAngle;  // angle between an half line and the dir
    float rotationAngle; // angle that rotates the axis to 0,0,1
  };

  /*! \brief Node for lightcuts tree*/
  class Node 
  {
  public:
    typedef BBox3f BB;

    /*! \brief Constructor.  */
    Node(vplVPL light[REP_LIGHT_NUM],BBox3f bb,Cone cone,std::vector<vplVPL>& vpls, bool right[REP_LIGHT_NUM], Node* lc=NULL, Node* rc=NULL);

    /*! \brief Destructor. */
    ~Node();

    /*! \brief Distance between two nodes for agglomerative clustering. Custom metric from Lightcuts.  */
    __forceinline float dist( Node* node, float sceneradiusSqr)
    {
      // if the two nodes are the same for technical reasons we define them to be inf far (no merging can happen)
      if (this==node)
        return FLT_MAX;

      BBox3f bb=embree::merge(this->cell,node->cell);
      float diagSqr = squarelength(size(bb)); 
      float coeffSqr = sceneradiusSqr;                          
      Color I = this->repLights[0].I;
      I += node->repLights[0].I;
      Cone newCone =  this->cone.merge(node->cone);
      float dist = 1-cos(newCone.halfAngle);

      dist = reduce_add(I) * (diagSqr + coeffSqr * dist * dist); 
      return dist;
    }


     /*! \brief Merge two nodes  */
    Node* merge(Node* node, Random& rnd);
   

    /*! \brief  Calculate the ub of the luminosity received by a shaded point from the node.  */
    __forceinline Color upperBoundLuminance(const DifferentialGeometry& dg, CompositedBRDF* brdfs, const LinearSpace3f& rotation,const LinearSpace3f& rotationR){

      const Vector3f m = this->cell.lower;
      const Vector3f M = this->cell.upper;
      Vector3f bp = clamp(dg.P,m,M);
      float lenSqr = max( (float) (CLAMPING_DISTANCE*CLAMPING_DISTANCE), squarelength(dg.P - bp));

      float cosBound = 1.0f;
      float cosHalfAngle = this->cone.cosHalfAngle;
      if (cosHalfAngle <= 0.0f)
        cosBound = 1.0f;
      else
      {
        Vector3f vv = cross(Vector3f(0,0,1),this->cone.dir);
        if(vv == zero) vv = Vector3f(0,1,0); // proper rotation axis if crossprod was 0
        vv = normalize(vv);
        LinearSpace3f rotation = LinearSpace3f::rotate(vv, this->cone.rotationAngle);

        Vector3f corners[8];
        BBox3f xbox = BBox3f(EmptyTy());
        this->cell.GetCorners(corners);
        for (int i = 0; i < 8; i++)
          xbox.grow(xfmPoint(rotation,dg.P - corners[i]));

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
      return abs(repLights[0].I * upperBoundMaterialTerm(this->cell, rotation,rotationR, dg,brdfs) * (cosBound / lenSqr));
    }


     /*! \brief  Calculate upper bound of cos for given direction and point with the node bbox*/
    __forceinline float cosUpperBoundForBB(Vector3f& lower, Vector3f& upper){

      // calculating cos upper bound 
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


    /*! \brief  Calculate upper bound for lambertian and specular material*/
    __forceinline Color upperBoundMaterialTerm( const BBox3f &bbox, const LinearSpace3f& rotation,const LinearSpace3f& rotationR, const DifferentialGeometry & dp,  CompositedBRDF* brdfs)
    {
      // get material properties
      float n=0;
      Color kd(zero),ks(zero);
      brdfs->getPhong(kd,ks,n);

      Vector3f corners[8];
      BBox3f newbox = BBox3f(EmptyTy());
      this->cell.GetCorners(corners);
      // transforming the bbox
      for (int i = 0; i < 8; i++)
      {
        corners[i] -= dp.P;
        newbox.grow(xfmPoint(rotation,corners[i]));
      }

      // get cos upper bound
      float cosTheta = this->cosUpperBoundForBB(newbox.lower , newbox.upper);

      if (ks != zero && n != 0)
      {
        BBox3f newbox = BBox3f(EmptyTy());
        // transforming the bbox
        for (int i = 0; i < 8; i++)
          newbox.grow(xfmPoint(rotationR,corners[i]));

        float glossyCosBound = cosUpperBoundForBB(newbox.lower,newbox.upper);

        return (float(one_over_pi)*kd + (float(one_over_two_pi)* (n+2) * ks * pow(glossyCosBound, n))) * cosTheta;
    }

    return float(one_over_pi) * kd * cosTheta ;
  }


     /*! \brief Calculate luminance from replight. */
    __forceinline Color Luminance(const DifferentialGeometry& dg, CompositedBRDF* brdfs,const Vector3f& wo,size_t rid){

      Color Lestim;
      
      float dist; Sample3f wi;

#if REP_LIGHT_NUM > 1
      Color I = repLights[rid].sample(dg,wi,dist);
#else
      Color I = repLights[0].sample(dg,wi,dist);
#endif

      Lestim = I * brdfs->eval(wo,dg,wi.value,ALL) * rcp(wi.pdf); 
      return Lestim;
    }


    /*! \brief Print  */
    friend std::ostream & operator<<(std::ostream &os, const Node& n);


  public:   
    vplVPL repLights[REP_LIGHT_NUM];                   //!< representative lights
    BBox3f cell;                                       //!< bb
    Cone cone;                                         //!< cone ofdirections 
    Node* lc, *rc, *parent;                            //!< ...
    KdNode<Node>* kdnode;                              //!< kdnode for efficient clustering
    std::vector<vplVPL>& vpls;                         //!< vpls
    bool valid;                                        //!< used for efficient clustering
    bool repLightFromRightChild[REP_LIGHT_NUM];        //!< sharing replight with which child
  };


  class ClusterLightCuts 
  {
  public:


    /*! \brief Constructor.  */
    ClusterLightCuts(std::vector<vplVPL>& vpls):root(NULL),vpls(vpls){}


    /*! \brief Copy constructor. Not implemented.  */
    ClusterLightCuts(ClusterLightCuts& cp);


    /*! \brief Destructor.  */
    ~ClusterLightCuts(){ if(root)delete root; }


    /*! \brief Builds tree from leaves. Naive method  */
    Node* buildTreeFromLeavesNaive();


    /*! \brief Builds tree from leaves. Optimized method with heap  */
    Node* buildTreeFromLeavesHeap();


    /*! \brief Builds tree from leaves. Optimized method utilizing locality (eg nearby nodes are likely to be merged)  */
    Node* buildTreeFromLeavesLocal();


    /*! \brief Initialize the tree.  */
    void initClusterStructure(float sceneradius);

  public:
    Node* root;                          //!< root of the tree
    MutexSys mutex;                      //!< mutex
    std::vector<vplVPL>& vpls;           //!< vpls
    float sceneradiusSqr;                //!< square of the diagonal of the scene's bb 
  };

}

#pragma warning ( default : 4068 )   //  enable unknown pragma warnings

#endif // clusterLightCuts_h__
