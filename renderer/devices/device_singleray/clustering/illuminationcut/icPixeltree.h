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
#ifndef icPixeltree_h__
#define icPixeltree_h__


#include "devices/device_singleray/default.h"
#include "devices/device_singleray/clustering/icParameters.h"
#include "devices/device_singleray/shapes/differentialgeometry.h"
#include "devices/device_singleray/brdfs/compositedbrdf.h"
#include "common/math/vec6.h"
#include <vector>
#include <stack>
#include <algorithm>


namespace embree {

  typedef BBox<Vec6f> BBox6f;

  class icPixel
  {
  public:

    /*! \brief Default constructor.  */
    icPixel(){}


    /*! \brief Constructor.  */
    icPixel(DifferentialGeometry& dg,Vector3f wo,Color kd, Color ks, float exp, size_t x,size_t y,size_t s);


    /*! \brief Destructor.  */
    ~icPixel(){}

    
    /*! \brief Evaluate the brdfs for a certain light direction.  */
    Color evalbrdfs(const Vector3f& wi) const {
      Color spec(zero);
      if (ks!=zero)
      {
        float beta = dot(reflected_wo,wi);
        if (beta > 0) 
          spec = ks * (exp+2) * float(one_over_two_pi) * pow(beta,exp) ;
      }
      return (kd * float(one_over_pi) + spec ) * clamp(dot(wi,Ns));
    }


    /*! \brief Evaluate the lambertian brdf but without kd (color).  */
    Color evalLambertianAngle(const Vector3f& wi) const {
      return Color(1.0f) * float(one_over_pi)  * clamp(dot(wi,Ns));
    }


    /*! \brief Evaluate the specular brdf but without ks (color).  */
    Color evalSpecularAngle(const Vector3f& wi) const {
      Color spec(zero);
      float beta = dot(reflected_wo,wi);
      if (beta > 0) 
        return Color(1.0f) * (exp+2) * float(one_over_two_pi) * pow(beta,exp) * clamp(dot(wi,Ns));
      else
        return  Color(zero);
    }


    Vector3f P, Ns,reflected_wo;  //! position, shading normal, reflected view ray
    float error;                  //! error for raytracing
    Color kd,ks,L;                //! brdf coeffs, luminance
    float exp;                    //! spec exponenet
    size_t x,y,s;                 //! pixel coords and sample id
  };


  class icPixelNode
  {
  public:

    friend class icPixeltree;
    enum { nChildren = 2 };


    /*! \brief Constructor.  */
    icPixelNode();


    /*! \brief Destructor.  */
    ~icPixelNode();

    /*! Create bbox for the i-th child in the octree. */
    BBox6f createBboxForChild( size_t i);


    /*! Calculate the dimension and the position of the next subdividing plane */
    __forceinline void getNextDividingPlanePosition(float& pos, size_t& dim)
    {
      BBox6f& b = this->cell;
      Vec6f diag = b.size();
      float dirDiam = length(diag.dir());
      bool dirDivision = dirDiam > MAXSIZE_DIRECTIONAL_BBOX;

      if (!dirDivision) // pos subdivison
      {
        dim = this->level % 3 ;
        pos = b.lower[dim]  + diag[dim] * 0.5f;
      }else // dir subdivison
      {
        dim = 3 + this->level % 3 ;
        pos = b.lower[dim]  + diag[dim] * 0.5f;
      }
      return;
    }


    /*! Calculate the number of nodes in the subtree. */
    size_t subTreeSize();


    /*! Calculate the depth of the subtree. */
    size_t subTreeDepth();


    /*! Upper bound the cos for a bb from (0,0,1). */
    __forceinline float CosUpperBoundForBB(Vector3f& lower, Vector3f& upper) const {
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


    /*! Upper bound material term. */
    __forceinline Color upperBoundMaterial(const BBox6f& bb) const {
      Vector3f corners[8];
      BBox3f newbox = BBox3f(EmptyTy());
      BBox3f b = BBox3f(bb.lower.pos() - this->cell.upper.pos(),bb.upper.pos() - this->cell.lower.pos());
      b.GetCorners(corners);
      for (int i = 0; i < 8; i++)
      {
        newbox.grow(xfmPoint(this->rotationNs,corners[i]));
      }

      // get cos upper bound
      float cosTheta = this->CosUpperBoundForBB(newbox.lower , newbox.upper);
      if (conjoint(this->pixel->P,BBox3f(bb.lower.pos(),bb.upper.pos())))
        cosTheta = 1.0f;

      if (maxKs != Color(zero) && maxExp != 0)
      {
        BBox3f newbox = BBox3f(EmptyTy());
        for (int i = 0; i < 8; i++)
          newbox.grow(xfmPoint(this->rotationR,corners[i]));

        float glossyCosBound = CosUpperBoundForBB(newbox.lower,newbox.upper);
        if (conjoint(this->pixel->P,BBox3f(bb.lower.pos(),bb.upper.pos())))
          glossyCosBound = 1.0f;

        return (float(one_over_pi)*maxKd + (float(one_over_two_pi)* (maxExp+2) * maxKs * pow(glossyCosBound, minExp))) * cosTheta;
      }
      return float(one_over_pi) * maxKd * cosTheta ;
    }
  


    /*! Number of pixels the node represents. */
    __forceinline size_t size()
    {
      size_t siz = end - start;
      return siz;
    }

    icPixelNode* parent;                        //< pointer to parent node
    icPixelNode* children[nChildren];           //< pointer to children, null if there is no child
    BBox6f cell;                                 //< bbox of the pixels in this node
    float posRadius;                             //< radius of spatial bbox
    icPixel* pixel;                             //< for leaf node the pixel, otherwise random rep
    icPixel* start;                             //< pixel range start
    icPixel* end;                               //< pixel range end
    Color maxKd, maxKs;                          //< max material coeffs
    float maxExp, minExp;                        //< spec exp for ub (min is for power of a value < 1)
    Color color;                                 //< luminance used for storing approx image and ub
    int level;                                   //< node level
    LinearSpace3f rotationNs, rotationR;         //< rotations for the ub calculations
  };


  class icPixeltree
  {
  public:
    enum{nChildren=icPixelNode::nChildren};


    /*! \brief Constructor.  */
    icPixeltree();;


    /*! \brief Destructor.  */
    ~icPixeltree();


    /*! \brief Initialize pixel storage.  */
    void initPixelStorage(size_t pixels);;


    /*! \brief Resize pixel storage to actual value.  */
    void shrink();


    /*! \brief Add a pixel to the storage. Threadsafe.  */
    void addPixel(DifferentialGeometry& dg,Vector3f& wo,Color& kd, Color& ks, float exp, size_t x,size_t y,size_t s);


    /*! \brief Build pixeltree.  */
    void buildTree();


    /*! \brief Set properties: bb, posradius, brdf bounds.  */
    void setProperties(icPixelNode* pnode);


    /*! Build compressed octree. Does not delete existing children.*/
    void createOctreeNonRecursiveAndCompress();


    icPixelNode* root;                       //< root of the tree
    std::vector<icPixel> pixels;             //< pixel storage
    Atomic pixelnum;                          //< atomic var for indexing the pixel storage
  };
}

#endif // icPixeltree_h__
