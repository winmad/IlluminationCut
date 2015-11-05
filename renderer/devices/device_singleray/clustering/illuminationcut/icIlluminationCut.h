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

#ifndef icIlluminationCut_h__
#define icIlluminationCut_h__


#include "devices/device_singleray/default.h"
#include "devices/device_singleray/clustering/icStatistics.h"
#include "devices/device_singleray/clustering/illuminationcut/icLighttree.h"
#include "devices/device_singleray/clustering/illuminationcut/icPixeltree.h"
#include "devices/device_singleray/api/scene.h"
#include "common/sys/taskscheduler.h"
#include <vector>
#include <queue>


namespace embree{

  class BackendScene;

  class icIlluminationCut: public RefCount
  {
  public:
    typedef enum {icGood=0,refineLight=1,refineNode=2,icGoodNoCheck=3}  IC;


    /*! Constructor, a illuminationcut */
    icIlluminationCut(float error,bool sampling,size_t spp);


    /*! Initialize the structures and parameters */
    void init(Ref<BackendScene> scene,size_t siz,float eps);


    /*! Render the image using the illuminationcut clustering */
    void illuminationcut();


    /* Render the approximate image using MT geometric cut */
    void geomcutRenderMT();


    /* Render the image using MT illuminationcut*/
    void illumcutRenderMT();


    // MT helper functions
    TASK_RUN_FUNCTION(icIlluminationCut,illumcut);
    TASK_RUN_FUNCTION(icIlluminationCut,geomcut);
    TASK_COMPLETE_FUNCTION(icIlluminationCut,donothing);


    /*! Creates a deque of pairs for MT algorithms */
    void initPairsMT(std::deque<std::pair<icLightNode*, icPixelNode*> > & pairs, size_t siz);

    /*! Render the pair using geometric cut for the approx image*/
    void geomcutRender(icLightNode* lnode,icPixelNode* pnode,Random& rnd);


    /*! Distribute color towards the leafs and set the maximum error threshold */
    Color distributeMinColor(icPixelNode* pnode,Color parentColor);


    /*! Renders the pairs according to illuminationcut */
    void illumcutRender(std::vector<std::pair<icLightNode*, icPixelNode* > >& localpairs,Random& rnd );

    /*! Checks of two nodes are in a geometric cut */
    __forceinline IC checkGeomCut(icLightNode* lnode,icPixelNode* pnode)
    {
      Vector3f lc = center(lnode->cell).pos();
      float    lr = lnode->posRadius;
      Vector3f pc = center(pnode->cell).pos();
      float    pr = pnode->posRadius;
      Vector3f wi = lc-pc;
      float len = length(wi);
      float dist = len - lr - pr;

#ifdef GLOSSY_DIR_CUT
      float angle =  rcp(len) * dot(wi,pnode->pixel->reflected_wo);
      bool incut;
      //in the specular direction we are more detailed;
      float epsilonMultiplier = 0.3;
      float minAngle = 0.95;
      if ( angle > minAngle && (pnode->maxKs!=zero) )  // only if there is Ks component
        incut = ( (epsilonMultiplier * geomEpsilon * dist) > max(lr,pr) );
      else
        incut = ( (geomEpsilon * dist) > max(lr,pr) );
#else
      bool  incut = ( (geomEpsilon * dist) > max(lr,pr) );
#endif

      if (!incut)
      {
        if (pnode->children[0] && lnode->children[0])
          return  ( lr > pr) ? refineLight : refineNode;
        if (lnode->children[0])
          return refineLight;
        if (pnode->children[0])
          return refineNode;
        return icGood;
      }
      else
#ifdef AGGLOMERATIVE_LIGHTTREE
        if (lnode->cone.cosHalfAngle > 0.99)
          return icGood;
        else
          if (lnode->children[0])
            return refineLight;
          else
            return icGood;
#else
        return icGood;
#endif
    }

    /*! Check illumination based well-separatedness */
    __forceinline  IC checkIlluminationCut(icLightNode* lnode,icPixelNode* pnode)
    {
      Vector3f lc = center(lnode->cell).pos();
      float    lr = lnode->posRadius;
      Vector3f pc = center(pnode->cell).pos();
      float    pr = pnode->posRadius;
      float realDist = length(lc-pc) ;
      float dist = realDist - lr - pr;

      float rcpDist = rcp( max(dist,(float) CLAMPING_DISTANCE) );

      TIMESTART(stats.atomic.timeUB,ub0);
      Color ubB = pnode->upperBoundMaterial(lnode->cell);
      Color ubI = lnode->upperBoundIntensity(pnode->cell);
      TIMESTOP(stats.atomic.timeUB,ub0);
      Color ub = ubB*ubI*rcpDist*rcpDist;

      bool incut = (pnode->color.r  >= ub.r) && (pnode->color.g  >= ub.g) && (pnode->color.b  >= ub.b);

      if (ub==Color(zero))
        return icGoodNoCheck;

      if (!incut)
      {
        bool pinner = pnode->children[0] ;
        bool linner = lnode->children[0] ;
        if (pinner && linner)
        {
          return  ( lr > pr) ? refineLight : refineNode;
        }
        if (linner)
          return refineLight;
        return icGood;
      }
      else
        return icGood;
    }


    /*! Shade only the random pixel from rep light */
    void shadeRepPixel(icLightNode* lnode,icPixelNode* pnode, Random& rnd);


    /*! Calculate for the node with upper bound visualization. Amortized cost of clustering. */
    void shadeSubtreeUB(icLightNode* lnode,icPixelNode* pnode);


    /*! Shade all pixels in the node from the replight*/
    void shadeSubtreeFull(icLightNode* lnode,icPixelNode* pnode);

    
    /*! Shade all pixels in the node from the replight but sample only a few visibilities*/
    void shadeSubtreeWithSamplingDescend(icLightNode* lnode,icPixelNode* pnode,Random& rnd,size_t depth);


    /*! Shade all pixels in the node from the replight but sample only a few visibilities*/
    void shadeSubArrayWithSamplingDescend(icLightNode* lnode,icPixel* start,icPixel* end,Random& rnd,size_t depth);


    /*! Checks the visibility to a few pixels and also shades these at the same time. */
    void checkVisibAndRender(icLightNode* lnode, icPixel** px, size_t siz,bool& allvisib, bool& alloccl);


    Ref<icLighttree> lighttree;                                       //< lightree
    icPixeltree pixeltree;                                            //< pixeltree
    std::deque<std::pair<icLightNode*, icPixelNode*> > pairs;         //< pairs for multithreaded algo
    float geomEpsilon;                                                //< controls the approximation image
    float error;                                                      //< controls the error for the final image
    const size_t exhaustiveRenderSize;                                //< groups of pixels less than this are not sampled for visibility
    Ref<BackendScene> scene;                                          //< scene
    MutexActive mutex;                                                //< multi threading
    bool sampling;                                                    //< control whether the visibility is sampled
  };
}
#endif 
