// ======================================================================== //
// Copyright 2009-2013 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

/*
*  Release date: 20/05/2015
*  Authors: Norbert Bus & Nabil Mustafa*  
*  If you find this code useful, please cite it by referencing the following paper:
*  @article{BMB15, title={Illumination{Cut}}, author={Bus, Norbert and Mustafa, Nabil H. and Biri, Venceslas},
*           journal = {Computer Graphics Forum (Proceedings of Eurographics 2015)},year={2015}}
*/

#ifndef __EMBREE_CT_INTEGRATOR_H__
#define __EMBREE_CT_INTEGRATOR_H__

#include "devices/device_singleray/integrators/integrator.h"
#include "devices/device_singleray/renderers//renderer.h"
#include "common/image/image.h"

namespace embree
{

  // class for lightcuts to represent a node in a cut
  struct CutTriple{

    typedef enum {rightchild=0,leftchild=1} Child;

    CutTriple()
      :estimatedIntesity(zero),errorUb(zero),ubMaxComponent(0),node(NULL){}
    CutTriple(Color e,Color m,Node* n)
      :estimatedIntesity(e),errorUb(m),ubMaxComponent(max(m.r,max(m.g,m.b))),node(n) {}

    __forceinline bool isErrorSmallEnough(Color& intensityToCompareTo, float errorThreshold){
      //return average(this->errorUb) <= errorThreshold * average(intensityToCompareTo);
      return  (this->errorUb.r <= errorThreshold * intensityToCompareTo.r &&
        this->errorUb.g <= errorThreshold * intensityToCompareTo.g &&
        this->errorUb.b <= errorThreshold * intensityToCompareTo.b );
    }

    __forceinline Color adjustIntesityForSameRepLight(Child whichChild, size_t rid){
      // changing the intensity if the rep lights only differ in intensity
      Node* child = (whichChild==rightchild) ? this->node->rc : this->node->lc ;

#if REP_LIGHT_NUM > 1 
      Color ratio = child->repLights[rid].I;
      Color c = this->node->repLights[rid].I ;
#else
      Color ratio = child->repLights[0].I;
      Color c = this->node->repLights[0].I ;
#endif

      if(c.r != 0.0f) ratio.r *= rcp(c.r);
      if(c.g != 0.0f) ratio.g *= rcp(c.g);
      if(c.b != 0.0f) ratio.b *= rcp(c.b);

      return this->estimatedIntesity * ratio; 
    }

    __forceinline Color calcIntesityForNewRepLight(Child whichChild,DifferentialGeometry& dg,const Vector3f& wo,CompositedBRDF& brdfs, const Ref<BackendScene>& scene,float& epsilon,size_t rid){

      LightSample ls;
      Node* child = (whichChild==rightchild) ? this->node->rc : this->node->lc ;

#if REP_LIGHT_NUM > 1 
      Color I = child->repLights[rid].sample(dg, ls.wi, ls.tMax);
#else
      Color I = child->repLights[0].sample(dg, ls.wi, ls.tMax);
#endif
      if ( I == Color(zero) || ls.wi.pdf == 0.0f ||  dot(dg.Ns,Vector3f(ls.wi)) <= 0.0f) return zero;

      Color br = brdfs.eval(wo,dg,ls.wi.value,ALL) ;
      if (br==Color(zero))
        return zero;

      Ray ray(dg.P,ls.wi,dg.error*epsilon, ls.tMax-dg.error*epsilon);
      TIMESTART(stats.atomic.timeRT,rt0);
      rtcOccluded(scene->scene,(RTCRay&)ray);
      TIMESTOP(stats.atomic.timeRT,rt0);
      STATS(stats.atomic.rays++);
      if(!ray)
      {
        Color Lestim = I * br * rcp(ls.wi.pdf);
        return Lestim;
      }
      else
        return zero;
    }

    // operator to be used in std algorithms for making a heap of this object
    friend bool operator< ( const CutTriple& lhs, const CutTriple& rhs);

    Color estimatedIntesity, errorUb;
    float ubMaxComponent;
    Node* node;
  };

  class ctIntegrator : public Integrator
  {
    /*! Tracks the state of the path. */
    class __align(16) LightPath
    {
    public:

      /*! Constructs a path. */
      __forceinline LightPath (const Ray& ray, const Medium& medium = Medium::Vacuum(), const int depth = 0,
                               const Color& throughput = one, const bool ignoreVisibleLights = false, const bool unbend = true)
        : lastRay(ray), lastMedium(medium), depth(depth), throughput(throughput), ignoreVisibleLights(ignoreVisibleLights), unbend(unbend) {}

      /*! Extends a light path. */
      __forceinline LightPath extended(const Ray& nextRay, const Medium& nextMedium, const Color& weight, const bool ignoreVL) const {
        return LightPath(nextRay, nextMedium, depth+1, throughput*weight, ignoreVL, unbend && (nextRay.dir == lastRay.dir));
      }

    public:
      Ray lastRay;                 /*! Last ray in the path. */
      Medium lastMedium;           /*! Medium the last ray travels inside. */
      uint32 depth;                /*! Recursion depth of path. */
      Color throughput;            /*! Determines the fraction of radiance reaches the pixel along the path. */
      bool ignoreVisibleLights;    /*! If the previous shade point used shadow rays we have to ignore the emission
                                       of geometrical lights to not double count them. */
      bool unbend;                 /*! True of the ray path is a straight line. */
    };

  public:

    /*! Construction of integrator from parameters. */
    ctIntegrator(const Parms& parms);

    /*! Registers samples we need tom the sampler. */
    void requestSamples(Ref<SamplerFactory>& samplerFactory, const Ref<BackendScene>& scene);

    /*! Computes the radiance arriving at the origin of the ray from the ray direction. */
    Color Li(Ray& ray, const Ref<BackendScene>& scene, IntegratorState& state);
    Color LiCut(Ray& oray, const Ref<BackendScene>& scene, IntegratorState& state,std::vector<CutTriple>& cut);

    /* Configuration. */
  private:
    size_t maxDepth;               //!< Maximal recursion depth (1=primary ray only)
    float minContribution;         //!< Minimal contribution of a path to the pixel.
    float epsilon;                 //!< Epsilon to avoid self intersections.
    Ref<Image> backplate;          //!< High resolution background.

    /*! Random variables. */
  private:
    int lightSampleID;            //!< 2D random variable to sample the light source.
    int firstScatterSampleID;     //!< 2D random variable to sample the BRDF.
    int firstScatterTypeSampleID; //!< 1D random variable to sample the BRDF type to choose.
    std::vector<int> precomputedLightSampleID;  //!< ID of precomputed light samples for lights that need precomputations.

    float error;                  //!< error for lightcuts
  };


 

}

#endif
