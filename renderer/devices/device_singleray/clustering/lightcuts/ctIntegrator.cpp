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

#include "devices/device_singleray/clustering/icStatistics.h"
#include "devices/device_singleray/clustering/lightcuts/ctIntegrator.h"

namespace embree
{
  ctIntegrator::ctIntegrator(const Parms& parms)
    : lightSampleID(-1), firstScatterSampleID(-1), firstScatterTypeSampleID(-1)
  {
    maxDepth        = parms.getInt  ("maxDepth"       ,10    );
    minContribution = parms.getFloat("minContribution",0.01f );
    epsilon         = parms.getFloat("epsilon"        ,32.0f)*float(ulp);
    backplate       = parms.getImage("backplate");
    error           = parms.getFloat("error",0.02f);
  }

  void ctIntegrator::requestSamples(Ref<SamplerFactory>& samplerFactory, const Ref<BackendScene>& scene)
  {
    precomputedLightSampleID.resize(scene->allLights.size());

    lightSampleID = samplerFactory->request2D();
    for (size_t i=0; i<scene->allLights.size(); i++) {
      precomputedLightSampleID[i] = -1;
      if (scene->allLights[i]->precompute())
        precomputedLightSampleID[i] = samplerFactory->requestLightSample(lightSampleID, scene->allLights[i]);
    }
    firstScatterSampleID = samplerFactory->request2D((int)maxDepth);
    firstScatterTypeSampleID = samplerFactory->request1D((int)maxDepth);
  }


  __forceinline bool operator< (const CutTriple& lhs, const CutTriple& rhs){
    return lhs.ubMaxComponent < rhs.ubMaxComponent;
  } 


  Color ctIntegrator::LiCut(Ray& oray, const Ref<BackendScene>& scene, IntegratorState& state,std::vector<CutTriple>& cut)
  {
    LightPath lightPath(oray); 
    size_t rid = state.spp_id;
    /*! Traverse ray. */
    DifferentialGeometry dg;
    rtcIntersect(scene->scene,(RTCRay&)lightPath.lastRay);
    scene->postIntersect(lightPath.lastRay,dg);
    state.numRays++;

    Color L = zero;
    const Vector3f wo = -lightPath.lastRay.dir;
    BRDFType directLightingBRDFTypes = (BRDFType)(ALL); 

    /*! Environment shading when nothing hit. */
    if (!lightPath.lastRay)
    {
      if (backplate && lightPath.unbend) {
        const int x = clamp(int(state.pixel.x * backplate->width ), 0, int(backplate->width )-1);
        const int y = clamp(int(state.pixel.y * backplate->height), 0, int(backplate->height)-1);
        L = backplate->get(x, y);
      }
      else {
        if (!lightPath.ignoreVisibleLights)
          for (size_t i=0; i<scene->envLights.size(); i++)
            L += scene->envLights[i]->Le(wo);
      }
#if defined(UB_IMAGE)
      float ubNumber = 0;
      float r,g,b;
      ColormapJet::mapToJet(ubNumber,r,g,b);
      return Color(r,g,b);
#else
      return L;
#endif
    }

    /*! face forward normals */
    bool backfacing = false;
    if (dot(dg.Ng, lightPath.lastRay.dir) > 0) {
      backfacing = true; dg.Ng = -dg.Ng; dg.Ns = -dg.Ns;
    }

    /*! Shade surface. */
    CompositedBRDF brdfs;
    if (dg.material) dg.material->shade(lightPath.lastRay, lightPath.lastMedium, dg, brdfs);

#ifndef UB_IMAGE
    /*! Add light emitted by hit area light source. */
    if (!lightPath.ignoreVisibleLights && dg.light && !backfacing)
      L += dg.light->Le(dg,wo);
#endif

    /*! Check if any BRDF component uses direct lighting. */
    bool useDirectLighting = false;
    for (size_t i=0; i<brdfs.size(); i++)
      useDirectLighting |= (brdfs[i]->type & directLightingBRDFTypes) != NONE;

    /*! Direct lighting. */
    if (useDirectLighting)
    {
      TIMESTART(stats.atomic.timeHeap,heap0);
      cut.clear();
      TIMESTOP(stats.atomic.timeHeap,heap0);

      // the total value for the current cut
      Color totalIntensity(zero);
      Color estimatedIntensity(zero);

      //calculation for root
      Node* root = scene->clusterlightcuts->root;
      LightSample ls;
      Color errorUbForRoot(zero);

      root->repLights[rid].sample(dg, ls.wi, ls.tMax);

      if (dot(dg.Ns,Vector3f(ls.wi)) > 0.0f )
      {
        estimatedIntensity = root->Luminance(dg,&brdfs,wo,rid);
        if (estimatedIntensity != Color(zero))
        {
          STATS(stats.atomic.rays++);
          Ray ray(dg.P, ls.wi, dg.error*epsilon, ls.tMax-dg.error*epsilon, lightPath.lastRay.time);
          TIMESTART(stats.atomic.timeRT,rt1);
          rtcOccluded(scene->scene,(RTCRay&)ray);
          TIMESTOP(stats.atomic.timeRT,rt1);
          if( ray )
            estimatedIntensity = Color(zero);
        }
      }
      else
        estimatedIntensity = Color(zero);

      // calculate rotations for the pixel
      STATS(stats.atomic.ub++);
      STATS(stats.atomic.clusters++);
      float angleN = -acosf(clamp(dot(Vector3f(0,0,1) , dg.Ns),-1.0f,1.0f));
      Vector3f rotationAxis = cross(Vector3f(0,0,1) , dg.Ns);
      if (rotationAxis == zero) rotationAxis = Vector3f(0,1,0);
      rotationAxis=normalize(rotationAxis);
      LinearSpace3f rotation = LinearSpace3f::rotate(rotationAxis, angleN);

      Vector3f R = - wo + dg.Ns*(2*(dot(wo,dg.Ns))); 
      if (dg.Ns == wo)
        R = wo;
      float angleR = -acosf(clamp(dot(Vector3f(0,0,1) , R),-1.0f,1.0f));
      rotationAxis = cross(Vector3f(0,0,1) , R);
      if (rotationAxis == zero) rotationAxis = Vector3f(0,1,0);
      rotationAxis=normalize(rotationAxis);
      LinearSpace3f rotationR = LinearSpace3f::rotate(rotationAxis, angleR);
      TIMESTART(stats.atomic.timeUB,ub1);
      errorUbForRoot = root->upperBoundLuminance(dg,&brdfs,rotation,rotationR);
      TIMESTOP(stats.atomic.timeUB,ub1);
      totalIntensity = estimatedIntensity;
      cut.push_back(CutTriple(estimatedIntensity,errorUbForRoot,root));

      size_t counter = 0;
      while (cut.size()) // while the first element of the heap has too big error bound or everything is popped
      {
        CutTriple cutToRefine = cut.front();

        if ( cutToRefine.isErrorSmallEnough(totalIntensity,this->error) || ( (cut.size() + counter) > CT_MAXCUT ) )
          break;        // the top element is ok hence we can break

        if( !cutToRefine.node->rc && !cutToRefine.node->lc )
        {   // this leaf node contributes a lot but we cant split it so we just remove it from the heap
          // its contribution will remain in totalIntensity
          TIMESTART(stats.atomic.timeHeap,heap1);
          std::pop_heap (cut.begin(),cut.end()); cut.pop_back();
          TIMESTOP(stats.atomic.timeHeap,heap1);
          counter++;
          continue;
        }
        else
        {
          // remove the current nodes contribution and delete it from the heap
          totalIntensity -= cutToRefine.estimatedIntesity;
          STATS(stats.atomic.clusters--);
          TIMESTART(stats.atomic.timeHeap,heap2);
          std::pop_heap (cut.begin(),cut.end()); 
          cut.pop_back();
          TIMESTOP(stats.atomic.timeHeap,heap2);

          //add the children and their contribution
          if(cutToRefine.node->rc)
          {
            Color estimatedIntensity = zero;

            if ( cutToRefine.node->repLightFromRightChild[rid])
              estimatedIntensity = cutToRefine.adjustIntesityForSameRepLight(CutTriple::rightchild,rid);
            else
              estimatedIntensity = cutToRefine.calcIntesityForNewRepLight(CutTriple::rightchild,dg,wo,brdfs,scene,epsilon,rid);

            // add it to the cut
            STATS(stats.atomic.ub++);
            STATS(stats.atomic.clusters++);
            totalIntensity += estimatedIntensity; 
            TIMESTART(stats.atomic.timeUB,ub2);
            Color ub = cutToRefine.node->rc->upperBoundLuminance(dg,&brdfs,rotation,rotationR);
            TIMESTOP(stats.atomic.timeUB,ub2);
            TIMESTART(stats.atomic.timeHeap,heap3);
            cut.push_back( CutTriple(estimatedIntensity,ub,cutToRefine.node->rc) ); 
            std::push_heap (cut.begin(),cut.end());
            TIMESTOP(stats.atomic.timeHeap,heap3);
            STATS(if (!cutToRefine.node->rc->rc && !cutToRefine.node->rc->lc)
              stats.atomic.singletonpair++);
          }

          if(cutToRefine.node->lc)
          {
            Color estimatedIntensity = zero;

            if ( !cutToRefine.node->repLightFromRightChild[rid] )
              estimatedIntensity = cutToRefine.adjustIntesityForSameRepLight(CutTriple::leftchild,rid);
            else
              estimatedIntensity = cutToRefine.calcIntesityForNewRepLight(CutTriple::leftchild,dg,wo,brdfs,scene,epsilon,rid);

            // add it to the cut
            STATS(stats.atomic.ub++);
            STATS(stats.atomic.clusters++);
            totalIntensity += estimatedIntensity;
            TIMESTART(stats.atomic.timeUB,ub3);
            Color ub = cutToRefine.node->lc->upperBoundLuminance(dg,&brdfs,rotation,rotationR);
            TIMESTOP(stats.atomic.timeUB,ub3);
            TIMESTART(stats.atomic.timeHeap,heap4);
            cut.push_back( CutTriple(estimatedIntensity,ub,cutToRefine.node->lc)  ); 
            std::push_heap (cut.begin(),cut.end());
            TIMESTOP(stats.atomic.timeHeap,heap4);
            STATS(if (!cutToRefine.node->lc->rc && !cutToRefine.node->lc->lc)
              stats.atomic.singletonpair++);
          }
        }
      }

#if defined(UB_IMAGE)
      float ubNumber = 2 * cut.size() + counter;
      float r,g,b;
      ColormapJet::mapToJet(ubNumber,r,g,b);
      L = Color(r,g,b);
#else
      // add the indirect illumination to the color of the pixel
      L += totalIntensity;
#endif
    }

    return L;
  }

  Color ctIntegrator::Li(Ray& ray, const Ref<BackendScene>& scene, IntegratorState& state) {
    throw std::runtime_error("Not implemented");
  }
}

