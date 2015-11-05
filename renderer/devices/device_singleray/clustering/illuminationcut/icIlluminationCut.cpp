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

#include "devices/device_singleray/clustering/illuminationcut/icIlluminationCut.h"
#include "devices/device_singleray/api/scene.h"


namespace embree{

  void icIlluminationCut::init(Ref<BackendScene> scene,size_t siz,float eps)
  {
    this->scene = scene;
    this->lighttree = scene->lighttree;
    this->pixeltree.initPixelStorage(siz);
    this->geomEpsilon = eps;
  }


  void icIlluminationCut::shadeSubtreeUB( icLightNode* lnode,icPixelNode* pnode )
  {
    STATS(stats.atomic.clusters += pnode->end - pnode->start);
    STATS(stats.atomic.gpsize += pnode->end - pnode->start);
    icPixel* pix = pnode->start;
    for ( ; pix != pnode->end ; ++pix)
    {
      pix->L += Color(rcp( float(pnode->end - pnode->start)));
    }
    return;
  }

  void icIlluminationCut::shadeSubtreeFull( icLightNode* lnode,icPixelNode* pnode )
  {
    STATS(stats.atomic.clusters += pnode->end - pnode->start);
    STATS(stats.atomic.gpsize += pnode->end - pnode->start);
    icPixel* pix = pnode->start;
    for ( ; pix != pnode->end ; ++pix)
    {
      Sample3f wi; float tMax;
#if REP_LIGHT_NUM > 1 
      Color I = lnode->repLights[pix->s].sample(pix->P,wi,tMax);
#else
      Color I = lnode->repLights[0].sample(pix->P,wi,tMax);
#endif
      if ( I == Color(zero) || wi.pdf == 0.0f) continue;

      /*! Evaluate BRDF */
      Color brdf = pix->evalbrdfs(wi);
      if (brdf == Color(zero)) continue;

      float err = pix->error * 32.0f * float(ulp);
      Ray shadowRay(pix->P, wi, err, tMax - err);
      TIMESTART(stats.atomic.timeRT,rt0);
      rtcOccluded(scene->scene,(RTCRay&)shadowRay);
      TIMESTOP(stats.atomic.timeRT,rt0);
      STATS(stats.atomic.rays++);
      if (shadowRay) continue;

      /*! Evaluate BRDF. */
      pix->L += I * brdf * rcp(wi.pdf);
    }
    return;
  }

  void icIlluminationCut::shadeSubArrayWithSamplingDescend( icLightNode* lnode,icPixel* start, icPixel* end,Random& rnd,size_t depth)
  {
    const size_t range = end - start;
    // if light is singleton we have no error bound therefore render it in full detail
    if (range > exhaustiveRenderSize &&  lnode->children[0] && depth--)
    {
      // get a few sampled pixels
      icPixel* px[16];
      size_t siz;
      if (bool deterministicSample = true)
      {
        if(range > 32)
          siz = 16;
        else if (range > 16) siz = 8;
        else siz = 4;
        size_t step = range/siz;
        for (size_t i = 1 ; i < siz-2 ; ++i )
          px[i] = start + i*step;
        px[0] = start;
        px[siz-2] = end-1;
        px[siz-1] = end;
      }
      else 
      {
        if(range > 32) siz = 16;
        else if (range > 16) siz = 8;
        else siz = 4;
        px[0] = start;
        px[siz-1] = end;
        size_t step = range/siz;
        for (size_t i = 1 ; i < siz-2 ; ++i )
          px[i] = start + rnd.getInt(range);
        //px[i] = start + i*step;
        px[siz-2] = end-1;
        std::sort(px+1, px+siz-2);
      }
      //check the visibility conditions
      bool allVisib = true;
      bool allOccl = true;
      this->checkVisibAndRender(lnode,px,siz,allVisib,allOccl);
      if (allVisib)
      {
        //render everything with visibility true
        for(size_t i = 0; i < (siz-1) ; ++i)
        {
          //if(px[i]==px[i+1])continue;
          icPixel* pix = px[i]+1;
          for ( ; pix < px[i+1] ; ++pix)
          {
            Sample3f wi; float tMax;
#if REP_LIGHT_NUM > 1 
            Color I = lnode->repLights[pix->s].sample(pix->P,wi,tMax);
#else
            Color I = lnode->repLights[0].sample(pix->P,wi,tMax);
#endif
            if ( I == Color(zero) || wi.pdf == 0.0f) continue;

            /*! Evaluate BRDF */
            Color brdf = pix->evalbrdfs(wi);
            if (brdf == Color(zero)) continue;

            /*! Evaluate BRDF. */
            pix->L += I * brdf * rcp(wi.pdf);
          }
        }
      } 
      else if (! allOccl)
      {
        for(size_t i = 0; i < (siz-1) ; ++i)
        {
          if ( (px[i+1]-px[i]) > 1)
            shadeSubArrayWithSamplingDescend(lnode,px[i]+1,px[i+1],rnd,depth);
        }
      }
    }
    else
    {
      // shade every pixel
      icPixel* pix = start;
      for ( ; pix < end ; ++pix)
      {
        Sample3f wi; float tMax;
#if REP_LIGHT_NUM > 1 
        Color I = lnode->repLights[pix->s].sample(pix->P,wi,tMax);
#else
        Color I = lnode->repLights[0].sample(pix->P,wi,tMax);
#endif
        if ( I == Color(zero) || wi.pdf == 0.0f) continue;

        /*! Evaluate BRDF */
        Color brdf = pix->evalbrdfs(wi);
        if (brdf == Color(zero)) continue;

        float err = pix->error * 32.0f * float(ulp);
        Ray shadowRay(pix->P, wi, err, tMax - err);
        TIMESTART(stats.atomic.timeRT,rt0);
        rtcOccluded(scene->scene,(RTCRay&)shadowRay);
        TIMESTOP(stats.atomic.timeRT,rt0);
        STATS(stats.atomic.rays++);
        if (shadowRay) continue;

        /*! Evaluate BRDF. */
        pix->L += I * brdf * rcp(wi.pdf);
      }
    }
  }
  void icIlluminationCut::shadeSubtreeWithSamplingDescend( icLightNode* lnode,icPixelNode* pnode,Random& rnd,size_t depth)
  {
    STATS(stats.atomic.clusters += pnode->end - pnode->start);
    STATS(stats.atomic.gpsize += pnode->end - pnode->start);
    this->shadeSubArrayWithSamplingDescend(lnode,pnode->start,pnode->end,rnd,depth);

  }

  void icIlluminationCut::illumcutRenderMT()
  {
    pairs.clear();
    pairs.push_back(std::pair<icLightNode*, icPixelNode*>(lighttree->root,pixeltree.root) );
    // add 2000 pairs to the global work queue
    this->initPairsMT(pairs,2000);
    auto comp = [](const std::pair<icLightNode*,icPixelNode*> p1,const std::pair<icLightNode*,icPixelNode*> p2)
    {
      return reduce_add(p1.second->color) < reduce_add(p2.second->color);
    };
    std::sort(pairs.begin(),pairs.end(),comp);

    // start multithreaded illuminationcut
    TaskScheduler::EventSync event;
    TaskScheduler::Task task(&event,_illumcut,this,TaskScheduler::getNumThreads(),_donothing,this,"render::getPixels");
    TaskScheduler::addTask(-1,TaskScheduler::GLOBAL_BACK,&task);
    event.sync();
  }

  void icIlluminationCut::illumcut(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event)
  {
    Random rnd(RND_SEED);
    std::vector<std::pair<icLightNode*, icPixelNode* > > localpairs;
    localpairs.reserve(2000);
    /*! tile pick loop */
    while (true)
    {
      /*! pick a new tile */
      std::pair<icLightNode*,icPixelNode*> pair;
      {
        Lock<MutexActive> lock(this->mutex);
        if (this->pairs.empty())
          break;
        else
        {
          pair = pairs.front();
          pairs.pop_front();
        }
      }
      localpairs.clear();
      localpairs.push_back(pair);
      this->illumcutRender(localpairs,rnd);
    }
  }
  void icIlluminationCut::geomcut(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event)
  {
    Random rnd(RND_SEED);
    /*! tile pick loop */
    while (true)
    {
      /*! pick a new tile */
      std::pair<icLightNode*,icPixelNode*> pair;
      {
        Lock<MutexActive> lock(this->mutex);
        if (this->pairs.empty())
          break;
        else
        {
          pair = pairs.front();
          pairs.pop_front();
        }
      }
      this->geomcutRender(pair.first,pair.second,rnd);
    }
  }

  void icIlluminationCut::donothing(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event)
  {
  }

  void icIlluminationCut::geomcutRenderMT()
  {
    pairs.clear();
    pairs.push_back(std::pair<icLightNode*, icPixelNode*>(lighttree->root,pixeltree.root) );
    this->initPairsMT(pairs,500);

    // start multithreaded wsi
    TaskScheduler::EventSync event;
    TaskScheduler::Task task(&event,_geomcut,this,TaskScheduler::getNumThreads(),_donothing,this,"render::getPixels");
    TaskScheduler::addTask(-1,TaskScheduler::GLOBAL_BACK,&task);
    event.sync();
  }






  void icIlluminationCut::checkVisibAndRender( icLightNode* lnode, icPixel** px, size_t siz,bool& allvisib, bool& alloccl )
  {
    allvisib = alloccl = true;
    for( size_t i = 0; i < (siz-1) ; ++i )
    {
      if (px[i]==px[i+1]) continue;
      Sample3f wi; float tMax;
#if REP_LIGHT_NUM > 1 
      Color I = lnode->repLights[px[i]->s].sample(px[i]->P,wi,tMax);
#else
      Color I = lnode->repLights[0].sample(px[i]->P,wi,tMax);
#endif
      if ( I == Color(zero) || wi.pdf == 0.0f)
      {
        alloccl &= true;
        allvisib &= false;
        continue;
      }

      /*! Evaluate BRDF */
      Color brdf = px[i]->evalbrdfs(wi);
      if (brdf == Color(zero)){
        alloccl &= true;
        allvisib &= false;
        continue;
      }

      float err = px[i]->error * 32.0f * float(ulp);
      Ray shadowRay(px[i]->P, wi, err, tMax - err);
      TIMESTART(stats.atomic.timeRT,rt0);
      rtcOccluded(scene->scene,(RTCRay&)shadowRay);
      TIMESTOP(stats.atomic.timeRT,rt0);
      STATS(stats.atomic.rays++);
      if (shadowRay) 
      {
        alloccl &= true;
        allvisib &= false;
        continue;
      }
      else
      {
        alloccl &= false;
        allvisib &= true;
      }

      /*! Evaluate BRDF. */
      px[i]->L += I * brdf * rcp(wi.pdf);
    }
  }

  void icIlluminationCut::illuminationcut()
  {
    // build pixeltree
    double t0 = getSeconds();
    pixeltree.buildTree();
    double dt = getSeconds() - t0;
    std::cout<<"Clustering preprocess for Pixeltree took: "<< dt * 1000.0f <<" ms"<<std::endl;

    // build and render the approx image
    t0 = getSeconds();
    this->geomcutRenderMT();
    dt = getSeconds() - t0;
    std::cout<<"Geometric cut: "<< dt * 1000.0f <<" ms"<<std::endl;
    t0 = getSeconds();
    this->distributeMinColor(this->pixeltree.root,Color(zero));
    dt = getSeconds() - t0;
    std::cout<<"DistributeMinColor: "<< dt * 1000.0f <<" ms"<<std::endl;
#ifdef APPROX_IMAGE
    return;
#endif

    // build and render illumation based wspd
    t0 = getSeconds();
    this->illumcutRenderMT();
    dt = getSeconds() - t0;
    std::cout<<"IlluminationCut: "<< dt * 1000.0f <<" ms"<<std::endl;
  }


  void icIlluminationCut::illumcutRender( std::vector<std::pair<icLightNode*, icPixelNode* > >& localpairs,Random& rnd )
  {
    while( !localpairs.empty() )
    {
      // pop a pair
      std::pair<icLightNode*,icPixelNode*> p = localpairs.back();
      localpairs.pop_back();
      icLightNode* lnode = p.first;
      icPixelNode* pnode = p.second;

      // check illuminationcut criteria
      STATS(stats.atomic.ub++);
      icIlluminationCut::IC code = checkIlluminationCut(lnode,pnode);

      // render or refine
      if(code==icGoodNoCheck)
      {
#if defined(UB_IMAGE)
        shadeSubtreeUB(lnode,pnode);
#endif
        STATS(stats.atomic.icpairs++);
        STATS(stats.atomic.clusters += pnode->end - pnode->start);
        STATS(stats.atomic.gpsize += pnode->end - pnode->start);
        continue;
      }
      else if (code == icGood)
      {
        if (!lnode->children[0] || !pnode->children[0])
        {
          STATS(stats.atomic.singletonpair++);
        }
        STATS(stats.atomic.icpairs++);
#  if defined(UB_IMAGE)
        shadeSubtreeUB(lnode,pnode);
#  else
        if (this->sampling)
          shadeSubtreeWithSamplingDescend(lnode,pnode,rnd,8);
        else
          shadeSubtreeFull(lnode,pnode);
#  endif
      } 
      else if(code == refineLight)
      {
        if (lnode->children[0])
          for (size_t i = 0 ;  i < icLightNode::nChildren ; ++i)
          {
            std::pair<icLightNode*,icPixelNode*> pp(lnode->children[i],pnode) ;
            localpairs.push_back(pp);
          }
      }
      else if(code == refineNode)
      {
        if (pnode->children[0])
          for (size_t i = 0 ;  i < icPixelNode::nChildren ; ++i)
          {
            std::pair<icLightNode*,icPixelNode*> pp(lnode,pnode->children[i]) ;
            localpairs.push_back(pp);
          }
      }
    }
  }

  Color icIlluminationCut::distributeMinColor( icPixelNode* pnode,Color parentColor )
  {
    float err = this->error;
    if (pnode->children[0] && pnode->children[1])
    {
      // recurse
      parentColor += pnode->color;
      pnode->color = Color(FLT_MAX);
      for (size_t i = 0 ; i < icPixelNode::nChildren ; ++i)
        if (pnode->children[i])
          pnode->color = min(pnode->color,distributeMinColor(pnode->children[i],parentColor));
    }
    else
    {
      // for historical reasons
      pnode->pixel->L = parentColor + pnode->color;
      pnode->color = err * pnode->pixel->L;
#ifndef APPROX_IMAGE
      pnode->pixel->L = Color(zero);
#endif
    }
    return pnode->color;
  }

  void icIlluminationCut::shadeRepPixel(icLightNode* lnode,icPixelNode* pnode, Random& rnd)
  {
    Sample3f wi; float tMax;
    icPixel* px = pnode->start + rnd.getInt(pnode->size());

#if REP_LIGHT_NUM > 1 
    Color I = lnode->repLights[px->s].sample(px->P,wi,tMax);
#else
    Color I = lnode->repLights[0].sample(px->P,wi,tMax);
#endif
    if ( I == Color(zero) || wi.pdf == 0.0f) return;

    I *= rcp(wi.pdf);
    /*! Evaluate BRDF */
    Color kdAngle = px->evalLambertianAngle(wi);
    Color ksAngle = px->evalSpecularAngle(wi);
    Color brdfWoKd = (px->kd == Color(zero))?Color(zero):px->kd*kdAngle;
    brdfWoKd += (px->ks == Color(zero))?Color(zero):px->ks*ksAngle;
    I *= brdfWoKd;
    if(I != Color(zero))
    {
      float err = px->error * 32.0f * float(ulp);
      Ray shadowRay(px->P, wi, err, tMax - err);
      STATS(stats.atomic.rays++);
      rtcOccluded(scene->scene,(RTCRay&)shadowRay);
      if (shadowRay) return;
    }
    pnode->color += I ;
    return;
  }

  icIlluminationCut::icIlluminationCut( float error,bool sampling,size_t spp) 
    :lighttree(NULL),error(error),sampling(sampling),exhaustiveRenderSize(spp>1?1:8)
  {
  }

  void icIlluminationCut::initPairsMT( std::deque<std::pair<icLightNode*, icPixelNode*> > & pairs, size_t siz )
  {
    auto morepixels = [] ( const std::pair<icLightNode*,icPixelNode*>& a, 
      const std::pair<icLightNode*,icPixelNode*>& b) 
    { 
      return (a.second->end - a.second->start) < (b.second->end - b.second->start); 
    };
    while ( pairs.size() && pairs.size() < siz )
    {
      // get the pair with the most pixels
      std::pop_heap(pairs.begin(),pairs.end(),morepixels);
      std::pair<icLightNode* , icPixelNode* > p = pairs.back(); 
      pairs.pop_back();

      icLightNode* lnode = p.first;
      icPixelNode* pnode = p.second;
      // subdivide only the pixels because otherwise there might be race conditions when MT
      if (pnode->children[0])
        for (size_t i = 0 ;  i < icPixelNode::nChildren ; ++i)
        {
          pairs.push_back( std::pair<icLightNode*,icPixelNode*>(lnode,pnode->children[i]));
          std::push_heap (pairs.begin(),pairs.end(),morepixels);
        }
      else
      {
        throw std::runtime_error("Not enough pixel nodes to fill the queue for MT rendering.");
      }
    }
  }

  void icIlluminationCut::geomcutRender( icLightNode* lnode,icPixelNode* pnode,Random& rnd )
  {
    // check geom condition
    icIlluminationCut::IC code = checkGeomCut(lnode,pnode);

    // refine or render according to the result
    if (code == icGood)
    {
      STATS(stats.atomic.geompairs++);
      shadeRepPixel(lnode,pnode,rnd);
    } 
    else if(code == refineLight)
    {
      // there should be no problem with leafs because then they are in the cut
      if (lnode->children[0])
        for (size_t i = 0 ;  i < icLightNode::nChildren ; ++i)
          geomcutRender(lnode->children[i],pnode,rnd);
    }
    else if(code == refineNode)
    {
      // there should be no problem with leafs because then they are in the cut
      if (pnode->children[0])
        for (size_t i = 0 ;  i < icPixelNode::nChildren ; ++i)
          geomcutRender(lnode,pnode->children[i],rnd);
    }
  }

 

}

