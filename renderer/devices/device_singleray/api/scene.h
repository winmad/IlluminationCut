// ======================================================================== //
// Copyright 2009-2013 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.
//          //
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

#ifndef __EMBREE_BACKEND_SCENE_H__
#define __EMBREE_BACKEND_SCENE_H__

#include "handle.h"
#include "instance.h"

#include "../lights/light.h"
#include "devices/device_singleray/lights/trianglelight.h"
#include "devices/device_singleray/lights/spotlight.h"
#include "devices/device_singleray/lights/directionallight.h"
#include "devices/device_singleray/lights/hdrilight.h"
#include "devices/device_singleray/clustering/vpl/vplVPL.h"
#include "devices/device_singleray/clustering/lightcuts/ctLightcuts.h"
#include "devices/device_singleray/clustering/illuminationcut/icIlluminationCut.h"
#include "devices/device_singleray/clustering/lsLightcuts/lcLightTree.h"
#include "../shapes/differentialgeometry.h"
#include <fstream>
#include <time.h>

/*! include interface to ray tracing core */
#include <embree2/rtcore.h>

namespace embree
{
  /*! Scene holding all geometry and lights. */
  class BackendScene : public RefCount
  {
  public:

    class Handle : public InstanceHandle<BackendScene> {
      ALIGNED_CLASS;
    public:

      Handle () : accelTy("default"), builderTy("default"), traverserTy("default"),vplFile(""),vplNum(1000),epsilon(0.5),vpldepth(5),error(0.1) {}

      void set(const std::string& property, const Variant& data)
      {
        if      (property == "accel") accelTy = data.getString();
        else if (property == "builder") builderTy = data.getString();
        else if (property == "traverser") traverserTy = data.getString();
        else if (property == "vplfile") vplFile = data.getString();
        else if (property == "vplnumber") vplNum = data.getInt();
        else if (property == "epsilon") epsilon = data.getFloat();
        else if (property == "vpldepth") vpldepth = data.getInt();
        else if (property == "error") error = data.getFloat();
        else if (property == "clustering") clustering = data.getString();
      }

      virtual void setPrimitive(size_t slot, Ref<PrimitiveHandle> prim) = 0;

    public:
      std::string accelTy;
      std::string builderTy;
      std::string traverserTy;
      std::string vplFile;
      size_t vplNum;
      float epsilon;
      size_t vpldepth;
      float error;
      std::string clustering;
    };

  public:

    BackendScene (RTCScene scene)
      : scene(scene),clusterlightcuts(NULL),lcLighttree(NULL),mdLighttree(NULL) {}

    ~BackendScene () {
      if (scene) rtcDeleteScene(scene);
      if (clusterlightcuts) delete clusterlightcuts;
      if (lcLighttree) delete lcLighttree;
      if (mdLighttree) delete mdLighttree;
    }

    /*! Adds a light to the scene. */
    void add(const Ref<Light>& light) {
      allLights.push_back(light);
      if (Ref<EnvironmentLight> envlight = dynamic_cast<EnvironmentLight*>(light.ptr)) envLights.push_back(envlight);
    }

    /*! Helper to call the post intersector of the shape instance,
    *  which will call the post intersector of the shape. */
    virtual void postIntersect(const Ray& ray, DifferentialGeometry& dg) const = 0;

    /*! Generate VPLs for current light setup and/or read/write if the supplied VPL file is valid. */
    void generateVPL(size_t vplNum, const char* vplFile,size_t vpldepth);

    void preprocessLightcuts()
    {
      double t0 = getSeconds();
      clusterlightcuts = new ClusterLightCuts(this->vpls);
      clusterlightcuts->initClusterStructure(sceneradius);
      double dt = getSeconds() - t0;
      std::cout<<"Clustering preprocess for Lightcuts took: "<< dt * 1000.0f <<" ms"<<std::endl;
    }
    
    void preprocessWsiWSPD()
    {
      double t0 = getSeconds();
      lighttree = new icLighttree(this->vpls,this->sceneradius);
      lighttree->buildTree();
      double dt = getSeconds() - t0;
      std::cout<<"Clustering preprocess for Lighttree took: "<< dt * 1000.0f <<" ms"<<std::endl;
    }

    void preprocessLcLighttree();
    void preprocessMdLighttree();

  public:
    std::vector<Ref<Light> > allLights;              //!< All lights of the scene
    std::vector<Ref<EnvironmentLight> > envLights;   //!< Environment lights of the scene
    std::vector<vplVPL>   vpls;                      //!< All VPL lights
    RTCScene scene;
    ClusterLightCuts* clusterlightcuts;
    Ref<icLighttree> lighttree;
    lcLightTree* lcLighttree;
    MdLightTree* mdLighttree;
    float sceneradius;
    Vector3f scenecenter;
  };
}

#endif
