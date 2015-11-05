#include "devices/device_singleray/api/scene.h"
#include "devices/device_singleray/clustering/lsLightcuts/lcDivisiveLightTreeBuilder.h"



namespace embree{

  void BackendScene::preprocessLcLighttree()
  {
    double t0 = getSeconds();
    lcLighttree = DivisiveLightTreeBuilder().Build(this, this->vpls);
    double dt = getSeconds() - t0;
    std::cout<<"Clustering preprocess for lcLighttree took: "<< dt * 1000.0f <<" ms"<<std::endl;
  }

  void BackendScene::preprocessMdLighttree()
  {
    double t0 = getSeconds();
    mdLighttree = DivisiveMdLightTreeBuilder().Build(this, this->vpls);
    double dt = getSeconds() - t0;
    std::cout<<"Clustering preprocess for mdLighttree took: "<< dt * 1000.0f <<" ms"<<std::endl;
  }

  void BackendScene::generateVPL( size_t vplNum, const char* vplFile,size_t vpldepth )
  {
    this->vpls.clear();
    SET_CLAMPING_DISTANCE( 0.05 * sceneradius );

    //------------------------- READ FILE IF IT DOES EXIST--------------------------------------
    if( strcmp(vplFile,"")){
      FILE *ptr;
      ptr = fopen(vplFile,"rb");  
      if (ptr)
      {
        size_t siz;
        fread(&siz,sizeof(size_t),1,ptr);
        vpls.resize(siz);
        for (size_t i = 0 ; i<siz; ++i)
          vpls[i].read(ptr);

        fclose(ptr);
        std::cout << "Vpl file read. Total # of vpls:"<< vpls.size()<<std::endl;
        return;
      }
      else  std::cout << "Unable to read vpl file."<< std::endl;
    }

    //------------------------- OTHERWISE GENERATE VPLS--------------------------------------

    Random rnd(RND_SEED);

    // Place VPLs that substitute the original lights (direct VPLs)
    std::vector<Ref<Light> >::iterator it = allLights.begin();
    for ( ; it !=allLights.end(); ++it)
    {
      // Replace TriangleLight
      if( (*it).dynamicCast<TriangleLight>() )
      {
        Ref<TriangleLight> tr = (*it).dynamicCast<TriangleLight>();
        Color I = tr->L * length(cross(tr->e1,tr->e2)) * 0.5;
        int triangleSamples = 1024;
        // HACK to avoid too many VPLs for many quads...
        if (allLights.size() > 2000)
          triangleSamples = 100;
        for (int i=0 ; i < triangleSamples ; ++i)
        {
          Vector3f pos = uniformSampleTriangle(rnd.getFloat(),rnd.getFloat(),tr->v0,tr->v1,tr->v2);
          this->vpls.push_back(vplVPL(pos, -tr->Ng, I * rcp((float)triangleSamples) ) );
        }
      }
      // Replace SpotLight
      else if ((*it).dynamicCast<SpotLight>())
      {
        Ref<SpotLight> sl = (*it).dynamicCast<SpotLight>();
        if ( sl->cosAngleMax >= 0.05 || sl->cosAngleMin <= 0.95 )
          throw std::runtime_error("Not supporting spotlights with other than 90 degree falloff for VPL generation.");

        this->vpls.push_back(vplVPL(sl->P,sl->_D,sl->I));
      }
      // Replace environment map
      else if ((*it).dynamicCast<HDRILight>())
      {
        Ref<HDRILight> hd = (*it).dynamicCast<HDRILight>();
        int envmapSamples = 32768;
        float dirlightdistance = DIR_LIGHT_DISTANCE;
        for (int i=0 ; i < envmapSamples ; ++i)
        {
          Vector3f pos(uniformSampleSphere(rnd.getFloat(),rnd.getFloat()));
          Vector3f dir = -pos;
          pos *= dirlightdistance; pos += scenecenter;
          // CHECK below line later it is not completely precise
          Color I = hd->Le(-dir) * four_pi * rcp((float)envmapSamples) * dirlightdistance * dirlightdistance;               
          this->vpls.push_back(vplVPL(pos, -dir , I ));
        }
      }
      else if ((*it).dynamicCast<DirectionalLight>())
      {
        Ref<DirectionalLight> sl = (*it).dynamicCast<DirectionalLight>();
        float distance = 100000;
        this->vpls.push_back(vplVPL(sl->_wo*distance,sl->_wo,sl->E*distance*distance));
      }
      else
        throw std::runtime_error("Not supported light type for VPL generation.");
    }

    // Place VPLs  (indirect VPLs)
    size_t direct = this->vpls.size();
    size_t orig = this->allLights.size();
    if (vpldepth <=1) return;
    while (this->vpls.size() < (vplNum+direct))
    {
      // Sample orig lights for photons
      Ref<Light> vp = this->allLights[rnd.getInt(orig)];
      Vector3f pos; 
      Sample3f dir;
      Vec2f vr(rnd.getFloat(),rnd.getFloat());
      Vec2f vs(rnd.getFloat(),rnd.getFloat());

      Color I = vp->samplePhoton(dir,pos,vr,vs,sceneradius,scenecenter);
      if (dir.pdf==0.0f) continue;
      I *= (float) orig  *  rcp(dir.pdf);

      // Generate a path for the photon and add VPLs 
      int maxdepth = int(vpldepth) - 1;   // -1 because 1 for pt means only direct light
      while ( maxdepth-- && (I != zero ) )
      {
        Ray ray (pos,dir,32*std::numeric_limits<float>::epsilon(), inf);
        DifferentialGeometry dg;
        rtcIntersect(this->scene,(RTCRay&)ray);
        this->postIntersect(ray,dg);

        if (!ray)//nothing is hit
          break;
        /*! face forward normals */
        if (dot(dg.Ng, ray.dir) > 0)
        {
          dg.Ng = -dg.Ng; dg.Ns = -dg.Ns;
        }
        CompositedBRDF brdfs;
        if (dg.material) dg.material->shade(ray, Medium::Vacuum(), dg, brdfs);
        BRDFType diffuseBRDFTypes = (BRDFType)(DIFFUSE); 
        Vec2f s = Vec2f(rnd.getFloat(),rnd.getFloat());
        Color c = brdfs.sample(-ray.dir,dg, dir,diffuseBRDFTypes,s,rnd.getFloat(),diffuseBRDFTypes);
        if (c==zero || dir.pdf==0.0f) break;

        // the dot prod should not be zero since then c would be zero
        float dott = dot(dir,dg.Ns);
        if (dott==0.0f) break;
        Color contrib = I * c * rcp(dott) * rcp((float)vplNum) ;
        //if (std::isnan(contrib.r) || std::isnan(contrib.g) || std::isnan(contrib.b))
        //  break;

        this->vpls.push_back(vplVPL(dg.P,-dg.Ns,contrib));

        Color contribScale =  c * rcp(dir.pdf);
        float rrProb = min(1.0f, reduce_avg(contribScale));
        if(rrProb == 0.0f)break;
        if (rnd.getFloat() > rrProb)
          break;

        I *= contribScale * rcp(rrProb);
        pos = dg.P;
      }
    }

    //------------------------------PERTURB VPLS WITH THE SAME POSITION-----------------------------------------------------------
    auto vplComp = [] (const vplVPL& vp1,const vplVPL& vp2) { return vp1.P < vp2.P; };
    std::sort(this->vpls.begin(),this->vpls.end(),vplComp);
    std::vector<vplVPL>::iterator itt = this->vpls.begin();
    while(itt != ( vpls.end() - 1 ) && itt != vpls.end())
    {
      vplVPL& v1 = *itt, & v2 = *(itt+1);
      if(v1.P == v2.P)
      {
        v2.I += v1.I;
        itt = this->vpls.erase(itt);
      }
      else
        ++itt;
    }

    //------------------------------WRITE OUT VPLS INTO A FILE-----------------------------------------------------------

    // writing out the VPLs
    FILE *ptr;
    ptr = fopen(vplFile,"wb");  
    if (ptr)
    {
      size_t siz = this->vpls.size();
      fwrite(&siz,sizeof(size_t),1,ptr);
      if (vpls.size())
        for (size_t i = 0 ; i<siz; ++i)
          vpls[i].write(ptr);

      fclose(ptr);
      std::cout << "Vpl file written. Total # of vpls:"<< vpls.size()<<std::endl;
    }
    else std::cout << "Unable to write vpl file. Total # of vpls: "<<vpls.size() << std::endl;

    return;
  }

}