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

#include "devices/device_singleray/clustering/lightcuts/ctLightcuts.h"
#include <vector>
#include <list>
#include <queue>
#include <time.h>

namespace embree
{
  std::ostream & operator<<(std::ostream &os, const Node& n)
  {
    os<<">>>>"  <<  std::endl;
    os << "address: " << &n << std::endl;
    os << "cell: " << n.cell << std::endl;
    os << "kids: ";
    os<<" "<<n.rc;
    os<<" "<<n.lc  <<  std::endl;
    os << "light: " <<  std::endl;
    os << "dir: " << n.cone.dir << std::endl;
    os << "halfangle: " << n.cone.halfAngle << std::endl;
    os<<">>>>" <<  std::endl;
    
    os << std::endl;
    return os;
  }

  
Node* ClusterLightCuts::buildTreeFromLeavesNaive(  )
  {
    std::cout<<"Naive agglomerative clustering"<<std::endl;
    // create a list for the leaves
    std::list<Node*> leafs;

    // get the lights according to the vector and add them as leaves
    for(size_t i=0; i<vpls.size(); ++i)
    {
      vplVPL leaf[REP_LIGHT_NUM];
      bool fromright[REP_LIGHT_NUM];
      for (size_t j = 0 ; j < REP_LIGHT_NUM ; ++j)
      {
        leaf[j] = vpls[i];
        fromright[j] = false;
      }
      // create the corresponding node with the id of the light
      Node* n = new Node(leaf,leaf[0].P,Cone( - leaf[0]._D),vpls,fromright);
      leafs.push_back(n);
    }
    // random number provider
    Random rnd;
    rnd.setSeed((int)time(0));
    // no heap and other optimizations
    while(leafs.size()!=1)
    {
      // naive, O(n^3) clustering
      std::list<Node*>::iterator it= leafs.begin(),it2=leafs.begin();
      Node *na,*nb;
      float mindist=FLT_MAX;
      for( ; it!=leafs.end(); it++){
        for( ; it2!=leafs.end(); it2++)
        {   
          if(*it==*it2)continue;
          float dist = (*it)->dist((*it2),sceneradiusSqr);
          if(dist<mindist) { mindist=dist; na=*it; nb=*it2; }
        }
      }
      Node* m  = na->merge(nb,rnd);

      leafs.remove(na);
      leafs.remove(nb);
      leafs.push_back(m);   
    }

    return *leafs.begin();
  }

  
  struct Triple
  {
    Triple(float d,size_t a, size_t b):dist(d),cl1(a),cl2(b){}
    friend bool operator< ( const Triple& a, const Triple& b);
    float dist;
    size_t cl1, cl2;
  };
  

  bool operator< (const Triple& a, const Triple& b){return a.dist > b.dist;} 


  Node* ClusterLightCuts::buildTreeFromLeavesHeap(  )
  {
    std::cout<<"Heap based agglomerative clustering"<<std::endl;
    // create a list for the leaves
    // bool indicates that it is still a valid cluster (haven't been merged yet)
    std::vector<std::pair<Node*,bool> > clusters;
    clusters.reserve( 2 * vpls.size() );

    std::priority_queue< Triple, std::vector<Triple> > distances;

    // get the lights according to the vector and add them as leaves
    for(size_t i=0; i<vpls.size(); ++i)
    {
      vplVPL leaf[REP_LIGHT_NUM];
      bool fromright[REP_LIGHT_NUM];
      for (size_t j = 0 ; j < REP_LIGHT_NUM ; ++j)
      {
        leaf[j] = vpls[i];
        fromright[j] = false;
      }
      // create the corresponding node with the id of the light
      Node* n = new Node(leaf,leaf[0].P,Cone( - leaf[0]._D),vpls,fromright);
      clusters.push_back(std::pair<Node*,bool>(n,true));
    }
    //initialize the distances
    for (size_t i = 0; i < clusters.size(); i++){
      float mindist = FLT_MAX;
      size_t minIdx = (size_t)-1;
      for (size_t j = 0; j < clusters.size(); j++){
        float dd = clusters[i].first->dist(clusters[j].first,sceneradiusSqr);
        if( mindist > dd ){
          mindist = dd;  minIdx = j;
        }
      }
      distances.push(Triple(mindist,i,minIdx));   // not the best since nlogn heap building
    }
    
    // random number provider
    Random rnd;
    rnd.setSeed((int)time(0));
   
    while(distances.size() > 1)
    {
      Triple tr = distances.top();
      distances.pop();
      // ... do something
      if ( !clusters[tr.cl1].second )
      {
        // cl1 was already clustered with somebody else
        // we just pop it and continue
      } else if ( !clusters[tr.cl2].second )
      {
          float minDist = FLT_MAX;
          size_t minIdx = (size_t) -1;
          //search for best pair of cl1;
          for (size_t j = 0; j < clusters.size(); j++)
          {
            if (clusters[j].second)
            {
              float dist = clusters[tr.cl2].first->dist(clusters[j].first,sceneradiusSqr);
              if (minDist > dist)
              {
                minDist = dist; minIdx = j;
              }
            }
          }
          distances.push(Triple(minDist,tr.cl1,minIdx));
      } else
      {
        if (tr.cl1 != tr.cl2)
        {
          clusters[tr.cl1].second = false;
          clusters[tr.cl2].second = false;
          Node* newcluster = clusters[tr.cl1].first->merge(clusters[tr.cl2].first,rnd);
          clusters.push_back(std::pair<Node*,bool>(newcluster,true));
          clusters[tr.cl1].first->parent = newcluster;
          clusters[tr.cl2].first->parent = newcluster;
          //find best pair
          float minDist = FLT_MAX;
          size_t minIdx = (size_t) -1;
          //search for best pair of cl1;
          // NAIVE SEARCH, not the kdtree version, should change it now
          for (size_t j = 0; j < clusters.size(); j++)
          {
            if (clusters[j].second)
            {
              float dist = newcluster->dist(clusters[j].first,sceneradiusSqr);
              if (minDist > dist)
              {
                minDist = dist; minIdx = j;
              }
            }
          }
          distances.push(Triple(minDist,clusters.size()-1,minIdx));
        }
      }
    }
    //taking care of last pair;
    Triple tr = distances.top();
    clusters[tr.cl1].second = false;
    clusters[tr.cl2].second = false;
    Node* newcluster = clusters[tr.cl1].first->merge(clusters[tr.cl2].first,rnd);
    clusters.push_back(std::pair<Node*,bool>(newcluster,true));
    clusters[tr.cl1].first->parent = newcluster;
    clusters[tr.cl2].first->parent = newcluster;
       
    return clusters.back().first;
  }


  Node* ClusterLightCuts::buildTreeFromLeavesLocal()
  { 
    // Fast Agglomerative Clustering for Rendering, Walter et al, 2008
    // Locally-ordered agglomerative clustering
    std::cout<<"Local agglomerative clustering"<<std::endl;

    size_t numNodes = 2 * vpls.size() - 1 ;      
    size_t numCurrentClusters = 0; 
    std::vector<Node*> clusters ; clusters.resize(numNodes);

    // create the initial clusters, corresponding to single vpls
    for(size_t i=0; i<vpls.size(); ++i){
      vplVPL leaf[REP_LIGHT_NUM];
      bool fromright[REP_LIGHT_NUM];
      for (size_t j = 0 ; j < REP_LIGHT_NUM ; ++j)
      {
        leaf[j] = vpls[i];
        fromright[j] = false;
      }
      // create the corresponding node with the id of the light
      Node* n = new Node(leaf,leaf[0].P,Cone( - leaf[0]._D),vpls,fromright);
      clusters[i] = n; 
      numCurrentClusters++;
    }

    // construct kdtree
    KdTree<Node> kdtree(this->sceneradiusSqr);
    kdtree.constructTree(clusters,numCurrentClusters);

    Random rnd;
    rnd.setSeed((int)time(0));

    Node* na = clusters[0];
    Node* nb  = kdtree.queryNN(na);

    if(!nb) throw std::runtime_error("local agglomerative clustering error");
    while ( numCurrentClusters  <  numNodes)
    {
      Node* nc = kdtree.queryNN(nb);
      if(!nc) throw std::runtime_error("local agglomerative clustering error");

      if ( na == nc )
      {
        na->valid = false;
        nb->valid = false;
        Node* naOld = na;
        na = na->merge(nb,rnd);
        kdtree.invalidateNodesAndUpdate(naOld,nb,na);
        clusters[numCurrentClusters++] = na;
        nb = NULL;
        nb = kdtree.queryNN(na);
      } else {
        na = nb;
        nb = nc;
      }
    }
    return clusters[numNodes - 1];
  }

  void ClusterLightCuts::initClusterStructure( float sceneradius )
  {
    this->sceneradiusSqr = sceneradius*sceneradius;
    root = this->buildTreeFromLeavesLocal();
    return;
  }



  Node* Node::merge( Node* node, Random& rnd )
  {
    // setup representative light 
    //randomly pick one according to the intensity
    float sum =reduce_add(this->repLights[0].I + node->repLights[0].I);
    if (sum < 0.000001f)
      sum = 0.000001f;

    vplVPL lights[REP_LIGHT_NUM];
    bool fromRight[REP_LIGHT_NUM];
    for (size_t i = 0 ; i< REP_LIGHT_NUM ; ++i)
    {
      if (rnd.getFloat() < reduce_add(this->repLights[i].I)/sum ){
        lights[i] = this->repLights[i];
        fromRight[i] = false;
      }
      else{
        lights[i] = node->repLights[i];
        fromRight[i] = true;
      }
      lights[i].I = this->repLights[i].I + node->repLights[i].I;
    }
    Node* parent = new Node(lights,embree::merge(this->cell,node->cell),this->cone.merge(node->cone),vpls,fromRight,this,node);
    return parent;
  }

  Node::Node( vplVPL light[REP_LIGHT_NUM],BBox3f bb,Cone cone,std::vector<vplVPL>& vpls, bool right[REP_LIGHT_NUM], Node* lc/*=NULL*/, Node* rc/*=NULL*/ ) 
    :cell(bb),cone(cone),lc(lc),rc(rc),vpls(vpls),valid(true)
  {
    for (size_t i = 0 ; i < REP_LIGHT_NUM; ++i)
    {
      repLights[i] = light[i];
      repLightFromRightChild[i] = right[i];
    }
  }

  Node::~Node()
  {
    if(lc) {delete lc;lc=NULL;}
    if(rc) {delete rc;rc=NULL;}
  }

  Cone::Cone( Vector3f dir,float halfAngle/*=0*/ ) 
    :dir(normalize(dir)),halfAngle(halfAngle),cosHalfAngle(cos(halfAngle)),rotationAngle(-acosf(clamp(dot(Vector3f(0,0,1) , dir),-1.0f,1.0f)))
  {

  }

}
