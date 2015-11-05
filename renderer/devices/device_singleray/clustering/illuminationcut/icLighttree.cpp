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
#include "devices/device_singleray/clustering/illuminationcut/icLighttree.h"
#include <algorithm>
#include <queue>
#include <time.h>


namespace embree{

  void icLighttree::buildTreeAgglomeraticeLocal()
  { 
    // Fast Agglomerative Clustering for Rendering, Walter et al, 2008
    // Locally-ordered agglomerative clustering
    std::cout<<"Local agglomerative clustering"<<std::endl;

    // clusters
    size_t numNodes = 2 * vpls.size() - 1 ;      
    size_t numCurrentClusters = 0; 
    std::vector<icLightNode*> clusters ; clusters.resize(numNodes);

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
      icLightNode* n = new icLightNode(leaf,Vec6f(leaf[0].P,Vector3f(0,0,0)),icCone( - leaf[0]._D),fromright);
      clusters[i] = n; 
      numCurrentClusters++;
    }

    // construct kdtree
    KdTree<icLightNode> kdtree(this->sceneradiusSqr);
    kdtree.constructTree(clusters,numCurrentClusters);

    icLightNode* na = clusters[0];
    // find best match
    icLightNode* nb  = kdtree.queryNN(na);

    if(!nb) throw std::runtime_error("local agglomerative clustering error");
    while ( numCurrentClusters  <  numNodes)
    {
      icLightNode* nc = kdtree.queryNN(nb);
      if(!nc) throw std::runtime_error("local agglomerative clustering error");

      if ( na == nc )
      {
        na->valid = false;
        nb->valid = false;
        icLightNode* naOld = na;
        na = na->merge(nb,this->rnd);
        kdtree.invalidateNodesAndUpdate(naOld,nb,na);
        clusters[numCurrentClusters++] = na;
        nb = NULL;
        nb = kdtree.queryNN(na);
      } else {
        na = nb;
        nb = nc;
      }
    }
    this->root = clusters[numNodes - 1];
  }


  icLightNode* icLightNode::merge( icLightNode* node, Random& rnd )
  {
    // setup representative light 

    //randomly pick one according to the intensity
    float sum =max( 0.000001f, reduce_add(this->repLights[0].I + node->repLights[0].I) );

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
    icLightNode* parent = new icLightNode(lights,embree::merge(this->cell,node->cell),this->cone.mergeCones(node->cone),fromRight,this,node);
    return parent;
  }


  size_t icLightNode::subTreeSize()
  {
    size_t size=0;
    for (size_t i = 0; i < nChildren ; i++)
    {
      if(this->children[i])
        size += this->children[i]->subTreeSize();
    }
    return size + 1;
  }


  size_t icLightNode::subTreeDepth()
  {
    size_t depth=0;
    for (size_t i = 0; i < nChildren ; i++)
      if(this->children[i])
        depth = max(this->children[i]->subTreeDepth(), depth);
    return depth + 1;
  }


  icLightNode::icLightNode() :parent(NULL)
  {
    for (size_t i = 0 ; i < nChildren ; ++i)
      children[i] = NULL;
  }


  icLightNode::icLightNode( vplVPL light[REP_LIGHT_NUM],BBox6f bb,icCone cone, bool right[REP_LIGHT_NUM], 
                                  icLightNode* lc/*=NULL*/, icLightNode* rc/*=NULL*/ ) 
    : cell(bb),cone(cone),valid(true)
  {
    children[0]=lc;
    children[1]=rc;
    for (size_t i = 0 ; i < REP_LIGHT_NUM; ++i)
    {
      repLights[i] = light[i];
      repLightFromRightChild[i] = right[i];
    }
  }


  icLightNode::~icLightNode()
  {
    for (size_t i = 0 ; i < nChildren ; ++i)
      if (children[i])
        delete children[i];
  }


  void icLightNode::getNextDividingPlanePosition( float& pos, size_t& dim )
  {
    BBox6f& b = this->cell;
    Vec6f diag = b.size();
    float dirDiam = length(diag.dir());
    bool dirDivision = dirDiam > MAXSIZE_DIRECTIONAL_BBOX;

    // return the middle point of the bb along the next divided dimension
    if (!dirDivision)
    {
      dim = this->level % 3 ;
      pos = b.lower[dim]  + diag[dim] * 0.5f;
    }
    else 
    {
      dim = 3 + this->level % 3 ;
      pos = b.lower[dim]  + diag[dim] * 0.5f;
    }
    return;
  }


  embree::BBox6f icLightNode::createBboxForChild( size_t i )
  {
    BBox6f b = this->cell;
    Vec6f diag = b.size();
    size_t dim; float pos;
    this->getNextDividingPlanePosition(pos,dim);
    b.lower[dim]  = ( i & 0x1 ) ? pos : b.lower[dim] ;
    diag[dim] *= 0.5f;
    b.upper = b.lower + diag;
    return b;
  }


   template<bool spatial>
  class separatorPredicate{
  public:
    separatorPredicate(float sep,size_t dim):sep(sep),dim(dim){}
    bool operator()(vplVPL& vpl){return vpl.P[dim]<sep;}
    float sep;
    size_t dim;
  };
  template<>
  class separatorPredicate<false>{
  public:
    separatorPredicate(float sep,size_t dim):sep(sep),dim(dim % 3){}
    bool operator()(vplVPL& vpl){return vpl._D[dim]<sep;}
    float sep;
    size_t dim;
  };


  void icLighttree::createOctreeNonRecursiveAndCompress()
  {
    std::cout<<"Divisive spatial clustering"<<std::endl;
    // root containing all vpls with a big bb
    this->root = new icLightNode;
    this->root->level = 0;
    this->root->parent = NULL;
    this->root->cell = BBox6f(empty);
    for (size_t i = 0 ; i < this->vpls.size() ; ++i)
    {
      this->root->cell.grow(Vec6f(vpls[i].P,vpls[i]._D));
    }
    root->cell.upper += Vec6f(Vector3f(2.0001f),Vector3f(0.15f));
    root->cell.lower -= Vec6f(Vector3f(2.0001f),Vector3f(0.05f));
    for(size_t i = 0 ; i < REP_LIGHT_NUM ; ++i)
    {
      root->repLights[i] = vpls[0];
      root->repLights[i].I = Color(FLT_MAX);
    }

    // stack of vpl ranges to
    std::stack<NodeRange> nodesToProcess;
    nodesToProcess.push(NodeRange(this->root,&(vpls[0]),&(vpls[0]) + vpls.size()));
    while ( !nodesToProcess.empty() )
    {
      // get next node to process
      NodeRange noderange = nodesToProcess.top(); nodesToProcess.pop();
      icLightNode* node = noderange.node;
     
      // put the lights in the proper child (bb)
      size_t dim ;   float separator;
      node->getNextDividingPlanePosition(separator,dim);
      vplVPL* middle;
      if (dim < 3)
      {
        separatorPredicate<true> sp(separator,dim);
        middle = std::partition(noderange.start,noderange.end,sp);
      } else
      {
        separatorPredicate<false> sp(separator,dim);
        middle = std::partition(noderange.start,noderange.end,sp);
      } 
      // node has to be compressed because there is only one child
      if (middle==noderange.start)
      {
        node->cell = node->createBboxForChild(1);
        node->level++;
        nodesToProcess.push(noderange);
        continue;
      }
      else if (middle==noderange.end)
      {
        node->cell = node->createBboxForChild(0);
        node->level++;
        nodesToProcess.push(noderange);
        continue;
      }
      // there are two new children, create new node
      for (size_t i = 0; i < nChildren; i++ )
      {
        node->children[i] = new icLightNode;
        node->children[i]->level = node->level + 1;
        node->children[i]->parent = node;
        node->children[i]->cell = node->createBboxForChild(i);
      }

      NodeRange node1(node->children[0],noderange.start,middle);
      NodeRange node2(node->children[1],middle,noderange.end);
      NodeRange nrs[2] = {node1,node2};

      // for each children add it to the stack or create a leaf
      for (size_t i = 0; i < nChildren; i++ )
      {
        if(nrs[i].size() > 1 ) // not a leaf node, process it later
        {
          nodesToProcess.push(nrs[i]);
        }
        else if(nrs[i].size() == 1)
        {
          for (size_t k = 0; k < REP_LIGHT_NUM ; ++k)
            nrs[i].node->repLights[k] = *(nrs[i].start);
          nrs[i].node->cell = BBox6f(Vec6f(nrs[i].node->repLights[0].P,nrs[i].node->repLights[0]._D));
        }
      }
    }
  }


  void icLighttree::setProperties( icLightNode* lnode )
  {
    // decide if leaf
    int numChildren = 0;
    for (size_t i = 0 ; i < icLightNode::nChildren ; ++i)
      if (lnode->children[i])
        numChildren++;

    if (numChildren)
    {
      // recurse
      for (size_t i = 0 ; i < icLightNode::nChildren ; ++i)
        if (lnode->children[i])
          setProperties(lnode->children[i]);

      // set cell and cone
      lnode->cell = BBox6f(empty);
      lnode->cone = lnode->children[0]->cone;
      for (size_t i = 0 ; i < icLightNode::nChildren ; ++i)
        if (lnode->children[i])
        {
          lnode->cell.grow(lnode->children[i]->cell);
          lnode->cone = lnode->cone.mergeCones(lnode->children[i]->cone);
        }
      // set posradius
      lnode->posRadius = 0.5f * length(size(lnode->cell).pos());
      // set replights
      float ratio = reduce_add(lnode->children[0]->repLights[0].I);
      Color I = lnode->children[0]->repLights[0].I + lnode->children[1]->repLights[0].I;
      ratio *= rcp(reduce_add(I));
      for (size_t k = 0 ; k < REP_LIGHT_NUM ; ++k)
      {
        if (ratio > rnd.getFloat())
          lnode->repLights[k] = lnode->children[0]->repLights[k];
        else
          lnode->repLights[k] = lnode->children[1]->repLights[k];
        lnode->repLights[k].I = I;
      }
    }
    else
    {
      // setup leafs
      lnode->cell = BBox6f(Vec6f(lnode->repLights[0].P,lnode->repLights[0]._D));
      lnode->cone = icCone(-lnode->repLights[0]._D);
      lnode->posRadius = 0;
    }
  }


  void icLighttree::buildTree()
  {
    std::cout<<"Building lighttree..."<<std::endl;
    rnd.setSeed(RND_SEED);
#ifdef AGGLOMERATIVE_LIGHTTREE
    this->buildTreeAgglomeraticeLocal();
#else
    this->createOctreeNonRecursiveAndCompress();
#endif
    this->setProperties(this->root);
    std::cout<<"Lighttree size: "<<this->root->subTreeSize()<<std::endl;
  }



  icCone::icCone( Vector3f dir,float halfAngle/*=0*/ ) :
    axis(normalize(dir)),
    halfAngle(halfAngle),
    cosHalfAngle(cos(halfAngle))
  {
    // initialize the rotation
    Vector3f rotationAxis = cross(Vector3f(0,0,1) , this->axis);
    float rotationAngle = -acosf(clamp(dot(Vector3f(0,0,1) , dir),-1.0f,1.0f));
    if (rotationAxis == zero) rotationAxis = Vector3f(0,1,0);
    rotationAxis=normalize(rotationAxis);
    this->rotation = LinearSpace3f::rotate(rotationAxis,rotationAngle);
  }

}