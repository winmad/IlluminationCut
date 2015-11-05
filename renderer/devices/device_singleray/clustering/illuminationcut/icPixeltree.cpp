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

#include "devices/device_singleray/clustering/illuminationcut/icPixeltree.h"
#include <algorithm>
#include <time.h>


namespace embree {

  size_t icPixelNode::subTreeSize()
  {
    size_t size=0;
    for (size_t i = 0; i < nChildren ; i++)
    {
      if(this->children[i])
        size += this->children[i]->subTreeSize();
    }
    return size + 1;
  }


  size_t icPixelNode::subTreeDepth()
  {
    size_t depth=0;
    for (size_t i = 0; i < nChildren ; i++)
      if(this->children[i])
        depth = max(this->children[i]->subTreeDepth(), depth);
    return depth + 1;
  }


  icPixelNode::icPixelNode() :parent(NULL),maxKd(FLT_MAX),maxKs(FLT_MAX),maxExp(FLT_MAX),minExp(zero),color(zero)
  {
    for (size_t i = 0 ; i < nChildren ; ++i)
      children[i] = NULL;
  }


  icPixelNode::~icPixelNode()
  {
    for (size_t i = 0 ; i < nChildren ; ++i)
      if (children[i])
        delete children[i];
  }


  embree::BBox6f icPixelNode::createBboxForChild( size_t i )
  {
    BBox6f b = this->cell;
    Vec6f diag = b.size();
    size_t dim; float pos;
    this->getNextDividingPlanePosition(pos,dim);
    b.lower[dim]  = ( i & 0x1 ) ? pos : b.lower[dim];
    diag[dim] *= 0.5f;
    b.upper = b.lower + diag;

    return b;
   }


  template<bool spatial>
  class separatorPredicate{
  public:
    separatorPredicate(float sep,size_t dim):sep(sep),dim(dim){}
    bool operator()(icPixel& p){return p.P[dim]<sep;}
    float sep;
    size_t dim;
  };
  template<>
  class separatorPredicate<false>{
  public:
    separatorPredicate(float sep,size_t dim):sep(sep),dim(dim % 3){}
    bool operator()(icPixel& p){return p.Ns[dim]<sep;}
    float sep;
    size_t dim;
  };


  void icPixeltree::createOctreeNonRecursiveAndCompress()
  {
    // create root
    this->root = new icPixelNode;
    this->root->level = 0;
    this->root->parent = NULL;
    this->root->cell = BBox6f(empty);
    for (size_t i = 0 ; i < this->pixels.size() ; ++i)
    {
      this->root->cell.grow(Vec6f(pixels[i].P,pixels[i].Ns));
    }
    root->cell.upper += Vec6f(Vector3f(2.0001f),Vector3f(0.15f));
    root->cell.lower -= Vec6f(Vector3f(2.0001f),Vector3f(0.05f));
    this->root->pixel = &(pixels[0]);
    this->root->start = &(pixels[0]);
    this->root->end = this->root->start + this->pixels.size();
    std::stack<icPixelNode*> nodesToProcess;
    nodesToProcess.push(this->root);
    Random rnd(RND_SEED);
    //Random rnd(int(time(0)));
    while ( !nodesToProcess.empty() )
    {
      // get next node to process
      icPixelNode* node = nodesToProcess.top(); nodesToProcess.pop();

      // put the lights in the proper child (bb)
      size_t dim ;
      float separator;
      node->getNextDividingPlanePosition(separator,dim);
      icPixel* middle;
      if (dim < 3)
      {
        separatorPredicate<true> sp(separator,dim);
        middle = std::partition(node->start,node->end,sp);
      } else
      {
        separatorPredicate<false> sp(separator,dim);
        middle = std::partition(node->start,node->end,sp);
      } 
      // node has to be compressed because there is only one child
      if (middle==node->start)
      {
        node->cell = node->createBboxForChild(1);
        node->level++;
        nodesToProcess.push(node);
        continue;
      }
      else if (middle==node->end)
      {
        node->cell = node->createBboxForChild(0);
        node->level++;
        nodesToProcess.push(node);
        continue;
      }
      // there are two new children
      for (size_t i = 0; i < nChildren; i++ )
      {
        node->children[i] = new icPixelNode;
        node->children[i]->level = node->level + 1;
        node->children[i]->parent = node;
        node->children[i]->cell = node->createBboxForChild(i);
      }

      node->children[0]->start = node->start;
      node->children[0]->end = middle;
      node->children[1]->start = middle;
      node->children[1]->end = node->end;
      icPixelNode* nrs[2] = {node->children[0],node->children[1]};

      // for each children create the subtree or stop or delete it
      for (size_t i = 0; i < nChildren; i++ )
      {
        if(nrs[i]->size() > 1 ) // not a leaf node, process it later
        {
          nodesToProcess.push(nrs[i]);
          nrs[i]->pixel = nrs[i]->start + rnd.getInt(nrs[i]->size());
        }
        else if(nrs[i]->size() == 1)
        {
          nrs[i]->pixel = nrs[i]->start;
          nrs[i]->cell = BBox6f(Vec6f(nrs[i]->pixel->P,nrs[i]->pixel->Ns));
        }
      }
    }
  }


  icPixeltree::~icPixeltree()
  {
    if (root)
      delete root;
  }


  icPixeltree::icPixeltree() :pixelnum(0),root(NULL)
  {

  }


  void icPixeltree::initPixelStorage( size_t pixels )
  {
    this->pixels.clear();this->pixels.resize(pixels);
  }


  void icPixeltree::shrink()
  {
    this->pixels.resize(pixelnum);
  }


  void icPixeltree::addPixel( DifferentialGeometry& dg,Vector3f& wo,Color& kd, Color& ks, float exp, size_t x,size_t y,size_t s )
  {
    this->pixels[pixelnum++] = icPixel(dg,wo,kd,ks,exp,x,y,s);
  }


  void icPixeltree::buildTree()
  {
    std::cout<<"Building pixeltree..."<<std::endl;
    this->createOctreeNonRecursiveAndCompress();
    this->setProperties(this->root);
    std::cout<<"Pixeltree size: "<<this->root->subTreeSize()<<std::endl;
  }


  void icPixeltree::setProperties( icPixelNode* pnode )
  {
    int counter = 0;
    for (size_t i = 0 ; i < icPixelNode::nChildren ; ++i)
      if (pnode->children[i])
        counter++;
    if (counter)
    {
      // recurse
      for (size_t i = 0 ; i < icPixelNode::nChildren ; ++i)
        if (pnode->children[i])
          setProperties(pnode->children[i]);
      pnode->maxKd = Color(zero);
      pnode->maxKs = Color(zero);
      pnode->maxExp = 0.0f;
      pnode->minExp = FLT_MAX;
      pnode->cell = BBox6f(empty);
      for (size_t i = 0 ; i < icPixelNode::nChildren ; ++i)
        if (pnode->children[i])
        {
          pnode->maxKd  = max(pnode->maxKd,pnode->children[i]->maxKd);
          pnode->maxKs  = max(pnode->maxKs,pnode->children[i]->maxKs);
          pnode->maxExp = max(pnode->maxExp,pnode->children[i]->maxExp);
          pnode->minExp = min(pnode->minExp,pnode->children[i]->minExp);
          pnode->cell.grow(pnode->children[i]->cell);
          pnode->posRadius = 0.5 * length(size(pnode->cell).pos());
        }
    }
    else
    {
      pnode->maxKd = pnode->pixel->kd;
      pnode->maxKs = pnode->pixel->ks;
      pnode->maxExp = pnode->pixel->exp;
      pnode->minExp = max(pnode->pixel->exp, 1.0f);
      pnode->cell = BBox6f(Vec6f(pnode->pixel->P,pnode->pixel->Ns));
      pnode->posRadius = 0.0f;
    }
    float rotationAngleNs = -acosf(clamp(dot(Vector3f(0,0,1) , pnode->pixel->Ns),-1.0f,1.0f));
    Vector3f rotationAxis = cross(Vector3f(0,0,1) , pnode->pixel->Ns);
    if (rotationAxis == zero) rotationAxis = Vector3f(0,1,0);
    rotationAxis=normalize(rotationAxis);
    pnode->rotationNs = LinearSpace3f::rotate(rotationAxis, rotationAngleNs);

    float rotationAngleR = -acosf(clamp(dot(Vector3f(0,0,1) , pnode->pixel->reflected_wo),-1.0f,1.0f));
    rotationAxis = cross(Vector3f(0,0,1) , pnode->pixel->reflected_wo);
    if (rotationAxis == zero) rotationAxis = Vector3f(0,1,0);
    rotationAxis=normalize(rotationAxis);
    pnode->rotationR = LinearSpace3f::rotate(rotationAxis, rotationAngleR);
  }


  icPixel::icPixel( DifferentialGeometry& dg,Vector3f wo,Color kd, Color ks, float exp, size_t x,size_t y,size_t s ) :
  P(dg.P),Ns(dg.Ns),error(dg.error),reflected_wo(reflect(wo,dg.Ns)),kd(kd),ks(ks),exp(exp),x(x),y(y),s(s),L(0,0,0)
  {

  }

}