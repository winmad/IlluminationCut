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

#ifndef kdtree_h__
#define kdtree_h__

#include "devices/device_singleray/default.h"
#include "common/math/bbox.h"
#include <vector>
#include <algorithm>

namespace embree{

  template<class Node>
  class KdNode{
  public:
    KdNode():depth(0),cluster(NULL),intensity(zero),bbox(zero),left(NULL),right(NULL),valid(false),parent(NULL){}
    ~KdNode(){
      if (left)
        delete left;
      if (right)
        delete right;
    }

    int depth;
    Node* cluster;
    Color intensity;
    typename Node::BB bbox;

    KdNode* left;
    KdNode* right;
    bool valid;
    KdNode* parent;
  };


  template<class Node>
  class KdTree{
  public:

    KdNode<Node>* root;
    float sceneradiusSqr;

    KdTree(float sr):root(NULL),sceneradiusSqr(sr){

    }

    ~KdTree(){
      if (root)
        delete root;
    }

    void constructTree(std::vector<Node*>& leafs,size_t numLeafs);
    KdNode<Node>* constructSubtree(std::vector<Node*>& points,int level);
    Node* queryNN(Node* node); 
    void invalidateNodesAndUpdate(Node* a, Node* b,Node* newNode);
    bool deleteInvalidLeaf(KdNode<Node>* leaf);
    void deleteSubtree(KdNode<Node>* node);
    void setParents(KdNode<Node>* node);
    void queryParent( KdNode<Node>* parent,KdNode<Node>* child, float& upperbound, float& currentdist, KdNode<Node>*& currentbest, KdNode<Node>* queryNode, int& explored );
    void queryChild( KdNode<Node>* parent, bool left , float& upperbound, float& currentdist, KdNode<Node>*& currentbest, KdNode<Node>* queryNode, int& explored );
  };



  ///////////////////////////////////////////// IMPLEMENTATION
  template<class Node>
  struct levComp {
    int level;
    static float dist(KdNode<Node>* split, KdNode<Node>* current, int level){
      switch ( level % 3 ){
      case 0:
        return (current->bbox.lower.x - split->bbox.lower.x);
      case 1:
        return (current->bbox.lower.y - split->bbox.lower.y);
      case 2:
        return (current->bbox.lower.z - split->bbox.lower.z);
      }

      throw std::runtime_error("compare error in kdtree construction");
      return false;
    }
    levComp(int level):level(level){}
    bool operator() (Node* i,Node* j) {
      switch ( level % 3 ){
      case 0:
        return (i->cell.lower.x < j->cell.lower.x);
      case 1:
        return (i->cell.lower.y < j->cell.lower.y);
      case 2:
        return (i->cell.lower.z < j->cell.lower.z);
      }

      throw std::runtime_error("compare error in kdtree construction");
      return false;
    }
  };

  template<class Node>
  void KdTree<Node>::constructTree(std::vector<Node*>& leafs,size_t numLeafs){
    std::vector<Node*> points;
    for (size_t i = 0; i< numLeafs; ++i)
      points.push_back(leafs[i]);

    this->root = constructSubtree(points,0);
    this->setParents(this->root);
    return;
  }


  template<class Node>
  KdNode<Node>* KdTree<Node>::constructSubtree(std::vector<Node*>& points,int level){
    if (points.size() <= 1)
    {
      if (points.size()==0)
        return NULL;
      KdNode<Node>* node = new KdNode<Node>;
      node->valid = true;
      node->depth = level;

      node->cluster = points[0];
      node->intensity = node->cluster->repLights[0].I;
      node->bbox = node->cluster->cell;
      node->cluster->kdnode = node;
      return node;
    }
    KdNode<Node>* node = new KdNode<Node>;
    node->valid = true;
    node->depth = level;

    size_t medianPlace = points.size()/2 ;
    std::sort(points.begin(),points.end(),levComp<Node>(level));

    // FIXME there is no need to create a new vector, should just pass the iterators
    std::vector<Node*> left(points.begin(), points.begin()+medianPlace);
    std::vector<Node*> right(points.begin()+medianPlace+1, points.end());

    node->cluster = points[medianPlace];
    node->intensity = node->cluster->repLights[0].I;
    node->bbox = node->cluster->cell;
    node->cluster->kdnode = node;
    node->left = constructSubtree(left,level + 1);
    node->right = constructSubtree(right,level + 1);

    return node;
  }

  template<class Node>
  void KdTree<Node>::setParents( KdNode<Node>* node )
  {
    if (node->left)
    {
      node->left->parent = node;
      this->setParents(node->left);
    }
    if (node->right)
    {
      node->right->parent = node;
      this->setParents(node->right);
    }
  }


  template<class Node>
  Node* KdTree<Node>::queryNN( Node* node)
  {
    // get the corresponding kdnode
    KdNode<Node>* queryNode = node->kdnode;
    //we have an upper bound on the diagonal alpha<sqrt(distmin/I)
    //initial guess root, and check other if equality
    KdNode<Node>* currentbest = queryNode;
    int explored = 0;

    // modified values for search start
    float currentdist = std::numeric_limits<float>::infinity();
    float upperbound = std::numeric_limits<float>::infinity();
    //check children
    this->queryChild( queryNode, true, upperbound, currentdist, currentbest, queryNode, explored );
    this->queryChild( queryNode, false, upperbound, currentdist, currentbest, queryNode, explored);
    // check upwards
    this->queryParent( queryNode->parent, queryNode, upperbound, currentdist, currentbest, queryNode ,explored);
    // return the best one;
    if (currentbest == queryNode)
    {
      return NULL;
    }
    return currentbest->cluster;

  }

  template<class Node>
  void KdTree<Node>::invalidateNodesAndUpdate( Node* a, Node* b,Node* newNode )
  {
    // invalidate the lower one, update intensity, assign new cluster
    if (a->kdnode->depth > b->kdnode->depth)
    {
      //invalidate kdnode
      a->kdnode->valid = false;
      b->kdnode->intensity += a->kdnode->intensity;
      this->deleteSubtree(a->kdnode);

      newNode->kdnode = b->kdnode;
      b->kdnode->cluster = newNode;
    }
    else
    {
      b->kdnode->valid = false;
      a->kdnode->intensity += b->kdnode->intensity;
      this->deleteSubtree(b->kdnode);

      newNode->kdnode = a->kdnode;
      a->kdnode->cluster = newNode;
    }
    return;
  }


  template<class Node>
  void KdTree<Node>::queryParent( KdNode<Node>* parent,KdNode<Node>* child, float& upperbound, float& currentdist, KdNode<Node>*& currentbest, KdNode<Node>* queryNode, int& explored )
  {
    // this the root (hopefully)
    if (!parent)
      return;
    //recurse into other child
    bool otherchild = (child != parent->left); 
    this->queryChild(parent,otherchild,upperbound,currentdist,currentbest,queryNode,explored);
    // recurse upwards
    this->queryParent(parent->parent,parent,upperbound,currentdist,currentbest,queryNode,explored);
    // end
    return;
  }
  template<class Node>
  void KdTree<Node>::queryChild( KdNode<Node>* parent, bool left , float& upperbound, float& currentdist, KdNode<Node>*& currentbest, KdNode<Node>* queryNode, int& explored )
  {

    if (!parent)
      return;

    explored++;
    //update mindist with parent
    if (parent->valid && parent!=queryNode)
    {
      float tmpDist = queryNode->cluster->dist(parent->cluster,this->sceneradiusSqr);
      if (currentdist > tmpDist)
      {
        currentdist = tmpDist;
        currentbest = parent;
        upperbound  = sqrt(currentdist/reduce_add(queryNode->intensity));
      }
    }

    float levelDistParent = levComp<Node>::dist(parent,queryNode, parent->depth); // negative means it is on the left
    if (left)
    {
      if (levelDistParent<=0)
      {
        this->queryChild(parent->left,true,upperbound,currentdist,currentbest,queryNode,explored);
        this->queryChild(parent->left,false,upperbound,currentdist,currentbest,queryNode,explored);
      }else
      {
        if (levelDistParent > upperbound)
        {
          return;
        }else
        {
          //recurse
          this->queryChild(parent->left,true,upperbound,currentdist,currentbest,queryNode,explored);
          this->queryChild(parent->left,false,upperbound,currentdist,currentbest,queryNode,explored);
        }
      }

    }
    if (!left)
    {
      if (levelDistParent>=0)
      {
        //recurse
        this->queryChild(parent->right,true,upperbound,currentdist,currentbest,queryNode,explored);
        this->queryChild(parent->right,false,upperbound,currentdist,currentbest,queryNode,explored);
      }else
      {
        if (-levelDistParent > upperbound)
        {
          return;
        }else
        {
          //recurse
          this->queryChild(parent->right,true,upperbound,currentdist,currentbest,queryNode,explored);
          this->queryChild(parent->right,false,upperbound,currentdist,currentbest,queryNode,explored);
        }
      }

    }
    return;
  }

  template<class Node>
  bool KdTree<Node>::deleteInvalidLeaf( KdNode<Node>* leaf )
  {
    if (!leaf)
      return true;

    bool deletedLeaf = false;
    // if invalid leaf 
    if (!leaf->valid && !leaf->right && !leaf->left)
    {
      bool left = (leaf->parent->left == leaf);
      if (left)
        leaf->parent->left =NULL;
      else
        leaf->parent->right = NULL;

      deletedLeaf = true;
      delete leaf;
    }
    return deletedLeaf;
  }

  template<class Node>
  void KdTree<Node>::deleteSubtree( KdNode<Node>* node )
  {
    if (node == this->root)
      return;

    KdNode<Node>* parent = node->parent;
    KdNode<Node>* otherChild = (parent->left==node)? parent->right:parent->left;
    if (this->deleteInvalidLeaf(node))
      if (this->deleteInvalidLeaf(otherChild))
        this->deleteSubtree(parent);
  }






}



#endif // kdtree_h__
