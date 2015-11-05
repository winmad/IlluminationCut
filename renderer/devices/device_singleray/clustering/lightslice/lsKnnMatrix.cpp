//#include <arma/armadillo>
#include "lsKnnMatrix.h"
#include "devices/device_singleray/clustering/lightslice/lsArrays.h"
#include "devices/device_singleray/clustering/lightslice/lsKdtree.h"
#include "devices/device_singleray/clustering/lightslice/lsLightEval.h"
#include "devices/device_singleray/clustering/lightslice/lsDistribution1D.h"
#include "devices/device_singleray/clustering/lightslice/lsKnnMatrixImpl.h"

#include "devices/device_singleray/clustering/lightslice/gsl/include/gsl/gsl_matrix_double.h"
#include "devices/device_singleray/clustering/lightslice/gsl/include/gsl/gsl_cblas.h"
#include "devices/device_singleray/clustering/lightslice/gsl/include/gsl/gsl_blas.h"

#include "devices/device_singleray/clustering/icStatistics.h"

#include "common/sys/taskscheduler.h"

#include <queue>
#include <map>

namespace embree{

  void KnnMatrix::init(const Ref<BackendScene> &scene,Ref<Camera>& camera) 
  {
    this->_scene = scene; 
    this->_camera = camera;
    float radius = 1000;//(length(_scene->accel->bounds().size()) / 2.0f) * 0.05f;
    //this->_clamp = 4 * radius * radius;
    this->_clamp = CLAMPING_DISTANCE*CLAMPING_DISTANCE;
    _normScale = _scene->sceneradius / 8.0f;
  }

void KnnMatrix::_SetBackground(Ref<SwapChain>& image)
{
  double t = getSeconds();
  for (uint32_t i = 0; i < _bkPixels.size(); i++)
  {
    BackgroundPixel &bkPixel = _bkPixels[i];
    Vec2i &pixel = bkPixel.pixel;
    image->accumulate(pixel.x,pixel.y, bkPixel.background , 1.0f);
  }
  double dt = getSeconds()-t;
  std::cout << "_SetBackground: " << dt*1000.0f << " ms " << std::endl;
}

struct ComputeColNorm
{
  ComputeColNorm(carray2<Color> &matrix, carray2<Color> &norms, const vector<GatherGroup> &groups) 
    : _matrix(matrix), _norms(norms), _groups(groups) {
      this->counter = 0;
      this->maxNum = (uint32_t) groups.size();
  }
  TASK_RUN_FUNCTION(ComputeColNorm,compute);
  TASK_COMPLETE_FUNCTION(ComputeColNorm,donothing);
  void compute()
  {
    TaskScheduler::EventSync event;
    TaskScheduler::Task task(&event,_compute,this,TaskScheduler::getNumThreads(),_donothing,this,"render::getPixels");
    TaskScheduler::addTask(-1,TaskScheduler::GLOBAL_BACK,&task);
    event.sync();
  }
 
  carray2<Color>                      &_norms;
  carray2<Color>				        &_matrix;
  const vector<GatherGroup>			&_groups;
  uint32_t maxNum;
  Atomic counter;
};
void ComputeColNorm::compute(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event)
{
  while(true)
  {
    uint32_t g = (uint32_t) counter++;
    if(g>=this->maxNum)break;
    const GatherGroup &group = _groups[g];
    for (uint32_t i = 0; i < _matrix.width(); i++)
    {
      Color n(zero);
      for (int j= 0; j < group.neighbors.size(); j++)
      {
        uint32_t gidx = group.neighbors[j];
        const Color &v = _matrix.at(i, gidx);
        n += v * v;
      }
      _norms.at(i, g) = sqrt(n);
    }
  }
}
void ComputeColNorm::donothing(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event)
{
}


struct FinalRenderThread
{
  typedef uint32_t argument_type;
  FinalRenderThread(KnnMatrix *knnMat, const vector<uint32_t> &indices, const vector<ScaleLight> &scaleLights, Ref<SwapChain>& image) 
    : _knnMat(knnMat), _indices(indices), _scaleLights(scaleLights), _image(image) {
      this->counter = 0;
      this->maxNumber = (uint32_t) _indices.size();
  };
  void compute()
  {
    while (true)
    {
      uint32_t g = (uint32_t) counter++;
      if(g>=this->maxNumber)break;
      _RenderGatherPoint(_indices[g]); 
    }
  }
  void _RenderGatherPoint( uint32_t g ) {
    GatherPoint &gp = _knnMat->_gatherPoints[g];
    Color L(zero);
    for (uint32_t c = 0; c < _scaleLights.size(); c++)
    {
      const ScaleLight &lc = _scaleLights[c];
      uint32_t idx = lc.indices[gp.index];
      const Color& weight = lc.weights[gp.index];
      //STATS(stats.atomic.nbClusters++;)
      if(zero != weight)
      {
        L += _knnMat->RenderCell(LightEvalUtil::EvalLight(_knnMat->_clamp), idx, gp) * weight;
      }
    }
    {
      STATS(stats.atomic.pixels++);
      Vec2i pos = gp.pixel ;
      //Color cc = gp.emission + L * gp.strength * gp.weight;
      Color cc = gp.emission + L * gp.weight;
      _image->accumulate(pos.x,pos.y,cc,1.0f);
    }
  }
private:
  KnnMatrix                   *_knnMat;
  const vector<ScaleLight>    &_scaleLights;
  const vector<uint32_t>      &_indices;
  Ref<SwapChain>            _image;
  Atomic counter;
  uint32_t maxNumber;
};


struct FinalGroupRenderThread
{
  typedef uint32_t argument_type;
  FinalGroupRenderThread(KnnMatrix *knnMat, Ref<SwapChain>& image)
    : _knnMat(knnMat), _groups(knnMat->_gpGroups), _lights(knnMat->_scaledLights), _image(image) {
  }
  void operator()(uint32_t gg)  { 
    GatherGroup &group = _groups[gg];
    vector<ScaleLight> &scaleLight = _lights[gg];
    for (int i = 0; i < group.indices.size(); i++)
    {
      int g = group.indices[i];
      GatherPoint &gp = _knnMat->_gatherPoints[g];
      Color L(zero);
      for (uint32_t c = 0; c < scaleLight.size(); c++)
      {
        const ScaleLight &lc = scaleLight[c];
        uint32_t idx = lc.indices[gp.index];
        const Color& weight = lc.weights[gp.index];
        if(zero != weight)
        {
          L += _knnMat->RenderCell(LightEvalUtil::EvalLight(_knnMat->_clamp), idx, gp) * weight;
        }
      }
      {
        Vec2i pos = gp.pixel ;
        Color cc = gp.emission + L * gp.strength * gp.weight;
        _image->accumulate(pos.x,pos.y,cc,1.0f);
      }
    }
  }
  KnnMatrix                         *_knnMat;
  std::vector<GatherGroup>          &_groups;
  std::vector<vector<ScaleLight> >  &_lights;
  Ref<SwapChain>                 _image;
};


void KnnMatrix::Render(Ref<SwapChain>& image, uint32_t spp, const std::vector<vplVPL>& indirect, uint32_t seedNum, uint32_t budget)
{
  _GenerateLights(indirect);
  _ShootGatherPoints(image->getWidth(), image->getHeight(), spp,image);
  _GroupGatherPoints(seedNum);
  _FindGatherGroupNeighbors();

  carray2<Color> matrix(_lightList.GetSize(), (uint32_t)(_gpGroups.size()));
  _RenderReducedMatrix(matrix);

   _SetBackground(image);

  vector<vector<uint32_t> > clusters;
  _InitialClusters(clusters, matrix, budget * 0.3f, false);
  _RefineClusters(clusters, matrix, budget, spp, image);
}


struct RefineClusterThread {
  typedef uint32_t argument_type;
  RefineClusterThread(KnnMatrix *knnMat, const vector<vector<uint32_t> > &clusters, vector<vector<ScaleLight> > &scaleLights, 
    carray2<Color> &matrix, carray<Color> &fullNorms, carray2<Color> &colNorms, vector<GatherGroup> &gpGroups, const carray<uint64_t> &randSeeds, uint32_t budget, uint32_t samples, Ref<SwapChain>& image) 
    : _knnMat(knnMat), _clusters(clusters), _scaleLights(scaleLights), _matrix(matrix), _norms(colNorms), _fullNorms(fullNorms), _gpGroups(gpGroups), _budget(budget), _samples(samples), _seeds(randSeeds), _image(image) {
      numGpGroups = _gpGroups.size();
      counterGroups = 0;
  } 
  TASK_RUN_FUNCTION(RefineClusterThread,compute);
  TASK_COMPLETE_FUNCTION(RefineClusterThread,donothing);
  void compute()
  {
    TaskScheduler::EventSync event;
    TaskScheduler::Task task(&event,_compute,this,TaskScheduler::getNumThreads(),_donothing,this,"render::getPixels");
    TaskScheduler::addTask(-1,TaskScheduler::GLOBAL_BACK,&task);
    event.sync();
  }


  float cost(const carray<std::pair<uint32_t, float> > &projection, const BBox<Vec1i> &range, const std::vector<uint32_t> &nsamples, const Color *norms) const {
    if (range.size().x <= 1)
      return 0.0f;

    Color nsum(zero);
    Color vsum(zero);
    for (int32_t i = range.lower.x; i < range.upper.x; i++)
    {
      uint32_t col = projection[i].first;
      const Color &n = norms[col];
      nsum += n;
    }
    for (uint32_t j = 0; j < nsamples.size(); j++)
    {
      Color sum(zero);
      for (int32_t i = range.lower.x; i < range.upper.x; i++)
      {
        const Color &v = _matrix.at(projection[i].first, nsamples[j]); 
        sum += v;
      }
      vsum += sum * sum;
    }
    float cost = reduce_add((nsum * nsum) - vsum) / 3.0f;
    return max(0.0f, cost);
  }
  KnnMatrix                       *_knnMat;
  uint32_t						_budget;
  uint32_t						_samples;
  vector<vector<ScaleLight> >		&_scaleLights;
  const carray<uint64_t>			&_seeds;
  const carray2<Color>			&_matrix;
  const carray2<Color>			&_norms;
  const carray<Color>				&_fullNorms;
  const vector<GatherGroup>		&_gpGroups;
  const vector<vector<uint32_t> > &_clusters;
  Ref<SwapChain> _image;
  size_t numGpGroups; 
  Atomic counterGroups ;
};

void RefineClusterThread::compute(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event)  
{
  while (true)
  {
    uint32_t g = (uint32_t) counterGroups++;
    if(g>=this->numGpGroups)break;

    const GatherGroup &gpGroup = _gpGroups[g];

    const Color *norms = _norms.row(g);
    const vector<uint32_t> &samples = gpGroup.neighbors;
    carray<std::pair<uint32_t, float> > projection(_matrix.width());
    carray<Color> line((uint32_t)samples.size());

    std::vector<std::pair<BBox<Vec1i>, float> > cluster_queue;
    cluster_queue.reserve(_budget);

    uint32_t off = 0;
    for (uint32_t c = 0; c < _clusters.size(); c++) {
      const vector<uint32_t> &cluster = _clusters[c];
      if (cluster.size() == 0)
        continue;

      BBox<Vec1i> range;
      range.lower = (Vec1i) off;
      for (int32_t i = 0; i < cluster.size(); i++)
        projection.at(off++).first = cluster[i];
      range.upper = (Vec1i) off;

      float cst = cost(projection, range, samples, norms);
      cluster_queue.push_back(std::make_pair(range, cst));

      std::push_heap(cluster_queue.begin(), cluster_queue.end(), [](std::pair<BBox<Vec1i>, float> &d1, std::pair<BBox<Vec1i>, float> &d2)->bool { return d1.second < d2.second;});
    }

    while(cluster_queue.size() < _budget)
    {
      float error = cluster_queue.front().second;
      if (error <= 0.0f)
        break;

      BBox<Vec1i> range = cluster_queue.front().first;
      assert(range.size().x >= 2);

      std::pop_heap(cluster_queue.begin(), cluster_queue.end(), [](std::pair<BBox<Vec1i>, float> &d1, std::pair<BBox<Vec1i>, float> &d2)->bool { return d1.second < d2.second;});
      cluster_queue.pop_back();

      float largest = -FLT_MAX, second = -FLT_MAX;
      uint32_t largestIdx = -1, secondIdx = -1;
      for (int32_t i = range.lower; i < range.upper; i++) {
        uint32_t col = projection[i].first;
        float v = reduce_add(norms[col])/3.0f;
        if (v > largest) {
          second = largest;
          secondIdx = largestIdx;
          largest = v;
          largestIdx = col;
        } 
        else if (v > second && v < largest) {
          second = v;
          secondIdx = col;
        }
      }
      assert(largestIdx != -1 && secondIdx != -1);


      for (uint32_t i = 0; i < samples.size(); i++) {
        uint32_t j = samples.at(i);
        line.at(i) = _matrix.at(largestIdx, j) - _matrix.at(secondIdx, j);
      }

      // Why is this not working?
      //BBox<Vec1i> r(EmptyTy());
      BBox<Vec1i> r; r.lower = pos_inf; r.upper = neg_inf;
      uint32_t offset = 0;
      for (int32_t i = range.lower; i < range.upper; i++)
      {
        Color dot(zero); 
        uint32_t col = projection[i].first;
        for (uint32_t j = 0; j < samples.size(); j++)
          dot += _matrix.at(col, samples.at(j)) * line[j];
        float d = reduce_add(dot)/3.0f;
        r.grow( Vec1i( (int) d ) );
        projection[i].second = d;
      }

      float pmid = (float) center(r).x;  // INTEGER DIVISON RESULT TO float?
      std::pair<uint32_t, float> *start = projection.data() + range.lower;
      std::pair<uint32_t, float> *end = projection.data() + range.upper;
      std::pair<uint32_t, float> *middle = 0;
      middle = std::partition(start, end, [pmid](const std::pair<uint32_t, float> &a) { return a.second <= pmid; });
      if (middle == start || middle == end )
        middle = start + (end - start) / 2;

      assert(middle != start && middle != end);

      int rangMid = (int)(middle - start) + range.lower;
      BBox<Vec1i> range1(range.lower, (Vec1i) rangMid);
      BBox<Vec1i> range2( (Vec1i) rangMid, range.upper);
      float cost1 = cost(projection, range1, samples, norms);
      float cost2 = cost(projection, range2, samples, norms);

      cluster_queue.push_back(std::make_pair(range1, cost1));
      std::push_heap(cluster_queue.begin(), cluster_queue.end(), [](std::pair<BBox<Vec1i>, float> &d1, std::pair<BBox<Vec1i>, float> &d2)->bool {
        return d1.second < d2.second;
      });
      cluster_queue.push_back(std::make_pair(range2, cost2));
      std::push_heap(cluster_queue.begin(), cluster_queue.end(), [](std::pair<BBox<Vec1i>, float> &d1, std::pair<BBox<Vec1i>, float> &d2)->bool {
        return d1.second < d2.second;
      });
    }

    RandomPathSamplerStd sampler;
    std::vector<ScaleLight> scaleLight;
    for(uint32_t k = 0; k < cluster_queue.size(); k++)
    {
      BBox<Vec1i> &range = cluster_queue[k].first;
      if ( range.empty() )
        continue;

      std::vector<Color> cnorms;
      Color cnormSum(zero);
      float maxL = -1;
      uint32_t index = -1;
      Vector3f maxNorm(zero);
      for (int32_t i = range.lower; i < range.upper; i++)
      {
        const Color &v = _fullNorms[projection[i].first];
        cnorms.push_back(v);
        cnormSum += v;
      }

      if (zero != cnormSum)
      {
        lsDistribution1D<float, Color> dist(&cnorms[0], (uint32_t)cnorms.size());
        float pdf;
        ScaleLight light;
        for (uint32_t s = 0; s < _samples; s++)
        {
          uint32_t idx = dist.SampleDiscrete(sampler.Next1D(), &pdf);
          light.indices.push_back(projection[range.lower + idx].first);
          //Vector3f weight = cnormSum | cnorms[idx];
          Color weight;
          weight.r = cnorms[idx].r!=0.0f ? cnormSum.r / cnorms[idx].r : 0 ;
          weight.g = cnorms[idx].g!=0.0f ? cnormSum.g / cnorms[idx].g : 0 ;
          weight.b = cnorms[idx].b!=0.0f ? cnormSum.b / cnorms[idx].b : 0 ;
          light.weights.push_back(weight);
        }
        scaleLight.push_back(light);
      }
    }
    FinalRenderThread thread(_knnMat, gpGroup.indices, scaleLight, _image);
    //scheduler->start();
    //scheduler->addTask(TaskScheduler::ThreadInfo(),TaskScheduler::GLOBAL_FRONT,
    //  (TaskScheduler::runFunction)&FinalRenderThread::run_thread,&thread,scheduler->getNumThreads(),NULL,NULL,"render::tile");
    //scheduler->stop();
    // Threading in threading is not applicable for embree so this will be singlethreaded
    thread.compute();
  }
}

void RefineClusterThread::donothing(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event)
{
}

void KnnMatrix::_RefineClusters(vector<vector<uint32_t> > &clusters, carray2<Color> &matrix, uint32_t budget, uint32_t samples, Ref<SwapChain>& image)
{
  double t1 = getSeconds();
  double t = getSeconds();
  
  carray<Color> fullNorms(matrix.width());
  ComputeFullColNorm computeFullNormThread(matrix, fullNorms);
  computeFullNormThread.compute();
 
  carray2<Color> colNorms(matrix.width(), _gpGroups.size());
  ComputeColNorm computeColNormThread(matrix, colNorms, _gpGroups);
  computeColNormThread.compute();

  double dt = getSeconds()-t;
  std::cout << " -- pre-compute norms: " << dt*1000.0f << " ms " << std::endl;

  _scaledLights.resize(_gpGroups.size());
  carray<uint64_t> randSeeds((uint32_t)_gpGroups.size());

  t = getSeconds();
  RefineClusterThread thread(this, clusters, _scaledLights, matrix, fullNorms, colNorms, _gpGroups, randSeeds, budget, samples, image);
  thread.compute();
  dt = getSeconds()-t;
  std::cout << " -- refine clusters: " << dt*1000.0f << " ms " << std::endl;
  dt = getSeconds()-t1;
  std::cout << "_RefineClusters: " << dt*1000.0f << " ms " << std::endl;
}


void KnnMatrix::_ShootGatherPoints(uint32_t width, uint32_t height, uint32_t samples,Ref<SwapChain>& image)
{
  double t = getSeconds();
  GatherPointShooter::Shoot(_scene,_camera, width, height, samples, _gatherPoints, _bkPixels, image);
  double dt = getSeconds()-t;
  std::cout << "_ShootGatherPoints: " << dt*1000.0f << " ms " << std::endl;
}

void KnnMatrix::_GenerateLights( const std::vector<vplVPL>& indirect )
{
  // convert our lights into their format
  //ListVirtualLightCache cache(_lightList);
  //_generator->Generate(indirect, &cache);
  double t  = getSeconds();
  for (size_t i = 0 ; i < indirect.size() ; ++i) 
  {
    // VLight is always a spotlight
    const vplVPL& vl = indirect[i];
    _lightList._otrLights.push_back(lsOrientedLight(vl.P,-vl._D,vl.I));
  }
  double dt = getSeconds()-t;
  std::cout << "_GenerateLights: " << dt*1000.0f << " ms " << std::endl;
}
/*
void KnnMatrix::RenderGatherGroup( Image<Vector3f> *gpImage )
{
  RandomPathSamplerStd sampler;
  if (_report) _report->beginActivity("rendering gather point cluster");
  for (uint32_t i = 0; i < _gpGroups.size(); i++)
  {
    const GatherGroup &gpGroup = _gpGroups[i];
    Vector3f color(sampler.Next1D(),sampler.Next1D(),sampler.Next1D());

    for (uint32_t j = 0; j < gpGroup.indices.size(); j++)
    {
      GatherPoint& gp = _gatherPoints[gpGroup.indices[j]];
      gpImage->ElementAt(gp.pixel.x, gpImage->Height()- gp.pixel.y - 1) = color;
    }
    if (_report) _report->progress(i / (float)_gpGroups.size(), 1);
  }
  if (_report) _report->endActivity();
}

*/
void KnnMatrix::_KdGatherGroup( vector<GatherKdItem>::iterator start, vector<GatherKdItem>::iterator end)
{
  BBox<Vec6f> bbox = BBox<Vec6f>(EmptyTy());
  for (vector<GatherKdItem>::iterator it = start; it != end; it++)
    bbox.grow(it->p);
  if(end - start < _maxGatherGroupSize)
  {
    _gpGroups.push_back(GatherGroup());
    GatherGroup &group = _gpGroups.back();
    for (vector<GatherKdItem>::iterator it = start; it != end; it++)
    {
      group.indices.push_back(it->idx);
      GatherPoint &gp = _gatherPoints[it->idx];
      group.bbox.grow(gp.isect.P);
      group.normal += gp.isect.Ns;
    }
    group.normal = normalize(group.normal);
    assert(zero != group.normal);
    uint32_t sidx = min((uint32_t)(_sampler.Next1D() * group.indices.size()), (uint32_t)(group.indices.size() - 1));
    group.seed = group.indices[sidx];
  }
  else
  {
    uint32_t dim = maxDim(bbox.size());

    float pmid = center(bbox)[dim];
    vector<GatherKdItem>::iterator mid = 
      std::partition(start, end, [dim, pmid](const GatherKdItem &a) { return a.p[dim] < pmid; });

    _KdGatherGroup(start, mid);
    _KdGatherGroup(mid, end);
  }
}

void KnnMatrix::_GroupGatherPoints( uint32_t seedNum )
{
  double t = getSeconds();
  _maxGatherGroupSize = max<uint32_t>(1, (uint32_t)(_gatherPoints.size() / seedNum));
  // Compute Gather Point Bounding Box
  BBox3f gpBBox = BBox3f(EmptyTy());
  for (uint32_t i = 0; i < _gatherPoints.size(); i++)
  {
    GatherPoint &gp = _gatherPoints[i];
    Vector3f &P = gp.isect.P;
    gpBBox.grow(P);
  }
  _normScale = reduce_add(gpBBox.size()) / 3.0f / 8.0f;
  _diagonal = length(gpBBox.size()) / 32.0f;

  vector<GatherKdItem> items(_gatherPoints.size());
  for (uint32_t i = 0; i < _gatherPoints.size(); i++)
  {
    GatherPoint &gp = _gatherPoints[i];
    GatherKdItem &item = items[i];
    item.idx = i;
    item.p = Vec6f(gp.isect.P, gp.isect.Ns * _normScale);
  }

  _KdGatherGroup(items.begin(), items.end());

  std::cout << "Gather group numbers: " << _gpGroups.size()<<std::endl;
  double dt = getSeconds()-t;
  std::cout << "_GroupGatherPoints: " << dt*1000.0f << " ms " << std::endl;
}

void KnnMatrix::_FindGatherGroupNeighbors()
{
  double t = getSeconds();
  vector<GatherKdItem> data;
  for (uint32_t i = 0; i < _gpGroups.size(); i++)
  {
    GatherGroup &gpGroup = _gpGroups[i];
    Vector3f P = center(gpGroup.bbox);
    Vector3f &N = gpGroup.normal;
    data.push_back(GatherKdItem(i, Vec6f(P, N * _normScale)));
  }

  lsKdTree<GatherKdItem> *kdTree = new lsKdTree<GatherKdItem>(data);
  for (uint32_t g = 0; g < _gpGroups.size(); g++)
  {
    GatherGroup &gpGroup = _gpGroups[g];
    Vec6f P(center(gpGroup.bbox), gpGroup.normal * _normScale);

    // find mini matrix
    SeedProcess proc(8);
    float r2 = FLT_MAX;
    kdTree->Lookup(P, proc, r2);

    for (uint32_t i = 0; i < proc.foundSeeds; i++)
    {
      CloseSeed &cSeed = proc.closeSeeds[i];
      uint32_t index = cSeed.seedKdItem->idx;
      gpGroup.neighbors.push_back(index);
    }
  }
  delete kdTree;
  double dt = getSeconds()-t;
  std::cout << "_FindGatherGroupNeighbors: " << dt*1000.0f << " ms " << std::endl;
}


class RenderReducedMatrixThread
{
public:
  typedef uint32_t argument_type;
  RenderReducedMatrixThread(carray2<Color> &matrix, KnnMatrix *knnMat) 
    : _matrix(matrix), _knnMat(knnMat),numGroup(knnMat->_gpGroups.size()) 
  {
    this->groupCounter = 0;
  };
  TASK_RUN_FUNCTION(RenderReducedMatrixThread,compute);
  TASK_COMPLETE_FUNCTION(RenderReducedMatrixThread,donothing);
  void compute()
  {
    TaskScheduler::EventSync event;
    TaskScheduler::Task task(&event,_compute,this,TaskScheduler::getNumThreads(),_donothing,this,"render::getPixels");
    TaskScheduler::addTask(-1,TaskScheduler::GLOBAL_BACK,&task);
    event.sync();
  }
private:
  carray2<Color>				&_matrix;
  KnnMatrix	                *_knnMat;
  Atomic                 groupCounter;
  const unsigned int     numGroup;
};

void RenderReducedMatrixThread::compute(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event)
{
  while (true)
  {
    unsigned int g = this->groupCounter++;
    if(g >= this->numGroup) break;
    const GatherGroup &gpGroup = _knnMat->_gpGroups[g];
    uint32_t gpIdx = gpGroup.seed;
    const GatherPoint &gp = _knnMat->_gatherPoints[gpIdx];
    for (uint32_t i = 0; i < _matrix.width(); i++)
      _matrix.at(i, g) = _knnMat->RenderCell(LightEvalUtil::EvalL(_knnMat->_clamp), i, gp);
  }
}
void RenderReducedMatrixThread::donothing(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event)
{
}

void KnnMatrix::_RenderReducedMatrix(carray2<Color> &matrix)
{
  double t = getSeconds();

  RenderReducedMatrixThread thread(matrix, this);
  thread.compute();

  double dt = getSeconds()-t;
  std::cout << "_RenderReducedMatrix: " << dt*1000.0f << " ms " << std::endl;
}


/*
void KnnMatrix::_NewInitialClusters(vector<vector<uint32_t> > &clusters, carray2<Vector3f> &matrix, uint32_t budget)
{
  RandomPathSamplerStd sampler;
  if (_report) _report->beginActivity("initial clustering");

  carray2<float> normalizedMatrix(matrix.width(), matrix.height());
  carray<float> norms(matrix.width());

  for (uint32_t i = 0; i < matrix.width(); i++)
  {
    float sum = 0;
    for (uint32_t j = 0; j < matrix.height(); j++)
    {
      float v = matrix.at(i, j).GetLength();
      sum += v * v;
    }
    float n = sqrt(sum);
    norms.at(i) = max(n, 0.0f);
    for (uint32_t j = 0; j < matrix.height(); j++)
    {
      float v = n <= 0.0f ? 0.0f : matrix.at(i, j).GetLength() / n;
      assert(v >= 0.0f);
      normalizedMatrix.at(i, j) = v;
    }
  }

  set<uint32_t> center_set;
  Distribution1Df dist(norms.data(), norms.size());
  if (!dist.IsValid()) {
    while (center_set.size() < budget)	{
      uint32_t idx = static_cast<uint32_t>(sampler.Next1D() * norms.size());
      center_set.insert(idx);
    }
  } else {
    uint32_t s = budget;
    while (s--) {
      float pdf; 
      uint32_t idx = static_cast<uint32_t>(dist.SampleDiscrete(sampler.Next1D(), &pdf));
      center_set.insert(idx);
    }
  }

  vector<uint32_t> centers(center_set.begin(), center_set.end());
  clusters.resize(centers.size());
  for (uint32_t i = 0; i < normalizedMatrix.width(); i++)
  {
    int mini = -1;
    double mind = FLT_MAX;
    for (uint32_t c = 0; c < centers.size(); c++)
    {
      uint32_t cindex = centers[c];
      float d = 0;
      for (uint32_t j = 0; j < normalizedMatrix.height(); j++)
      {
        float v = normalizedMatrix.at(i, j) - normalizedMatrix.at(cindex, j);
        d += v * v;
        assert(d >= 0.0f);
      }

      if (mind >= d)
      {
        mind = d;
        mini = c;
        assert(mini < clusters.size());
      }
    }
    clusters[mini].push_back(i);
  }

  if (_report) _report->endActivity();
}

*/

float col3fLength(const Color& c){ return sqrt( c.r*c.r + c.g*c.g + c.b*c.b); }

void KnnMatrix::_InitialClusters(vector<vector<uint32_t> > &clusters, carray2<Color> &matrix, uint32_t budget, bool randProj )
{
  double t = getSeconds();

  RandomPathSamplerStd sampler;

  gsl_matrix* origin = gsl_matrix_alloc(matrix.height(), matrix.width());
  for(uint32_t j = 0; j < matrix.height(); j++)
    for (uint32_t i = 0; i < matrix.width(); i++)
      gsl_matrix_set(origin, j, i, col3fLength(matrix.at(i, j)) );

  gsl_matrix* input;
  // random projection?
  if(randProj)
  {
    gsl_matrix* rand_mat = gsl_matrix_alloc(50, matrix.height());
    for(uint32_t j = 0; j < 50; j++)
      for (uint32_t i = 0; i < matrix.height(); i++)
        gsl_matrix_set(rand_mat, j, i, sampler.Next1D());
    input = gsl_matrix_alloc(50, matrix.width());
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, rand_mat, origin, 0.0, input);
    gsl_matrix_free(origin);
  }
  else
    input = origin;


  gsl_vector* norms = gsl_vector_alloc(input->size2);

  vector<uint32_t> nzIndices;
  vector<uint32_t> zIndices;

  vector<double> nzNorms;
  nzNorms.reserve(norms->size);
  for (uint32_t i = 0; i < input->size2; i++)
  {
    gsl_vector_view column = gsl_matrix_column(input, i);
    double d = gsl_blas_dnrm2(&column.vector);
    gsl_vector_set(norms, i, d);
    if (d > 0.0)
    {
      nzNorms.push_back(d);
      nzIndices.push_back(i);
    }
    else
      zIndices.push_back(i);
  }

  gsl_matrix* nzinput = gsl_matrix_alloc(input->size1, nzNorms.size());
  gsl_vector* nznorms = gsl_vector_alloc(nzNorms.size());
  for (uint32_t i = 0; i < nzNorms.size(); i++)
  {
    gsl_vector_set(nznorms, i, nzNorms[i]);
    gsl_vector_view column = gsl_matrix_column(input, nzIndices[i]);
    gsl_matrix_set_col(nzinput, i, &column.vector);
  }

  gsl_vector* sinput = gsl_vector_alloc(nzinput->size1);
  for (uint32_t j = 0; j < nzinput->size1; j++)
  {
    gsl_vector_view row = gsl_matrix_row(nzinput, j);
    double s = gsl_blas_dasum(&row.vector);
    gsl_vector_set(sinput, j, s);
  }

  double snorm = gsl_blas_dasum(nznorms);
  gsl_vector* alphas = gsl_vector_alloc(nzinput->size2);
  {
    gsl_vector *y = gsl_vector_alloc(nzinput->size2);
    gsl_blas_dgemv(CblasTrans, 1.0, nzinput, sinput, 0.0, y);

    gsl_vector_memcpy(alphas, nznorms);
    gsl_vector_scale(alphas, snorm);

    gsl_vector_sub(alphas, y);
    gsl_vector_free(y);
  }

  std::map<uint32_t, double> centers;
  lsDistribution1Dd alphaDist(alphas->data, (uint32_t)alphas->size);
  sampler.BeginPixel(budget);
  if (!alphaDist.IsValid())
  {
    while (centers.size() < budget)
    {
      uint32_t idx = static_cast<uint32_t>(sampler.Time() * alphas->size);
      centers[nzIndices[idx]] += alphas->size;
      sampler.NextPixelSample();
    }
  }
  else
  {
    uint32_t s = budget;
    while (s--)
    {
      double pdf; 
      uint32_t idx = static_cast<uint32_t>(alphaDist.SampleDiscrete(sampler.Time(), &pdf));
      centers[nzIndices[idx]] += pdf == 0.0 ? 0.0 : 1.0 / pdf;
      sampler.NextPixelSample();
    }
  }
  sampler.EndPixel();
  assert(centers.size() > 0);

  gsl_vector* knorms = gsl_vector_alloc((uint32_t)centers.size());
  gsl_matrix* kcolumns = gsl_matrix_alloc(input->size1, (uint32_t)centers.size());
  std::map<uint32_t, double>::iterator it = centers.begin();
  for (uint32_t idx = 0; it != centers.end(); idx++, it++)
  {
    uint32_t i = (uint32_t)it->first;
    assert(i >= 0 && i < norms->size);
    double n = gsl_vector_get(norms, i);
    gsl_vector_set(knorms, idx, n * it->second);

    gsl_vector_view column = gsl_matrix_column(input, i);
    gsl_matrix_set_col(kcolumns, idx, &column.vector);
    gsl_vector_view kcolumn = gsl_matrix_column(kcolumns, idx);
    gsl_vector_scale(&kcolumn.vector, it->second);
  }

  gsl_matrix *mdist = gsl_matrix_alloc(knorms->size, nznorms->size);
  {
    gsl_matrix_view x = gsl_matrix_view_vector(knorms, knorms->size, 1);
    gsl_matrix_view y = gsl_matrix_view_vector(nznorms, nznorms->size, 1);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &x.matrix, &y.matrix, 0.0, mdist);

    gsl_matrix *col2 = gsl_matrix_alloc(kcolumns->size2, nzinput->size2);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, kcolumns, nzinput, 0.0, col2);

    gsl_matrix_sub(mdist, col2);
    gsl_matrix_free(col2);
  }

  clusters.resize(kcolumns->size2);
  for (uint32_t i = 0; i < nzinput->size2; i++)
  {
    int mini = -1;
    double mind = FLT_MAX;
    for (uint32_t j = 0; j < kcolumns->size2; j++)
    {
      double dist = max(0.0, gsl_matrix_get(mdist, j, i));
      if (mind >= dist)
      {
        mind = dist;
        mini = j;
      }
    }
    clusters[mini].push_back(nzIndices[i]);
  }
  gsl_matrix_free(mdist);
  gsl_matrix_free(nzinput);
  gsl_vector_free(nznorms);
  gsl_vector_free(alphas);
  gsl_vector_free(norms);
  gsl_matrix_free(input);

  // fill up lights and ranges

  if (zIndices.size())
    clusters.push_back(zIndices);

  double dt = getSeconds()-t;
  std::cout << "_InitialClusters: " << dt*1000.0f << " ms " << std::endl;
}

/*
void KnnMatrix::_KmeanGatherGroup( vector<GatherKdItem> items, uint32_t seedNum )
{
  vector<Vec6f> seeds;
  StratifiedPathSamplerStd sampler;
  seedNum = sampler.RoundSamples(seedNum);
  sampler.BeginPixel(seedNum);
  for (uint32_t i = 0; i < seedNum; i++)
  {
    Ray ray = _scene->MainCamera()->GenerateRay(sampler.Pixel(), sampler.Lens(), 0.0f);
    Intersection isect;
    if(_engine->Intersect(ray, &isect))
    {
      BxdfUnion msu;
      isect.m->SampleReflectance(isect.dp, msu);
      if (msu.HasSmooth())
        seeds.push_back(Vec6f(isect.dp.P, isect.dp.N * _normScale));
    }
    sampler.NextPixelSample();
  }
  sampler.EndPixel();

  vector<concurrent_vector<uint32_t> > clusters(seeds.size());

  GatherPointAssign assign(seeds, items, clusters);
  GatherPointUpdate update(seeds, items, clusters);
  for (uint32_t i = 0; i < 4; i++)
  {
    for (uint32_t i = 0; i < clusters.size(); i++)
      clusters[i].clear();
    parallel_for(blocked_range<uint32_t>(0, (uint32_t)items.size()), assign);
    parallel_for(blocked_range<uint32_t>(0, (uint32_t)seeds.size()), update);
  }

  for (uint32_t i = 0; i < clusters.size(); i++)
  {
    if (clusters[i].size() == 0)
      continue;
    _gpGroups.push_back(GatherGroup());
    GatherGroup &group = _gpGroups.back();
    for (concurrent_vector<uint32_t>::iterator it = clusters[i].begin(); it != clusters[i].end(); it++)
    {
      group.indices.push_back(items[*it].idx);
      GatherPoint &gp = _gatherPoints[items[*it].idx];
      group.bbox.Grow(gp.isect.dp.P);
      group.normal += (gp.isect.dp.N);
    }
    group.normal.Normalize();
  }
}

*/
}
