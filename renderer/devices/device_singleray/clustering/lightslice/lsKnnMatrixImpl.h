#ifndef _CO_CLUSTER_MATRIX_UTILS_H_
#define _CO_CLUSTER_MATRIX_UTILS_H_

#include "devices/device_singleray/default.h"
#include "common/sys/taskscheduler.h"
#include "devices/device_singleray/clustering/lightslice/lsArrays.h"

namespace embree{

  class KnnMatrix;

  struct ComputeFullColNorm
  {
    ComputeFullColNorm(carray2<Color> &matrix, carray<Color> &norms) 
      : _matrix(matrix), _norms(norms)
    {
      this->counter = 0 ;
      this->maxNum = _norms.size();
    }
    /*! start functon */
    TASK_RUN_FUNCTION(ComputeFullColNorm,compute);
    TASK_COMPLETE_FUNCTION(ComputeFullColNorm,donothing);
    void compute()
    {
      TaskScheduler::EventSync event;
      TaskScheduler::Task task(&event,_compute,this,TaskScheduler::getNumThreads(),_donothing,this,"render::getPixels");
      TaskScheduler::addTask(-1,TaskScheduler::GLOBAL_BACK,&task);
      event.sync();
    }
   
    carray2<Color>		        &_matrix;
    carray<Color>                       &_norms;
    Atomic counter;
    uint32_t maxNum;
  };

  /*
  struct ComputeRowSum
  {
    ComputeRowSum(Image<Vec3f> &matrix, const vector<uint32_t> &cluster, const vector<uint32_t> &nsamples, vector<Vec3f> &norms) 
      : _matrix(matrix), _rsums(norms), _nsamples(nsamples), _cluster(cluster) { }
    void operator() (const blocked_range<uint32_t> &r) const;
    vector<Vec3f>                       &_rsums;
    Image<Vec3f>				        &_matrix;
    const vector<uint32_t>              &_nsamples;
    const vector<uint32_t>              &_cluster;
  };
  */
}
#endif // _CO_CLUSTER_MATRIX_UTILS_H_
