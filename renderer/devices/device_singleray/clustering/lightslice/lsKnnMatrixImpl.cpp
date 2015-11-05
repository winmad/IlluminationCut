//#include <arma/armadillo>
#include "lsKnnMatrixImpl.h"
#include "lsKnnMatrix.h"

namespace embree{
  void ComputeFullColNorm::compute(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event)
  {
    while (true)
    {
      uint32_t i = (uint32_t) counter++;
      if(i>=this->maxNum) break;
      Color n(zero);
      for (uint32_t j = 0; j < _matrix.height(); j++)
      {
        const Color &v = _matrix.at(i, j);
        n += v * v;
      }
      _norms.at(i) = sqrt(n);
    }
  }

  void ComputeFullColNorm::donothing(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event)
  {
  }
}
