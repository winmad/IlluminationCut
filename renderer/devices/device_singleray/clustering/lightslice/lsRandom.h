#ifndef _RANDOM_H_
#define _RANDOM_H_

#include "devices/device_singleray/default.h"
#include <random>

namespace embree{

#define EMBREERANDOM
#ifdef EMBREERANDOM

  template<typename real>
  inline real uniform1D01(Random& rnd) {
    return (real) rnd.getFloat();
  }

  template<typename real>
  inline Vec2<real> uniform2D01(Random& rnd) {
    return Vec2<real>(rnd.getFloat(),rnd.getFloat());
  }

  template<typename integer>
  inline integer uniform1D(Random& rnd, integer m) {
    return (integer) rnd.getInt(m);
  }
#else

#ifdef WIN32
  using std::tr1::minstd_rand;
  using std::tr1::mt19937;
  using std::tr1::uniform_real_distribution;
  using std::tr1::uniform_int_distribution;
#else
  using std::minstd_rand;
  using std::mt19937;
  using std::uniform_real_distribution;
  using std::uniform_int_distribution;
#endif

  typedef minstd_rand random_std_engine;
  typedef mt19937 random_mt_engine;

  template<typename real, class Engine>
  inline real uniform1D01(Engine& engine) {
    static uniform_real_distribution<real> dist(0,1);
    return dist(engine);
  }

  template<typename real, class Engine>
  inline Vec2<real> uniform2D01(Engine& engine) {
    static uniform_real_distribution<real> dist(0,1);
    return Vec2<real>(dist(engine),dist(engine));
  }

  template<typename integer, class Engine>
  inline integer uniform1D(Engine& engine, integer m) {
    static uniform_real_distribution<double> dist(0,1);
    double f = dist(engine);
    return min<integer>((integer)(f * m), (integer)(m - 1));
  }
#endif
}
#endif
