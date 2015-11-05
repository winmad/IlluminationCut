#ifndef _PIXELSAMPLER_H_
#define _PIXELSAMPLER_H_

#include "devices/device_singleray/default.h"
#include "devices/device_singleray/clustering/lightslice/lsRandom.h"
#include "devices/device_singleray/clustering/lightslice/lsArrays.h"
#include <iostream>

  using std::cerr;
  using std::endl;

namespace embree {

  // HACK: should use std::swap, but VS10 does not compile with it
#define _SWAP_HACK
#ifdef _SWAP_HACK
  template<typename T>
  inline void __swap(T& a, T& b) {
    T c = a;
    a = b;
    b = c;
  }
#else
  using std::swap;
#endif

  class PathSampler {
  public:
    // round the samples to the supported number
    virtual int RoundSamples(int ns) { return ns; }

    // ns is the minimum estimated number of samples
    virtual void BeginPixel(int ns) { }
    virtual void EndPixel() { }

    // next pixel sample (override if special per pixel behaviour)
    virtual void NextPixelSample() { }

    // image sampling (override if better)
    virtual Vec2f Pixel() { return Next2D(); }
    virtual float Time() { return Next1D(); }
    virtual Vec2f Lens() { return Next2D(); }

    // other dimensions (all in 0..1 domain)
    virtual float Next1D() = 0;
    virtual Vec2f Next2D() = 0;

    virtual PathSampler* Fork() = 0;

    // scrable samples
    template<class T>
    void RandomShuffle(T* d, int l) {
      for(int i = 0; i < l; i ++) { 
#ifdef _SWAP_HACK
        int j = _RandomShuffleNext(l);
        __swap(d[i], d[j]); // HACK: should use std::swap, but VS10 does not compile with it
#else
        swap<T>(d[i], d[_RandomShuffleNext(l)]);
#endif
      }
    }
    template<class T>
    void RandomShuffle(carray<T>& s) {
      RandomShuffle(s.data(), static_cast<int>(s.size()));
    }

  private:
    virtual int _RandomShuffleNext(int m) = 0;
  };

  class RandomPathSamplerStd : public PathSampler {
  public:
    RandomPathSamplerStd() {}
    RandomPathSamplerStd(unsigned int m) : seed(m){rnd.setSeed(m);}
    // other dimensions (all in 0..1 domain)
    virtual float Next1D() { return uniform1D01<float>(rnd); }
    virtual Vec2f Next2D() { return Vec2f(uniform1D01<float>(rnd),uniform1D01<float>(rnd)); }
    virtual PathSampler* Fork() { return new RandomPathSamplerStd(seed); }
  protected:
    Random rnd;
    unsigned int seed;

    virtual int _RandomShuffleNext(int m) { return uniform1D(rnd,m); }
  };


  template<bool Jitter>
  class StratifiedPathSampler : public PathSampler {
  public:
    StratifiedPathSampler() {}
    StratifiedPathSampler(unsigned int m) : seed(m) {rnd.setSeed(m);}
    // round the samples to the supported number
    virtual int RoundSamples(int ns) { int nss = (int)sqrt((double)ns+0.5); return nss*nss; }

    // ns is the estimated number of samples
    virtual void BeginPixel(int ns) { _InitSamples(ns); curSample = 0; }
    virtual void EndPixel() { }

    // next pixel sample (override if special per pixel behaviour)
    virtual void NextPixelSample() { curSample ++; }

    // image sampling (override if better)
    virtual Vec2f Pixel() { if(curSample < pixelSamples.size()) return pixelSamples[curSample]; else return Next2D(); }
    virtual float Time() { if(curSample < timeSamples.size()) return timeSamples[curSample]; else return Next1D(); }
    virtual Vec2f Lens() { if(curSample < lensSamples.size()) return lensSamples[curSample]; else return Next2D(); }

    // other dimensions (all in 0..1 domain)
    virtual float Next1D() { return uniform1D01<float>(rnd); }
    virtual Vec2f Next2D() { return Vec2f(uniform1D01<float>(rnd),uniform1D01<float>(rnd)); }

    virtual PathSampler* Fork() { return new StratifiedPathSampler<Jitter>(seed); }
  protected:
    uint32_t curSample;
    carray<Vec2f> pixelSamples;
    carray<float> timeSamples;
    carray<Vec2f> lensSamples;

    Random rnd;
    unsigned int seed;

    void _InitSamples(unsigned int ns) {
      unsigned int nss = (int)sqrt((double)ns+0.5);
      if(nss*nss != ns) { cerr << "samples must be k^2, instead is " << ns << endl; }

      if(pixelSamples.size() != ns) {
        pixelSamples.resize(ns);
        timeSamples.resize(ns);
        lensSamples.resize(ns);
      }

      // stratify samples
      _Stratify2D(pixelSamples, nss, nss);
      _Stratify1D(timeSamples, ns);
      _Stratify2D(lensSamples, nss, nss);

      // shuffle time and lens
      RandomShuffle(timeSamples);
      RandomShuffle(lensSamples);
    }

    void _Stratify1D(carray<float>& s, int l) {
      for(int i = 0; i < l; i ++) s[i] = (((Jitter)?uniform1D01<float>(rnd):0.5f) + i) / l;
    }
    void _Stratify2D(carray<Vec2f>& s, int w, int h) {
      for(int i = 0, c = 0; i < w; i ++) {
        for(int j = 0; j < h; j ++, c++) {
          s[c] = Vec2f( (((Jitter)?uniform1D01<float>(rnd):0.5f) + i) / w, (((Jitter)?uniform1D01<float>(rnd):0.5f) + j) / h);
        }
      }
    }

    virtual int _RandomShuffleNext(int m) { return uniform1D(rnd,m); }
  };

  typedef StratifiedPathSampler<false> StratifiedPathSamplerStd;
  typedef StratifiedPathSampler<true> JitteredPathSamplerStd;

}
#endif
