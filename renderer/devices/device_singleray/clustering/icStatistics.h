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

#ifndef statistics_h__
#define statistics_h__

#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "devices/device_singleray/default.h"
#include "devices/device_singleray/clustering/icParameters.h"


#ifdef STATISTICS_ENABLED
#  define STATS(x) x
#  define _EASY_STATS
#  define _COMPLEX_STATS
#else
#  define STATS(x) (void) 0
#endif

#ifdef TIMING_ENABLED
#  define TIMESTART(x,name) volatile double name = getSeconds()
      
#  define TIMESTOP(x,name)  x += 1000*(getSeconds() - name)
#else
#  define TIMESTART(x,name) (void) 0
#  define TIMESTOP(x,name) (void) 0
#endif

// Two levels of stats (basically locking ones and not locking ones)
#ifdef _EASY_STATS
#  define EASY_STATS(x) x
#else
#  define EASY_STATS(x) 
#endif
#ifdef _COMPLEX_STATS
#  define COMPLEX_STATS(x) x
#else
#  define COMPLEX_STATS(x) 
#endif


namespace embree 
{
  /*! \brief Class for writing something to a file.*/
  class SimpleFileWriter
  {
  public:

    SimpleFileWriter(){};

    SimpleFileWriter(SimpleFileWriter const&);  

    void operator=(SimpleFileWriter const&); 

    ~SimpleFileWriter(){EASY_STATS(
      if (!strcmp(fileName.c_str(),""))
      {
        std::cout << "Not used SimpleFilewriter" << std::endl;
      }else
      {
        std::ofstream myfile (fileName.c_str());
        if (myfile.is_open())
        {
          myfile << buffer.str()<<std::endl;
          myfile.close();
        }
      }
      )}

    void setFilename(std::string fname){EASY_STATS(
      this->fileName = fname; 
    return;
    )}

    template <typename T>
    SimpleFileWriter& operator<<(T simple){EASY_STATS(
      buffer<< simple; 
    )
      return *this;
    }

    SimpleFileWriter& newline(){EASY_STATS(
      buffer<< std::endl; 
    )
      return *this;
    }

  private:
    std::string fileName;
    std::stringstream buffer;
  };


  /*! \brief Thread-safe class for stats (atomic).*/
  class AtomicCounters{
  public:

    AtomicCounters():
        pixels(0),gatherpoints(0),clusters(0),ub(0),rays(0),
        icpairs(0),geompairs(0),gpsize(0),singletonpair(0),
        timeRT(0),timeUB(0),timeHeap(0),timeTotal(0){}

    ~AtomicCounters(){
      std::cout<< "Avg clusters/gatherpoints: " << clusters / float(gatherpoints) << std::endl;
      std::cout<< "Avg ub/gatherpoints: " << ub / float(gatherpoints) << std::endl;
      std::cout<< "Avg rays/gatherpoints: " << rays / float(gatherpoints) << std::endl;
      std::cout<< "Avg ic/gatherpoints: " << icpairs / float(gatherpoints) << std::endl;
      std::cout<< "Avg geom/gatherpoints: " << geompairs / float(gatherpoints) << std::endl;
      std::cout<< "Avg gpsize/wspairs: " << gpsize / float(icpairs) << std::endl;
      std::cout<< "Avg singletonpair/gatherpoints: " << singletonpair / float(gatherpoints) << std::endl;
      std::cout<< "RT time: "<< timeRT * 0.001f <<" s"<<std::endl;
      std::cout<< "UB time: "<< timeUB * 0.001f <<" s"<<std::endl;
      std::cout<< "Heap time: "<< timeHeap * 0.001f <<" s"<<std::endl;
      std::cout<< "Total time: "<< timeTotal * 0.001f <<" s"<<std::endl;
    }

  public:
    Atomic pixels;
    Atomic gatherpoints;
    Atomic clusters;
    Atomic ub;
    Atomic rays;
    Atomic icpairs;
    Atomic geompairs;
    Atomic gpsize;
    Atomic singletonpair;
    Atomic timeRT;
    Atomic timeUB;
    Atomic timeHeap;
    Atomic timeTotal;
  };

  /*! \brief Thread-safe class for stats (mutex).*/
  class SafeCounters{
  public:

    SafeCounters():somenumber(0){}

    ~SafeCounters(){}

    __forceinline void adjust(int g)
    {
      Lock<MutexActive> lock(mutex);
      somenumber += g;
    }

  public:
    int somenumber;
    MutexActive mutex;
  };


  /*! \brief Global object for stats.*/
  class Stats{
  public:
    SafeCounters safe;
    AtomicCounters atomic;
  };
  extern Stats stats;


  class ColormapJet
  {
  public:
    /*! \brief Map values to a colormap.*/
    static void mapToJet(float value, float& r, float& g, float& b);
  };
}


#endif // statistics_h__
