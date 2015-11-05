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

#ifndef icparameters_h__
#define icparameters_h__

// multiple representative lights (if > 1 it has to be bigger than the samples per pixel)
#define REP_LIGHT_NUM 1

// dir light distance: directional lights are simulated with spotlights sufficiently far away
#define DIR_LIGHT_DISTANCE 40000

// fixed clamping distance or scene dependent (this is the lower threshold for distance to avoid singularities)
#if 1
#  define CLAMPING_DISTANCE 200
#  define SET_CLAMPING_DISTANCE(d) void(0)
#else
  // helper to avoid having a cpp file for static member variable
  template< class T >
  struct init_helper {  static float dist; };
  template< class T >
  float init_helper<T>::dist =0.0f;
  struct clamping: public init_helper<void>{};

# define CLAMPING_DISTANCE clamping::dist
# define SET_CLAMPING_DISTANCE(d) clamping::dist=d
#endif


// maximum size of directional bbox: for the 6dim trees the first levels are subdivided according 
// to the directional space, as long as the bounding box of the normals have smaller diameter 
// than this value (e.g., 0.2 means at most 4 such levels since: 2*sqrt(2)/2^4=0.18)
// after the bb of the directions is small enough normal spatial subdivisions continue
#define MAXSIZE_DIRECTIONAL_BBOX 0.6f

// define fixed randomness for debugging
#define RND_SEED (int)time(int(0))
//#define RND_SEED           int(0)

// control lighttree type: lightcuts style agglomerative clustering or divise octree based
#define AGGLOMERATIVE_LIGHTTREE 

// CT descend until fixed cut size (not used for complex scenes)
#define CT_MAXCUT 1000000

// statistics macro
//#define STATISTICS_ENABLED

// output the approximate image instead of the final one
//#define APPROX_IMAGE

// output the approximate image instead of the final one
//#define UB_IMAGE

// Timing for performance influencing timers (e.g., timers that measure small repeated portions of the code)
//#define TIMING_ENABLED

// enable more fine glossy direction WSPD for approximation image
//#define GLOSSY_DIR_CUT

#endif // icparameters_h__


