## ======================================================================== ##
## Copyright 2009-2013 Intel Corporation                                    ##
##                                                                          ##
## Licensed under the Apache License, Version 2.0 (the "License");          ##
## you may not use this file except in compliance with the License.         ##
## You may obtain a copy of the License at                                  ##
##                                                                          ##
##     http://www.apache.org/licenses/LICENSE-2.0                           ##
##                                                                          ##
## Unless required by applicable law or agreed to in writing, software      ##
## distributed under the License is distributed on an "AS IS" BASIS,        ##
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. ##
## See the License for the specific language governing permissions and      ##
## limitations under the License.                                           ##
## ======================================================================== ##

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/clustering/lightslice/gsl/include )

SET (SOURCES
    api/singleray_device.cpp
    api/scene.cpp
    lights/hdrilight.cpp   
    shapes/trianglemesh_normals.cpp   
    shapes/trianglemesh_full.cpp       
    samplers/sampler.cpp
    samplers/distribution1d.cpp
    samplers/distribution2d.cpp
    integrators/pathtraceintegrator.cpp
    filters/filter.cpp
    renderers/debugrenderer.cpp
    renderers/integratorrenderer.cpp
    renderers/progress.cpp
    clustering/icStatistics.cpp
    clustering/lightcuts/ctRenderer.cpp
    clustering/lightcuts/ctLightcuts.cpp
    clustering/lightcuts/ctKdTree.cpp
    clustering/lightcuts/ctIntegrator.cpp
    clustering/vpl/vplRenderer.cpp
    clustering/vpl/vplIntegrator.cpp
    clustering/illuminationcut/icIlluminationCut.cpp
    clustering/illuminationcut/icRenderer.cpp
    clustering/illuminationcut/icPixeltree.cpp
    clustering/illuminationcut/icLighttree.cpp
    clustering/illuminationcut/icIntegrator.cpp
    clustering/lightslice/lsGatherPointShooter.cpp
    clustering/lightslice/lsIntegrator.cpp
    clustering/lightslice/lsKnnMatrix.cpp
    clustering/lightslice/lsKnnMatrixImpl.cpp
    clustering/lightslice/lsLightData.cpp
    clustering/lightslice/lsLightEval.cpp
    clustering/lightslice/lsRenderer.cpp
    clustering/lsLightcuts/lcDivisiveLightTreeBuilder.cpp
    clustering/lsLightcuts/lcIntegrator.cpp
    clustering/lsLightcuts/lcLightcutter.cpp
    clustering/lsLightcuts/lcLightTree.cpp
    clustering/lsLightcuts/lcMTLightcutter.cpp
    clustering/lsLightcuts/lcRenderer.cpp
    clustering/mdLightcuts/mdGatherTreeBuilder.cpp
    clustering/mdLightcuts/mdIntegrator.cpp
    clustering/mdLightcuts/mdLightcutter.cpp
    clustering/mdLightcuts/mdMTLightcutter.cpp
    clustering/mdLightcuts/mdRenderer.cpp

    )
SET (CMAKE_CXX_FLAGS " -fPIC -std=c++0x ")

#static libraries for the server
#SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#SET(GSL_LIBRARY "/usr/lib64/libgsl.a")
#SET(GSL_LIBRARY "/usr/lib64/libgslcblas.a")
SET(GSL_LIBRARY "/usr/lib/libgsl.so")
SET(CBLAS_LIBRARY "/usr/lib/libgslcblas.so")

IF (__XEON__)
  ADD_LIBRARY(device_singleray SHARED ${SOURCES})
  SET (CMAKE_SHARED_LINKER_FLAGS " -lgsl -lgslcblas -lm")
  TARGET_LINK_LIBRARIES(device_singleray sys ${EMBREE_LIBRARY} image ${GSL_LIBRARY} ${CBLAS_LIBRARY})
ELSE ()
  INCLUDE(icc_mic)
  ADD_EXECUTABLE(device_singleray_knc ${SOURCES})
  TARGET_LINK_LIBRARIES(device_singleray_knc sys_knc ${EMBREE_LIBRARY_MIC} coi_server_knc)
ENDIF ()
