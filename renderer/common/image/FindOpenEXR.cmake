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

SET(LIBRARY_PATHS
  /usr/lib /usr/local/lib /opt/local/lib)

FIND_PATH(OPENEXR_INCLUDE_PATH ImfHeader.h
  /usr/include       /usr/include/OpenEXR 
  /usr/local/include /usr/local/include/OpenEXR 
  /opt/local/include /opt/local/include/OpenEXR)

FIND_LIBRARY(OPENEXR_HALF_LIBRARY      NAMES Half      PATHS ${LIBRARY_PATHS})
FIND_LIBRARY(OPENEXR_IEX_LIBRARY       NAMES Iex       PATHS ${LIBRARY_PATHS})
FIND_LIBRARY(OPENEXR_IMATH_LIBRARY     NAMES Imath     PATHS ${LIBRARY_PATHS})
FIND_LIBRARY(OPENEXR_ILMIMF_LIBRARY    NAMES IlmImf    PATHS ${LIBRARY_PATHS})
FIND_LIBRARY(OPENEXR_ILMTHREAD_LIBRARY NAMES IlmThread PATHS ${LIBRARY_PATHS})

IF (OPENEXR_INCLUDE_PATH AND OPENEXR_HALF_LIBRARY AND OPENEXR_IEX_LIBRARY AND OPENEXR_IMATH_LIBRARY AND OPENEXR_ILMIMF_LIBRARY AND OPENEXR_ILMTHREAD_LIBRARY)
   SET(OPENEXR_FOUND TRUE)
   SET(OPENEXR_INCLUDE_PATHS ${OPENEXR_INCLUDE_PATH} CACHE STRING "OpenEXR include paths")
   SET(OPENEXR_LIBRARIES ${OPENEXR_HALF_LIBRARY} ${OPENEXR_IEX_LIBRARY} ${OPENEXR_IMATH_LIBRARY} ${OPENEXR_ILMIMF_LIBRARY} ${OPENEXR_ILMTHREAD_LIBRARY} CACHE STRING "OpenEXR Libraries")
ENDIF ()

MARK_AS_ADVANCED(
  OPENEXR_FOUND
  OPENEXR_INCLUDE_PATHS
  OPENEXR_LIBRARIES
  OPENEXR_INCLUDE_PATH
  OPENEXR_HALF_LIBRARY
  OPENEXR_IEX_LIBRARY
  OPENEXR_IMATH_LIBRARY
  OPENEXR_ILMIMF_LIBRARY
  OPENEXR_ILMTHREAD_LIBRARY
)
