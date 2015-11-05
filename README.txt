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
*  Authors: Norbert Bus & Nabil Mustafa
*  
*  If you find this code useful, please cite it by referencing the following paper:
*  @article{BMB15, title={Illumination{Cut}}, author={Bus, Norbert and Mustafa, Nabil H. and Biri, Venceslas},
*           journal = {Computer Graphics Forum (Proceedings of Eurographics 2015)},year={2015}}
*/

===================================  IlluminationCut v1.0 ===================================

Brief: 
Implementation of a few rendering algorithms in the many-lights framework with the Embree renderer and ray tracing kernels. For detailed descriptions use the corresponding papers.
The algorithms are the following (with the abbreviations used for them):
-ct: Lightcuts with agglomerative clustering for the lightree              (implemented by the authors of IlluminationCut)
-ic: IlluminationCut with agglomerative clustering for the lightree        (implemented by the authors of IlluminationCut)
-lc: Lightcuts with divisive clustering for the lightree                   (implemented by the authors of LightSlice)
-md: Multidimensional Lightcuts with divisive clustering for the lightree  (implemented by the authors of LightSlice)
-ls: LightSlice                                                            (implemented by the authors of LightSlice)
The relevant source files are in the folder renderer\devices\device_singleray\clustering.

Disclaimer: 
This is not the exact same version of the source code that is used to measure performance for the corresponding paper. Performance might have been affected during refactoring.

Compilation on linux:
Tested on Ubuntu 14 x64 with gcc 4.7. Use the following commands to compile and run the rendering engine.
You will need the glut and gsl libraries, hence run the  following commands:
 1. sudo apt-get install freeglut3 freeglut3-dev
 2. sudo apt-get install libxmu-dev libxi-dev
 3. sudo apt-get install libgsl0ldbl
 4. sudo apt-get install gsl-bin libgsl0-dev
Uncompress the project and cd into it. Then execute the following commands:
 5. cd renderer
 6. mkdir build
 7. cd build
 8. ccmake ..
 9. make
Run the application with the provided configuration file (models/conference/conference.ecs). 
10. ./renderer -c ../models/conference/conference.ecs

 
Compilation on Windows x64 with Visual Studio 2010:
Assuming you have the source located at yourpath
1. add the folder yourpath/embree_rtc/lib/x64 to your PATH
2. set the EMBREE_INSTALL_DIR environment variable to the main folder of Embree i.e., yourpath/embree_rtc.
3. use the provided Microsoft Visual Studio 2010 solution file (in the folder yourpath/renderer).
4. run the application in VS2010 with the provided configuration file (yourpath/renderer/models/conference/conference.ecs) by setting the command line arguments to:
 -c yourpath/renderer/models/conference/conference.ecs

Usage Details:
Detailed readme files for Embree can be found in the folders: embree_rtc (embree_rtc_linux) and renderer. The new command line parameters are all listed and explained in the renderer/models/conference/conference.ecs file. 

For the sake of completeness we list them here as well:
Location of the model file
-i conference.obj
Setting the used renderer. Different renderers are: vplrenderer, ctrenderer, icrenderer, lsrenderer, lcrenderer, mdrenderer. Simply substitute the one to run into the below option.
-renderer ctrenderer { depth = 9 }
Rendered image size
-size 256 256
The number of samples per pixel
-spp 1
The number of virtual point lights
-vplnumber 1000
The file to read/write the vpl configuration (if the file is present it overwrites the vplnumber setting)
-vplfile conference.vpl
The maximum error for the methods in the lightcuts family
-error 0.01
Option to swith sampling on or off for illuminationcut: full, sampling
-sampling sampling
Parameter that controls the slice number for lightslice
-seed 400
Parameter that controls the columns number for lightslice
-budget 400
The output filename (currently only tga but other formats can be supported as well)
-o conference.tga 
The setup of lightsources. Complete list of different types (e.g. spotlight) can be found in the embree renderer readme or other ecs files. For quadlight the parameters are the position of a corner (3) two vectors pointing to the neighboring vertices (6) and intensity (3)
-quadlight -120 100 -150  0 0 80 50 0 0 50 50 50	
Camera position
-vp 1478 553 -1000 
Camera view direction
-vi 478 173 000 
Camera up
-vu 0 1 0
Camera field of view 
-fov 50 
  
 