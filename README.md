# Integrated Localization Analyzer for Nanoscale Distributions (ILAND)

ILAND is a software package for Matlab that enables quantitative 2D and 3D analysis of single molecule localization microscopy (SMLM) data. The package includes density-based clustering algorithms like DBSCAN and algorithms used in spatial statistics like the radial density function or Ripley's function. In addition, it contains algorithms for quantifying the conformation and texture of the nuclear nanostructure. ILAND has been specifically designed for the evaluation of large sample sizes and data with high emitter densities.

## Getting Started

The follwong sections describe how to get a copy of the software and how to install it on your local machine. Detailed instructions on how to use the software including examples can be found in the manual [ILAND_manual.pdf](ILAND_manual.pdf).

### Requirements

* Matlab R2014b or newer
	* Statistics and Machine Learning Toolbox 
	* Image Processing Toolbox

* *(optional, but highly recommended for a much faster computation)* [Andrea Tagliasacchi's kdtree](https://github.com/ataiya/kdtree)

* [multiWaitbar](https://de.mathworks.com/matlabcentral/fileexchange/26589-multiwaitbar-label-varargin) (a copy is included in this distribution)

At least 8 GByte RAM are recommended.

### Installing

* download the software package from https://github.com/Jan-NM/ILAND
* extract **ILAND-master.zip**
* copy the generated **ILAND-master** directory into your local Matlab working directory (on Windows machines, usually C:\Users\\*user name*\Documents\MATLAB)

To use Andrea Tagliasacchi's kdtree, download the following files from https://github.com/ataiya/kdtree/tree/master/toolbox:
* KDTree.h
* kdtree_ball_query.cpp
* kdtree_ball_query.m
* kdtree_build.cpp
* kdtree_build.m
* kdtree_delete.cpp
* kdtree_delete.m
* MyHeaps.h

These files must now be compiled.
* create a folder named *kdtree* in the directory *...\ILAND-master\utilities\...** 
* configure the [MEX environment](https://de.mathworks.com/help/matlab/matlab_external/what-you-need-to-build-mex-files.html)
* download and install a Matlab-supported C++ compiler
The following description is for Windows machines:
* type *mex -setup C++* in Matlab's command window to setup the mex environment
* open **...\ILAND-master\utilities\ataiya_kdtree\...** in Matlab's current folder panel
* type *mex kdtree_build* in Matlab's command window
* type *mex kdtree_ball_query* in Matlab's command window
* type *mex kdtree_delete* in Matlab's command window

* to use ILAND, right click on ILAND-master in Matlab's current folder panel, go to *Add to Path* and click on *Selected Folders and Subfolders*

ILAND can be used via the command window or by opening the user interface. To open the user interface type *startClusterAnalysis* in the command window. Detailed instructions on how to use the software including examples can be found in the manual [ILAND_manual.pdf](ILAND_manual.pdf).

## License

ILAND is licensed under the GNU GPL - see the [License.md](License.md) file for details. ILAND includes [multiWaitbar](https://de.mathworks.com/matlabcentral/fileexchange/26589-multiwaitbar-label-varargin), which comes with a separate license.

## Notes

DBSCAN is based on the paper: Ester at al., "A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases with Noise", (1996).
