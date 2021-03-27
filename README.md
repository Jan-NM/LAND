# Localization Analyzer for Nanoscale Distributions (LAND)

![GitHub tag (latest SemVer)](https://img.shields.io/github/v/tag/Jan-NM/LAND?sort=semver)
![GitHub](https://img.shields.io/github/license/Jan-NM/LAND)

LAND is a software package written in MATLAB that enables quantitative 2D and 3D analysis of single molecule localization microscopy (SMLM) data. The package includes density-based clustering algorithms like DBSCAN and algorithms used in spatial statistics like the radial density function or Ripley's function. In addition, it contains algorithms for quantifying the conformation and texture of the nuclear nanostructure. LAND has been specifically designed for the evaluation of large sample sizes and data with high emitter densities.

## Getting Started

The follwong sections describe how to get a copy of the software and how to install it on your local machine. Detailed instructions on how to use the software including examples can be found in the [manual](help/LAND_manual.pdf).

### Requirements

* MATLAB R2014b or newer
	* Statistics and Machine Learning Toolbox 
	* Image Processing Toolbox


* *(optional, but highly recommended for a much faster computation)* [ataiya/kdtree](https://github.com/ataiya/kdtree)

* [multiWaitbar](https://de.mathworks.com/matlabcentral/fileexchange/26589-multiwaitbar-label-varargin) (a copy is included in this distribution)

At least 8 GByte RAM are recommended.

### Installation

* download the software package from https://github.com/Jan-NM/LAND/releases
* extract `LAND-master.zip`
* copy the generated `LAND-master` directory into your local MATLAB working directory
* to use LAND, right click on `LAND-master` in MATLAB's current folder panel, go to `Add to Path` and click on `Selected Folders and Subfolders`

LAND can be used via the command window or by opening a user interface. To open the user interface type `startClusterAnalysis` in the command window. Detailed instructions on how to use the software including examples can be found in the [manual](help/LAND_manual.pdf).

**Recommended Installation: ataiya/kdtree**

* configure a [MEX environment](https://de.mathworks.com/support/requirements/supported-compilers.html) by downloading and installing a Matlab-supported C++ compiler
* open "...\LAND-master\utilities\ataiya_kdtree\..." in Matlab's current folder panel
* execute the following commands in MATLAB's command window to setup the mex environment
```
mex -setup C++
mex kdtree_build.cpp
mex kdtree_ball_query.cpp
mex kdtree_delete.cpp
```

### Input Data Format

Input data should be in contained in a numeric nx12 `.mat` file. Each row (n) should correspond to a single molecule signal. The columns should be ordered in the following way:
- *column 2/3/11 = x/y/z-position*
- *column 4/5/12 = x/y/z-localization precision*
- *column 9 = frame number*

The remaining columns can be filled with zeros. If the data does not contain any localization precision and z-position (i.e. 2D data), these columns should be also filled with zeros.

## Contributing and Support

You are encouraged to contribute to this project. Feel free to open Issues for feature contributions or bug reports.

## License

LAND is licensed under the GNU GPL - see the [LICENSE](LICENSE) file for details. LAND includes [multiWaitbar](https://de.mathworks.com/matlabcentral/fileexchange/26589-multiwaitbar-label-varargin), which comes with a separate license.

## Included Third-Party Software Packages

* [multiWaitbar](https://de.mathworks.com/matlabcentral/fileexchange/26589-multiwaitbar-label-varargin)

* ataiya/kdtree: Andrea Tagliasacchi (2020). ataiya/kdtree (https://www.github.com/ataiya/kdtree), GitHub. Retrieved April 5, 2020. 

## Cite As

If you are using this code in one of your publications, please cite this [paper](https://doi.org/10.1039/C9NR00943D):

Neumann et al., "Nanoscale distribution of TLR4 on primary human macrophages stimulated with LPS and ATI", (2019).

## Notes

DBSCAN is based on the paper: Ester at al., "A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases with Noise", (1996).
