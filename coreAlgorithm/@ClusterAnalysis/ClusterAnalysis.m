classdef ClusterAnalysis < handle
%ClusterAnalysis Object class for performing cluster analysis of SMLM data
% 
% Class constructor extracts positions and localization precision of SMLM
% data and creates automatically a random position table. Input
% localization table should be in Orte file format i.e. 
% column 2/3/11 = x/y/z-position, column 4/5/12 = x/y/z-localization precision,
% column 9 = frame number
% If your localization table does not contain the x-, y- or z localization
% precisions, a warning will be displayed. One option to avoid possible
% errors due to missing localization precisions could be to use an average
% localization precision as calculated by nearest neighbor based analysis
% (Ulrike Endesfelder, Sebastian Malkusch, Franziska Fricke, and Mike Heilemann (June
% 2014). “A simple method to estimate the average localization precision of a single-molecule
% localization microscopy experiment.” In: Histochem. Cell Biol. 141.6,
% pp. 629–38. ISSN: 1432-119X. DOI: 10.1007/s00418- 014- 1192- 3) for example.
%
% To create an instance of this class call:
%
%   variable name = ClusterAnalysis(localization data, varargin);
%
% varargin comprise the following additional arguments that can be handed over
%   varargin{1}: name or number of the sample (char or numeric value)
%   varargin{2}: dimensionality of the sample (2 for 2D or 3 for 3D) -
%                note: some function may not fully support 3D data
%   varargin{3}: algorithm used to calculate local density. allows the following
%                parameters: 'averageDensity', 'kernelDensity', 'specificValue'
%                if 'kernelDensity' is selected and varargin{4} is empty
%                default paramters are used: sampling distance = 100 nm and
%                grid spacing = 10 nm
%   varargin{4}: additional parameter for selected algorithm to create local density: varargin{3}
%                if 'kernelDensity' is selected the input should be a
%                1x2 matrix, where the first argument defines the sampling
%                distance (default: 100 nm) and the second argument the grid
%                spacing (default: 10 nm)
%   varargin{5}: algorithm used to create random data from calculated local
%                density (true = complete spatial randomness (poisson point process);
%                false = same number of points as in experimental data with calculated local density))
%
% default values:
%   varargin{1}: '001'
%   varargin{2}: 2
%   varargin{3}: 'averageDensity'
%   varargin{4}: empty array
%   varargin{5}: true
%
% The created position table / random table has the following order:
% column 1/2/3 = x/y/z-position
% column 4/5/6 = x/y/z-localization precision
% column 7 = frame number (not available in random table)
%
% The following analysis algorithms are currently implemented:
%
%   DBSCAN (A Density-Based Algorithm for Discovering Clusters in Large
%   Spatial Databases with Noise. Ester at al. 1996) - 2D + 3D
%
%	Voronoi - 2D + 3D
%
%   k-nearest neighbor distance - 2D + 3D
%
%   distance analysis - 2D + 3D
%
%   radial density function - 2D + 3D
%
%   ripley's function - 2D + 3D
%
%   compaction parameter - 2D
%
%   parameter estimation for DBSCAN - 2D
%
%   plotting of data (scatterplot) - 2D + 3D
%
% To perform cluster analysis type:
%   e.g.  cell01.function name(parameter values);
%
%
% Requirements 
%   Matlab version: 2014b or newer recommended (some commands in some functions
%   like boundary or histogram do not work on older versions).
%   Statistics and Machine Learning Toolbox 
%   Image Processing Toolbox
%   multiWaitbar (https://de.mathworks.com/matlabcentral/fileexchange/26589-multiwaitbar-label-varargin,
%   from Matlab File Exchange, included in this distribution)
%   Andrea Tagliasacchi's C++ kdtree implementation(https://github.com/ataiya/kdtree,
%   version is not included in this distribution)
%   At least 8 GByte RAM are recommended.
%
%
% Jan Neumann, 27.06.17
%
%     Copyright (C) <2017>  <Jan Neumann>
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (SetAccess = protected)
        clusterStruct
        randomClusterStruct
        positionTable
        randomTable
        queryData
        randomQueryData
        dimension
        NoOfPoints
        randomNoOfPoints
        localDensity
        physicalDimension
        randomPhysicalDimension
        sampleID
    end
    
    properties (Access = private, Hidden = true)
       flagSearchTree
       flagSearchTreeRandomData
    end
        
    methods
        % constructor
        function obj = ClusterAnalysis(localizationFile, varargin)
            %% input parser           
            % check number of input arguments
            narginchk(1, 6);
                        
            % default values
            obj.dimension = 2;
            obj.sampleID = '001';
            % set default values if random algorithm is specified
            if nargin >= 4
                randomAlgorithm = varargin{3};
                    randomAlgorithm = validatestring(randomAlgorithm, ...
                        {'averageDensity', 'kernelDensity', 'specificValue'});
                    switch randomAlgorithm                           
                        case 'averageDensity'
                            % does not need any parameters
                            randomAlgorithmValue = [];
                        case 'kernelDensity'
                            % first value is sampling distance, second value is grid spacing
                            randomAlgorithmValue = [100 10]; % in nm
                        case 'specificValue'
                            if nargin < 5
                                error('If you want to use a certain value for the density of the signals, you have to add it as an argument to the function!')
                            end
                    end    
            else            
                randomAlgorithm = 'averageDensity';
                randomAlgorithmValue = [];
            end
            assumeCSR = true;
            % check input arguments            
            switch nargin
                case 1
                    % nothing to do, default values are already assigned
                case 2
                    validateattributes(varargin{1}, {'char', 'numeric'}, {'nonempty'})
                    obj.sampleID = varargin{1};
                case 3
                    validateattributes(varargin{1}, {'char', 'numeric'}, {'nonempty'})
                    obj.sampleID = varargin{1};
                    switch varargin{2}
                        case 2
                            obj.dimension = 2;
                        case 3
                            obj.dimension = 3;
                        otherwise
                            error('Unknown dimension!')
                    end
                case 4
                    validateattributes(varargin{1}, {'char', 'numeric'}, {'nonempty'})
                    obj.sampleID = varargin{1};
                    switch varargin{2}
                        case 2
                            obj.dimension = 2;
                        case 3
                            obj.dimension = 3;
                        otherwise
                            error('Unknown dimension!')
                    end
                    % random algortihm already set with default values for
                    % randomAlgorithmValue
                case 5
                    validateattributes(varargin{1}, {'char', 'numeric'}, {'nonempty'})
                    obj.sampleID = varargin{1};
                    switch varargin{2}
                        case 2
                            obj.dimension = 2;
                        case 3
                            obj.dimension = 3;
                        otherwise
                            error('Unknown dimension!')
                    end
                    switch randomAlgorithm                       
                        case 'averageDensity'
                            % does not need any parameters
                            randomAlgorithmValue = [];
                        case 'kernelDensity'
                            validateattributes(varargin{4}, {'numeric'}, {'nonnegative', 'size', [1, 2]})
                            if isempty(varargin{4}) 
                                warning('No value for random algorithm specified. Using default values!')
                                % first value is sampling distance, second value is grid spacing
                                randomAlgorithmValue = [100 10]; % values are in nm
                            else    
                                randomAlgorithmValue = varargin{4};
                            end
                        case {'specificValue'}     
                            validateattributes(varargin{4}, {'numeric'}, {'nonnegative'})
                            randomAlgorithmValue = varargin{4};
                    end    
                case 6
                    validateattributes(varargin{1}, {'char', 'numeric'}, {'nonempty'})
                    obj.sampleID = varargin{1};
                    switch varargin{2}
                        case 2
                            obj.dimension = 2;
                        case 3
                            obj.dimension = 3;
                        otherwise
                            error('Unknown dimension!')
                    end
                    switch randomAlgorithm                       
                        case 'averageDensity'
                            % does not need any parameters
                            randomAlgorithmValue = [];
                        case 'kernelDensity'
                            validateattributes(varargin{4}, {'numeric'}, {'nonnegative', 'size', [1, 2]})
                            if isempty(varargin{4}) 
                                warning('No value for random algorithm specified. Using default values!')
                                % first value is sampling distance, second value is grid spacing
                                randomAlgorithmValue = [100 10]; % values are in nm
                            else    
                                randomAlgorithmValue = varargin{4};
                            end
                        case {'specificValue'}     
                            validateattributes(varargin{4}, {'numeric'}, {'nonnegative'})
                            randomAlgorithmValue = varargin{4};
                    end
                    assumeCSR = varargin{5};
            end
            %% create cluster class with specified parameters
            % if cooordinates in list of localizations contains negativ values,
            % they are shifted towards positive values
            % checks also for NaN values and empty fields
            if  min(localizationFile(:, 2)) < 0
                localizationFile(:, 2) = localizationFile(:, 2) - min(localizationFile(:, 2)) + 1;
                warning('Negative values detected in list of localizations. Shifting coordinates towards positive values!')
            end
            if min(localizationFile(:, 3)) < 0
                localizationFile(:, 3) = localizationFile(:, 3) - min(localizationFile(:, 3)) + 1;
                warning('Negative values detected in list of localizations. Shifting coordinates towards positive values!')
            end
            if any(isnan(localizationFile(:, 2))) || any(isnan(localizationFile(:, 3))) ||...
                    any(isempty(localizationFile(:, 2))) || any(isempty(localizationFile(:, 3)))
                warning('List of localizations contains NaN values or empty fields!')
            end
            if obj.dimension == 3
                try
                    if min(localizationFile(:, 11)) < 0
                        localizationFile(:, 11) = localizationFile(:, 11) - min(localizationFile(:, 11)) + 1;
                        warning('Negative values detected in list of localizations. Shifting coordinates towards positive values!')
                    end
                    if any(isnan(localizationFile(:, 11))) || any(isempty(localizationFile(:, 11)))
                        warning('List of localizations contains NaN values or empty fields!')
                    end
                catch
                    error('No z-coordinate found. Are you sure your data contains a z-coordinate and the data is in the correct format?');
                end
            end
            % extract positions and localization precision from
            % localization table in Orte format (i.e. column 2/3/11 = x/y/z-position, column 4/5/12 = x/y/z-localization precision, column 9 = frame number)
            obj.NoOfPoints = size(localizationFile, 1);
            obj.physicalDimension = [max(localizationFile(:, 2)); max(localizationFile(:, 3)); 0];
            % stay with double data type, otherwise the kdtree (cTree)
            % might not work
            obj.positionTable = zeros(obj.NoOfPoints, 7, 'double');
            obj.positionTable(:, 1) = localizationFile(:, 2);
            obj.positionTable(:, 2) = localizationFile(:, 3);
            % add frame numbering
            obj.positionTable(:, 7) = localizationFile(:, 9);
            try
                obj.positionTable(:, 4) = localizationFile(:, 4);
                obj.positionTable(:, 5) = localizationFile(:, 5);
            catch
                warning('Localization precision is missing!');
            end
            if obj.dimension == 3
                obj.positionTable(:, 3) = localizationFile(:, 11);
                obj.positionTable(:, 6) = localizationFile(:, 12);
                try
                    x = obj.positionTable(:, 1);
                    y = obj.positionTable(:, 2);
                    z = obj.positionTable(:, 3);
                    obj.queryData = kdtree_build([x y z]);
                    obj.flagSearchTree = 'cTree';
                catch
                    warning('Unable to create c++ based kdtree, using Matlab''s internal implementation.');
                    obj.queryData = createns(obj.positionTable(:, 1:3));
                    obj.flagSearchTree = 'matlabTree';
                end
                obj.physicalDimension(3, 1) = max(localizationFile(:, 11));
            else
                obj.positionTable(:, 3) = NaN;
                obj.positionTable(:, 6) = NaN;
                try
                    x = obj.positionTable(:, 1);
                    y = obj.positionTable(:, 2);
                    obj.queryData = kdtree_build([x y]);
                    obj.flagSearchTree = 'cTree';
                catch
                    warning('Unable to create c++ based kdtree, using Matlab''s internal implementation.');
                    obj.queryData = createns(obj.positionTable(:, 1:2));
                    obj.flagSearchTree = 'matlabTree';
                end
            end
            obj.clusterStruct = struct([]);
            % choose algorithm for random point distribution
            switch randomAlgorithm
                case 'averageDensity'
                    if obj.dimension == 3
                        obj.localDensity = obj.NoOfPoints / (obj.physicalDimension(1, 1)/1000 * obj.physicalDimension(2, 1)/1000 * obj.physicalDimension(3, 1)/1000); % in µm^3
                    else
                        obj.localDensity = obj.NoOfPoints / (obj.physicalDimension(1, 1)/1000 * obj.physicalDimension(2, 1)/1000); % in µm^2
                    end
                case 'kernelDensity'
                    multiWaitbar('calculating kernel density...', 'Busy');  
                    if obj.dimension == 3
                        xGridSpace = round((max(obj.positionTable(:, 1)) - min(obj.positionTable(:, 1))) / randomAlgorithmValue(2));
                        yGridSpace = round((max(obj.positionTable(:, 2)) - min(obj.positionTable(:, 2))) / randomAlgorithmValue(2));
                        zGridSpace = round((max(obj.positionTable(:, 3)) - min(obj.positionTable(:, 3))) / randomAlgorithmValue(2));
                        [x1, x2, x3] = meshgrid(min(obj.positionTable(:, 1)):randomAlgorithmValue(1):max(obj.positionTable(:, 1)),...
                            min(obj.positionTable(:, 2)):randomAlgorithmValue(1):max(obj.positionTable(:, 2)),...
                            min(obj.positionTable(:, 3)):randomAlgorithmValue(1):max(obj.positionTable(:, 3)));
                        x1 = x1(:);
                        x2 = x2(:);
                        x3 = x3(:);
                        xi = [x1 x2 x3];
                        bw = obj.NoOfPoints * (-1 / (obj.dimension + 4)); % bandwidth is estimated by Scott's rule
                        KDEestimate = mvksdensity(obj.positionTable(:, 1:3), xi, 'Bandwidth', bw);
                        % transform KDEestimate into density per µm²
                        x = linspace(min(x1), max(x1), xGridSpace);
                        y = linspace(min(x2), max(x2), yGridSpace);
                        z = linspace(min(x3), max(x3), zGridSpace);
                        KDEestimate = KDEestimate * randomAlgorithmValue(1)^3;
                        [xq, yq, zq] = meshgrid(x, y, z);
                        z = griddata(x1, x2, x3, KDEestimate, xq, yq, zq);
                        z = z ./ (max(x1) * max(x2)* max(x3)/(1000*randomAlgorithmValue(2))^3);
                        obj.localDensity = sum(sum(sum(z)));
                    else
                        xGridSpace = round((max(obj.positionTable(:, 1)) - min(obj.positionTable(:, 1))) / randomAlgorithmValue(2));
                        yGridSpace = round((max(obj.positionTable(:, 2)) - min(obj.positionTable(:, 2))) / randomAlgorithmValue(2));
                        [x1, x2] = meshgrid(min(obj.positionTable(:, 1)):randomAlgorithmValue(1):max(obj.positionTable(:, 1)), min(obj.positionTable(:, 2)):randomAlgorithmValue(1):max(obj.positionTable(:, 2)));
                        x1 = x1(:);
                        x2 = x2(:);
                        xi = [x1 x2];
                        KDEestimate = ksdensity(obj.positionTable(:, 1:2), xi);
                        % transform z from nm into density per µm²
                        x = linspace(min(x1), max(x1), xGridSpace);
                        y = linspace(min(x2), max(x2), yGridSpace);
                        KDEestimate = KDEestimate * randomAlgorithmValue(1).^2;
                        [xq, yq] = meshgrid(x, y);
                        z = griddata(x1, x2, KDEestimate, xq, yq);
                        z = z ./ (max(x1) * max(x2)/(1000*randomAlgorithmValue(2))^2);
                        obj.localDensity = sum(sum(z));
                    end
                    multiWaitbar('calculating kernel density...', 'Close');
                case 'specificValue'
                     obj.localDensity = randomAlgorithmValue;
            end
            if assumeCSR == true
                % create CSR with estimated localDensity
                if obj.dimension == 3
                    nPoints = poissrnd(obj.localDensity*obj.physicalDimension(1, 1)/1000*obj.physicalDimension(2, 1)/1000*obj.physicalDimension(3, 1)/1000); % 1000 because localDensity is in µm^3 
                else
                    nPoints = poissrnd(obj.localDensity*obj.physicalDimension(1, 1)/1000*obj.physicalDimension(2, 1)/1000); % 1000 because localDensity is in µm^2 
                end
                obj.randomTable = zeros(nPoints, 6, 'double');
                obj.randomNoOfPoints = nPoints;
                obj.randomTable(:, 1) = (obj.physicalDimension(1, 1)).*rand(nPoints, 1);
                obj.randomTable(:, 2) = (obj.physicalDimension(2, 1)).*rand(nPoints, 1);
                % loc. prec. of random table stays empty - if loc. prec. for
                % random data is desired, insert lines here
                if obj.dimension == 3
                    obj.randomTable(:, 3) = (obj.physicalDimension(3, 1)).*rand(nPoints, 1);
                end
            else
                % place measured points randomly, uniformly distributed in
                % an area calculated from localDensity
                randomImageSize = obj.NoOfPoints / obj.localDensity; % image size to incoporate points with local density
                obj.randomTable = zeros(obj.NoOfPoints, 6, 'double');
                obj.randomTable(:, 1) = (sqrt(randomImageSize)*1000).*rand(obj.NoOfPoints, 1); % image is squared
                obj.randomTable(:, 2) = (sqrt(randomImageSize)*1000).*rand(obj.NoOfPoints, 1);
                obj.randomNoOfPoints = obj.NoOfPoints;
                % loc. prec. of random table stays empty - if loc. prec. for
                % random data is desired, insert lines here
                if obj.dimension == 3
                    obj.randomTable(:, 3) = (sqrt(randomImageSize)*1000).*rand(obj.NoOfPoints, 1);
                end
            end
            % create kdtree for random data
            if obj.dimension == 3
                try
                    x = obj.randomTable(:, 1);
                    y = obj.randomTable(:, 2);
                    z = obj.randomTable(:, 3);
                    obj.randomQueryData = kdtree_build([x y z]);
                    obj.flagSearchTreeRandomData = 'cTree';
                catch
                    warning('Unable to create c++ based kdtree, using Matlab''s internal implementation.');
                    obj.randomQueryData = createns(obj.randomTable(:, 1:3));
                    obj.flagSearchTreeRandomData = 'matlabTree';
                end
                obj.randomPhysicalDimension = [max(obj.randomTable(:, 1)); max(obj.randomTable(:, 2));  max(obj.randomTable(:, 3))];
            else
                try
                    x = obj.randomTable(:, 1);
                    y = obj.randomTable(:, 2);
                    obj.randomQueryData = kdtree_build([x y]);
                    obj.flagSearchTreeRandomData = 'cTree';
                catch
                    warning('Unable to create c++ based kdtree, using Matlab''s internal implementation.');
                    obj.randomQueryData = createns(obj.randomTable(:, 1:2));
                    obj.flagSearchTreeRandomData = 'matlabTree';
               end
                obj.randomPhysicalDimension = [max(obj.randomTable(:, 1)); max(obj.randomTable(:, 2)); 0];
            end
            obj.randomClusterStruct = struct([]);
            if isempty(obj.randomNoOfPoints) || obj.randomNoOfPoints < 2
                warning('Number of random points is lower than 2. This could lead to errors, if random data is used for calculations.')
            end
        end
        % destructor
        function delete(obj)
            try
                kdtree_delete(obj.queryData);
                kdtree_delete(obj.randomQueryData);
                obj.queryData = NaN;
                obj.randomQueryData = NaN;
            catch
            end
        end
        % implemented algorithms
        [varargout] = DBSCAN(obj, radius, minPoints, isRandom, maxDiameter, showImage)
        
        [varargout] = kNNDistance(obj, k, isRandom, maxDistance, showPlot)
        
        [varargout] = radialDensityFunction(obj, binSize,  maxRadius, isRandom, showImage)
        
        [varargout] = ripley(obj, samplingDistance,  maxRadius, isRandom, showImage)
        
        [varargout] = distanceAnalysis(obj, maxDistance, isRandom, showPlot)
		
		[ varargout ] = voronoiCluster( obj, isRandom, showImage )
        
        [varargout] = compactionParameter(obj, varargin)
        
        [varargout] = parameterEstimation(obj, method)
        
        [varargout] = scatterPlot(obj)        
    end
end