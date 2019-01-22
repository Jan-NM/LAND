classdef ClusterAnalysis < handle
%ClusterAnalysis Object class for performing cluster analysis of SMLM data
% 
% Class constructor extracts positions and localization precision of SMLM
% data and creates automatically a random position table. Input
% localization table should be in Orte file format i.e. 
% column 2/3 = x/y-position, column 3/5 = x/y-localization precision,
% column 9 = frame number
% If your localization table does not contain the x- and y-localization
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
%   e.g. cell01 = ClusterAnalysis(localization data, varargin);
%
% varargin comprise the following additional arguments that can be hand over
%   varargin{1}: name or number of the sample (char or numeric value)
%   varargin{2}: dimensionality of the sample (2 for 2D or 3 for 3D) - note 3D is not completely implemented 
%   varargin{3}: algorithm used to calculate local density. allows the following
%               parameters: 'nearestNeighbor', 'averageDensity',
%               'convolution', 'kernelDensity', 'specificValue'
%   varargin{4}: additional parameter for selected algorithm to create local density: varargin{3}
%   varargin{5}: algorithm used to create random data from calculated local
%               density (true = complete spatial randomness (poisson point process);
%               false = same number of points as in experimental data with calculated local density))
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
% column 7 = frame number (not in random table)
%
% The following cluster algorithms are currently implemented:
%   DBSCAN (A Density-Based Algorithm for Discovering Clusters in Large
%   Spatial Databases with Noise. Ester at al. 1996)
%
%   k-nearest neighbor distance
%
%   distance analysis
%
%   pair correlation algortihm based on distance measure
%
%   currently in development:
%   Voronoi cluster algorithm (SR-Tesseler: a method to segment and 
%   quantify localization-based super-resolution microscopy data.
%   Levet et al. 2015; ClusterViSu, a method for clustering of protein 
%   complexes by Voronoi tessellation in super- resolution microscopy.
%   Andronov et al. 2016)
%
%   parameterEstimation(see DBSCAN)
%
%   scatterPlot
%
%   To perform cluster analysis type:
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
%   version is included in this distribution - but should be compiled on system used)
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
                        {'nearestNeighbor', 'averageDensity', 'convolution', 'kernelDensity', 'specificValue'});
                    switch randomAlgorithm
                        case 'nearestNeighbor'
                            % does not need any parameters
                            randomAlgorithmValue = [];                             
                        case 'averageDensity'
                            % does not need any parameters
                            randomAlgorithmValue = [];
                        case {'convolution', 'kernelDensity'}
                            % use half of the estimated spatial resolution
                            % as filter radius (sigma)
                            [~, nnDistance] = knnsearch(localizationFile(:, 2:3), localizationFile(:, 2:3), 'K', 2);
                            randomAlgorithmValue = 0.5 * sqrt( (mean(nnDistance(:, 2), 1))^2 + mean(mean(localizationFile(:, 4:5), 2))^2 );
                        case 'specificValue'
                            if nargin < 5
                                error('If you want to use a specific value, you must add it as argument to the function!')
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
                        case 'nearestNeighbor'
                            % does not need any parameters
                            randomAlgorithmValue = [];                             
                        case 'averageDensity'
                            % does not need any parameters
                            randomAlgorithmValue = [];
                        case {'convolution', 'kernelDensity'}
                            validateattributes(varargin{4}, {'numeric'}, {'nonnegative'})
                            if isempty(varargin{4}) 
                                warning('No value for random algorithm specified. Using default values!')
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
                        case 'nearestNeighbor'
                            % does not need any parameters
                            randomAlgorithmValue = [];                             
                        case 'averageDensity'
                            % does not need any parameters
                            randomAlgorithmValue = [];
                        case {'convolution', 'kernelDensity'}
                            validateattributes(varargin{4}, {'numeric'}, {'nonnegative'})
                            if isempty(varargin{4}) 
                                warning('No value for random algorithm specified. Using default values!')
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
            % extract positions and localization precision from
            % localization table in Orte format (i.e. column 2/3 = x/y-position, column 3/5 = x/y-localization precision, column 9 = frame number)
            obj.NoOfPoints = size(localizationFile, 1);
            obj.physicalDimension = [max(localizationFile(:, 2)); max(localizationFile(:, 3)); 0];
            % change to double if higher precision is needed, to reduce
            % memory consumption stay with single data type
            obj.positionTable = zeros(obj.NoOfPoints, 7, 'double');
            obj.positionTable(:, 1) = localizationFile(:, 2);
            obj.positionTable(:, 2) = localizationFile(:, 3);
            % add frame numbering
            obj.positionTable(:, 7) = localizationFile(:, 9);
            try
                obj.positionTable(:, 4) = localizationFile(:, 4);
                obj.positionTable(:, 5) = localizationFile(:, 5);
            catch
                warning('Localization precision is missing?');
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
                    warning('Unable to create Ataiya''s kdtree, using Matlab''s internal implementation.');
                    obj.queryData = createns(obj.positionTable(:, 1:3));
                    obj.flagSearchTree = 'matlabTree';
                end
                obj.physicalDimension(3, 1) = max(localizationFile(:, 3));
            else
                obj.positionTable(:, 3) = NaN;
                obj.positionTable(:, 6) = NaN;
                try
                    x = obj.positionTable(:, 1);
                    y = obj.positionTable(:, 2);
                    obj.queryData = kdtree_build([x y]);
                    obj.flagSearchTree = 'cTree';
                catch
                    warning('Unable to create Ataiya''s kdtree, using Matlab''s internal implementation.');
                    obj.queryData = createns(obj.positionTable(:, 1:2));
                    obj.flagSearchTree = 'matlabTree';
                end
            end
            obj.clusterStruct = struct([]);
            % choose algorithm for random point distribution
            switch randomAlgorithm
                case 'nearestNeighbor'
                    % assumes that one point accumulates a cicular area
                    % with half of the mean 2-NN-distance
                    [~, nnDistance] = knnsearch(obj.positionTable(:, 1:2), obj.positionTable(:, 1:2), 'K', 2);
                    obj.localDensity = 1 / (pi*((mean(nnDistance(:, 2), 1)) /1000 /2)^2); % in 1/ µm^2
                case 'averageDensity'
                    obj.localDensity = obj.NoOfPoints / (obj.physicalDimension(1, 1)/1000 * obj.physicalDimension(2, 1)/1000); % in µm^2
                case 'convolution' 
                    % estimation of local density using a disk shaped
                    % kernel and convolution
                    pixelSize = ceil(randomAlgorithmValue / 2.3); % pixel size calculated from estimated resolution with nyquist criteria
                    radius = randomAlgorithmValue; % radius of disk shaped kernel in nm
                    binnedImage = double(visModuleCluster(obj.positionTable, 'histogramBinning', pixelSize));
                    % create disk shaped kernel
                    kernel = fspecial('disk', radius);
                    imageConv = conv2(binnedImage, kernel);
                    localArea = sum(sum(imageConv > 0)) * (pixelSize/1000)^2; % local area in µm^2
                    obj.localDensity = obj.NoOfPoints / localArea; 
                 case 'kernelDensity'
                     kernelSize = randomAlgorithmValue; % in nm
                     radius = sqrt(-2*kernelSize^2 *  log(10^(-6)) ); % in nm, calculate neglectable distance
                     pointDensity = zeros(obj.NoOfPoints, 1, 'single');
                     % this part is slow and should be improved in feature
                     % versions
                     % for waitbar
                     warning('Kernel density algorithm is slow for large localization tables. This part will be improved in future versions!')
                     multiWaitbar('computing kernel density...', 0);
                     prevPercent = 0;
                     counter = 0;
                     for ii = 1:obj.NoOfPoints
                         [~, distances] = rangesearch(obj.positionTable(ii, 1:2), obj.queryData.X, radius);
                         distances = cell2mat(distances.');
                         kernel = 1 / (2*pi*kernelSize^2) .* exp(- (distances).^2 ./ (2*kernelSize^2));
                         pointDensity(ii, 1) = sum(kernel) / size(obj.NoOfPoints, 2) * 1000^2;
                         % waitbar
                         currentPercent = fix(100*counter/obj.NoOfPoints);
                         if currentPercent > prevPercent
                             multiWaitbar( 'computing kernel density...', 'Value', counter/obj.NoOfPoints);
                             prevPercent = currentPercent;
                         end
                         counter = counter + 1;
                     end
                     obj.localDensity = mean(pointDensity);
                case 'specificValue'
                     obj.localDensity = randomAlgorithmValue;
            end
            if assumeCSR == true
                % create CSR with estiamted localDensity
                nPoints = poissrnd(obj.localDensity*obj.physicalDimension(1, 1)/1000*obj.physicalDimension(2, 1)/1000); % 1000 because localDensity is in µm^2      
                obj.randomTable = zeros(nPoints, 6, 'double');
                obj.randomNoOfPoints = nPoints;
                obj.randomTable(:, 1) = (obj.physicalDimension(1, 1)).*rand(nPoints, 1);
                obj.randomTable(:, 2) = (obj.physicalDimension(2, 1)).*rand(nPoints, 1);
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
                    obj.flagSearchTree = 'cTree';
                catch
                    warning('Unable to create Ataiya''s kdtree, using Matlab''s internal implementation.');
                    obj.randomQueryData = createns(obj.randomTable(:, 1:3));
                    obj.flagSearchTree = 'matlabTree';
                end
                obj.randomPhysicalDimension = [max(obj.randomTable(:, 1)); max(obj.randomTable(:, 2));  max(obj.randomTable(:, 3))];
            else
                try
                    x = obj.randomTable(:, 1);
                    y = obj.randomTable(:, 2);
                    obj.randomQueryData = kdtree_build([x y]);
                    obj.flagSearchTree = 'cTree';
                catch
                    warning('Unable to create Ataiya''s kdtree, using Matlab''s internal implementation.');
                    obj.randomQueryData = createns(obj.randomTable(:, 1:2));
                    obj.flagSearchTree = 'matlabTree';
               end
                obj.randomPhysicalDimension = [max(obj.randomTable(:, 1)); max(obj.randomTable(:, 2)); 0];
            end
            obj.randomClusterStruct = struct([]);
            if isempty(obj.randomNoOfPoints) || obj.randomNoOfPoints < 2
                warning('Number of random points is lower than 2. This could lead to errors, if random data is used for calculation. Maybe parameters specified in calculation of random data are not choosen appropriate.')
            end
        end
        % destructor
        function delete(obj)
            try
                kdtree_delete(obj.queryData);
                kdtree_delete(obj.randomQueryData);
            catch
            end
        end
        % implemented cluster algorithms
        [varargout] = DBSCAN(obj, radius, minPoints, isRandom, maxDiameter, showImage)
        
        [varargout] = DBSCANdistanceMatrix(obj, radius, minPoints, isRandom, maxDiameter, showImage)
        
        % [varargout] = kNNDistance(obj, k, isRandom, maxDistance, showPlot)
        
        [varargout] = radialDensityFunction(obj, binSize,  maxRadius, isRandom, showImage)
        
        [varargout] = radialDensityFunctionDistanceMatrix(obj, binSize,  maxRadius, isRandom, showImage)
        
        [varargout] = ripley(obj, samplingDistance,  maxRadius, isRandom, showImage)
        
        [varargout] = ripleyDistanceMatrix(obj, samplingDistance,  maxRadius, isRandom, showImage)
        
        % [varargout] = voronoiCluster(obj, isRandom, showImage)
        
        [varargout] = distanceAnalysis(obj, maxDistance, isRandom, showPlot)
        
        % [varargout] = gridAnalysis(obj, gridSize, showPlot)
        
        % [varargout] = parameterEstimation(obj, method)
        
        [varargout] = scatterPlot(obj)

    end
end