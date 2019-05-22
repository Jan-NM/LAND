function [] = ripley( obj, samplingDistance, maxRadius, isRandom, showPlot)
%ripley computes Ripley's function using the point
%       coordinates
%
% Input:
%   samplingDistance:   step size for the radius, in nm
%   maxRadius:          maximum distance up to which the ripley function is calculated
%                       (only points with this minimum distance to the image boundary 
%                       are considered), in nm
%   isRandom:           true if random positions should be used
%   showPlot:           plot ripley finction
%   
% Output:
%   (1) L(r)-r values for each bin
%   (2) center of each bin, in nm
%
% Jan Neumann, 02.06.2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% init
multiWaitbar('computing Ripley''s H-function...', 0);
switch nargin  
    case 3
        isRandom = 0;
        showPlot = 0;
    case 4
        showPlot = 0;
    case 5

    otherwise
        error('Wrong number of input arguments!')
end
if isRandom == true
    positions = obj.randomTable;
    dataMat = obj.randomClusterStruct;
    nPoints = obj.randomNoOfPoints;
    treeData = obj.randomQueryData;
    flagSearchTree = obj.flagSearchTreeRandomData;
    physicalDimension(1, :) = obj.randomPhysicalDimension(1, :);
    physicalDimension(2, :) = obj.randomPhysicalDimension(2, :);
    if obj.dimension == 3
        physicalDimension(3, :) = obj.randomPhysicalDimension(3, :);
    end
else
    positions = obj.positionTable;
    dataMat = obj.clusterStruct;
    nPoints = obj.NoOfPoints;
    treeData = obj.queryData;
    flagSearchTree = obj.flagSearchTree;
    physicalDimension(1, :) = obj.physicalDimension(1, :);
    physicalDimension(2, :) = obj.physicalDimension(2, :);
    if obj.dimension == 3
        physicalDimension(3, :) = obj.physicalDimension(3, :);
    end
end
% initialize output fields of dataMat
dataMat(1).ripley = [];
dataMat(2).ripley = [];
%% compute Ripley's function 
% create bins for Ripley histogram
nbins = ceil(maxRadius/samplingDistance);
[binArray, edges] = histcounts([0 maxRadius], nbins);
binArray = zeros(size(binArray, 2), 1);

% to avoid edge effects only points within a minimum distance of
% maxDiameter to the image boundaries are considered
if obj.dimension == 3
    mapSize = ceil(max(positions(:, 1:3)));
    currentPositions = positions(positions(:, 1) > maxRadius & positions(:, 1) < mapSize(1) - maxRadius...
        & positions(:, 2) > maxRadius & positions(:, 2) < mapSize(2) - maxRadius...
        & positions(:, 3) > maxRadius & positions(:, 3) < mapSize(3) - maxRadius, :);
else
    mapSize = ceil(max(positions(:, 1:2)));
    currentPositions = positions(positions(:, 1) > maxRadius & positions(:, 1) < mapSize(1)-maxRadius...
        & positions(:, 2) > maxRadius & positions(:, 2) < mapSize(2)-maxRadius, :);
end
% check of currentPositions is empty
if isempty(currentPositions)
    error('Image size is to small or maxRadius is to big! Can not find any points which have maxRadius distance from the image borders! Try to decrease maxRadius or increase image region.');
end
nPointsFinal = size(currentPositions, 1);

% for waitbar
prevPercent = 0;
counter = 1;
% loop through all points
for ii = 1:size(currentPositions, 1)
    % loop through all radii
    for ll = 1:size(binArray, 1)
        % get number of points within radii
        switch flagSearchTree
            case 'cTree'
                radDistances = kdtree_ball_query(treeData, currentPositions(ii, 1:obj.dimension), edges(ll + 1));
            case 'matlabTree'
                distance = rangesearch(currentPositions(ii, 1:obj.dimension), treeData.X, edges(ll + 1));
                radDistances = (find(~cellfun('isempty', distance))).';
        end
        binArray(ll, 1) = binArray(ll, 1) + numel(radDistances) - 1; % sum up points within r, ignore current point
    end
    % waitbar
    currentPercent = fix(100*counter/nPointsFinal);
    if currentPercent > prevPercent
        multiWaitbar( 'computing Ripley''s H-function...', 'Value', counter/nPointsFinal);
        prevPercent = currentPercent;
    end
    counter = counter + 1;
end

%% normalize data
if obj.dimension == 3
    density = nPoints / (physicalDimension(1, :) * physicalDimension(2, :) * physicalDimension(3, :));
else
    density = nPoints / (physicalDimension(1, :) * physicalDimension(2, :));
end
binArray = binArray ./ nPointsFinal ./density;
for ll = 1:size(binArray, 1)
    if obj.dimension == 3
        binArray(ll) = nthroot(binArray(ll, 1) / pi * 3/4, 3) - (edges(1, ll) + samplingDistance);
    else
        binArray(ll) = sqrt(binArray(ll, 1) / pi) - (edges(1, ll) + samplingDistance);
    end
end
multiWaitbar( 'computing Ripley''s H-function...', 'Close');
%% visualization
if showPlot == true
    figure( 'Name', num2str(obj.sampleID) );
    plot(edges(1 ,1:end-1) + samplingDistance/2, binArray);
    grid on;
    title('Ripleys function');
    xlabel('distance [nm]');
    ylabel('L(r)-r');
    xlim([min(edges) max(edges)]);
end
%% output
dataMat(1).ripley = binArray;
dataMat(2).ripley = (edges(1, 1:end-1) + samplingDistance/2).';
if isRandom == true
    obj.randomClusterStruct = dataMat;
else
    obj.clusterStruct = dataMat;
end
end