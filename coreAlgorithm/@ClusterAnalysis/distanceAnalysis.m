function [] = distanceAnalysis(obj, maxDistance, isRandom, showPlot)
%distanceAnalyses calculates pairwise distances
%
% Image is divided into a grid with a cell size specified by maxDistance.
% Distances from points in the outer cells of the grid are not caluclated to
% avoid edge effects.
%
% Image based waitbar displays progress of computation. This could be
% replaced in future versions using a conventional waitbar.
% Variable distances changes size in every loop iteration. Future version
% should implement a cell array.
%
% Input:
%   maxDistance: maximum allowed distance to which distances are calculated in nm (default = 200)
%   isRandom: true if random positons should be used
%   showPlot: show histogram of NN distribution (default = true)
% Output:
%   (1) distances of all points
%
%   Note: if clustering with DBSCAN was done before, the function
%   automatically calculates distances of all the points, of points
%   within cluster and of the points outside of the cluster
%   Additional output:
%   (2) nearest neighbor distances of points within clusters
%   (3) nearest neighbor distances of points outside of clusters
%
%   Matlab 2014b or newer and Statistics and Machine Learning Toolbox
%
% Important when using histogram:
% data are normalized to all distances from 0 to maxDistance
%
% Jan Neumann, 13.07.2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% init
isDBSCAN = false; % remembers if DBSCAN data already exist
switch nargin
    case 1
        maxDistance = 200;
        isRandom = 0;
        showPlot = 0;
    case 2
        isRandom = 0;
        showPlot = 0;
    case 3
        showPlot = 0;
    case 4
        
    otherwise
        error('Wrong number of input arguments!')
end
if isRandom == true
    positions = obj.randomTable;
    dataMat = obj.randomClusterStruct;
    if isfield(obj.randomClusterStruct, 'clusterDBSCAN')
        clusterAssignment = obj.randomClusterStruct(3).clusterDBSCAN;
        numberOfSample = obj.clusterStruct(1).clusterDBSCAN;
        isDBSCAN = true;
    end
else
    positions = obj.positionTable;
    dataMat = obj.clusterStruct;
    if isfield(obj.clusterStruct, 'clusterDBSCAN')
        clusterAssignment = obj.clusterStruct(3).clusterDBSCAN;
        isDBSCAN = true;
    end
end
if isRandom == true
    if obj.dimension == 3
        dataSize = min(obj.randomPhysicalDimension(1:3));
    else
        dataSize = min(obj.randomPhysicalDimension(1:2));
    end
else
    if obj.dimension == 3
        dataSize = min(obj.physicalDimension(1:3));
    else
        dataSize = min(obj.physicalDimension(1:2));
    end
end
if dataSize < 3*maxDistance
    disp('Image size is smaller then 3 times the cut-off distance! Try to decrease cut-off distance or increase image region.');
    return
end
%% distance calculation
% re-sort points in such a way that boundary points are at the end of the
% list
if obj.dimension == 3
    mapSize = ceil(max(positions(:, 1:3)));
    currentPositions = positions(positions(:, 1) > maxDistance & positions(:, 1) < mapSize(1)-maxDistance...
        & positions(:, 2) > maxDistance & positions(:, 2) < mapSize(2)-maxDistance...
        & positions(:, 3) > maxDistance & positions(:, 3) < mapSize(3) - maxDistance, :);
    outsidePositions = positions( (positions(:, 1) <= maxDistance | positions(:, 1) >= mapSize(1)-maxDistance)...
        | (positions(:, 2) <= maxDistance | positions(:, 2) >= mapSize(2)-maxDistance)...
        | positions(:, 3) <= maxDistance | positions(:, 3) >= mapSize(3) - maxDistance, :);
    positions = cat(1, currentPositions, outsidePositions);
else
    mapSize = ceil(max(positions(:, 1:2)));
    currentPositions = positions(positions(:, 1) > maxDistance & positions(:, 1) < mapSize(1)-maxDistance...
        & positions(:, 2) > maxDistance & positions(:, 2) < mapSize(2)-maxDistance, :);
    outsidePositions = positions( (positions(:, 1) <= maxDistance | positions(:, 1) >= mapSize(1)-maxDistance)...
        | (positions(:, 2) <= maxDistance | positions(:, 2) >= mapSize(2)-maxDistance), :);
    positions = cat(1, currentPositions, outsidePositions);
end
% cell array for saving distances
distanceCell = cell(1, size(currentPositions, 1));
nDistanceCell = 1;
multiWaitbar('calculating distances...', 0);
% for waitbar
prevPercent = 0;
counter = 1;

for ii = 1 : size(currentPositions, 1)
    if obj.dimension == 3
        % check if point is within the image border
        if (positions(ii, 1) > maxDistance && positions(ii, 1) < mapSize(1)-maxDistance...
                && positions(ii, 2) > maxDistance && positions(ii, 2) < mapSize(2)-maxDistance...
                && positions(ii, 3) > maxDistance && positions(ii, 3) < mapSize(3) - maxDistance)
            % calculate distance of current point to all other points that have
            % not yet been visited
            segment = positions(ii:end, :);
            currentDistance = single(pdist2(segment(:, 1:3), positions(ii, 1:3)));
            % removes distances larger then maxDistance
            currentIdx = (currentDistance <= maxDistance & currentDistance > 0);
            currentDistance = currentDistance(currentIdx);
            distanceCell{nDistanceCell} = reshape(currentDistance, [], 1);
            nDistanceCell = nDistanceCell + 1;
        end
    else
        % check if point is within the image border
        if (positions(ii, 1) > maxDistance && positions(ii, 1) < mapSize(1)-maxDistance...
                && positions(ii, 2) > maxDistance && positions(ii, 2) < mapSize(2)-maxDistance)
            % calculate distance of current point to all other points that have
            % not yet been visited
            segment = positions(ii:end, :);
            currentDistance = single(pdist2(segment(:, 1:2), positions(ii, 1:2)));
            % removes distances larger then maxDistance
            currentIdx = (currentDistance <= maxDistance & currentDistance > 0);
            currentDistance = currentDistance(currentIdx);
            distanceCell{nDistanceCell} = reshape(currentDistance, [], 1);
            nDistanceCell = nDistanceCell + 1;
        end
    end
    % waitbar
    currentPercent = fix(100*counter/size(currentPositions, 1));
    if currentPercent > prevPercent
        multiWaitbar( 'calculating distances...', 'Value', counter/size(currentPositions, 1));
        prevPercent = currentPercent;
    end
    counter = counter + 1;
end
distances = (cell2mat(distanceCell.'));
%% visualization
if showPlot == true
    figure( 'Name', num2str(obj.sampleID) );
    histogram(distances, maxDistance, 'Normalization', 'probability'); % data are normalized to all distances with 0 and maxDistance
    grid on;
    title('Distance Analysis - all points');
    xlabel('distance [nm]');
    ylabel('normalized frequency');
    ax = gca;
    ax.XLim = [0 maxDistance];
end
%% output
dataMat(1).allDistancesZ = distances;
if isDBSCAN == true
    dataMat(2).allDistances = innerDistance; % points inside of clusters
    dataMat(3).allDistances = outerDistance; % points outside of clusters
end
if isRandom == true
    obj.randomClusterStruct = dataMat;
else
    obj.clusterStruct = dataMat;
end
multiWaitbar('CLOSEALL');
end