function [] = distanceAnalysis(obj, maxDistance, isRandom, showPlot)
%distanceAnalyses calculates pairwise distances
%
% Calculates pairwise distances between detected signals. To
% avoid edge effects, distances from points that are located within 
% maxDistance from the image border are not caluclated .
%
% Input:
%   maxDistance: maximum allowed distance to which distances are calculated in nm (default = 200)
%   isRandom: true if random positons should be used (default = false)
%   showPlot: show histogram of distance distribution (default = false)
% Output:
%   (1) distances of all points
%
%   Note: If computations have already been carried out using DBSCAN, the function
%   automatically calculates the distances of all the points, of points
%   within the clusters and of the points outside of the clusters.
%   Additional output:
%   (2) pairwise distances of points within clusters
%   (3) pairwise distances of points outside of clusters
%
%   Matlab 2014b or newer and Statistics and Machine Learning Toolbox
%
% Note:
% Histograms are normalized (according to all distances measured from zero to maxDistance).
%
% Jan Neumann, 13.07.2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% init
isDBSCAN = false; % flag, DBSCAN was carried out before
supressWaitbar = false; % flag for waitbar
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
    dataSize = min(obj.randomPhysicalDimension(1:obj.dimension));
    if isfield(obj.randomClusterStruct, 'clusterDBSCAN')
        clusterAssignment = obj.randomClusterStruct(3).clusterDBSCAN;
        isDBSCAN = true;
    end
else
    positions = obj.positionTable;
    dataMat = obj.clusterStruct;
    dataSize = min(obj.physicalDimension(1:obj.dimension));
    if isfield(obj.clusterStruct, 'clusterDBSCAN')
        clusterAssignment = obj.clusterStruct(3).clusterDBSCAN;
        isDBSCAN = true;
    end
end
if dataSize < 3*maxDistance
    disp('Image size is smaller then 3 times the cut-off distance! Try to decrease cut-off distance or increase image region.');
    return
end
%% perform distance calculation
imageSize = ceil(max(positions(:, 1:obj.dimension)));
distances = distanceCalculation(positions, obj.dimension, imageSize, maxDistance, supressWaitbar);
% distance calculation for points inside and outside of clusters
if isDBSCAN == true
    % for outside points clusterAssignment is 0
    tempCluster = clusterAssignment(:, 1) == 0;
    outerPositions = positions(tempCluster, 1:obj.dimension);
    outerDistances = distanceCalculation(outerPositions, obj.dimension, imageSize, maxDistance, supressWaitbar);
    % for points within clusters
    if any(clusterAssignment)
        supressWaitbar = true;
        innerDistances = cell(1, max(clusterAssignment));
        multiWaitbar('calculating distances...', 0);
        % for waitbar
        prevPercent = 0;
        counter = 1;
        for kk = 1:max(clusterAssignment)
            tempCluster = clusterAssignment(:, 1) == kk;
            innerPositions = positions(tempCluster, 1:obj.dimension);
            innerDistances{kk} = distanceCalculation(innerPositions, obj.dimension, imageSize, maxDistance, supressWaitbar);
            % waitbar
            currentPercent = fix(100*counter/max(clusterAssignment));
            if currentPercent > prevPercent
                multiWaitbar( 'calculating distances...', 'Value', counter/double(max(clusterAssignment)));
                prevPercent = currentPercent;
            end
            counter = counter + 1;
        end
        % remove empty cell arrays -> cluster had been outside the
        % boundary
        innerDistances = innerDistances(~cellfun('isempty',innerDistances));
        innerDistances = (cell2mat(innerDistances.'));
        multiWaitbar('calculating distances...', 'Close');
    else
        innerDistances = 0;
    end
end
%% visualization
if showPlot == true
    figure( 'Name', num2str(obj.sampleID) );
    if isDBSCAN == true && any(clusterAssignment)
        subplot(1, 3, 1)
        h = histogram(distances, ceil(max(distances(:, 1))), 'Normalization', 'probability');
        grid on;
        title('Distance Analysis - all points');
        xlabel('distance [nm]');
        ylabel('normalized frequency');
        ax = gca;
        ax.XLim = [0 maxDistance];
        subplot(1, 3, 2)
        h = histogram(innerDistances, ceil(max(innerDistances(:, 1))), 'Normalization', 'probability');
        grid on;
        title('Distance Analysis - points inside of clusters');
        xlabel('distance [nm]');
        ylabel('normalized frequency');
        ax = gca;
        ax.XLim = [0 maxDistance];
        subplot(1, 3, 3)
        h = histogram(outerDistances, ceil(max(outerDistances(:, 1))), 'Normalization', 'probability');
        grid on;
        title('Distance Analysis - points outside of clusters');
        xlabel('distance [nm]');
        ylabel('normalized frequency');
        ax = gca;
        ax.XLim = [0 maxDistance];
    else
        histogram(distances, maxDistance, 'Normalization', 'probability'); % data are normalized to all distances with 0 and maxDistance
        grid on;
        title('Distance Analysis - all points');
        xlabel('distance [nm]');
        ylabel('normalized frequency');
        ax = gca;
        ax.XLim = [0 maxDistance];
    end
end
%% output
dataMat(1).allDistances = distances;
if isDBSCAN == true
    dataMat(2).allDistances = innerDistances; % points inside of clusters
    dataMat(3).allDistances = outerDistances; % points outside of clusters
end
if isRandom == true
    obj.randomClusterStruct = dataMat;
else
    obj.clusterStruct = dataMat;
end
multiWaitbar('CLOSEALL');
end


function dist = distanceCalculation(positions, dim, mapSize, maxDistance, supressWaitbar)
% re-sort points in such a way that boundary points are at the end of the
% list
if dim == 3
    currentPositions = positions(positions(:, 1) > maxDistance & positions(:, 1) < mapSize(1)-maxDistance...
        & positions(:, 2) > maxDistance & positions(:, 2) < mapSize(2)-maxDistance...
        & positions(:, 3) > maxDistance & positions(:, 3) < mapSize(3) - maxDistance, :);
    outsidePositions = positions( (positions(:, 1) <= maxDistance | positions(:, 1) >= mapSize(1)-maxDistance)...
        | (positions(:, 2) <= maxDistance | positions(:, 2) >= mapSize(2)-maxDistance)...
        | positions(:, 3) <= maxDistance | positions(:, 3) >= mapSize(3) - maxDistance, :);
    positions = cat(1, currentPositions, outsidePositions);
else
    currentPositions = positions(positions(:, 1) > maxDistance & positions(:, 1) < mapSize(1)-maxDistance...
        & positions(:, 2) > maxDistance & positions(:, 2) < mapSize(2)-maxDistance, :);
    outsidePositions = positions( (positions(:, 1) <= maxDistance | positions(:, 1) >= mapSize(1)-maxDistance)...
        | (positions(:, 2) <= maxDistance | positions(:, 2) >= mapSize(2)-maxDistance), :);
    positions = cat(1, currentPositions, outsidePositions);
end
% cell array for saving distances
distanceCell = cell(1, size(currentPositions, 1));
nDistanceCell = 1;
if supressWaitbar == false
    multiWaitbar('calculating distances...', 0);
    % for waitbar
    prevPercent = 0;
    counter = 1;
end

for ii = 1 : size(currentPositions, 1)
    % calculate distance of current point to all other points that have
    % not yet been visited
    currentDistance = single(pdist2(positions(ii:end, 1:dim), positions(ii, 1:dim)));
    % removes distances larger then maxDistance
    currentIdx = (currentDistance <= maxDistance & currentDistance > 0);
    currentDistance = currentDistance(currentIdx);
    distanceCell{nDistanceCell} = reshape(currentDistance, [], 1);
    nDistanceCell = nDistanceCell + 1;
    if supressWaitbar == false
        % waitbar
        currentPercent = fix(100*counter/size(currentPositions, 1));
        if currentPercent > prevPercent
            multiWaitbar( 'calculating distances...', 'Value', counter/size(currentPositions, 1));
            prevPercent = currentPercent;
        end
        counter = counter + 1;
    end
end
dist = (cell2mat(distanceCell.'));
if supressWaitbar == false
    multiWaitbar('calculating distances...', 'Close');
end
end