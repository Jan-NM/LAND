function [] = kNNDistance(obj, k, isRandom, maxDistance, showPlot)
%kNNDistance calculates distances to k-th nearest neighbor
% 
% k = 2 calculates distance to next neighbor
% 
% Input:
%   k:              number of next neigbor to calculate distance (default: 2)
%   isRandom:       true if random positons should be used
%   maxDistance:    maximum allowed distance for nearest neighbor in nm (default = size of image)
%   showPlot:       show histogram of NN distribution (default = true)
% Output:
%   (1) k-value
%   (2) nearest neighbor distances of all points
%
%   Note: If clustering with DBSCAN was done before, the function
%   automatically calculates the NN-distances of all points, of points
%   within clusters and of the points outside of the clusters.
%   Additional output:
%   (3) k-nearest neighbor distances of points within clusters
%   (4) k-nearest neighbor distances of points outside of clusters
%
%   Matlab 2014b or newer and Statistics and Machine Learning Toolbox
%
% Jan Neumann, 27.06.2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% init
isDBSCAN = false; % remembers if DBSCAN data already exist
switch nargin
    case 1
        k = 2;
        isRandom = 0;
        if isRandom == true
            maxDistance = max(obj.randomPhysicalDimension);
        else
            maxDistance = max(obj.physicalDimension);
        end
        showPlot = 0;
    case 2
        isRandom = 0;
        if isRandom == true
            maxDistance = max(obj.randomPhysicalDimension);
        else
            maxDistance = max(obj.physicalDimension);
        end
        showPlot = 0;
    case 3
        if isRandom == true
            maxDistance = max(obj.randomPhysicalDimension);
        else
            maxDistance = max(obj.physicalDimension);
        end
        showPlot = 0;
    case 4
        showPlot = 0;
    case 5
        
    otherwise
        error('Wrong number of input arguments!')
end
if maxDistance <= 0
    error('maxDistance must be larger than zero!')
end
if isRandom == true
    positions = obj.randomTable;
    dataMat = obj.randomClusterStruct;
    treeData = obj.randomQueryData;
    if isfield(obj.randomClusterStruct, 'clusterDBSCAN')
        clusterAssignment = obj.randomClusterStruct(3).clusterDBSCAN;
        isDBSCAN = true;
    end
else
    positions = obj.positionTable;
    dataMat = obj.clusterStruct;
    treeData = obj.queryData;
    if isfield(obj.clusterStruct, 'clusterDBSCAN')
        clusterAssignment = obj.clusterStruct(3).clusterDBSCAN;
        isDBSCAN = true;
    end
end
if k == 1 || k > obj.NoOfPoints
    error('k must be greater than 1 and less than the maximum number of points within the data set!')
end
% initialize output fields of dataMat
dataMat(1).kNNDistance = [];
dataMat(2).kNNDistance = [];
if isDBSCAN == true
    dataMat(3).kNNDistance = [];
    dataMat(4).kNNDistance = [];
end
%% k-NN-calculation for all distances
switch obj.flagSearchTree
    case 'cTree' % for knn search, use Matlab's internal function
        if obj.dimension == 3
            [~, nnDistance] = knnsearch(positions(:, 1:3), positions(:, 1:3), 'K', k);
        else
            [~, nnDistance] = knnsearch(positions(:, 1:2), positions(:, 1:2), 'K', k);
        end
    case 'matlabTree'
        [~, nnDistance] = knnsearch(treeData.X, treeData.X, 'K', k);
end
% mean over NN distance
nnDistances = mean(nnDistance(:, 2:end), 2);
nnDistances = nnDistances(nnDistances(:, 1) <= maxDistance, :);
%% k-NN-calculation for points inside and outside of clusters
if isDBSCAN == true
    % for points of clusters the clusterAssignment equals 0
    tempCluster = clusterAssignment(:, 1) == 0; % calculate NN distance outside of clusters
    if obj.dimension == 3
        outerPoints = positions(tempCluster, 1:3);
        [~, outerNNDistance] = knnsearch(outerPoints(:, 1:3), outerPoints(:, 1:3), 'K', k);
    else
        outerPoints = positions(tempCluster, 1:2);
        [~, outerNNDistance] = knnsearch(outerPoints(:, 1:2), outerPoints(:, 1:2), 'K', k);
    end
    outerNNDistance = mean(outerNNDistance(:, 2:end), 2);
    outerNNDistance = outerNNDistance(outerNNDistance(:, 1) <= maxDistance, :);
    innerNNDistance = []; % calculate NN distance within clusters
    if any(clusterAssignment)
        for kk = 1:max(clusterAssignment)
            tempCluster = clusterAssignment(:, 1) == kk;
            if obj.dimension == 3
                innerPoints = positions(tempCluster, 1:3);
                [~, tempInnerNNDistance] = knnsearch(innerPoints(:, 1:3), innerPoints(:, 1:3), 'K', k);
            else
                innerPoints = positions(tempCluster, 1:2);
                [~, tempInnerNNDistance] = knnsearch(innerPoints(:, 1:2), innerPoints(:, 1:2), 'K', k);
            end
            currentInnerNNDistance = mean(tempInnerNNDistance(:, 2:end), 2);
            currentInnerNNDistance = currentInnerNNDistance(currentInnerNNDistance(:, 1) <= maxDistance, :);
            innerNNDistance = cat(1, innerNNDistance, currentInnerNNDistance);
        end
    end
end
%% visualization
if showPlot == true
    figure('Name', [num2str(k), '-Nearest Neighbor Distances: ',num2str(obj.sampleID)] );
    if isDBSCAN == true
            plotCounter = 3; % check if outer or inner k-NN distances exist
            if isempty(innerNNDistance)
                plotCounter = plotCounter - 1;
            end
            if isempty(outerNNDistance)
                plotCounter = plotCounter - 1;
            end
            subplot(1, plotCounter, 1) 
            histogram(nnDistances, ceil(max(nnDistances(:, 1))), 'Normalization', 'probability');
            grid on;
            title('all points');
            xlabel('distance [nm]');
            ylabel('normalized frequency');
            ax = gca;
            ax.XLim = [0 maxDistance];
            if ~isempty(innerNNDistance)
                subplot(1, plotCounter, 2)
                histogram(innerNNDistance, ceil(max(innerNNDistance(:, 1))), 'Normalization', 'probability');
                grid on;
                title('points inside of clusters');
                xlabel('distance [nm]');
                ylabel('normalized frequency');
                ax = gca;
                ax.XLim = [0 maxDistance];
            end
            if ~isempty(outerNNDistance)
                subplot(1, plotCounter, 3)
                histogram(outerNNDistance, ceil(max(outerNNDistance(:, 1))), 'Normalization', 'probability');
                grid on;
                title('points outside of clusters');
                xlabel('distance [nm]');
                ylabel('normalized frequency');
                ax = gca;
                ax.XLim = [0 maxDistance];
            end
    else
        histogram(nnDistances, ceil(max(nnDistances(:, 1))), 'Normalization', 'probability');
        grid on;
        title('all points');
        xlabel('distance [nm]');
        ylabel('normalized frequency');
        ax = gca;
        ax.XLim = [0 maxDistance];
    end
end
%% output
dataMat(1).kNNDistance = k;
dataMat(2).kNNDistance = nnDistances; % all points
if isDBSCAN == true
    dataMat(3).kNNDistance = innerNNDistance; % points inside of clusters
    dataMat(4).kNNDistance = outerNNDistance; % points outside of clusters
end
if isRandom == true
    obj.randomClusterStruct = dataMat;
else
    obj.clusterStruct = dataMat;
end
end