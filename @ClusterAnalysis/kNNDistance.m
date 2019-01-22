function [] = kNNDistance(obj, k, isRandom, maxDistance, showPlot)
%kNNDistance calculates distance to k-th nearest neighbor
% 
% k = 2 calculates distance to next neighbor
% 
% Input:
%   k: number of next neigbor to calculate distance (default: 2)
%   isRandom: true if random positons should be used
%   maxDistance: maximum allowed distance for nearest neighbor in nm (default = size of image)
%   showPlot: show histogram of NN distribution (default = true)
% Output:
%   (1) k-value
%   (2) nearest neighbor distances of all points
%
%   Note: if clustering with DBSCAN was done before, the function
%   automatically calculates NN-distance of all the points, of points
%   within cluster and of the points outside of the cluster
%   Additional output:
%   (3) nearest neighbor distances of points within clusters
%   (2) nearest neighbor distances of points outside of clusters
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
    error('k must be larger then 1 and smaller then the maximum number of points within the data set!')
end
%% k-NN-calculation for all distances
[~, nnDistance] = knnsearch(treeData.X, treeData.X, 'K', k);
% mean over NN distance
nnDistances = mean(nnDistance(:, 2:end), 2);
nnDistances = nnDistances(nnDistances(:, 1)<= maxDistance, :);
%% k-NN-calculation for points inside and outside of clusters
if isDBSCAN == true
    % for outside points clusterAssignment is 0
    tempCluster = clusterAssignment(:, 1) == 0;
    outerPoints = positions(tempCluster, 1:2);
    [~, outerNNDistance] = knnsearch(outerPoints(:, 1:2), outerPoints(:, 1:2), 'K', k);
    % mean over outer NN distance
    outerNNDistance = mean(outerNNDistance(:, 2:end), 2);
    outerNNDistance = outerNNDistance(outerNNDistance(:, 1)<= maxDistance, :);
    innerNNDistance = [];
    if any(clusterAssignment)
        for kk = 1:max(clusterAssignment)
            tempCluster = clusterAssignment(:, 1) == kk;
            innerPoints = positions(tempCluster, 1:2);
            % calculate NN distance within cluster
            [~, tempInnerNNDistance] = knnsearch(innerPoints(:, 1:2), innerPoints(:, 1:2), 'K', k);
            currentInnerNNDistance = mean(tempInnerNNDistance(:, 2:end), 2);
            currentInnerNNDistance = currentInnerNNDistance(currentInnerNNDistance(:, 1) <= maxDistance, :);
            innerNNDistance = cat(1, innerNNDistance, currentInnerNNDistance);
        end
    end
end
%% visualization
if showPlot == true
    figure( 'Name', num2str(obj.sampleID) );
    if isDBSCAN == true
        subplot(1, 3, 1)
        h = histogram(nnDistances, ceil(max(nnDistances(:, 1))), 'Normalization', 'probability');
        grid on;
        title('Nearest Neighbor Distance - all points');
        xlabel('distance [nm]');
        ylabel('normalized frequency');
        ax = gca;
        ax.XLim = [0 maxDistance];
        subplot(1, 3, 2)
        h = histogram(innerNNDistance, ceil(max(innerNNDistance(:, 1))), 'Normalization', 'probability');
        grid on;
        title('Nearest Neighbor Distance - points inside of clusters');
        xlabel('distance [nm]');
        ylabel('normalized frequency');
        ax = gca;
        ax.XLim = [0 maxDistance];
        subplot(1, 3, 3)
        h = histogram(outerNNDistance, ceil(max(outerNNDistance(:, 1))), 'Normalization', 'probability');
        grid on;
        title('Nearest Neighbor Distance - points outside of clusters');
        xlabel('distance [nm]');
        ylabel('normalized frequency');
        ax = gca;
        ax.XLim = [0 maxDistance];
    else
        h = histogram(nnDistances, ceil(max(nnDistances(:, 1))), 'Normalization', 'probability');
        grid on;
        title('Nearest Neighbor Distance');
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