function [] = DBSCANdistanceMatrix(obj, radius, minPoints, isRandom, maxDiameter, showImage)
%DBSCAN for details see A Density-Based Algorithm for Discovering Clusters in Large
%   Spatial Databases with Noise. Ester at al. 1996
%
% The current version is implemented for 2D data sets.
%
% Input:
%   radius: radius within a certian number of points (minPoints) must be
%           found to accept point as part of a cluster
%   minPoints: minimum number of points within specified distance (radius)
%           to fullfill cluster criteria
%   isRandom: true = use random data; false = use experimental data
%   maxDiameter: maximum diameter of a cluster, cluster larger than
%           maxDiamter will removed (default: maximum dimension of the image)
%   showImage: true = show image with noise points in black and coloured cluster
%       overlaid to histogram binned image of the data. White cross
%       indicates center of mass of the cluster
%
% Output
%   (1) number of points in cluster
%   (2) obj.clusterStruct or obj.randomclusterStruct (depending on value of isRandom)
%   containing the cluster and thier properties
%   (3) clusterAssignment - every point in the localization data set is
%   assigned to noise (value = 0) or to a cluster (removed cluster - see maxDiameter
%   will be set to 0 as well)
%
% Note:
% The current version removes clusters, which consists of two or less
% points. Cluster touching the border of an image will be not removed.
%
% Matlab version: 2014b or newer (because of boundary command).
%
% Cluster areas are calculated using Matlab's boundary and polyarea command,
% alternativly convhull could be used, but this does not allow for
% shrinikng for irregular shaped clusters.
% Diameter property assumes more or less circular clusters.
% Large localization tables will split in chunks. Next, a distance matrix
% is pre-computed and is saved as sparse matrix to save memory.
%
% Jan Neumann, 27.06.2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% default parameters
memoryThresholdParam = 0.8; % maximum of total memory that will be used 
visPixelsize = 10;
%% init
multiWaitbar('computing distance matrix...', 0);
switch nargin
    case 3
        isRandom = 0;
        if isRandom == true
            maxDiameter = max(obj.randomPhysicalDimension);
        else
            maxDiameter = max(obj.physicalDimension);
        end
    case 4
        if isRandom == true
            maxDiameter = max(obj.randomPhysicalDimension);
        else
            maxDiameter = max(obj.physicalDimension);
        end
        showImage = 0;
    case 5
        showImage = 0;
    case 6
        
    otherwise
        error('Wrong number of input arguments!')
end
if isRandom == true
    positions = obj.randomTable;
    dataMat = obj.randomClusterStruct;
    nPoints = obj.randomNoOfPoints;
else
    positions = obj.positionTable;
    dataMat = obj.clusterStruct;
    nPoints = obj.NoOfPoints;
end
%% pre-compute distance matrix
% check available memory
spaceOnRam = memory;
% limit memory consumption to memoryThresholdParam times total available memory
memoryThreshold = memoryThresholdParam*spaceOnRam.MaxPossibleArrayBytes;
maxNoOfDistances = ceil(memoryThreshold / 8 / nPoints); % /8 for data type double
if maxNoOfDistances > nPoints
    maxNoOfDistances = nPoints;
end
% split Orte matrix
chunkSize = diff(fix(linspace(0, nPoints, ceil(nPoints/maxNoOfDistances) + 1))); % + 1 should stay
splittedPositions = mat2cell(positions, chunkSize, size(positions, 2));
offset = 0;
ivec = cell(1, size(splittedPositions, 1)); % pre-allocate column and row vectors for sparse matrix
jvec = cell(1, size(splittedPositions, 1));
% for waitbar
prevPercent = 0;
counter = 1;
for ii = 1:size(splittedPositions, 1)
    distances = pdist2(splittedPositions{ii}(:, 1:2), positions(:, 1:2));
    [row, column] = find(distances <= radius); % keep distances which are smaller than radius
    % also includes distance to each point itself (which is zero) - will be
    % considered in DBSCAN algorithm (minPoints + 1)
    ivec{ii} = (row + offset).';
    jvec{ii} = column.';
    offset = offset + size(splittedPositions{ii}, 1);
    % waitbar
    currentPercent = fix(100*counter/size(splittedPositions, 1));
    if currentPercent > prevPercent
        multiWaitbar( 'computing distance matrix...', 'Value', counter/size(splittedPositions, 1));
        prevPercent = currentPercent;
    end
    counter = counter + 1;
end
% create sparse matrix
jvec = (cell2mat(jvec)).';
ivec = (cell2mat(ivec)).';
xvec = true(length(ivec), 1);
distances = sparse(ivec, jvec, xvec, nPoints, nPoints);
multiWaitbar( 'computing distance matrix...', 'Close' );
%% DBSCAN algorithm
clusterID = 0; % number labeling of cluster
% initialize Cluster matrix - all points are at beginning
% unassigned []
clusterAssignment = zeros(nPoints, 1);
for ii = 1:nPoints
    clusterAssignment(ii, 1) = -2;
end
multiWaitbar('DBSCAN algorithm running...', 'Busy');
for ii = 1:nPoints
    if ~(clusterAssignment(ii, 1) == -2)
        % point is already assigned
        continue
    end
    clusterAssignment(ii, 1) = -1; % point is visited
    if sum(distances(:, ii), 1) < minPoints + 1 % + 1 ignore current point
        clusterAssignment(ii, 1) = 0; % point is noise
    else
        clusterID = clusterID + 1;
        clusterAssignment(ii, 1) = clusterID;
        indice = 1;
        distanceAppend = find(distances(:, ii));
        while indice < length(distanceAppend)
            ll = distanceAppend(indice);
            indice = indice + 1;
            if clusterAssignment(ll, 1) == -2 % if point is not visited
                clusterAssignment(ll, 1) = -1;  
                if sum(distances(:, ll), 1) >= minPoints + 1 % + 1
                    newAppend = find(distances(:, ll));
                    distanceAppend = cat(1, distanceAppend, newAppend); % do not implement union or unique - this may change the order of the array
                end
            end
            if clusterAssignment(ll, 1) == -1 % not member of cluster
                clusterAssignment(ll, 1) = clusterID;
            end
        end
    end
end
%% calculate cluster properties
for kk = 1:max(clusterAssignment)
    dataMat(2).clusterDBSCAN(kk).PointsCluster = find(clusterAssignment(:, 1) == kk);
    tempCluster = [dataMat(2).clusterDBSCAN(kk).PointsCluster];
    % extract position of cluster points
    x = double(positions(tempCluster, 1));
    y = double(positions(tempCluster, 2));
    [~, dataMat(2).clusterDBSCAN(kk).Area] = boundary(x,y);
    dataMat(2).clusterDBSCAN(kk).Molecules = size(tempCluster, 1);
    dataMat(2).clusterDBSCAN(kk).Diameter = 2.*sqrt(dataMat(2).clusterDBSCAN(kk).Area./pi); % in nm
    dataMat(2).clusterDBSCAN(kk).Area = dataMat(2).clusterDBSCAN(kk).Area*(1/1000).^2; % µm^2
    % calculate center of mass
    dataMat(2).clusterDBSCAN(kk).xCoM = 1/length(tempCluster)*sum(positions(tempCluster, 1));
    dataMat(2).clusterDBSCAN(kk).yCoM = 1/length(tempCluster)*sum(positions(tempCluster, 2));
end
% filter cluster according to diamter and calculate with remaining cluster additional properties

%implement area
if ~isempty(dataMat)
    % filter data
    deleteInd = [dataMat(2).clusterDBSCAN.Diameter] > maxDiameter;
    for ii = 1:size(dataMat(2).clusterDBSCAN, 2)
        if deleteInd(ii) == 1
            tempCluster = clusterAssignment(:, 1) == ii;
            clusterAssignment(tempCluster, 1) = 0;
        end
    end
    dataMat(2).clusterDBSCAN([dataMat(2).clusterDBSCAN.Diameter] > maxDiameter) = [];
    % filter clusters which consists of two or less points
    deleteInd = [dataMat(2).clusterDBSCAN.Diameter] == 0;
    for ii = 1:size(dataMat(2).clusterDBSCAN, 2)
        if deleteInd(ii) == 1
            tempCluster = clusterAssignment(:, 1) == ii;
            clusterAssignment(tempCluster, 1) = 0;
        end
    end
    dataMat(2).clusterDBSCAN([dataMat(2).clusterDBSCAN.Diameter] == 0) = [];
    % shift cluster number in clusterAssignment, so starting from one in
    % ascending order
    counter = 0;
    for ii = 1:max(clusterAssignment)
        tempVal = clusterAssignment(:, 1) == ii;
        if any(tempVal)
            clusterAssignment(tempVal, 1) = counter + 1;
            counter = counter + 1;   
        end
    end
    % calculate distance to next neighboring cluster
    [~, distance] = knnsearch(cat(2,vertcat(dataMat(2).clusterDBSCAN.xCoM),vertcat(dataMat(2).clusterDBSCAN.yCoM)),...
        cat(2,vertcat(dataMat(2).clusterDBSCAN.xCoM),vertcat(dataMat(2).clusterDBSCAN.yCoM)), 'k', 2);
    % calculate total number of molecules in cluster
    dataMat(1).clusterDBSCAN = sum([dataMat(2).clusterDBSCAN.Molecules]);
    for ii = 1:size(dataMat(2).clusterDBSCAN, 2)
        if ~isempty(distance) && ~isequal(distance, 0)
            dataMat(2).clusterDBSCAN(ii).NNCluster = distance(ii, 2);
        end
        % check for two points and if Area = zero
        if ~isempty(dataMat(2).clusterDBSCAN)
            dataMat(2).clusterDBSCAN(ii).ClusterDensity = dataMat(2).clusterDBSCAN(ii).Molecules./dataMat(2).clusterDBSCAN(ii).Area;
        end
    end
else
    % no points in cluster
    dataMat(1).clusterDBSCAN = 0;
end
% save cluster assignment for distance calculation (to calculate distances from points inside and outside of clusters)
dataMat(3).clusterDBSCAN = uint32(clusterAssignment);
multiWaitbar( 'DBSCAN algorithm running...', 'Close' );
%% visualization
% plots all clusters, filtered clusters are considered
if showImage == true && ~isempty(clusterAssignment)
    % create histogram binned image and overlay clustered points
    clusterImage = visModuleCluster(positions, 'histogramBinning', visPixelsize);
    figure( 'Name', num2str(obj.sampleID) );
    imagesc(clusterImage);
    hold on
    multiWaitbar('creating DBSCAN image...', 0);
    % for waitbar
    prevPercent = 0;
    counter = 1;
    for kk = 1:max(clusterAssignment)
        tempCluster = find(clusterAssignment(:, 1) == kk);
        % plot points inside a cluster - each cluster with different color
        scatter(double((positions(tempCluster, 2)) - min(positions(:, 2)))./visPixelsize + 0.5, double((positions(tempCluster, 1)) - min(positions(:, 2)))./visPixelsize + 0.5, 100, '.');
        hold on
        % plot center of mass
        scatter(double((dataMat(2).clusterDBSCAN(kk).yCoM) - min(positions(:, 2)))./visPixelsize + 0.5, double((dataMat(2).clusterDBSCAN(kk).xCoM) - min(positions(:, 2)))./visPixelsize + 0.5, 100, 'wx');
        % waitbar
        currentPercent = fix(100*counter/max(clusterAssignment));
        if currentPercent > prevPercent
            multiWaitbar( 'creating DBSCAN image...', 'Value', counter/max(clusterAssignment));
            prevPercent = currentPercent;
        end
        counter = counter + 1;
    end
    % add points outside of cluster
    tempCluster = find(clusterAssignment(:, 1) == 0);
    scatter(double((positions(tempCluster, 2)) - min(positions(:, 2)) )./visPixelsize + 0.5, double((positions(tempCluster, 1)) - min(positions(:, 2)) )./visPixelsize + 0.5, 50, 'k.');
    title('Histogram binned image with points overlaid - points outside of clusters and filtered points are black')
    axis equal
    axis off
    hold off
end
%% output
if isRandom == true
    obj.randomClusterStruct = dataMat;
else
    obj.clusterStruct = dataMat;
end
multiWaitbar('CLOSEALL');
end