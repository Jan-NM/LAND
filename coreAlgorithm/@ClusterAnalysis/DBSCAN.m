function [] = DBSCAN(obj, radius, minPoints, isRandom, maxDiameter, showImage)
%DBSCAN for details see A Density-Based Algorithm for Discovering Clusters in Large
%   Spatial Databases with Noise. Ester at al. 1996
%
% The current version is implemented for 2D + 3D data sets.
%
% Input:
%   radius: radius by which a certain number of points (minPoints) must be,
%           to accept the point as part of a cluster
%   minPoints: minimum number of points within specified distance (radius)
%           to fullfill cluster criteria
%   isRandom: true = use random data; false = use experimental data
%   maxDiameter: maximum diameter of a cluster, clusters larger than
%           maxDiameter will be removed (default: maximum dimension of the image)
%   showImage: true = show scatterplot, noise points and filtered clusters 
%           are displayed in black, clusters are displayed with a random
%           colour. The black cross indicates the center of mass of the
%           clusters.
%
% Output
%   (1) total number of points that are part of a cluster
%   (2) obj.clusterStruct or obj.randomClusterStruct (depending on value of isRandom)
%   containing the clusters and thier properties
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
% Cluster areas are calculated using Matlab's boundary command.
% Diameter/Volume property assumes more or less circular clusters.
%
% Jan Neumann, 27.06.2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% default parameters
visMethod = 'gaussianBlur'; % alternative: 'histogramBinning'
visPixelsize = 10;
%% init
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
    treeData = obj.randomQueryData;
    flagSearchTree = obj.flagSearchTreeRandomData;
else
    positions = obj.positionTable;
    dataMat = obj.clusterStruct;
    nPoints = obj.NoOfPoints;
    treeData = obj.queryData;
    flagSearchTree = obj.flagSearchTree;
end
% initialize output fields of dataMat
dataMat(1).clusterDBSCAN = [];
dataMat(2).clusterDBSCAN = [];
dataMat(3).clusterDBSCAN = [];
%% DBSCAN algorithm
clusterID = 0; % number labeling of cluster
% initialize cluster matrix - all points are at the beginning
% unassigned []
clusterAssignment = zeros(nPoints, 1);
for ii = 1:nPoints
    clusterAssignment(ii, 1) = -2;
end
multiWaitbar('DBSCAN algorithm running...', 0);
% for waitbar
prevPercent = 0;
for ii = 1:nPoints
    if ~(clusterAssignment(ii, 1) == -2)
        % point is already assigned
        continue
    end
    clusterAssignment(ii, 1) = -1; % point is visited
    switch flagSearchTree
        case 'cTree'
            distanceAppend = kdtree_ball_query(treeData, positions(ii, 1:obj.dimension), radius);
        case 'matlabTree'
            distance = rangesearch(positions(ii, 1:obj.dimension), treeData.X, radius);
            distanceAppend = (find(~cellfun('isempty', distance))).';
    end
    if numel(distanceAppend) < minPoints + 1 % + 1 ignore current point
        clusterAssignment(ii, 1) = 0; % point is noise
    else
        clusterID = clusterID + 1;
        clusterAssignment(ii, 1) = clusterID;
        indice = 1;
        while indice < numel(distanceAppend)
            ll = distanceAppend(indice);
            indice = indice + 1;
            if clusterAssignment(ll, 1) == -2 % if point is not visited
                clusterAssignment(ll, 1) = -1;
                switch flagSearchTree
                    case 'cTree'
                        newAppend = kdtree_ball_query(treeData, positions(ll, 1:obj.dimension), radius);
                    case 'matlabTree'
                        newDistance = rangesearch(positions(ll, 1:obj.dimension), treeData.X, radius);
                        newAppend = (find(~cellfun('isempty', newDistance))).';
                end
                if numel(newAppend) >= minPoints + 1 % + 1
                    distanceAppend = cat(1, distanceAppend, newAppend); % do not implement union or unique - this may change the order of the array
                end
            end
            if clusterAssignment(ll, 1) == -1 % not member of a cluster
                clusterAssignment(ll, 1) = clusterID;
            end
        end
    end
    currentPercent = fix(100*sum(clusterAssignment~=-2)/nPoints);
    if currentPercent > prevPercent
        multiWaitbar( 'DBSCAN algorithm running...', 'Value', sum(clusterAssignment~=-2)/nPoints);
        prevPercent = currentPercent;
    end
end
%% calculate cluster properties
for kk = 1:max(clusterAssignment)
    dataMat(2).clusterDBSCAN(kk).PointsCluster = find(clusterAssignment(:, 1) == kk);
    tempCluster = [dataMat(2).clusterDBSCAN(kk).PointsCluster];
    % number of molecules per cluster
    dataMat(2).clusterDBSCAN(kk).Molecules = size(tempCluster, 1);
    % extract positions of clustered points
    x = double(positions(tempCluster, 1));
    y = double(positions(tempCluster, 2));
    % calculate center of mass
    dataMat(2).clusterDBSCAN(kk).xCoM = 1/length(tempCluster)*sum(positions(tempCluster, 1));
    dataMat(2).clusterDBSCAN(kk).yCoM = 1/length(tempCluster)*sum(positions(tempCluster, 2));
    if obj.dimension == 2
        [~, dataMat(2).clusterDBSCAN(kk).Area] = boundary(x, y);
         % calculate diameter of clusters (assume that clusters are circular)
        dataMat(2).clusterDBSCAN(kk).Diameter = 2.*sqrt(dataMat(2).clusterDBSCAN(kk).Area./pi); % in nm
        % calculate area in µm^2
        dataMat(2).clusterDBSCAN(kk).Area = dataMat(2).clusterDBSCAN(kk).Area*(1/1000).^2; % in µm^2
    else
        % extract z position
        z = double(positions(tempCluster, 3));
        % calculate center of mass in z direction
        dataMat(2).clusterDBSCAN(kk).zCoM = 1/length(tempCluster)*sum(positions(tempCluster, 3));
        [~, dataMat(2).clusterDBSCAN(kk).Volume] = boundary(x, y, z);
        % calculate diameter of clusters (assume that clusters are spherical)
        dataMat(2).clusterDBSCAN(kk).Diameter = nthroot(6.* dataMat(2).clusterDBSCAN(kk).Volume./pi , 3); % in nm
        % calculate volume in µm^3
        dataMat(2).clusterDBSCAN(kk).Volume = dataMat(2).clusterDBSCAN(kk).Volume*(1/1000).^3; % in µm^3
    end
end

% filter cluster according to diameter and calculate with remaining cluster additional properties
if ~isempty(dataMat(2).clusterDBSCAN)
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
    if obj.dimension == 2
        [~, distance] = knnsearch(cat(2, vertcat(dataMat(2).clusterDBSCAN.xCoM), vertcat(dataMat(2).clusterDBSCAN.yCoM)),...
            cat(2, vertcat(dataMat(2).clusterDBSCAN.xCoM), vertcat(dataMat(2).clusterDBSCAN.yCoM)), 'k', 2);
    else
        [~, distance] = knnsearch(cat(2, vertcat(dataMat(2).clusterDBSCAN.xCoM), vertcat(dataMat(2).clusterDBSCAN.yCoM), vertcat(dataMat(2).clusterDBSCAN.zCoM)),...
            cat(2, vertcat(dataMat(2).clusterDBSCAN.xCoM), vertcat(dataMat(2).clusterDBSCAN.yCoM), vertcat(dataMat(2).clusterDBSCAN.zCoM)), 'k', 2);
    end
    % calculate total number of molecules in cluster
    dataMat(1).clusterDBSCAN = sum([dataMat(2).clusterDBSCAN.Molecules]);
    for ii = 1:size(dataMat(2).clusterDBSCAN, 2)
        if ~isempty(distance) && ~isequal(distance, 0)
            dataMat(2).clusterDBSCAN(ii).NNCluster = distance(ii, 2);
        end
        % check for two points and if Area = zero
        if ~isempty(dataMat(2).clusterDBSCAN)
            if obj.dimension == 2
                dataMat(2).clusterDBSCAN(ii).ClusterDensity = dataMat(2).clusterDBSCAN(ii).Molecules./dataMat(2).clusterDBSCAN(ii).Area; % molecules per µm^2
            else
                dataMat(2).clusterDBSCAN(ii).ClusterDensity = dataMat(2).clusterDBSCAN(ii).Molecules./dataMat(2).clusterDBSCAN(ii).Volume; % molecules per µm^3
            end
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
    figure( 'Name', num2str(obj.sampleID) );
    multiWaitbar('creating DBSCAN image...', 0);
    % for 2D data overlay, create histogram binned image and overlay clustered points
    if obj.dimension == 2
        clusterImage = visModuleCluster(positions, visMethod, visPixelsize);
        subplot(1, 2, 1);
        imshow(clusterImage, 'Colormap', hot);
        hbar = colorbar;
        ylabel(hbar, 'intensity [arb. unit]')
        axis image
        xlabel(['x', num2str(visPixelsize), ' nm'])
        ylabel(['x', num2str(visPixelsize), ' nm'])
        title('Gaussian blurred image')
    end
    % for waitbar
    prevPercent = 0;
    counter = 1;
    for kk = 1:max(clusterAssignment)
        tempCluster = find(clusterAssignment(:, 1) == kk);
        % plot points inside a cluster - each cluster with different color
        if obj.dimension == 2
            subplot(1, 2, 2);
            hold on
            scatter(double((positions(tempCluster, 2)) - min(positions(:, 2)))./visPixelsize, double((positions(tempCluster, 1)) - min(positions(:, 1)))./visPixelsize, 100, '.');
            % plot center of mass
            scatter(double((dataMat(2).clusterDBSCAN(kk).yCoM) - min(positions(:, 2)))./visPixelsize, double((dataMat(2).clusterDBSCAN(kk).xCoM) - min(positions(:, 1)))./visPixelsize, 150, 'kx');
        else
            scatter3(double((positions(tempCluster, 2)) - min(positions(:, 2)))./visPixelsize, double((positions(tempCluster, 1)) - min(positions(:, 1)))./visPixelsize, double((positions(tempCluster, 3)) - min(positions(:, 3)))./visPixelsize, 100,'.');
            hold on
            % plot center of mass
            scatter3(double((dataMat(2).clusterDBSCAN(kk).yCoM) - min(positions(:, 2)))./visPixelsize, double((dataMat(2).clusterDBSCAN(kk).xCoM) - min(positions(:, 1)))./visPixelsize, double((dataMat(2).clusterDBSCAN(kk).zCoM) - min(positions(:, 3)))./visPixelsize, 150, 'kx');
        end
        % waitbar
        currentPercent = fix(100*counter/max(clusterAssignment));
        if currentPercent > prevPercent
            multiWaitbar( 'creating DBSCAN image...', 'Value', counter/max(clusterAssignment));
            prevPercent = currentPercent;
        end
        counter = counter + 1;
    end
    % add points outside of clusters
    tempCluster = find(clusterAssignment(:, 1) == 0);
    if obj.dimension == 2
        scatter(double((positions(tempCluster, 2)) - min(positions(:, 2)) )./visPixelsize, double((positions(tempCluster, 1)) - min(positions(:, 1)) )./visPixelsize, 50, 'k.');
    else
        scatter3(double((positions(tempCluster, 2)) - min(positions(:, 2)) )./visPixelsize, double((positions(tempCluster, 1)) - min(positions(:, 1)) )./visPixelsize, double((positions(tempCluster, 3)) - min(positions(:, 3)) )./visPixelsize, 50, 'k.');
    end
    title('Detected clusters and noise points (filtered points are assigned as noise)')
    set(gca, 'Ydir', 'reverse')
    axis image
    xlabel(['x', num2str(visPixelsize), ' nm'])
    ylabel(['x', num2str(visPixelsize), ' nm'])
    if obj.dimension == 3
        zlabel(['x', num2str(visPixelsize), ' nm'])
    end
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