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
%% distance calculation
% calculate distance using cell grid based algorithm
% check size of localization table and split localization table
if isRandom == true
    dataSize = min(obj.randomPhysicalDimension(1), obj.randomPhysicalDimension(2));
else
    dataSize = min(obj.physicalDimension(1), obj.physicalDimension(2));
end
if dataSize < 3*maxDistance
    disp('Image size is smaller then 3 times the cut-off distance! Try to decrease cut-off distance or increase image region.');
    return
end
% init distance matrix
distances = [];

mapSize = ceil(max(positions(:, 1:2)));
xStepSize = ceil(mapSize(1)/maxDistance);
yStepSize = ceil(mapSize(2)/maxDistance);
cellVerification = zeros(xStepSize, yStepSize, 9, 'uint16'); % for verification which pair of cells had been already visited
multiWaitbar('calculating distances...', 0);
% for waitbar
prevPercent = 0;
counter = 1;
for stripe = 2:xStepSize-1
    for line = 2:yStepSize-1
        % compute boundaries of current cell grid
        upperStripe = max(stripe-1, 1); % take only positive values
        lowerStripe = min(xStepSize, stripe+1); % prevent exceeding image boundaries
        leftStripe = max(line-1, 1); % take only positive values
        rightStripe = min(yStepSize, line+1);% prevent exceeding image boundaries
        % get signals from current position
        currentPositions = positions(positions(:, 1) >= upperStripe*maxDistance & positions(:, 1) < (lowerStripe-1)*maxDistance...
            & positions(:, 2) >= leftStripe*maxDistance & positions(:, 2) < (rightStripe-1)*maxDistance, :);
        % do something with first stripe...check first if already calculated
        currentDistance = single(pdist(currentPositions(:, 1:2)))'; % distance within cell
        distances = [distances; currentDistance];
        % set verifiction of current stripe
        horizontalPosition = line;
        verticalPosition = stripe;
        currentIdx = strcat(num2str(verticalPosition),num2str(horizontalPosition));
        cellVerification(stripe, line, 5) = str2double(currentIdx);
        % position 1
        % do somthing with position 1 upper left
        if ~(cellVerification(upperStripe, leftStripe, :) == str2double(currentIdx))
            pointerPositions = positions(positions(:, 1) >= (upperStripe-1)*maxDistance & positions(:, 1) < (lowerStripe-2)*maxDistance...
                & positions(:, 2) >= (leftStripe-1)*maxDistance & positions(:, 2) < (rightStripe-2)*maxDistance, :);
            % do something with first stripe...check first if already calculated
            currentDistance = single(pdist2(currentPositions(:, 1:2), pointerPositions(:, 1:2)));
            distances = [distances; reshape(currentDistance,[],1)];
            % set verifiction of current stripe
            horizontalPosition = leftStripe;
            verticalPosition = upperStripe;
            upperIdx = strcat(num2str(verticalPosition),num2str(horizontalPosition));
            cellVerification(stripe, line, 1) = str2double(upperIdx);
        end
        % position 2
        % do somthing with position 2 upper
        if ~(cellVerification(upperStripe, line, :) == str2double(currentIdx))
            pointerPositions = positions(positions(:, 1) >= (upperStripe-1)*maxDistance & positions(:, 1) < (lowerStripe-2)*maxDistance...
                & positions(:, 2) >= (leftStripe)*maxDistance & positions(:, 2) < (rightStripe-1)*maxDistance, :);
            % do something with first stripe...check first if already calculated
            currentDistance = single(pdist2(currentPositions(:, 1:2), pointerPositions(:, 1:2)))';
            distances = [distances; reshape(currentDistance,[],1)];
            % set verifiction of current stripe
            horizontalPosition = line;
            verticalPosition = upperStripe;
            upperIdx = strcat(num2str(verticalPosition),num2str(horizontalPosition));
            cellVerification(stripe, line, 2) = str2double(upperIdx);
        end
        % position 3
        % do somthing with position 3 upper right
        if ~(cellVerification(upperStripe, rightStripe, :) == str2double(currentIdx))
            pointerPositions = positions(positions(:, 1) >= (upperStripe-1)*maxDistance & positions(:, 1) < (lowerStripe-2)*maxDistance...
                & positions(:, 2) >= (leftStripe+1)*maxDistance & positions(:, 2) < (rightStripe)*maxDistance, :);
            % do something with first stripe...check first if already calculated
            currentDistance = single(pdist2(currentPositions(:, 1:2), pointerPositions(:, 1:2)))';
            distances = [distances; reshape(currentDistance,[],1)];
            % set verifiction of current stripe
            horizontalPosition = rightStripe;
            verticalPosition = upperStripe;
            upperIdx = strcat(num2str(verticalPosition),num2str(horizontalPosition));
            cellVerification(stripe, line, 3) = str2double(upperIdx);
        end
        % position 4
        % do somthing with position 4 left
        if ~(cellVerification(stripe, leftStripe, :) == str2double(currentIdx))
            pointerPositions = positions(positions(:, 1) >= (upperStripe)*maxDistance & positions(:, 1) < (lowerStripe-1)*maxDistance...
                & positions(:, 2) >= (leftStripe-1)*maxDistance & positions(:, 2) < (rightStripe-2)*maxDistance, :);
            % do something with first stripe...check first if already calculated
            currentDistance = single(pdist2(currentPositions(:, 1:2), pointerPositions(:, 1:2)))';
            distances = [distances; reshape(currentDistance,[],1)];
            % set verifiction of current stripe
            horizontalPosition = leftStripe;
            verticalPosition = stripe;
            upperIdx = strcat(num2str(verticalPosition),num2str(horizontalPosition));
            cellVerification(stripe, line, 4) = str2double(upperIdx);
        end
        % position 6
        % do somthing with position 6 right
        if ~(cellVerification(stripe, rightStripe, :) == str2double(currentIdx))
            pointerPositions = positions(positions(:, 1) >= (upperStripe)*maxDistance & positions(:, 1) < (lowerStripe-1)*maxDistance...
                & positions(:, 2) >= (leftStripe+1)*maxDistance & positions(:, 2) < (rightStripe)*maxDistance, :);
            % do something with first stripe...check first if already calculated
            currentDistance = single(pdist2(currentPositions(:, 1:2), pointerPositions(:, 1:2)))';
            distances = [distances; reshape(currentDistance,[],1)];
            % set verifiction of current stripe
            horizontalPosition = rightStripe;
            verticalPosition = stripe;
            upperIdx = strcat(num2str(verticalPosition),num2str(horizontalPosition));
            cellVerification(stripe, line, 6) = str2double(upperIdx);
        end
        % position 7
        % do somthing with position 7 lower left
        if ~(cellVerification(lowerStripe, leftStripe, :) == str2double(currentIdx))
            pointerPositions = positions(positions(:, 1) >= (upperStripe+1)*maxDistance & positions(:, 1) < (lowerStripe)*maxDistance...
                & positions(:, 2) >= (leftStripe-1)*maxDistance & positions(:, 2) < (rightStripe-2)*maxDistance, :);
            % do something with first stripe...check first if already calculated
            currentDistance = single(pdist2(currentPositions(:, 1:2), pointerPositions(:, 1:2)))';
            distances = [distances; reshape(currentDistance,[],1)];
            % set verifiction of current stripe
            horizontalPosition = leftStripe;
            verticalPosition = lowerStripe;
            upperIdx = strcat(num2str(verticalPosition),num2str(horizontalPosition));
            cellVerification(stripe, line, 7) = str2double(upperIdx);
        end
        % position 8
        % do somthing with position 8 lower
        if ~(cellVerification(lowerStripe, line, :) == str2double(currentIdx))
            pointerPositions = positions(positions(:, 1) >= (upperStripe+1)*maxDistance & positions(:, 1) < (lowerStripe)*maxDistance...
                & positions(:, 2) >= (leftStripe)*maxDistance & positions(:, 2) < (rightStripe-1)*maxDistance, :);
            % do something with first stripe...check first if already calculated
            currentDistance = single(pdist2(currentPositions(:, 1:2), pointerPositions(:, 1:2)))';
            distances = [distances; reshape(currentDistance,[],1)];
            % set verifiction of current stripe
            horizontalPosition = line;
            verticalPosition = lowerStripe;
            upperIdx = strcat(num2str(verticalPosition),num2str(horizontalPosition));
            cellVerification(stripe, line, 8) = str2double(upperIdx);
        end
        % position 9
        % do somthing with position 9 lower right
        if ~(cellVerification(lowerStripe, rightStripe, :) == str2double(currentIdx))
            pointerPositions = positions(positions(:, 1) >= (upperStripe+1)*maxDistance & positions(:, 1) < (lowerStripe)*maxDistance...
                & positions(:, 2) >= (leftStripe+1)*maxDistance & positions(:, 2) < (rightStripe)*maxDistance, :);
            % do something with first stripe...check first if already calculated
            currentDistance = single(pdist2(currentPositions(:, 1:2), pointerPositions(:, 1:2)))';
            distances = [distances; reshape(currentDistance,[],1)];
            % set verifiction of current stripe
            horizontalPosition = rightStripe;
            verticalPosition = lowerStripe;
            upperIdx = strcat(num2str(verticalPosition),num2str(horizontalPosition));
            cellVerification(stripe, line, 9) = str2double(upperIdx);
        end
    end
    % waitbar
    currentPercent = fix(100*counter/(xStepSize - 2));
    if currentPercent > prevPercent
        multiWaitbar( 'calculating distances...', 'Value', counter/(xStepSize - 2));
        prevPercent = currentPercent;
    end
    counter = counter + 1;
end
% removes distances larger then maxDistance
currentIdx = (distances <= maxDistance & distances > 0);
distances = distances(currentIdx);
%% distance calculation for points inside and outside of clusters
if isDBSCAN == true
    % for outside points clusterAssignment is 0
    tempCluster = clusterAssignment(:, 1) == 0;
    outerPoints = positions(tempCluster, 1:2);
    if isRandom == true
        [~, idx] = datasample(positions, numberOfSample);
        outerPoints = positions(idx, 1:2);
    end
    % calculate distance using cell grid based algorithm
    % check size of localization table and split localization table
    if min(obj.physicalDimension(1), obj.physicalDimension(2)) < 3*maxDistance
        disp('Image size is smaller then 3 times the cut-off distance! Try to decrease cut-off distance or increase image region.');
        return
    end
    % init distance matrix
    outerDistance = [];
    
    mapSize = ceil(max(outerPoints(:, 1:2)));
    xStepSize = ceil(mapSize(1)/maxDistance);
    yStepSize = ceil(mapSize(2)/maxDistance);
    cellVerification = zeros(xStepSize, yStepSize, 9, 'uint16'); % for verification which pair of cells had been already visited
    multiWaitbar('calculating distances outside of cluster...', 0);
    % for waitbar
    prevPercent = 0;
    counter = 1;
    for stripe = 2:xStepSize-1
        for line = 2:yStepSize-1
            % compute boundaries of current cell grid
            upperStripe = max(stripe-1, 1); % take only positive values
            lowerStripe = min(xStepSize, stripe+1); % prevent exceeding image boundaries
            leftStripe = max(line-1, 1); % take only positive values
            rightStripe = min(yStepSize, line+1);% prevent exceeding image boundaries
            % get signals from current position
            currentouterPoints = outerPoints(outerPoints(:, 1) >= upperStripe*maxDistance & outerPoints(:, 1) < (lowerStripe-1)*maxDistance...
                & outerPoints(:, 2) >= leftStripe*maxDistance & outerPoints(:, 2) < (rightStripe-1)*maxDistance, :);
            % do something with first stripe...check first if already calculated
            currentDistance = single(pdist(currentouterPoints(:, 1:2)))'; % distance within cell
            outerDistance = [outerDistance; currentDistance];
            % set verifiction of current stripe
            horizontalPosition = line;
            verticalPosition = stripe;
            currentIdx = strcat(num2str(verticalPosition),num2str(horizontalPosition));
            cellVerification(stripe, line, 5) = str2double(currentIdx);
            % position 1
            % do somthing with position 1 upper left
            if ~(cellVerification(upperStripe, leftStripe, :) == str2double(currentIdx))
                pointerouterPoints = outerPoints(outerPoints(:, 1) >= (upperStripe-1)*maxDistance & outerPoints(:, 1) < (lowerStripe-2)*maxDistance...
                    & outerPoints(:, 2) >= (leftStripe-1)*maxDistance & outerPoints(:, 2) < (rightStripe-2)*maxDistance, :);
                % do something with first stripe...check first if already calculated
                currentDistance = single(pdist2(currentouterPoints(:, 1:2), pointerouterPoints(:, 1:2)));
                outerDistance = [outerDistance; reshape(currentDistance,[],1)];
                % set verifiction of current stripe
                horizontalPosition = leftStripe;
                verticalPosition = upperStripe;
                upperIdx = strcat(num2str(verticalPosition),num2str(horizontalPosition));
                cellVerification(stripe, line, 1) = str2double(upperIdx);
            end
            % position 2
            % do somthing with position 2 upper
            if ~(cellVerification(upperStripe, line, :) == str2double(currentIdx))
                pointerouterPoints = outerPoints(outerPoints(:, 1) >= (upperStripe-1)*maxDistance & outerPoints(:, 1) < (lowerStripe-2)*maxDistance...
                    & outerPoints(:, 2) >= (leftStripe)*maxDistance & outerPoints(:, 2) < (rightStripe-1)*maxDistance, :);
                % do something with first stripe...check first if already calculated
                currentDistance = single(pdist2(currentouterPoints(:, 1:2), pointerouterPoints(:, 1:2)))';
                outerDistance = [outerDistance; reshape(currentDistance,[],1)];
                % set verifiction of current stripe
                horizontalPosition = line;
                verticalPosition = upperStripe;
                upperIdx = strcat(num2str(verticalPosition),num2str(horizontalPosition));
                cellVerification(stripe, line, 2) = str2double(upperIdx);
            end
            % position 3
            % do somthing with position 3 upper right
            if ~(cellVerification(upperStripe, rightStripe, :) == str2double(currentIdx))
                pointerouterPoints = outerPoints(outerPoints(:, 1) >= (upperStripe-1)*maxDistance & outerPoints(:, 1) < (lowerStripe-2)*maxDistance...
                    & outerPoints(:, 2) >= (leftStripe+1)*maxDistance & outerPoints(:, 2) < (rightStripe)*maxDistance, :);
                % do something with first stripe...check first if already calculated
                currentDistance = single(pdist2(currentouterPoints(:, 1:2), pointerouterPoints(:, 1:2)))';
                outerDistance = [outerDistance; reshape(currentDistance,[],1)];
                % set verifiction of current stripe
                horizontalPosition = rightStripe;
                verticalPosition = upperStripe;
                upperIdx = strcat(num2str(verticalPosition),num2str(horizontalPosition));
                cellVerification(stripe, line, 3) = str2double(upperIdx);
            end
            % position 4
            % do somthing with position 4 left
            if ~(cellVerification(stripe, leftStripe, :) == str2double(currentIdx))
                pointerouterPoints = outerPoints(outerPoints(:, 1) >= (upperStripe)*maxDistance & outerPoints(:, 1) < (lowerStripe-1)*maxDistance...
                    & outerPoints(:, 2) >= (leftStripe-1)*maxDistance & outerPoints(:, 2) < (rightStripe-2)*maxDistance, :);
                % do something with first stripe...check first if already calculated
                currentDistance = single(pdist2(currentouterPoints(:, 1:2), pointerouterPoints(:, 1:2)))';
                outerDistance = [outerDistance; reshape(currentDistance,[],1)];
                % set verifiction of current stripe
                horizontalPosition = leftStripe;
                verticalPosition = stripe;
                upperIdx = strcat(num2str(verticalPosition),num2str(horizontalPosition));
                cellVerification(stripe, line, 4) = str2double(upperIdx);
            end
            % position 6
            % do somthing with position 6 right
            if ~(cellVerification(stripe, rightStripe, :) == str2double(currentIdx))
                pointerouterPoints = outerPoints(outerPoints(:, 1) >= (upperStripe)*maxDistance & outerPoints(:, 1) < (lowerStripe-1)*maxDistance...
                    & outerPoints(:, 2) >= (leftStripe+1)*maxDistance & outerPoints(:, 2) < (rightStripe)*maxDistance, :);
                % do something with first stripe...check first if already calculated
                currentDistance = single(pdist2(currentouterPoints(:, 1:2), pointerouterPoints(:, 1:2)))';
                outerDistance = [outerDistance; reshape(currentDistance,[],1)];
                % set verifiction of current stripe
                horizontalPosition = rightStripe;
                verticalPosition = stripe;
                upperIdx = strcat(num2str(verticalPosition),num2str(horizontalPosition));
                cellVerification(stripe, line, 6) = str2double(upperIdx);
            end
            % position 7
            % do somthing with position 7 lower left
            if ~(cellVerification(lowerStripe, leftStripe, :) == str2double(currentIdx))
                pointerouterPoints = outerPoints(outerPoints(:, 1) >= (upperStripe+1)*maxDistance & outerPoints(:, 1) < (lowerStripe)*maxDistance...
                    & outerPoints(:, 2) >= (leftStripe-1)*maxDistance & outerPoints(:, 2) < (rightStripe-2)*maxDistance, :);
                % do something with first stripe...check first if already calculated
                currentDistance = single(pdist2(currentouterPoints(:, 1:2), pointerouterPoints(:, 1:2)))';
                outerDistance = [outerDistance; reshape(currentDistance,[],1)];
                % set verifiction of current stripe
                horizontalPosition = leftStripe;
                verticalPosition = lowerStripe;
                upperIdx = strcat(num2str(verticalPosition),num2str(horizontalPosition));
                cellVerification(stripe, line, 7) = str2double(upperIdx);
            end
            % position 8
            % do somthing with position 8 lower
            if ~(cellVerification(lowerStripe, line, :) == str2double(currentIdx))
                pointerouterPoints = outerPoints(outerPoints(:, 1) >= (upperStripe+1)*maxDistance & outerPoints(:, 1) < (lowerStripe)*maxDistance...
                    & outerPoints(:, 2) >= (leftStripe)*maxDistance & outerPoints(:, 2) < (rightStripe-1)*maxDistance, :);
                % do something with first stripe...check first if already calculated
                currentDistance = single(pdist2(currentouterPoints(:, 1:2), pointerouterPoints(:, 1:2)))';
                outerDistance = [outerDistance; reshape(currentDistance,[],1)];
                % set verifiction of current stripe
                horizontalPosition = line;
                verticalPosition = lowerStripe;
                upperIdx = strcat(num2str(verticalPosition),num2str(horizontalPosition));
                cellVerification(stripe, line, 8) = str2double(upperIdx);
            end
            % position 9
            % do somthing with position 9 lower right
            if ~(cellVerification(lowerStripe, rightStripe, :) == str2double(currentIdx))
                pointerouterPoints = outerPoints(outerPoints(:, 1) >= (upperStripe+1)*maxDistance & outerPoints(:, 1) < (lowerStripe)*maxDistance...
                    & outerPoints(:, 2) >= (leftStripe+1)*maxDistance & outerPoints(:, 2) < (rightStripe)*maxDistance, :);
                % do something with first stripe...check first if already calculated
                currentDistance = single(pdist2(currentouterPoints(:, 1:2), pointerouterPoints(:, 1:2)))';
                outerDistance = [outerDistance; reshape(currentDistance,[],1)];
                % set verifiction of current stripe
                horizontalPosition = rightStripe;
                verticalPosition = lowerStripe;
                upperIdx = strcat(num2str(verticalPosition),num2str(horizontalPosition));
                cellVerification(stripe, line, 9) = str2double(upperIdx);
            end
        end
        % waitbar
        currentPercent = fix(100*counter/(xStepSize - 2));
        if currentPercent > prevPercent
            multiWaitbar( 'calculating distances outside of cluster...', 'Value', counter/(xStepSize - 2));
            prevPercent = currentPercent;
        end
        counter = counter + 1;
    end
    currentIdx = (outerDistance <= maxDistance & outerDistance > 0);
    outerDistance = outerDistance(currentIdx);
    innerDistance = []; % assumes that size of the cluster and number of points within a cluster is small enough to be able to use pdist function
    if any(clusterAssignment)
        for kk = 1:max(clusterAssignment)
            tempCluster = clusterAssignment(:, 1) == kk;
            innerPoints = positions(tempCluster, 1:2);
            % calculates pairwise distances using pdist function
            tempInnerDistance = single(pdist(innerPoints(:, 1:2)))';
            % removes distances larger then maxDistance
            currentIdx = (tempInnerDistance <= maxDistance & tempInnerDistance > 0);
            innerDistance = cat(1, innerDistance, tempInnerDistance(currentIdx));
        end
    end
end
%% visualization
if showPlot == true
    figure( 'Name', num2str(obj.sampleID) );
    if isDBSCAN == true
        subplot(1, 3, 1)
        h = histogram(distances, ceil(max(distances(:, 1))), 'Normalization', 'probability');
        grid on;
        title('Distance Analysis - all points');
        xlabel('distance [nm]');
        ylabel('normalized frequency');
        ax = gca;
        ax.XLim = [0 maxDistance];
        subplot(1, 3, 2)
        h = histogram(innerDistance, ceil(max(innerDistance(:, 1))), 'Normalization', 'probability');
        grid on;
        title('Distance Analysis - points inside of clusters');
        xlabel('distance [nm]');
        ylabel('normalized frequency');
        ax = gca;
        ax.XLim = [0 maxDistance];
        subplot(1, 3, 3)
        h = histogram(outerDistance, ceil(max(outerDistance(:, 1))), 'Normalization', 'probability');
        grid on;
        title('Distance Analysis - points outside of clusters');
        xlabel('distance [nm]');
        ylabel('normalized frequency');
        ax = gca;
        ax.XLim = [0 maxDistance];
    else
        h = histogram(distances, maxDistance, 'Normalization', 'probability'); % data are normalized to all distances with 0 and maxDistance
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