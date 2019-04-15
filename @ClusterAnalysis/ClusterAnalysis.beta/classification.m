%% classification of radial density function, requires c++ kdtree
%% evaluation parameters
maxRadius = 1000; % in nm, maximum radius to evaluate radial density function
binSize = 5; % distance between rings
gridSize = 25; % in nm, maximum size of a bin - smaller values need much more computational time
nCluster = 7; % number of clusters that should be detected
%% init
multiWaitbar('computing radial density function...', 0);
% create kdtree
positions(:, 1) = Orte(:, 2);
positions(:, 2) = Orte(:, 3);
treeData = kdtree_build([positions(:, 1) positions(:, 2)]);
nPoints = numel(Orte(:, 1));
physicalDimension(1, :) = max(positions(:, 1));
physicalDimension(2, :) = max(positions(:, 2));
%% extract positions
% to avoid edge effects only points within a minimum distance of
% maxDiameter to the image boundaries are considered
mapSize = ceil(max(positions(:, 1:2)));
currentPositions = positions(positions(:, 1) > maxRadius & positions(:, 1) <= mapSize(1)-maxRadius...
            & positions(:, 2) > maxRadius & positions(:, 2) <= mapSize(2)-maxRadius, :);
% check if currentPositions is empty
if isempty(currentPositions)
   error('Image size is to small or maxRadius is to big! Can not find any points which have maxRadius distance from the image borders! Try to decrease maxRadius or increase image region.'); 
end
%% create subgrid
mapSizeInnerBound = floor( (max(currentPositions(:, 1:2)) - min(currentPositions(:, 1:2)))./gridSize );
gridPositions = cell(mapSizeInnerBound(1), mapSizeInnerBound(2));
xPos = min(currentPositions(:, 2));
% yPos = (min(currentPositions(:,
% 2))):gridSize:(mapSizeInnerBound(1)*gridSize); %for parfor loop
yPos = min(currentPositions(:, 1));
for ii = 1:mapSizeInnerBound(1)
    for ll = 1:mapSizeInnerBound(2)
        gridPositions{ii, ll} = currentPositions(currentPositions(:, 1) > xPos & currentPositions(:, 1) <= (xPos + gridSize)...
            & currentPositions(:, 2) > yPos & currentPositions(:, 2) <= (yPos + gridSize), :);
        yPos = yPos + gridSize;
    end
     % set back y for next grid line
    yPos = min(currentPositions(:, 1));
    xPos = xPos + gridSize;
end
%% radial density function algorithm
% create bins for radial density function
nbins = ceil(maxRadius/binSize);
[binArray, edges] = histcounts([0 maxRadius], nbins);
binArray = zeros(size(binArray, 2), mapSizeInnerBound(1)*mapSizeInnerBound(2));

assignedGrid = zeros(mapSizeInnerBound(1)*mapSizeInnerBound(2), 1);

% for waitbar
prevPercent = 0;
counter = 1;
% loop through all grids
nGrids = 1;
for ii = 1:size(gridPositions, 1)
    for kk = 1:size(gridPositions, 2)
        % loop through all radii
        if isempty(gridPositions{ii, kk})
            binArray(:, nGrids) = 0;
            assignedGrid(nGrids, 1) = numel(gridPositions{ii ,kk}(:,1));
            nGrids = nGrids + 1;
            % waitbar
            currentPercent = fix(100*counter/(mapSizeInnerBound(1)*mapSizeInnerBound(2)));
            if currentPercent > prevPercent
                multiWaitbar( 'computing radial density function...', 'Value', counter/(mapSizeInnerBound(1)*mapSizeInnerBound(2)));
                prevPercent = currentPercent;
            end
            counter = counter + 1;
            continue
        else
            for hh = 1:numel(gridPositions{ii ,kk}(:,1))
                for ll = 1:size(binArray, 1)
                    % get number of points within radii
                    radInDistances = kdtree_ball_query(treeData, [ gridPositions{ii, kk}(hh, 1) gridPositions{ii, kk}(hh ,2)], edges(ll)); % get number of points with in r
                    radOutDistances = kdtree_ball_query(treeData, [ gridPositions{ii ,kk}(hh ,1) gridPositions{ii, kk}(hh ,2)], edges(ll + 1)); % get number of points with in r + dr
                    binArray(ll, nGrids) = binArray(ll, nGrids) + numel(radOutDistances) - numel(radInDistances); % number of points within shell dr
                end
            end
            assignedGrid(nGrids, 1) = numel(gridPositions{ii ,kk}(:,1));
            nGrids = nGrids + 1;
        end
        % waitbar
        currentPercent = fix(100*counter/(mapSizeInnerBound(1)*mapSizeInnerBound(2)));
        if currentPercent > prevPercent
            multiWaitbar( 'computing radial density function...', 'Value', counter/(mapSizeInnerBound(1)*mapSizeInnerBound(2)));
            prevPercent = currentPercent;
        end
        counter = counter + 1;
    end
end
%% normalize data (g(r) = binArray(ii) / (binArea*density*nPoints)) - center of each bin is taken as r
density = nPoints / (physicalDimension(1, :) * physicalDimension(2, :)); % take thresholded area? thats why nuclear area is increasing
for ff = 1: size(binArray, 2) % second dimension determines total number of grids
    for ii = 1 : nbins
        innerPart = edges(1, ii) + binSize/2;
        outerPart = innerPart + binSize;
        binArray(ii, ff) = binArray(ii, ff) / (pi*outerPart^2 - pi*innerPart^2);
        binArray(ii, ff) = binArray(ii, ff) ./ assignedGrid(ff, 1);
    end
end
binArray = binArray ./ density ;
binArray(isnan(binArray))=0;
multiWaitbar( 'computing radial density function...', 'Close');
%% kmeans clustering
[idx, C] = kmedoids(binArray.', nCluster); % choose different distance
%% plotting 2
% plot radial density functions
figure
for gg = 1:nCluster
    p1 = mean(binArray(:, find((idx==gg)).'), 2);
    plot(edges(1, 1:end-1) + binSize/2, p1);
    hold on
end
hold off
% plot color coded grids
B = reshape(idx, mapSizeInnerBound(2), mapSizeInnerBound(1)).'; % switch 2 and 1
figure
imagesc(B)
axis equal
hold on
%scatter(positions(:, 2)./gridSize-maxRadius./gridSize, positions(:, 1)./gridSize-maxRadius./gridSize) % offset is wrong
hold off

% save results before 