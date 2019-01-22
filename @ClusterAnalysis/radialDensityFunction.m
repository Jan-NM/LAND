function [] = radialDensityFunction( obj, binSize,  maxRadius, isRandom, showPlot )
%radialDensityFunction computes radial density function function using the point
%coordinates
%
% Input:
%   binSize: size of the shell (dr) in nm
%   maxRadius: maximum radius for investigation in nm (only points with this 
%           minimum distance to the image boundary are considered)
%   isRandom: true if random positions should be used
%   showPlot: plot pair correlation finction
%   
% Output:
%   (1) g(r) values for each bin
%   (2) center of each bin
%
% Jan Neumann, 27.06.2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% init
multiWaitbar('computing radial density function...', 0);
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
    physicalDimension(1, :) = obj.randomPhysicalDimension(1, :);
    physicalDimension(2, :) = obj.randomPhysicalDimension(2, :);
else
    positions = obj.positionTable;
    dataMat = obj.clusterStruct;
    nPoints = obj.NoOfPoints;
    treeData = obj.queryData;
    physicalDimension(1, :) = obj.physicalDimension(1, :);
    physicalDimension(2, :) = obj.physicalDimension(2, :);
end
%% radial density function algorithm
% create bins for radial density function
nbins = ceil(maxRadius/binSize);
[binArray, edges] = histcounts([0 maxRadius], nbins);
binArray = zeros(size(binArray, 2), 1);

% to avoid edge effects only points within a minimum distance of
% maxDiameter to the image boundaries are considered
mapSize = ceil(max(positions(:, 1:2)));
currentPositions = positions(positions(:, 1) > maxRadius & positions(:, 1) < mapSize(1)-maxRadius...
            & positions(:, 2) > maxRadius & positions(:, 2) < mapSize(2)-maxRadius, :);
% check of currentPositions is empty
if isempty(currentPositions)
   error('Image size is to small or maxRadius is to big! Can not find any points which have maxRadius distance from the image borders! Try to decrease maxRadius or increase image region.'); 
end

% for waitbar
prevPercent = 0;
counter = 1;
% loop through all points
for ii = 1:size(currentPositions, 1)
    % loop through all radii
    for ll = 1:size(binArray, 1)
        % get number of points within radii
        switch obj.flagSearchTree
            case 'cTree'
                radInDistances = kdtree_ball_query(treeData, currentPositions(ii, 1:2), edges(ll)); % get number of points with in r 
                radOutDistances = kdtree_ball_query(treeData, currentPositions(ii, 1:2), edges(ll + 1)); % get number of points with in r + dr
            case 'matlabTree'
                radInDistances = rangesearch(currentPositions(ii, 1:2), treeData.X, edges(ll));
                radOutDistances = rangesearch(currentPositions(ii, 1:2), treeData.X, edges(ll + 1));
                radInDistances = (find(~cellfun('isempty', radInDistances))).';
                radOutDistances = (find(~cellfun('isempty', radOutDistances))).';
        end       
        binArray(ll, 1) = binArray(ll, 1) + numel(radOutDistances) - numel(radInDistances); % number of points within shell dr
    end
    % waitbar
    currentPercent = fix(100*counter/nPoints);
    if currentPercent > prevPercent
        multiWaitbar( 'computing radial density function...', 'Value', counter/nPoints);
        prevPercent = currentPercent;
    end
    counter = counter + 1;
end

binArray = full(binArray);
nPointsFinal = size(currentPositions, 1);
%% normalize data (g(r) = binArray(ii) / (binArea*density*nPoints)) - center of each bin is taken as r
density = nPoints / (physicalDimension(1, :) * physicalDimension(2, :));
for ii = 1 : nbins
   innerPart = edges(1, ii) + binSize/2;
   outerPart = innerPart + binSize;
   binArray(ii) = binArray(ii) / (pi*outerPart^2 - pi*innerPart^2);
end
binArray = binArray ./ nPointsFinal ./ density ;
multiWaitbar( 'computing radial density function...', 'Close');
%% visualization
if showPlot == true
    figure( 'Name', num2str(obj.sampleID) );
    h =  plot(edges(1 ,1:end-1) + binSize/2, binArray);
    grid on;
    title('radial density function');
    xlabel('distance [nm]');
    ylabel('g(r)');
    xlim([min(edges) max(edges)]);
end
%% output
dataMat(1).RDF = binArray;
dataMat(2).RDF = (edges(1, 1:end-1) + binSize/2).';
if isRandom == true
    obj.randomClusterStruct = dataMat;
else
    obj.clusterStruct = dataMat;
end
end