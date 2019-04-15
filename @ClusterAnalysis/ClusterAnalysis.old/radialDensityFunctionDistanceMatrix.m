function [] = radialDensityFunctionDistanceMatrix( obj, binSize,  maxRadius, isRandom, showPlot )
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
%% default parameters
memoryThresholdParam = 0.8; % maximum of total memory that will be used 
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
    physicalDimension(1, :) = obj.randomPhysicalDimension(1, :);
    physicalDimension(2, :) = obj.randomPhysicalDimension(2, :);
else
    positions = obj.positionTable;
    dataMat = obj.clusterStruct;
    nPoints = obj.NoOfPoints;
    physicalDimension(1, :) = obj.physicalDimension(1, :);
    physicalDimension(2, :) = obj.physicalDimension(2, :);
end
%% radial density function algorithm
% check available memory
[~, spaceOnRam] = memory;
% limit memory consumption to memoryThresholdParam times total available memory
memoryThreshold = memoryThresholdParam*spaceOnRam.PhysicalMemory.Available;
maxNoOfDistances = ceil(memoryThreshold / 8 / nPoints); % /8 for data type double
if maxNoOfDistances > nPoints
    maxNoOfDistances = nPoints;
end
% to avoid edge effects only points within a minimum distance of
% maxDiameter to the image boundaries are considered
mapSize = ceil(max(positions(:, 1:2)));
currentPositions = positions(positions(:, 1) > maxRadius & positions(:, 1) < mapSize(1)-maxRadius...
            & positions(:, 2) > maxRadius & positions(:, 2) < mapSize(2)-maxRadius, :);
% check of currentPositions is empty
if isempty(currentPositions)
   error('Image size is to small or maxRadius is to big! Can not find any points which have maxRadius distance from the image borders! Try to decrease maxRadius or increase image region.'); 
end
% split Orte matrix       
chunkSize = diff(fix(linspace(0, size(currentPositions, 1), ceil(size(currentPositions, 1)/maxNoOfDistances) + 1))); % + 1 should stay
splittedPositions = mat2cell(currentPositions, chunkSize, size(currentPositions, 2));

% create bins for radial density histogram
nbins = ceil(maxRadius/binSize);
[binArray, edges] = histcounts([0 maxRadius], nbins);
binArray = zeros(size(binArray, 2), 1);

% for waitbar
prevPercent = 0;
counter = 1;
for ii = 1:size(splittedPositions, 1)
    distances = pdist2(splittedPositions{ii}(:, 1:2), positions(:, 1:2));
    % keep distances which are smaller than maxRadius and exclude
    % self-distances
    [ivec, jvec] = find(distances <= maxRadius & distances > 0);
    xvec = double(distances( sub2ind([size(distances, 1) size(distances, 2)], ivec, jvec) ).');
        
    ivec = ivec.';
    jvec = jvec.';
    % create sparse matrix
    distances = sparse(ivec, jvec, xvec, chunkSize(ii), nPoints);

    for ll = 1:size(binArray, 1)
        binArray(ll, 1) = binArray(ll, 1) + sum(nonzeros(distances(:, :)) > edges(ll) ) - sum(nonzeros(distances(:, :)) >= edges(ll + 1) );
    end
    % waitbar
    currentPercent = fix(100*counter/size(splittedPositions, 1));
    if currentPercent > prevPercent
        multiWaitbar( 'computing radial density function...', 'Value', counter/size(splittedPositions, 1));
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