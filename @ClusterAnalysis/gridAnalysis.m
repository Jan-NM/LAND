function [] = gridAnalysis( obj, gridSize, showPlot )
%gridAnalysis Summary of this function goes here
%   Detailed explanation goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get experimental and random data
% check that both have the same dimension
% bin data according to gridSize
% get number of elements in each bin
% compare to random data (divide by random data)
% apply color code
% check if random data exist - if not algorithm not possible
%% init
switch nargin
    case 2
        showPlot = 0;
    case 3
        
    otherwise
        error('Wrong number of input arguments!')
end
dataMat = obj.clusterStruct;
%% bin data
x1 = double(ceil(obj.positionTable(:, 1)./gridSize));
y1 = double(ceil(obj.positionTable(:, 2)./gridSize));
% random data
x1Random = double(ceil(obj.randomTable(:, 1)./gridSize));
y1Random = double(ceil(obj.randomTable(:, 2)./gridSize));

minX = min(min(x1), min(x1Random));
minY = min(min(y1), min(y1Random));
maxX = max(max(x1), max(x1Random));
maxY = max(max(y1), max(y1Random));

[xSize, ySize] = size(full(sparse(maxX-minX+1,maxY-minY+1,1)));
% bin (experimental) localization table
x1 = double(ceil(obj.positionTable(:, 1)./ gridSize));
y1 = double(ceil(obj.positionTable(:, 2)./ gridSize));
% generate histogram binned image
binnedExp = sparse(x1-min(x1) + 1, y1-min(y1) + 1, 1, xSize, ySize);
binnedExp = full(binnedExp) ./ (gridSize/1000 * gridSize/1000);
% binnedExp = full(binnedExp);
% bin (random) localization table
x1 = double(ceil(obj.randomTable(:, 1)./gridSize));
y1 = double(ceil(obj.randomTable(:, 2)./gridSize));
% generate histogram binned image
binnedRan = sparse(x1-min(x1)+1, y1-min(y1)+1, 1, xSize, ySize);
binnedRan = full(binnedRan) ./ (gridSize/1000 * gridSize/1000);


% calculate mean of binned random image
meanBinnedRan = mean(mean(binnedRan));
meanBinnedExp = mean(mean(binnedExp));
% calculate percent deviation from random
compositeImage = ((binnedExp ./ meanBinnedExp));
if showPlot == true
     figure
     % compositeImage(compositeImage > 20 * (meanBinnedExp/meanBinnedExp)) = 2;
     % compositeImage(compositeImage~=2) = compositeImage(compositeImage~=2) ./10;
     % compositeImage = log(compositeImage);
     imagesc(compositeImage);
     axis off
     axis image
end

%% calculate Variance Mean Ratio (VMR)
% var = (sum(sum(binnedExp.^2)) - (sum(sum(binnedExp)).^2/ (size(binnedExp, 1) * size(binnedExp, 2))) )/ ((size(binnedExp, 1) * size(binnedExp, 2)) - 1);
% meanVal = sum(sum(binnedExp)) / (size(binnedExp, 1) * size(binnedExp, 2));
% dataMat(1).VMR = var/meanVal
% 
% obj.clusterStruct = dataMat;
% dof = (size(binnedExp, 1) * size(binnedExp, 2));
% value = (sum(sum(binnedExp.^2)) - (sum(sum(binnedExp)).^2/ (size(binnedExp, 1) * size(binnedExp, 2))) ) / meanVal;
% 
% chi2cdf(value, dof)

% var = (sum(sum(binnedRan.^2)) - (sum(sum(binnedRan)).^2/ (size(binnedRan, 1) * size(binnedRan, 2))) )/ ((size(binnedRan, 1) * size(binnedRan, 2)) - 1);
% mean = sum(sum(binnedRan)) / (size(binnedRan, 1) * size(binnedRan, 2));
% var/mean
end