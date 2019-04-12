function [] = gridAnalysis( obj, gridSize, isRandom, showPlot )
%gridAnalysis calculates the density of signals per bin area
%   Calculates the density of signals within a specified bin area and plots
%   the data color-coded.
%
% Input: 
%   gridSize: size of bin in nm
%   isRandom: true = use random data; false = use experimental data
%   showPlot: true = show color-coded image
%
% Output:
%   (1) returns pixel(bin) values
%
% Jan Neumann, 29.01.2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% init
switch nargin
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
else
    positions = obj.positionTable;
    dataMat = obj.clusterStruct;
end
%% bin data
x1 = double(ceil(positions(:, 1)./gridSize));
y1 = double(ceil(positions(:, 2)./gridSize));
% generate histogram binned image
binnedExp = sparse(x1 - min(x1) + 1, y1 - min(y1) + 1, 1);
dataMat(1).gridAnalysis = full(binnedExp) ./ (gridSize/1000 * gridSize/1000); % conversion to 1/µm²

if showPlot == true
     figure;
     imagesc(binnedExp);
     cbar = colorbar;
     ylabel(cbar, 'signals per µm²')
     colormap hot;
     axis off
     axis image
     ax = gca;
     ax.FontSize = 30;
end
%% output
if isRandom == true
    obj.randomClusterStruct = dataMat;
else
    obj.clusterStruct = dataMat;
end
end