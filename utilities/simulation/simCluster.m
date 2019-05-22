function [Orte] = simCluster(parameters, dim, varargin)
% simCluster simulates clusters (2D/3D)
%   to evaluate localization microscopy clustering algorithms
%   
%   Input:  parameters      options: 'default' uses default parameters
%                                    'non-default' uses user-supplied
%                                    parameters
%           dim             2 or 3 for 2D or 3D simulation
%           If 'non-default' is specified, the user must supply the
%           following parameters: (1) image size in nm, (2) density of points
%           in 1/µm² over entire image, (3) density of cluster in 1/µm²
%           over entire image, (4) fraction of points within clusters,
%           (5) diameter of cluster in nm, (6) mean localization precision
%           in nm, (6) standard deviation of localization precision in nm
%           example: simulatedTable = simClusters('non-default', 1000, 50000, 5, 0.8, 100, 20, 10)           
%
%   Output: list of localizations. Localization table should be in Orte file format i.e. 
%           column 2/3/11 = x/y/z-position, column 3/5/12 = x/y/z-localization precision,
%           column 9 = frame number.
%
%   Note: Points of clusters that are located outside the image border are
%   removed.
%
% Jan Neumann, 19.11.2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% init parameters
narginchk(2, 9);
switch parameters
    case 'default'
        imageSize = 20000; % image size in nm
        density = 100; % density of points in 1/µm² over entire image
        densityCluster = 5; % density of cluster in 1/µm² over entire image
        clusterPoints = 0.7; % fraction of points within clusters
        clusterSize = 60; % diameter of cluster in nm;
        meanLoc = 9;
        stdLoc = 5;
    case 'non-default'
        imageSize = varargin{1};
        density = varargin{2};
        densityCluster = varargin{3};
        clusterPoints = varargin{4};
        clusterSize = varargin{5};
        meanLoc = varargin{6};
        stdLoc = varargin{7};
    otherwise
        disp('Wrong number of input arguments!')
end
%% calculate parameters
if dim == 2
    imageDim = imageSize/1000*imageSize/1000; % 1000 because localDensity is in µm^2
else
    imageDim = imageSize/1000*imageSize/1000*imageSize/1000;
end    
% calculate number of points from density
nPoints = poissrnd(density*imageDim);
% get number of clustered points
clusteredPoints = ceil(clusterPoints * nPoints);
% number of non-clustered points
nonClusterPoints = nPoints - clusteredPoints;
% get number of cluster, assume cluster are randomly distributed over the
% entire image
nCluster = ceil(poissrnd(densityCluster*imageDim));
% spread clustered points evenly across the cluster
pointsPerCluster = ceil(clusteredPoints / nCluster);
% due to rounding reassign points
nonClusterPoints = nonClusterPoints - (pointsPerCluster * nCluster - clusteredPoints);
%% simulate non-clustered points
% pre-allocate Orte matrix for non-clustered points
nonClusterOrte = zeros(nonClusterPoints, 9);
% create noise points, random distributed, representing monomer fraction
% create loclization precision with mean of meanLoc nm and standard deviation of stdLoc nm
nonClusterOrte(:, 4) = abs(normrnd(meanLoc, stdLoc, [nonClusterPoints 1]));
nonClusterOrte(:, 5) = nonClusterOrte(:, 4);
nonClusterOrte(:, 12) = nonClusterOrte(:, 4);
nonClusterOrte(:, 2) = rand(nonClusterPoints, 1).*imageSize;
nonClusterOrte(:, 3) = rand(nonClusterPoints, 1).*imageSize;
nonClusterOrte(:, 11) = rand(nonClusterPoints, 1).*imageSize;
%% simulate clustered points
xPoints = [];
yPoints = [];
if dim == 3
    zPoints = [];
end
for ii = 1:nCluster
    % positon of cluster
    xPos = rand(1, 1).*imageSize;
    yPos = rand(1, 1).*imageSize;
    % currentXPoints = (rand(pointsPerCluster, 1).* clusterSize) - clusterSize/2 + xPos;
    % currentYPoints = (rand(pointsPerCluster, 1).* clusterSize) - clusterSize/2 + yPos;
    if dim == 2
        angle = 2*pi*rand(pointsPerCluster, 1);
        % uniform distribution of points
        % radius = (clusterSize/2)*(rand(pointsPerCluster, 1));
    
        % gaussian distribution of points
        radius = abs(normrnd(0, clusterSize/2, [pointsPerCluster, 1]));
        currentXPoints = radius.*cos(angle)+ xPos;
        currentYPoints = radius.*sin(angle)+ yPos;
    else
        zPos = rand(1, 1).*imageSize;
        azimuthAngle = 2*pi*rand(pointsPerCluster, 1);
        polarAngle = acos(2*rand(pointsPerCluster, 1)-1);
        radius = abs(normrnd(0, clusterSize/2, [pointsPerCluster, 1]));
        currentXPoints = radius.*sin(polarAngle).*cos(azimuthAngle)+ xPos;
        currentYPoints = radius.*sin(polarAngle).*sin(azimuthAngle)+ yPos;
        currentZPoints = radius.*cos(polarAngle)+ zPos;
        zPoints = cat(1, currentZPoints, zPoints);
    end
    xPoints = cat(1, currentXPoints, xPoints);
    yPoints = cat(1, currentYPoints, yPoints);
end
% cut-off clustered points that are outside of the image border
if dim == 2
    coordinates = [xPoints, yPoints];
    coordinates = coordinates(coordinates(:, 1) > 0 & coordinates(:, 1) < imageSize & coordinates(:, 2) > 0 & coordinates(:, 2) < imageSize, :);
else
    coordinates = [xPoints, yPoints, zPoints];
    coordinates = coordinates(coordinates(:, 1) > 0 & coordinates(:, 1) < imageSize &...
        coordinates(:, 2) > 0 & coordinates(:, 2) < imageSize & coordinates(:, 3) > 0 & coordinates(:, 3) < imageSize, :);
end
    % pre-allocate Orte matrix for non-clustered points
clusteredPoints = size(coordinates, 1);
clusterOrte = zeros(clusteredPoints, 12);

clusterOrte(:, 2) = coordinates(:, 1);
clusterOrte(:, 3) = coordinates(:, 2);
clusterOrte(:, 4) = abs(normrnd(meanLoc, stdLoc, [clusteredPoints 1]));
clusterOrte(:, 5) = clusterOrte(:, 4);
if dim == 3
    clusterOrte(:, 11) = coordinates(:, 3);
    clusterOrte(:, 12) = clusterOrte(:, 4);
end
%% combine both Orte files
Orte = cat(1, nonClusterOrte, clusterOrte);
end