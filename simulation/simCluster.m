function [Orte] = simCluster(parameters, varargin)
% simCluster simulates membrane receptor clusters
%   to evaluate localization microscopy clustering algorithms
%
%   Output: list of localizations. Localization table should be in Orte file format i.e. 
%           column 2/3 = x/y-position, column 3/5 = x/y-localization precision,
%           column 9 = frame number.
%
% Jan Neumann, 19.11.2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% init parameters
narginchk(1, 7);
switch parameters
    case 'default'
        imageSize = 20000; % image size in nm
        density = 100; % density of points in 1/µm² over entire image
        densityCluster = 5; % density of cluster in 1/µm² over entire image
        clusterPoints = 0.3; % percentage of points within cluster
        clusterSize = 60; % diameter of cluster in nm;
    case 'non-default'
        imageSize = varargin{1};
        density = varargin{2};
        densityCluster = varargin{3};
        clusterPoints = varargin{4};
        clusterSize = varargin{5};        
    otherwise
        disp('Wrong number of input arguments!')
end
%% calculate parameters
% calculate number of points from density
nPoints = poissrnd(density*imageSize/1000*imageSize/1000); % 1000 because localDensity is in µm^2
% get number of clustered points
clusteredPoints = ceil(clusterPoints * nPoints);
% number of non-clustered points
nonClusterPoints = nPoints - clusteredPoints;
% get number of cluster, assume cluster are randomly distributed over the
% entire image
nCluster = ceil(poissrnd(densityCluster*imageSize/1000*imageSize/1000)); % 1000 because localDensity is in µm^2
% spread clustered points evenly across the cluster
pointsPerCluster = ceil(clusteredPoints / nCluster);
% due to rounding reassign points
nonClusterPoints = nonClusterPoints - (pointsPerCluster * nCluster - clusteredPoints);
clusteredPoints = clusteredPoints + (pointsPerCluster * nCluster - clusteredPoints);
%% simulate non-clustered points
% pre-allocate Orte matrix for non-clustered points
nonClusterOrte = zeros(nonClusterPoints, 9);
% create noise points, random distributed, representing monomer fraction
% create loclization precision with mean of 9 nm and standard deviation of 5 nm
nonClusterOrte(:, 4) = abs(normrnd(9 ,5 , [nonClusterPoints 1]));
nonClusterOrte(:, 5) = nonClusterOrte(:, 4);
nonClusterOrte(:, 2) = rand(nonClusterPoints, 1).*imageSize;
nonClusterOrte(:, 3) = rand(nonClusterPoints, 1).*imageSize;
%% simulate clustered points
% pre-allocate Orte matrix for non-clustered points
clusterOrte = zeros(clusteredPoints, 9);
xPoints = [];
yPoints = [];
for ii = 1:nCluster
    % positon of cluster
    xPos = rand(1, 1).*imageSize;
    yPos = rand(1, 1).*imageSize;
    % currentXPoints = (rand(pointsPerCluster, 1).* clusterSize) - clusterSize/2 + xPos;
    % currentYPoints = (rand(pointsPerCluster, 1).* clusterSize) - clusterSize/2 + yPos;
    angle = 2*pi*rand(pointsPerCluster, 1);
    % uniform distribution of points
    % radius = (clusterSize/2)*(rand(pointsPerCluster, 1));
    
    % gaussian distribution of points
    radius = abs(normrnd(0, clusterSize/2, [pointsPerCluster, 1]));
    currentXPoints = radius.*cos(angle)+ xPos;
    currentYPoints = radius.*sin(angle)+ yPos;
    xPoints = cat(1, currentXPoints, xPoints);
    yPoints = cat(1, currentYPoints, yPoints);
end

clusterOrte(:, 2) = xPoints;
clusterOrte(:, 3) = yPoints;
clusterOrte(:, 4) = abs(normrnd(9 ,5 , [clusteredPoints 1]));
clusterOrte(:, 5) = clusterOrte(:, 4);
%% combine both Orte files
Orte = cat(1, nonClusterOrte, clusterOrte);
end