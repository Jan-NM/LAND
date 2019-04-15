%% for 3D
imageSize = 3000; % image size in nm
density = 150; % density of points in 1/µm² over entire image
densityCluster = 5; % density of cluster in 1/µm² over entire image
clusterPoints = 0.6; % percentage of points within cluster
clusterSize = 30; % diameter of cluster in nm;

%% calculate parameters
% calculate number of points from density
nPoints = poissrnd(density*imageSize/1000*imageSize/1000*imageSize/1000); % 1000 because localDensity is in µm^2
% get number of clustered points
clusteredPoints = ceil(clusterPoints * nPoints);
% number of non-clustered points
nonClusterPoints = nPoints - clusteredPoints;
% get number of cluster, assume cluster are randomly distributed over the
% entire image
nCluster = ceil(poissrnd(densityCluster*imageSize/1000*imageSize/1000*imageSize/1000)); % 1000 because localDensity is in µm^2
% spread clustered points evenly across the cluster
pointsPerCluster = ceil(clusteredPoints / nCluster);
% due to rounding reassign points
nonClusterPoints = nonClusterPoints - (pointsPerCluster * nCluster - clusteredPoints);
clusteredPoints = clusteredPoints + (pointsPerCluster * nCluster - clusteredPoints);
%% simulate non-clustered points
% pre-allocate Orte matrix for non-clustered points
nonClusterOrte = zeros(nonClusterPoints, 12);
% create noise points, random distributed, representing monomer fraction
% create loclization precision with mean of 9 nm and standard deviation of 5 nm
nonClusterOrte(:, 4) = abs(normrnd(9 ,5 , [nonClusterPoints 1]));
nonClusterOrte(:, 5) = nonClusterOrte(:, 4);
nonClusterOrte(:, 12) = nonClusterOrte(:, 4);
nonClusterOrte(:, 2) = rand(nonClusterPoints, 1).*imageSize;
nonClusterOrte(:, 3) = rand(nonClusterPoints, 1).*imageSize;
nonClusterOrte(:, 11) = rand(nonClusterPoints, 1).*imageSize;
%% simulate clustered points
% pre-allocate Orte matrix for non-clustered points
clusterOrte = zeros(clusteredPoints, 12);
xPoints = [];
yPoints = [];
zPoints = [];
for ii = 1:nCluster
    % positon of cluster
    xPos = rand(1, 1).*imageSize;
    yPos = rand(1, 1).*imageSize;
    zPos = rand(1, 1).*imageSize;
    % currentXPoints = (rand(pointsPerCluster, 1).* clusterSize) - clusterSize/2 + xPos;
    % currentYPoints = (rand(pointsPerCluster, 1).* clusterSize) - clusterSize/2 + yPos;
    azimuthAngle = 2*pi*rand(pointsPerCluster, 1);
    polarAngle = acos(2*rand(pointsPerCluster, 1)-1);
    % uniform distribution of points
    % radius = (clusterSize/2)*(rand(pointsPerCluster, 1));
    
    % gaussian distribution of points
    radius = abs(normrnd(0, clusterSize/2, [pointsPerCluster, 1]));
    currentXPoints = radius.*sin(polarAngle).*cos(azimuthAngle)+ xPos;
    currentYPoints = radius.*sin(polarAngle).*sin(azimuthAngle)+ yPos;
    currentZPoints = radius.*cos(polarAngle)+ zPos;
    xPoints = cat(1, currentXPoints, xPoints);
    yPoints = cat(1, currentYPoints, yPoints);
    zPoints = cat(1, currentZPoints, zPoints);
end

clusterOrte(:, 2) = xPoints;
clusterOrte(:, 3) = yPoints;
clusterOrte(:, 11) = zPoints;
clusterOrte(:, 4) = abs(normrnd(9 ,5 , [clusteredPoints 1]));
clusterOrte(:, 5) = clusterOrte(:, 4);
clusterOrte(:, 12) = clusterOrte(:, 4);
%% combine both Orte files
Orte = cat(1, nonClusterOrte, clusterOrte);

x = Orte(:, 2);
y = Orte(:, 3);
z = Orte(:, 11);
scatter3(x,y,z,'.')
axis equal vis3d