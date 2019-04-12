function [ varargout ] = voronoiCluster( obj, isRandom, showImage )
%voronoiCluster Summary of this function goes here
%   Detailed explanation goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    
    case 1
        isRandom = 0;
        showImage = 0;
    case 2
        showImage = 0;
    case 3
        
    otherwise
        error('Wrong number of input arguments!')
end
if isRandom == true
    positions = obj.randomTable;
    dataMat = obj.randomClusterStruct;
    treeData = obj.randomQueryData;
else
    positions = obj.positionTable;
    dataMat = obj.clusterStruct;
    treeData = obj.queryData;
end

%% Voronoi
[vx,vy]=voronoi(double(positions(:,1)), double(positions(:,2)));
figure
plot(positions(:, 2), positions(:, 1), 'r+', vy, vx, 'b-')
axis equal
xlim([min(positions(:, 2)) max((positions(:, 2)))])
ylim([min(positions(:, 1)) max(positions(:, 1))])
for i = 1:length(vy)
    if all(vy{i}~=1)
        h = patch(vx(vy{i},1),vx(vy{i},2),i);
       % get(h);
    end
end
% A = polyarea(x,y);
A = polyarea(vx(vy{1},1),vx(vy{1},2));
for i = 1:length(vy)
    if all(vy{i}~=1)
        A(i) = polyarea(vx(vy{i},1),vx(vy{i},2));
    end
end
A=A.';
histogram(A(:,1),1000000)
xlim([0 10000])

end
