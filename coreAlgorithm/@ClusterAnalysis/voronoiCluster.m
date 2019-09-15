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
dt = delaunayTriangulation(double(positions(:,1)),double(positions(:,2)));
%% voronoi vertices and regions
[V,R] = voronoiDiagram(dt);

A = 0;
smallest_area = Inf;
biggest_area = 0;
xmin = min(positions(:, 1));
xmax = max(positions(:, 1));
ymin = min(positions(:, 2));
ymax = max(positions(:, 2));
reasonables = false(1,length(R));
j = 1;
%% calculate biggest and smallest area if whole voronoi cluster is inside the cell area
%% for each region check if inside
for i = 1:length(R)
    reasonable = true;
    for jxy=V(R{i},:)'
        if jxy(1) <= xmin || jxy(1) >= xmax || jxy(2) <= ymin || jxy(2) >= ymax
            reasonable = false;
            break
        end
    end
    if reasonable == false
        continue
    end
    reasonables(i) = 1;
    area = polyarea(V(R{i}([1:end,1]),1),V(R{i}([1:end,1]),2));
    if isnan(area) == false && isinf(area) == false
     A(j) = area;
     j = j + 1;
     if area < smallest_area
         smallest_area = area;
     end
     if area > biggest_area
         biggest_area = area;
     end
    end
end
disp("Size A: ")
A = A';
disp(size(A))
densities = 1./A;
figure
histogram(densities,40)
set(gca,'yscale','log')
title("Density histogram")
xlabel("Density in nm^{-2}")
ylabel("Number of clusters")

disp("quantile 0.5%")
lowest_plot_density = quantile(densities, 0.005)
disp("quantile 99.5%")
highest_plot_density = quantile(densities, 0.995)

% cumulative plot
figure
[yvals,xvals] = ecdf(densities);
semilogx(xvals,yvals,'-r')
title("Cumulative density plot")
ylabel("Cumulative fraction")
xlabel("Voronoi density")


figure
histogram(A,40)
set(gca,'yscale','log')
xlabel("Area in nm^{2}")
ylabel("Number of clusters")

title("Area histogram")
disp("Smallest area: ")
disp(smallest_area)

disp("Biggest area: ")
disp(biggest_area)

figure 
plot(double(positions(:,1)),double(positions(:,2)),'.r')
axis equal
xlim([xmin,xmax])
ylim([ymin,ymax])
title("Localizations")
xlabel("x in nm")
ylabel("y in nm")

figure
hold on
j = 1;
cmap = parula(1000);
for i=1:length(R)
    if reasonables(i)
       cmap_row = floor(size(cmap, 1)*((1/A(j)-lowest_plot_density)/(highest_plot_density-lowest_plot_density)))+1;
       if cmap_row >= size(cmap, 1)+1
           cmap_row = size(cmap, 1);
       end
       color = [0 0 0];
       if cmap_row > 0
          color = cmap(cmap_row,:);
       end
       h = fill(V(R{i},1),V(R{i},2), color);
       set(h,'EdgeColor','none')
       j = j+1;
    else
        %plot(V(R{i}([1:end,1]),1),V(R{i}([1:end,1]),2),'-b')
    end
end
%plot(double(positions(:,1)),double(positions(:,2)),'.w');
hold off
%%plot(vx,vy,'-b',double(positions(:,1)),double(positions(:,2)),'.r');
axis equal
xlim([xmin,xmax])
ylim([ymin,ymax])
set(gca,'Color','k')
title({"Voronoi diagram"; "Blue -> yellow (range: 0.5% - 99.5%)"})
xlabel("x in nm")
ylabel("y in nm")

end
