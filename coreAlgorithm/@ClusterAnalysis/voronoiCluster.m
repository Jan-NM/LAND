function [ varargout ] = voronoiCluster( obj, isRandom, dim, showImage, saveImageFolder, saveImageFilename )
%voronoiCluster Summary of this function goes here
%   Detailed explanation goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    
    case 1
        isRandom = 0;
        dim = 2;
        showImage = 0;
        saveImageFolder = "";
        saveImageFilename = "";
    case 2
        dim = 2;
        showImage = 0;
        saveImageFolder = "";
        saveImageFilename = "";
    case 3
        showImage = 0;
        saveImageFolder = "";
        saveImageFilename = "";
    case 4
        saveImageFolder = "";
        saveImageFilename = "";
    case 5
        saveImageFilename = "";
    case 6
        
        
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
if dim == 2
    dt = delaunayTriangulation(double(positions(:,1)),double(positions(:,2)));
elseif dim == 3
    dt = delaunayTriangulation(double(positions(:,1)),double(positions(:,2)),double(positions(:,3)));
else
    error("Only dimension 2 and 3 are supported")
end
%% voronoi vertices and regions
[V,R] = voronoiDiagram(dt);

A = 0;
smallest_vol = Inf;
biggest_vol = 0;
xmin = min(positions(:, 1));
xmax = max(positions(:, 1));
ymin = min(positions(:, 2));
ymax = max(positions(:, 2));
if dim == 3
   zmin = min(positions(:, 3));
   zmax = max(positions(:, 3)); 
end
reasonables = false(1,length(R));
j = 1;

area_or_vol = "Area";
if dim == 3
    area_or_vol = "Volume";
end

saveImage = false;
if saveImageFolder ~= ""
   saveImage = true;
end

%% calculate biggest and smallest area/volume if whole voronoi cluster is inside the cell area
%% for each region check if inside
for i = 1:length(R)
    reasonable = true;
    for jxy=V(R{i},:)
        if jxy(1) <= xmin || jxy(1) >= xmax || jxy(2) <= ymin || jxy(2) >= ymax
            reasonable = false;
            break
        end
    end
    if reasonable == false
        continue
    end
    reasonables(i) = 1;
    if dim == 2
        volume = polyarea(V(R{i}([1:end,1]),1),V(R{i}([1:end,1]),2));
    else
         XR10 = V(R{i},:);
         [K,volume] = convhull(XR10);
    end
    if isnan(volume) == false && isinf(volume) == false
     A(j) = volume;
     j = j + 1;
     if volume < smallest_vol
         smallest_vol = volume;
     end
     if volume > biggest_vol
         biggest_vol = volume;
     end
    end
end
disp("Size A: ")
A = A';
disp(size(A))
densities = 1./A;
if showImage || saveImage
    if ~showImage
        f = figure('visible', 'off');
    else
        f = figure
    end
    histogram(densities,40)
    set(gca,'yscale','log')
    title("Density histogram")
    xlabel("Density in nm^{-"+dim+"}")
    ylabel("Number of clusters")
    if saveImage
        saveas(f, saveImageFolder+"/"+saveImageFilename+"_density_histo", "png");
    end
end

disp("quantile 0.5%")
lowest_plot_density = quantile(densities, 0.005)
disp("quantile 99.5%")
highest_plot_density = quantile(densities, 0.995)

% cumulative plot
if showImage || saveImage
    if ~showImage
        f = figure('visible', 'off');
    else
        f = figure
    end
    [yvals,xvals] = ecdf(densities);
    semilogx(xvals,yvals,'-r')
    title("Cumulative density plot")
    ylabel("Cumulative fraction")
    xlabel("Voronoi density")
    if saveImage
        saveas(f, saveImageFolder+"/"+saveImageFilename+"_cum_density", "png");
    end
    
    if ~showImage
        f = figure('visible', 'off');
    else
        f = figure
    end
    histogram(A,40)
    set(gca,'yscale','log')
    xlabel(area_or_vol+" in nm^{"+dim+"}")
    ylabel("Number of clusters")
    title(area_or_vol+" histogram")
    
    if saveImage
        saveas(f, saveImageFolder+"/"+saveImageFilename+"_"+area_or_vol+"_histo", "png");
    end
end
    

disp("Smallest "+area_or_vol+": ")
disp(smallest_vol)

disp("Biggest "+area_or_vol+": ")
disp(biggest_vol)

if showImage || saveImage
    if ~showImage
        f = figure('visible', 'off');
    else
        f = figure
    end
    if dim == 2
        scatter(double(positions(:,1)),double(positions(:,2)),'.r')
    else
        scatter3(double(positions(:,1)),double(positions(:,2)),double(positions(:,3)),'.r')
    end
    axis equal
    xlim([xmin,xmax])
    ylim([ymin,ymax])
    if dim == 3
        zlim([zmin,zmax])
        zlabel("z in nm")
    end
    title("Localizations")
    xlabel("x in nm")
    ylabel("y in nm")
    if saveImage
        saveas(f, saveImageFolder+"/"+saveImageFilename+"_locs", "png");
    end
 
    if ~showImage
        f = figure('visible', 'off');
    else
        f = figure
    end
    hold on
    j = 1;
    if dim == 3
        z_mean = mean(double(positions(:,3)));
    end
        
    cmap = parula(1000);
    caxis([lowest_plot_density highest_plot_density]);
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
           if dim == 2
                h = fill(V(R{i},1),V(R{i},2), color);
                set(h,'EdgeColor','none')
           else 
               % take only a slice for faster plotting
                if any(V(R{i},3) >= z_mean - 20) && any(V(R{i},3) <= z_mean + 20)
                    XR10 = V(R{i},:);
                    K = convhull(XR10);
                    trisurf(K, XR10(:,1) ,XR10(:,2) ,XR10(:,3) , 'FaceColor', color, 'FaceAlpha',0.9, 'EdgeColor', 'none')
                end
           end
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
    if dim == 3
        zlim([zmin,zmax])
        zlabel("z in nm")
    end
    set(gca,'Color','k')
    if dim == 2
        title({"Voronoi diagram"; "Blue -> Yellow (range: 0.5% - 99.5%)"})
    else 
        title({"Voronoi diagram (z in ["+(z_mean-20)+","+(z_mean+20)+"])"; "Blue -> Yellow (range: 0.5% - 99.5%)"})
    end
    xlabel("x in nm")
    ylabel("y in nm")
    c = colorbar;
    c.Label.String = "Density";
    c.Limits = [lowest_plot_density highest_plot_density];
    if saveImage
        saveas(f, saveImageFolder+"/"+saveImageFilename+"_voronoi", "png");
    end
end
varargout{1} = smallest_vol;
varargout{2} = biggest_vol;
varargout{3} = lowest_plot_density;
varargout{4} = highest_plot_density;

end
