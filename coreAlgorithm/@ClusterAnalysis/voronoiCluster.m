function [] = voronoiCluster( obj, varargin )
%voronoiCluster computes Voronoi cells around each data point and uses them
%for density estimation
%
% The current version is implemented for 2D.
%
% Input:
%   obj:		 The cell object (given by cell.voronoiCluster)
%	varargin     struct containing the following fields in the same order:
%                   isRandom (default: false): true = use random data; false = use experimental data
%                   showImage (default: false): true = shows Voronoi
%							clustering plot, where intensity scales from
%							lowestVoronoiCellDensity to highestVoronoiCellDensity
%					lowestVoronoiCellDensity (default: 0.005 quantile):
%							cumulative probability of lowest voronoi cell 
%							density to plot, values lower than the given value
%							are assigned to that value
%					highestVoronoiCellDensity (default: 0.995 quantile):
%							cumulative probability of highest voronoi cell
%							density to plot, values	higher than the given
%							value are assigned to that value
%
% example for usage (will use default values for lowestVoronoiCellDensity 
% and highestVoronoiCellDensity):
%
% input = struct('isRandom', false, 'showImage', true);
% cell01 = ClusterAnalysis(dataInOrteFormat);
% cell01.voronoiCluster(input);
%
% Output
%   (1) Voronoi vertices
%   (2) Voronoi regions (as indices of corresponding vertices in (1) )
%   (3) Voronoi area (2D: µm^2) / volume (3D: µm^3)
%   (4) Voronoi density (2D: 1/µm^2) or (3D: 1/µm^3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% init
if nargin > 2
	error('Too many input arguments!')
end
if obj.dimension == 3
	error('Evaluetion of voronoi cells is not yet fully supported for 3D images!');
end
% voronoi cell density will be set later on to the default values
optargs = {false, false, 0.005, 0.995};
if ~isempty(varargin)
    numvarargs = numel(fieldnames(varargin{1}));
    fields = fieldnames(varargin{1});
    for ii = 1:numvarargs
        optargs{ii} = varargin{1}.(fields{ii});
    end
end
[isRandom, showImage, lowestVoronoiCellDensity,...
	highestVoronoiCellDensity] = optargs{:};

if isRandom == true
	positions = obj.randomTable;
	dataMat = obj.randomClusterStruct;
else
	positions = obj.positionTable;
	dataMat = obj.clusterStruct;
end
dim = obj.dimension;
multiWaitbar('Calculating Voronoi cells. Please wait!', 'busy');
%% Calculate voronoi vertices and regions via delaunay triangulation
if dim == 2
	dt = delaunayTriangulation( positions(:,1), positions(:,2) );
elseif dim == 3
	dt = delaunayTriangulation( positions(:,1) , positions(:,2), positions(:,3));
else
	error("Error in voronoi calculation: Function supports only 2D or 3D data!")
end
% V = Vertices
% R = each row represents indices of vertices of a voronoi region
[V, R] = voronoiDiagram(dt);
%% calculate voronoi areas / volumes
% filter regions that are not valid
xmin = min(positions(:, 1));
xmax = obj.physicalDimension(1);
ymin = min(positions(:, 2));
ymax = obj.physicalDimension(2);
if dim == 3
	zmin = min(positions(:, 3));
	zmax = obj.physicalDimension(3);
end
reasonables = false( numel(R), 1 );

for i = 1:numel(R)
	voronoiRegion = V(R{i},:);
	if dim == 2
		isAnyVerticeOutsideBoundary = any ( (voronoiRegion(:, 1) <= xmin) | (voronoiRegion(:, 1) >= xmax) |...
			(voronoiRegion(:, 2) <= ymin) | (voronoiRegion(:, 2) >= ymax) );
	else
		isAnyVerticeOutsideBoundary = any ( (voronoiRegion(:, 1) <= xmin) | (voronoiRegion(:, 1) >= xmax) |...
			(voronoiRegion(:, 2) <= ymin) | (voronoiRegion(:, 2) >= ymax) |...
			(voronoiRegion(:, 3) <= zmin) | (voronoiRegion(:, 3) >= zmax) );
	end
	if isAnyVerticeOutsideBoundary == true
		continue
	end
	reasonables(i) = true;
end

% filter reasonable regions
R = R(reasonables);

% calculate volume of voronoi regions
voronoi_region_size = nan( nnz(reasonables), 1 );
for i = 1:numel(R)
	[~, volume] = convhull(V(R{i}, :));
	if isnan(volume) == false && isinf(volume) == false
		voronoi_region_size(i) = volume;
	end
end

if dim == 2
	micro_densities = 10^6.* ( 1./voronoi_region_size );
	micro_area = 10^(-6).*voronoi_region_size;
end
%% plotting
% set quantile densities to corresponding densities in 1/µm^2
lowestVoronoiCellDensity = quantile(micro_densities, lowestVoronoiCellDensity);
highestVoronoiCellDensity = quantile(micro_densities, highestVoronoiCellDensity);

% works only for 2D
if showImage == true
	% set color scheme
	ncolors = 256;
	cmap = parula(ncolors);
	
	% find all R's that fit into same cmap_row
	% loop over cmap_row
	spacing = linspace(lowestVoronoiCellDensity, highestVoronoiCellDensity, ncolors);
	nBins = numel(spacing) - 1;
	A = cell( nBins + 2, 1 ); % + 2 for black values and values > then highest_plot_density
	for i = 1:nBins % last two elements will get specific attention later on
		ind = micro_densities(:) >= spacing(i) &  micro_densities(:) < spacing(i+1);
		% extract cells
		tempR = R(ind);
		% find max length
		maxSize = max( cellfun('size', tempR, 2) );
		% allocate NaN array
		tempArray_x = NaN( size(tempR, 1), maxSize );
		% replace values with voronoi vertices
		for k = 1:size(tempR, 1)
			temp_row = tempR{k};
			temp_max = size(temp_row, 2);
			tempArray_x(k, 1:temp_max) = temp_row;
		end
		A{i} =  tempArray_x;
	end
	% handle voronoi cells smaller than lowestVoronoiCellDensity
	ind = micro_densities(:) < spacing(1);
	tempR = R(ind);
	maxSize = max( cellfun('size', tempR, 2) );
	tempArray_x = NaN( size(tempR, 1), maxSize );
	for k = 1:size(tempR, 1)
		temp_row = tempR{k};
		temp_max = size(temp_row, 2);
		tempArray_x(k, 1:temp_max) = temp_row;
	end
	A{nBins + 1} =  tempArray_x;
	% handle voronoi cells larger than highestVoronoiCellDensity
	ind = micro_densities(:) >= spacing(nBins + 1);
	tempR = R(ind);
	maxSize = max( cellfun('size', tempR, 2) );
	tempArray_x = NaN( size(tempR, 1), maxSize );
	for k = 1:size(tempR, 1)
		temp_row = tempR{k};
		temp_max = size(temp_row, 2);
		tempArray_x(k, 1:temp_max) = temp_row;
	end
	A{nBins + 2} =  tempArray_x;
	
	% set figure
	h = figure('Name', 'Voronoi diagram', 'Visible', 'off');
	ax = axes(h, 'Color', cmap(1, :) );
	axis equal
	ax.XLim = [xmin xmax];
	ax.YLim = [ymin ymax];
	ax.XLabel.String = 'x in nm';
	ax.YLabel.String = 'y in nm';
	ax.CLimMode = 'manual';
	ax.CLim = [lowestVoronoiCellDensity highestVoronoiCellDensity];
	c = colorbar(ax);
	c.Label.String = "Density in μm^{-"+dim+"}";
	c.Limits = [lowestVoronoiCellDensity highestVoronoiCellDensity];
	hold on
	for i=1:nBins
		ntempColors = repmat(cmap(i, :), size(A{i}, 1), 1);
		patch(ax, 'Faces', A{i}, 'Vertices', V, 'FaceVertexCData', ntempColors,...
			'EdgeColor','none','FaceColor', cmap(i, :));
	end
	% handle voronoi cells smaller than lowestVoronoiCellDensity
	ntempColors = repmat(cmap(1, :), size(A{nBins + 1}, 1), 1);
		patch(ax, 'Faces', A{nBins + 1}, 'Vertices', V, 'FaceVertexCData', ntempColors,...
			'EdgeColor','none','FaceColor', cmap(1, :));
	% handle voronoi cells smaller than highestVoronoiCellDensity
	ntempColors = repmat(cmap(nBins, :), size(A{nBins + 2}, 1), 1);
		patch(ax, 'Faces', A{nBins + 2}, 'Vertices', V, 'FaceVertexCData', ntempColors,...
			'EdgeColor','none','FaceColor', cmap(nBins, :));
	hold off
	h.Visible = 'on';
end
%% output
dataMat(1).voronoi = V; % vertices
dataMat(2).voronoi = R; % regions

if obj.dimension == 2
	dataMat(3).voronoi = 10^(-6).*voronoi_region_size; % area in µm^2
	dataMat(4).voronoi = 10^6.* ( 1./voronoi_region_size ); % density in µm^2
else
	dataMat(3).voronoi = 10^(-9).*voronoi_region_size; % area in µm^3
	dataMat(4).voronoi = 10^9.* ( 1./voronoi_region_size ); % density in µm^3
end
if isRandom == true
	obj.randomClusterStruct = dataMat;
else
	obj.clusterStruct = dataMat;
end
%% UNDER DEVELOPMENT 
% include information on percentage of vertices are saturated
% plotting of different histograms
mdensities_range_idx = lowestVoronoiCellDensity <= micro_densities &  micro_densities <= highestVoronoiCellDensity;
mdensities_range = micro_densities(mdensities_range_idx);

% density histo plots
if showImage == false
	figure('visible', 'on');
	histogram(micro_densities, 40);  % add log scale and bin selection
	set(gca, 'yscale' ,'log')
	title({"Density histogram"; " #voronoi cells: "+length(micro_densities)})
	xlabel("Density in μm^{-"+dim+"}")
	ylabel("Number of voronoi cells")
	
	figure('visible', 'on');
	histogram(mdensities_range, 40);
	set(gca,'yscale', 'log')
	title({"Density histogram 0.5% - 99.5%"; " #voronoi cells: "+length(mdensities_range)})
	xlabel("Density in μm^{-"+dim+"}")
	ylabel("Number of voronoi cells")
end

% area
lowest_plot_area = quantile(micro_area, 0.005);
highest_plot_area = quantile(micro_area, 0.995);
marea_range_idx = lowest_plot_area <= micro_area &  micro_area <= highest_plot_area;
marea_range = micro_area(marea_range_idx);
area_or_vol = "Area";
if dim == 3
	area_or_vol = "Volume";
end

% cumulative plots and area histo plots
if showImage == false
	% cumulative density
	figure('visible', 'on');
	[yvals, xvals] = ecdf(micro_densities);
	semilogx(xvals,yvals,'-r')
	title("Cumulative density plot")
	ylabel("Cumulative fraction")
	xlabel("Voronoi density in μm^{-"+dim+"}")
	
	% area histogram
	figure('visible', 'on');
	histogram(micro_area, 40)
	set(gca,'yscale','log')
	xlabel(area_or_vol+" in μm^{"+dim+"}")
	ylabel("Number of voronoi cells")
	title({area_or_vol+" histogram"; "Type: "+title_class+", #voronoi cells: "+length(micro_area)})
	
	% area histogram - qunatile
	figure('visible', 'on');
	histogram(marea_range,40)
	set(gca,'yscale','log')
	xlabel(area_or_vol+" in μm^{"+dim+"}")
	ylabel("Number of voronoi cells")
	title({area_or_vol+" histogram 0.5% - 99.5%"; "Type: "+title_class+", #voronoi cells: "+length(marea_range)})
end

new_positions = positions(reasonables, 1:2);
pixelSize = 1; % in nm
% griddata not working properly -> temp
if showImage == false
	voronoi_region_size = zeros( nnz(reasonables), 1 );
	for i = 1:numel(R)
		[~, volume] = convhull(V(R{i}, :));
		if isnan(volume) == false && isinf(volume) == false
			voronoi_region_size(i) = volume;
		end
	end
	micro_densities= 10^6.* ( 1./voronoi_region_size );
	[xq,yq] = meshgrid( 1:ceil(x_max/pixelSize), 1:ceil(y_max/pixelSize) );
	xq = xq * pixelSize;
	yq = yq * pixelSize;
	vq = griddata( new_positions(:, 1), new_positions(:, 2), micro_densities(:,1), xq, yq, 'natural');
	vq = flip(vq, 1);
	% set figure
	h = figure('Name', 'Voronoi diagram', 'Visible', 'on');
	ax = axes(h);
	imagesc(ax, vq);
	axis image
	ax.CLimMode = 'manual';
	ax.CLim = [lowestVoronoiCellDensity highestVoronoiCellDensity];
	c = colorbar(ax);
	c.Label.String = "Density in μm^{-"+dim+"}";
	c.Limits = [lowestVoronoiCellDensity highestVoronoiCellDensity];
end
%% UNDER DEVELOPMENT
multiWaitbar('CLOSEALL');
end