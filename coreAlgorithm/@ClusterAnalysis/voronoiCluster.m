function [] = voronoiCluster( obj, varargin )
%voronoiCluster computes Voronoi cells around each data point and uses them
%for density estimation
%
% The current version is implemented for 2D and 3D.
%
% Input:
%   obj:		 The cell object (given by cell.voronoiCluster)
%	varargin     struct containing the following fields in the same order:
%                   isRandom (default: false):
%							true = use random data; false = use experimental data
%                   showImage (default: false):
%							true = shows Voronoi clustering plot, where intensity scales from
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
% Voronoi cell density will be set later on to the default values
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
%% Calculate Voronoi vertices and regions via Delaunay triangulation
% vertices: array that contains the point coordinates of all vertices
% voronoiRegions = cell array where each row represents a Voronoi region,
% hence each row contains the indices of the vertices that can be found in
% the vertices array
if dim == 2
	[vertices, voronoiRegions] = voronoiDiagram( delaunayTriangulation(...
		positions(:, 2), positions(:, 1) ) );
elseif dim == 3
	[vertices, voronoiRegions] = voronoiDiagram( delaunayTriangulation(...
		positions(:, 2) , positions(:, 1), positions(:, 3) ) );
else
	error("Error in voronoi calculation: Function supports only 2D or 3D data!")
end
%% calculate Voronoi areas / volumes
% filter regions that are not valid
xmin = min(positions(:, 2));
xmax = obj.physicalDimension(2);
ymin = min(positions(:, 1));
ymax = obj.physicalDimension(1);
if dim == 3
	zmin = min(positions(:, 3));
	zmax = obj.physicalDimension(3);
end
regionsInsideBoundary = false(numel(voronoiRegions), 1);

for iRegion = 1:numel(voronoiRegions)
	voronoiRegionVertices = vertices(voronoiRegions{iRegion}, :);
	if dim == 2
		isRegionOutsideBoundary = any ( (voronoiRegionVertices(:, 1) <= xmin) | (voronoiRegionVertices(:, 1) >= xmax) |...
			(voronoiRegionVertices(:, 2) <= ymin) | (voronoiRegionVertices(:, 2) >= ymax) );
	else
		isRegionOutsideBoundary = any ( (voronoiRegionVertices(:, 1) <= xmin) | (voronoiRegionVertices(:, 1) >= xmax) |...
			(voronoiRegionVertices(:, 2) <= ymin) | (voronoiRegionVertices(:, 2) >= ymax) |...
			(voronoiRegionVertices(:, 3) <= zmin) | (voronoiRegionVertices(:, 3) >= zmax) );
	end
	if isRegionOutsideBoundary == true
		continue
	end
	regionsInsideBoundary(iRegion) = true;
end

% keep Voronoi regions inside the boundaries
voronoiRegions = voronoiRegions(regionsInsideBoundary);

% calculate size (area / volume ) of voronoi regions
voronoiRegionsSize = NaN( nnz(regionsInsideBoundary), 1 );
for iRegion = 1:numel(voronoiRegions)
	[~, regionSize] = convhull( vertices(voronoiRegions{iRegion}, : ) );
	if isnan(regionSize) == false && isinf(regionSize) == false
		voronoiRegionsSize(iRegion) = regionSize; % in nm^2 / nm^3
	end
end

% TODO work with same unit i.e. add standardized unit converter,
% convert input / output of units to internal format
if dim == 2
	voronoiRegionsDensity = 10^6.* ( 1./voronoiRegionsSize ); % in µm^2
else
	voronoiRegionsDensity = 10^9.* ( 1./voronoiRegionsSize ); % in µm^3
end
%% plotting
if obj.dimension == 3
	warning('Plotting of Voronoi cells is supported for 2D images only!');
end
% set quantile densities to corresponding densities in 1/µm^2
lowestVoronoiCellDensity = quantile(voronoiRegionsDensity, lowestVoronoiCellDensity);
highestVoronoiCellDensity = quantile(voronoiRegionsDensity, highestVoronoiCellDensity);

% works only for 2D
if showImage == true && obj.dimension == 2
	% set color scheme
	nColors = 256;
	colorMap = turbo(nColors);
	
	% find all Voronoi regions that fit into same row of the color map
	voronoiDensityMappedToColorMap = linspace(lowestVoronoiCellDensity, highestVoronoiCellDensity, nColors);
	nDensityBins = numel(voronoiDensityMappedToColorMap) - 1;
	% add two additional for cells that are below or above of
	% lowestVoronoiCellDensity or highestVoronoiCellDensity
	regionsMappedToVoronoiDensity = cell( nDensityBins + 2, 1 );
	for iDensityBin = 1:nDensityBins
		mappedRegions = voronoiRegions( voronoiRegionsDensity(:) >= voronoiDensityMappedToColorMap(iDensityBin) &...
										voronoiRegionsDensity(:) < voronoiDensityMappedToColorMap( iDensityBin + 1 ) );
		% To plot the Voronoi regions, all regions contained in the mapped
		% cell array must be collected in an conventional array first, where
		% places not filled by indices (due to the unequal number of vertices
		% for every Voronoi region) must contain NaN values.
		% These NaN values are later on ignored by the plotting function.
	 
		% get the size of the cell containing the highest number of indices
		% to allocate an NaN array
		maxNumIndices = max( cellfun('size', mappedRegions, 2) );
		tempIndiceArray = NaN( size(mappedRegions, 1), maxNumIndices );
		% replace NaN values with the indices of the Voronoi vertices
		for iMappedRegion = 1:size(mappedRegions, 1)
			mappedRegion = mappedRegions{iMappedRegion};
			numIndices = size(mappedRegion, 2);
			tempIndiceArray(iMappedRegion, 1:numIndices) = mappedRegion;
		end
		regionsMappedToVoronoiDensity{iDensityBin} =  tempIndiceArray;
	end
	% handle Voronoi cells smaller than the lowestVoronoiCellDensity
	mappedRegions = voronoiRegions( voronoiRegionsDensity(:) < voronoiDensityMappedToColorMap(1) );
	maxNumIndices = max( cellfun('size', mappedRegions, 2) );
	tempIndiceArray = NaN( size(mappedRegions, 1), maxNumIndices );
	for iMappedRegion = 1:size(mappedRegions, 1)
		mappedRegion = mappedRegions{iMappedRegion};
		numIndices = size(mappedRegion, 2);
		tempIndiceArray(iMappedRegion, 1:numIndices) = mappedRegion;
	end
	regionsMappedToVoronoiDensity{nDensityBins + 1} =  tempIndiceArray;
	% handle Voronoi cells larger than the highestVoronoiCellDensity
	mappedRegions = voronoiRegions( voronoiRegionsDensity(:) >= voronoiDensityMappedToColorMap(nDensityBins + 1) );
	maxNumIndices = max( cellfun('size', mappedRegions, 2) );
	tempIndiceArray = NaN( size(mappedRegions, 1), maxNumIndices );
	for iMappedRegion = 1:size(mappedRegions, 1)
		mappedRegion = mappedRegions{iMappedRegion};
		numIndices = size(mappedRegion, 2);
		tempIndiceArray(iMappedRegion, 1:numIndices) = mappedRegion;
	end
	regionsMappedToVoronoiDensity{nDensityBins + 2} =  tempIndiceArray;
	
	% set figure
	h = figure('Name', 'Voronoi diagram', 'Visible', 'off');
	ax = axes(h, 'Color', colorMap(1, :) );
	% reverse axis to be aligned with original image
	ax.YDir = 'reverse';
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
	for iRegion=1:nDensityBins
		ntempColors = repmat(colorMap(iRegion, :), size(regionsMappedToVoronoiDensity{iRegion}, 1), 1);
		patch(ax, 'Faces', regionsMappedToVoronoiDensity{iRegion}, 'Vertices', vertices, 'FaceVertexCData', ntempColors,...
			'EdgeColor','none','FaceColor', colorMap(iRegion, :));
	end
	% handle voronoi cells smaller than lowestVoronoiCellDensity
	ntempColors = repmat(colorMap(1, :), size(regionsMappedToVoronoiDensity{nDensityBins + 1}, 1), 1);
	patch(ax, 'Faces', regionsMappedToVoronoiDensity{nDensityBins + 1}, 'Vertices', vertices, 'FaceVertexCData', ntempColors,...
		'EdgeColor','none','FaceColor', colorMap(1, :));
	% handle voronoi cells smaller than highestVoronoiCellDensity
	ntempColors = repmat(colorMap(nDensityBins, :), size(regionsMappedToVoronoiDensity{nDensityBins + 2}, 1), 1);
	patch(ax, 'Faces', regionsMappedToVoronoiDensity{nDensityBins + 2}, 'Vertices', vertices, 'FaceVertexCData', ntempColors,...
		'EdgeColor','none','FaceColor', colorMap(nDensityBins, :));
	hold off
	h.Visible = 'on';
end
%% output
dataMat(1).voronoiVertices = vertices;
dataMat(2).voronoiRegions = voronoiRegions;

if obj.dimension == 2
	dataMat(3).voronoiArea = 10^(-6).*voronoiRegionsSize; % area in µm^2
	dataMat(4).voronoiDensity = 10^6.* ( 1./voronoiRegionsSize ); % density in µm^2
else
	dataMat(3).voronoiVolume = 10^(-9).*voronoiRegionsSize; % volume in µm^3
	dataMat(4).voronoiDensity = 10^9.* ( 1./voronoiRegionsSize ); % density in µm^3
end
if isRandom == true
	obj.randomClusterStruct = dataMat;
else
	obj.clusterStruct = dataMat;
end
multiWaitbar('CLOSEALL');
end
% %% UNDER DEVELOPMENT
% % include information on percentage of vertices are saturated
% % plotting of different histograms
% mdensities_range_idx = lowestVoronoiCellDensity <= micro_densities &  micro_densities <= highestVoronoiCellDensity;
% mdensities_range = micro_densities(mdensities_range_idx);
%
% % density histo plots
% if showImage == false
% 	figure('visible', 'on');
% 	histogram(micro_densities, 40);  % add log scale and bin selection
% 	set(gca, 'yscale' ,'log')
% 	title({"Density histogram"; " #voronoi cells: "+length(micro_densities)})
% 	xlabel("Density in μm^{-"+dim+"}")
% 	ylabel("Number of voronoi cells")
%
% 	figure('visible', 'on');
% 	histogram(mdensities_range, 40);
% 	set(gca,'yscale', 'log')
% 	title({"Density histogram 0.5% - 99.5%"; " #voronoi cells: "+length(mdensities_range)})
% 	xlabel("Density in μm^{-"+dim+"}")
% 	ylabel("Number of voronoi cells")
% end
%
% % area
% lowest_plot_area = quantile(micro_area, 0.005);
% highest_plot_area = quantile(micro_area, 0.995);
% marea_range_idx = lowest_plot_area <= micro_area &  micro_area <= highest_plot_area;
% marea_range = micro_area(marea_range_idx);
% area_or_vol = "Area";
% if dim == 3
% 	area_or_vol = "Volume";
% end
%
% % cumulative plots and area histo plots
% if showImage == false
% 	% cumulative density
% 	figure('visible', 'on');
% 	[yvals, xvals] = ecdf(micro_densities);
% 	semilogx(xvals,yvals,'-r')
% 	title("Cumulative density plot")
% 	ylabel("Cumulative fraction")
% 	xlabel("Voronoi density in μm^{-"+dim+"}")
%
% 	% area histogram
% 	figure('visible', 'on');
% 	histogram(micro_area, 40)
% 	set(gca,'yscale','log')
% 	xlabel(area_or_vol+" in μm^{"+dim+"}")
% 	ylabel("Number of voronoi cells")
% 	title({area_or_vol+" histogram"; "Type: "+title_class+", #voronoi cells: "+length(micro_area)})
%
% 	% area histogram - qunatile
% 	figure('visible', 'on');
% 	histogram(marea_range,40)
% 	set(gca,'yscale','log')
% 	xlabel(area_or_vol+" in μm^{"+dim+"}")
% 	ylabel("Number of voronoi cells")
% 	title({area_or_vol+" histogram 0.5% - 99.5%"; "Type: "+title_class+", #voronoi cells: "+length(marea_range)})
% end
%
% new_positions = positions(regionsInsideBoundary, 1:2);
% pixelSize = 1; % in nm
% % griddata not working properly -> temp
% if showImage == false
% 	voronoiRegionSize = zeros( nnz(regionsInsideBoundary), 1 );
% 	for iRegion = 1:numel(voronoiRegions)
% 		[~, size] = convhull(vertices(voronoiRegions{iRegion}, :));
% 		if isnan(size) == false && isinf(size) == false
% 			voronoiRegionSize(iRegion) = size;
% 		end
% 	end
% 	micro_densities= 10^6.* ( 1./voronoiRegionSize );
% 	[xq,yq] = meshgrid( 1:ceil(x_max/pixelSize), 1:ceil(y_max/pixelSize) );
% 	xq = xq * pixelSize;
% 	yq = yq * pixelSize;
% 	vq = griddata( new_positions(:, 1), new_positions(:, 2), micro_densities(:,1), xq, yq, 'natural');
% 	vq = flip(vq, 1);
% 	% set figure
% 	h = figure('Name', 'Voronoi diagram', 'Visible', 'on');
% 	ax = axes(h);
% 	imagesc(ax, vq);
% 	axis image
% 	ax.CLimMode = 'manual';
% 	ax.CLim = [lowestVoronoiCellDensity highestVoronoiCellDensity];
% 	c = colorbar(ax);
% 	c.Label.String = "Density in μm^{-"+dim+"}";
% 	c.Limits = [lowestVoronoiCellDensity highestVoronoiCellDensity];
% end
%% UNDER DEVELOPMENT