function [] = DBSCAN(obj, radius, minPoints, isRandom, maxDiameter, showImage)
%DBSCAN for details see A Density-Based Algorithm for Discovering Clusters in Large
%   Spatial Databases with Noise. Ester at al. 1996
%
% The current version is implemented for 2D + 3D data sets.
%
% Input:
%   radius: radius by which a certain number of points (minPoints) must be,
%           to accept the point as part of a cluster
%   minPoints: minimum number of points within specified distance (radius)
%           to fullfill cluster criteria
%   isRandom: true = use random data; false = use experimental data
%   maxDiameter: maximum diameter of a cluster, clusters larger than
%           maxDiameter will be removed (default: maximum dimension of the image)
%   showImage: true = show scatterplot, noise points and filtered clusters
%           are displayed in black, clusters are displayed with a random
%           colour. The black cross indicates the center of mass of the
%           clusters.
%
% Output
%   (1) total number of points that are part of a cluster
%   (2) obj.clusterStruct or obj.randomClusterStruct (depending on value of isRandom)
%   containing the clusters and thier properties
%   (3) clusterAssignment - every point in the localization data set is
%   assigned to noise (value = 0) or to a cluster (removed cluster - see maxDiameter
%   will be set to 0 as well)
%
% Note:
% The current version removes clusters, which consists of two or less
% points. Cluster touching the border of an image will be not removed.
%
% Matlab version: 2014b or newer (because of boundary command).
%
% Cluster areas are calculated using Matlab's boundary command.
% Diameter/Volume property assumes more or less circular clusters.
%
% Jan Neumann, 27.06.2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% default parameters
visMethod = 'gaussianBlur'; % alternative: 'histogramBinning'
visPixelsize = 10;
%% init
if nargin < 3
	warning("You have to specify at least the following parameters for DBSCAN: radius, minimum number of points.")
	return
end
switch nargin
	case 3
		isRandom = 0;
		if isRandom == true
			maxDiameter = max(obj.randomPhysicalDimension);
		else
			maxDiameter = max(obj.physicalDimension);
		end
	case 4
		if isRandom == true
			maxDiameter = max(obj.randomPhysicalDimension);
		else
			maxDiameter = max(obj.physicalDimension);
		end
		showImage = 0;
	case 5
		showImage = 0;
	case 6
		
	otherwise
		error('Wrong number of input arguments!')
end
if isRandom == true
	positions = obj.randomTable;
	dataMat = obj.randomClusterStruct;
	nPoints = obj.randomNoOfPoints;
	treeData = obj.randomQueryData;
	flagSearchTree = obj.flagSearchTreeRandomData;
else
	positions = obj.positionTable;
	dataMat = obj.clusterStruct;
	nPoints = obj.NoOfPoints;
	treeData = obj.queryData;
	flagSearchTree = obj.flagSearchTree;
end
% initialize output fields of dataMat
dataMat(1).clusterDBSCAN = [];
dataMat(2).clusterDBSCAN = [];
dataMat(3).clusterDBSCAN = [];
%% DBSCAN algorithm
clusterID = 0; % number labeling of cluster
% initialize cluster matrix - all points are at the beginning
% unassigned []
clusterAssignment = zeros(nPoints, 1);
for ii = 1:nPoints
	clusterAssignment(ii, 1) = -2;
end
multiWaitbar('DBSCAN algorithm running...', 0);
% for waitbar
prevPercent = 0;
for ii = 1:nPoints
	if ~(clusterAssignment(ii, 1) == -2)
		% point is already assigned
		continue
	end
	clusterAssignment(ii, 1) = -1; % point is visited
	switch flagSearchTree
		case 'cTree'
			distanceAppend = kdtree_ball_query(treeData, positions(ii, 1:obj.dimension), radius);
		case 'matlabTree'
			distance = rangesearch(positions(ii, 1:obj.dimension), treeData.X, radius);
			distanceAppend = (find(~cellfun('isempty', distance))).';
	end
	if numel(distanceAppend) < minPoints + 1 % + 1 ignore current point
		clusterAssignment(ii, 1) = 0; % point is noise
	else
		clusterID = clusterID + 1;
		clusterAssignment(ii, 1) = clusterID;
		indice = 1;
		while indice < numel(distanceAppend)
			ll = distanceAppend(indice);
			indice = indice + 1;
			if clusterAssignment(ll, 1) == -2 % if point is not visited
				clusterAssignment(ll, 1) = -1;
				switch flagSearchTree
					case 'cTree'
						newAppend = kdtree_ball_query(treeData, positions(ll, 1:obj.dimension), radius);
					case 'matlabTree'
						newDistance = rangesearch(positions(ll, 1:obj.dimension), treeData.X, radius);
						newAppend = (find(~cellfun('isempty', newDistance))).';
				end
				if numel(newAppend) >= minPoints + 1 % + 1
					distanceAppend = cat(1, distanceAppend, newAppend); % do not implement union or unique - this may change the order of the array
				end
			end
			if clusterAssignment(ll, 1) == -1 % not member of a cluster
				clusterAssignment(ll, 1) = clusterID;
			end
		end
	end
	currentPercent = fix(100*sum(clusterAssignment~=-2)/nPoints);
	if currentPercent > prevPercent
		multiWaitbar( 'DBSCAN algorithm running...', 'Value', sum(clusterAssignment~=-2)/nPoints);
		prevPercent = currentPercent;
	end
end
%% calculate cluster properties
% for storing boundaries
% TODO check if clusters corresponding to a boundary are removed -> remove
% boundary as well
k_boundary = cell(1, max(clusterAssignment));
for kk = 1:max(clusterAssignment)
	dataMat(2).clusterDBSCAN(kk).PointsCluster = find(clusterAssignment(:, 1) == kk);
	temp_cluster = [dataMat(2).clusterDBSCAN(kk).PointsCluster];
	% number of molecules per cluster
	dataMat(2).clusterDBSCAN(kk).Molecules = size(temp_cluster, 1);
	% extract positions of clustered points
	x = double(positions(temp_cluster, 1));
	y = double(positions(temp_cluster, 2));
	% calculate center of mass
	dataMat(2).clusterDBSCAN(kk).xCoM = 1/length(temp_cluster)*sum(positions(temp_cluster, 1));
	dataMat(2).clusterDBSCAN(kk).yCoM = 1/length(temp_cluster)*sum(positions(temp_cluster, 2));
	if obj.dimension == 2
		[k_boundary{kk}, dataMat(2).clusterDBSCAN(kk).Area] = boundary(x, y);
		% calculate diameter of clusters (assume that clusters are circular)
		dataMat(2).clusterDBSCAN(kk).Diameter = 2.*sqrt(dataMat(2).clusterDBSCAN(kk).Area./pi); % in nm
		% calculate area in µm^2
		dataMat(2).clusterDBSCAN(kk).Area = dataMat(2).clusterDBSCAN(kk).Area*(1/1000).^2; % in µm^2
	else
		% extract z position
		z = double(positions(temp_cluster, 3));
		% calculate center of mass in z direction
		dataMat(2).clusterDBSCAN(kk).zCoM = 1/length(temp_cluster)*sum(positions(temp_cluster, 3));
		[k_boundary{kk}, dataMat(2).clusterDBSCAN(kk).Volume] = boundary(x, y, z);
		% calculate diameter of clusters (assume that clusters are spherical)
		dataMat(2).clusterDBSCAN(kk).Diameter = nthroot(6.* dataMat(2).clusterDBSCAN(kk).Volume./pi , 3); % in nm
		% calculate volume in µm^3
		dataMat(2).clusterDBSCAN(kk).Volume = dataMat(2).clusterDBSCAN(kk).Volume*(1/1000).^3; % in µm^3
	end
end

% filter cluster according to diameter and calculate with remaining cluster additional properties
if ~isempty(dataMat(2).clusterDBSCAN)
	% filter data
	deleteInd = [dataMat(2).clusterDBSCAN.Diameter] > maxDiameter;
	for ii = 1:size(dataMat(2).clusterDBSCAN, 2)
		if deleteInd(ii)
			temp_cluster = clusterAssignment(:, 1) == ii;
			clusterAssignment(temp_cluster, 1) = 0;
		end
	end
	k_boundary([dataMat(2).clusterDBSCAN.Diameter] > maxDiameter) = [];
	dataMat(2).clusterDBSCAN([dataMat(2).clusterDBSCAN.Diameter] > maxDiameter) = [];
	% filter clusters which consists of two or less points
	deleteInd = [dataMat(2).clusterDBSCAN.Diameter] == 0;
	for ii = 1:size(dataMat(2).clusterDBSCAN, 2)
		if deleteInd(ii)
			temp_cluster = clusterAssignment(:, 1) == ii;
			clusterAssignment(temp_cluster, 1) = 0;
		end
	end
	k_boundary([dataMat(2).clusterDBSCAN.Diameter] == 0) = [];
	dataMat(2).clusterDBSCAN([dataMat(2).clusterDBSCAN.Diameter] == 0) = [];
	% shift cluster number in clusterAssignment, so starting from one in
	% ascending order
	counter = 0;
	for ii = 1:max(clusterAssignment)
		tempVal = clusterAssignment(:, 1) == ii;
		if any(tempVal)
			clusterAssignment(tempVal, 1) = counter + 1;
			counter = counter + 1;
		end
	end
	
	% calculate distance to next neighboring cluster
	if obj.dimension == 2
		[~, distance] = knnsearch(cat(2, vertcat(dataMat(2).clusterDBSCAN.xCoM), vertcat(dataMat(2).clusterDBSCAN.yCoM)),...
			cat(2, vertcat(dataMat(2).clusterDBSCAN.xCoM), vertcat(dataMat(2).clusterDBSCAN.yCoM)), 'k', 2);
	else
		[~, distance] = knnsearch(cat(2, vertcat(dataMat(2).clusterDBSCAN.xCoM), vertcat(dataMat(2).clusterDBSCAN.yCoM), vertcat(dataMat(2).clusterDBSCAN.zCoM)),...
			cat(2, vertcat(dataMat(2).clusterDBSCAN.xCoM), vertcat(dataMat(2).clusterDBSCAN.yCoM), vertcat(dataMat(2).clusterDBSCAN.zCoM)), 'k', 2);
	end
	% calculate total number of molecules in cluster
	dataMat(1).clusterDBSCAN = sum([dataMat(2).clusterDBSCAN.Molecules]);
	for ii = 1:size(dataMat(2).clusterDBSCAN, 2)
		if ~isempty(distance) && ~isequal(distance, 0)
			dataMat(2).clusterDBSCAN(ii).NNCluster = distance(ii, 2);
		end
		% check for two points and if Area = zero
		if ~isempty(dataMat(2).clusterDBSCAN)
			if obj.dimension == 2
				dataMat(2).clusterDBSCAN(ii).ClusterDensity = dataMat(2).clusterDBSCAN(ii).Molecules./dataMat(2).clusterDBSCAN(ii).Area; % molecules per µm^2
			else
				dataMat(2).clusterDBSCAN(ii).ClusterDensity = dataMat(2).clusterDBSCAN(ii).Molecules./dataMat(2).clusterDBSCAN(ii).Volume; % molecules per µm^3
			end
		end
	end
else
	% no points in cluster
	dataMat(1).clusterDBSCAN = 0;
end
% save cluster assignment for distance calculation (to calculate distances from points inside and outside of clusters)
dataMat(3).clusterDBSCAN = uint32(clusterAssignment);
multiWaitbar( 'DBSCAN algorithm running...', 'Close' );
%% visualization
% plots all clusters, filtered clusters are considered as noise
if showImage == true && any(clusterAssignment)
	% for waitbar
	multiWaitbar('Generating DBSCAN image. Please wait!', 'busy');
	% perform calculations before plotting
	% for 2D data overlay, create histogram binned image and overlay clustered points
	if obj.dimension == 2
		clusterImage = visModuleCluster(positions, visMethod, visPixelsize);
	end	
	% do plotting
	fig = figure( 'Name', num2str(obj.sampleID), 'visible', 'off');
	if obj.dimension == 2
		ax1 = subplot(1, 2, 1); % tiledlayout(m,n)
		imshow(clusterImage, 'Colormap', hot);
		hbar = colorbar;
		ylabel(hbar, 'intensity [arb. unit]')
		axis image
		xlabel(ax1, ['x', num2str(visPixelsize), ' nm'])
		ylabel(ax1, ['x', num2str(visPixelsize), ' nm'])
		title('Gaussian blurred image')
		subplot(1, 2, 2);
	end
	hold on
	n_cluster = max(clusterAssignment);
	color_array = rand(n_cluster, 3);
	% plot cluster points and draw boundary
	pl_handle_boundary = gobjects(n_cluster, 1);
	pl_handle_cluster_points = gobjects(n_cluster, 1);
	for i = 1 : n_cluster
		temp_cluster = dataMat(2).clusterDBSCAN(i).PointsCluster;
		if obj.dimension == 2
			% extract positions of clustered points for plotting
			x_pos = ( positions(temp_cluster, 1) - min(positions(:, 1)) )./visPixelsize;
			y_pos = ( positions(temp_cluster, 2) - min(positions(:, 2)) )./visPixelsize;
			pl_handle_cluster_points(i, 1) = line( y_pos, x_pos, 'Color', color_array(i,:), 'LineStyle', 'None', 'Marker', '.' );
			% extract border points
			x_pos = x_pos(k_boundary{i});
			y_pos = y_pos(k_boundary{i});
			pl_handle_boundary(i, 1) = line( y_pos, x_pos, 'Color', color_array(i,:) );
			% pl_handle_boundary(i, 1) = fill(y_pos, x_pos, color_array(i,:), 'EdgeColor', color_array(i,:));
		else
			view(3);
			pan off
			% extract positions of clustered points for plotting
			x_pos = ( positions(temp_cluster, 1) - min(positions(:, 1)) )./visPixelsize;
			y_pos = ( positions(temp_cluster, 2) - min(positions(:, 2)) )./visPixelsize;
			z_pos = ( positions(temp_cluster, 3) - min(positions(:, 3)) )./visPixelsize;
			pl_handle_cluster_points(i, 1) = line( y_pos, x_pos, z_pos, 'Color', color_array(i,:), 'LineStyle', 'None', 'Marker', '.' );
			% extract border points
			temp_k = k_boundary{i}(:, 1);
			k_boundary{i}(:, 1) = k_boundary{i}(:, 2);
			k_boundary{i}(:, 2) = temp_k;
			pl_handle_boundary(i, 1) = trisurf(k_boundary{i}, y_pos, x_pos, z_pos, 'Facecolor', color_array(i,:),...
				'FaceAlpha', 0.4, 'EdgeColor', color_array(i,:), 'EdgeAlpha', 0.4);
		end
	end
	temp_cluster = find(clusterAssignment(:, 1) == 0);
	if obj.dimension == 2
		% add noise points
		pl_noise_points = scatter(double((positions(temp_cluster, 2)) - min(positions(:, 2)) )./visPixelsize,...
			double((positions(temp_cluster, 1)) - min(positions(:, 1)) )./visPixelsize, 50, 'k.');
		%add center of mass
		pl_center_mass = scatter(double(([dataMat(2).clusterDBSCAN.yCoM]) - min(positions(:, 2)))./visPixelsize,...
			double(([dataMat(2).clusterDBSCAN.xCoM]) - min(positions(:, 1)))./visPixelsize, 150, 'kx');
	else
		% add noise points
		pl_noise_points = scatter3(double((positions(temp_cluster, 2)) - min(positions(:, 2)) )./visPixelsize,...
			double((positions(temp_cluster, 1)) - min(positions(:, 1)) )./visPixelsize,...
			double((positions(temp_cluster, 3)) - min(positions(:, 3)) )./visPixelsize, 50, 'k.');
		% plot center of mass
		pl_center_mass = scatter3(double(([dataMat(2).clusterDBSCAN.yCoM]) - min(positions(:, 2)))./visPixelsize,...
			double(([dataMat(2).clusterDBSCAN.xCoM]) - min(positions(:, 1)))./visPixelsize,...
			double(([dataMat(2).clusterDBSCAN.zCoM]) - min(positions(:, 3)))./visPixelsize, 150, 'kx');
	end
	title('Detected clusters and noise points (filtered points are assigned as noise)')
	set(gca, 'Ydir', 'reverse')
	axis image
	xlim([min(positions(:, 2))./visPixelsize, max(positions(:, 2))./visPixelsize]); 
	ylim([min(positions(:, 1))./visPixelsize, max(positions(:, 1))./visPixelsize]);
	xlabel(['x', num2str(visPixelsize), ' nm'])
	ylabel(['x', num2str(visPixelsize), ' nm'])
	if obj.dimension == 3
		zlabel(['x', num2str(visPixelsize), ' nm'])
	end
	hold off
	set(fig, 'Visible', 'on');
	
	% for interactive plotting - enter datacursor mode
	dcm_obj = datacursormode(fig);
	pos_list = ( (positions(:, 1:2)- min(positions(:, 1:2)))./visPixelsize );
	posCoM_list(:, 1) = (( [dataMat(2).clusterDBSCAN.xCoM].' - min(positions(:, 1)) )./visPixelsize);
	posCoM_list(:, 2) = (( [dataMat(2).clusterDBSCAN.yCoM].' - min(positions(:, 2)) )./visPixelsize);
	set(dcm_obj, 'UpdateFcn', {@myupdatefcn, dataMat, clusterAssignment, pos_list, posCoM_list, obj.dimension});
	hold off
end
%% output
if isRandom == true
	obj.randomClusterStruct = dataMat;
else
	obj.clusterStruct = dataMat;
end
multiWaitbar('CLOSEALL');
end

function txt = myupdatefcn(~, event_obj, dataMat, clusterAssignment, pos_list, posCoM_list, dim)
% Customizes text of data tips
pos = get(event_obj, 'Position');
idx_1 = find( pos_list(:, 1) == pos(2));
idx_2 = find( pos_list(:, 2) == pos(1));
if idx_1 == idx_2
	if clusterAssignment(idx_1) == 0
		txt = {'Noise Point',...
			['xPos: ', num2str(pos(1))],...
			['yPos: ', num2str(pos(2))]};
	else
		if dim == 2
			txt = {['nMolecules: ', num2str( dataMat(2).clusterDBSCAN(clusterAssignment(idx_1)).Molecules )],...
				['Diameter [nm]: ', num2str( dataMat(2).clusterDBSCAN(clusterAssignment(idx_1)).Diameter )],...
				['Area [µm^2]: ', num2str(dataMat(2).clusterDBSCAN(clusterAssignment(idx_1)).Area)],...
				['Next-Cluster-Distance [nm]: ', num2str(dataMat(2).clusterDBSCAN(clusterAssignment(idx_1)).NNCluster)],...
				['Cluster Density [µm^2]: ', num2str(dataMat(2).clusterDBSCAN(clusterAssignment(idx_1)).ClusterDensity)]};
		else
			txt = {['nMolecules: ', num2str( dataMat(2).clusterDBSCAN(clusterAssignment(idx_1)).Molecules )],...
				['Diameter [nm]: ', num2str( dataMat(2).clusterDBSCAN(clusterAssignment(idx_1)).Diameter )],...
				['Volume [µm^3]: ', num2str(dataMat(2).clusterDBSCAN(clusterAssignment(idx_1)).Volume)],...
				['Next-Cluster-Distance [nm]: ', num2str(dataMat(2).clusterDBSCAN(clusterAssignment(idx_1)).NNCluster)],...
				['Cluster Density [µm^3]: ', num2str(dataMat(2).clusterDBSCAN(clusterAssignment(idx_1)).ClusterDensity)]};
		end
	end
elseif pos(2) == posCoM_list(posCoM_list(:, 1) == pos(2))
	idx_1 = find(posCoM_list(:, 1) == pos(2));
	if dim == 2
		txt = {['nMolecules: ', num2str( dataMat(2).clusterDBSCAN(idx_1).Molecules )],...
			['Diameter [nm]: ', num2str( dataMat(2).clusterDBSCAN(idx_1).Diameter )],...
			['Area [µm^2]: ', num2str(dataMat(2).clusterDBSCAN(idx_1).Area)],...
			['Next-Cluster-Distance [nm]: ', num2str(dataMat(2).clusterDBSCAN(idx_1).NNCluster)],...
			['Cluster Density [µm^2]: ', num2str(dataMat(2).clusterDBSCAN(idx_1).ClusterDensity)]};
	else
		txt = {['nMolecules: ', num2str( dataMat(2).clusterDBSCAN(idx_1).Molecules )],...
			['Diameter [nm]: ', num2str( dataMat(2).clusterDBSCAN(idx_1).Diameter )],...
			['Volume [µm^3]: ', num2str(dataMat(2).clusterDBSCAN(idx_1).Volume)],...
			['Next-Cluster-Distance [nm]: ', num2str(dataMat(2).clusterDBSCAN(idx_1).NNCluster)],...
			['Cluster Density [µm^3]: ', num2str(dataMat(2).clusterDBSCAN(idx_1).ClusterDensity)]};
	end
else
	txt = {'Undefined Error'};
end

end
