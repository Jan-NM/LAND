function [] = clusterBatchProcessing(Clusterparamstruct)
%clusterBatchProcessing batch processes multiple localization files for
%cluster analysis
%
%
% Jan Neumann, 23.07.17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% init parameters
global files
% file and sampleID handling
[~, FILESmax] = size(files);
if FILESmax < 10
    jm = ['0' num2str(FILESmax)];
else
    jm = num2str(FILESmax);
end
pathname = [Clusterparamstruct.DIR_input,filesep];
%% start batch processing
clusterData = struct([]);
% loop through files
for ii = 1:FILESmax
    if ii < 10
        js = ['0' num2str(ii)];    
    else
        js = num2str(ii);
    end
    Clusterparamstruct.displayStatus(['running Evaluation of File No. ' js ' of ' jm '...']);
    filename = char(files(ii));
    sampleID = filename(1:end-4);
    if isequal(filename,0)||isequal(pathname,0)
        Clusterparamstruct.displayStatus('Error - Evaluation not successful');
        Clusterparamstruct.displayError('No File Selected');
        Clusterparamstruct.enableStartButton();
        error('No File Selected');
    end
    % load file for evaluation
    localizationFile = load([pathname filename]);
    oldField = char(fieldnames(localizationFile));
    if ~(strcmp(oldField, 'Orte'))
        newField = 'Orte';
        [localizationFile.(newField)] = localizationFile.(oldField);
        localizationFile = rmfield(localizationFile, oldField);   
    end
    failedFiles = []; % for error handling
    try
        clusterData(ii).sampleID = sampleID;
        clusterData(ii).Analysis = ClusterAnalysis(localizationFile.Orte, sampleID,...
            Clusterparamstruct.dimension, Clusterparamstruct.densityAlgorithm, ...
            Clusterparamstruct.randomAlgorithmValue, Clusterparamstruct.randomAlgorithm);       
        % process corresponding analysis programs, always run DBSCAN before
        % distance calculations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DBSCAN
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Clusterparamstruct.DBSCAN.algorithm == true
            clusterData(ii).Analysis.DBSCAN(Clusterparamstruct.DBSCAN.radius, Clusterparamstruct.DBSCAN.minPoint,...
                0, Clusterparamstruct.DBSCAN.maxDiameter, Clusterparamstruct.showPlots);
            if Clusterparamstruct.compareRandomData == true
                clusterData(ii).Analysis.DBSCAN(Clusterparamstruct.DBSCAN.radius, Clusterparamstruct.DBSCAN.minPoint,...
                    1, Clusterparamstruct.DBSCAN.maxDiameter, Clusterparamstruct.showPlots);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Nearest Neighbor Analysis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Clusterparamstruct.kNNDistance.algorithm == true
            clusterData(ii).Analysis.kNNDistance(Clusterparamstruct.kNNDistance.k, 0,...
                Clusterparamstruct.kNNDistance.maxDistance, Clusterparamstruct.showPlots);
            if Clusterparamstruct.compareRandomData == true
                clusterData(ii).Analysis.kNNDistance(Clusterparamstruct.kNNDistance.k, 1,...
                    Clusterparamstruct.kNNDistance.maxDistance, Clusterparamstruct.showPlots);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Distance Analysis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Clusterparamstruct.Distance.algorithm == true
            clusterData(ii).Analysis.distanceAnalysis(Clusterparamstruct.Distance.maxDistance, 0,...
                Clusterparamstruct.showPlots);
            if Clusterparamstruct.compareRandomData == true
                clusterData(ii).Analysis.distanceAnalysis(Clusterparamstruct.Distance.maxDistance, 1,...
                    Clusterparamstruct.showPlots);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Radial Density Function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Clusterparamstruct.RDF.algorithm == true
            clusterData(ii).Analysis.radialDensityFunction(Clusterparamstruct.RDF.binSize,...
                Clusterparamstruct.RDF.maxDistance, 0, Clusterparamstruct.showPlots);
            if Clusterparamstruct.compareRandomData == true
                clusterData(ii).Analysis.radialDensityFunction(Clusterparamstruct.RDF.binSize,...
                Clusterparamstruct.RDF.maxDistance, 1, Clusterparamstruct.showPlots);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Ripley's Function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Clusterparamstruct.ripley.algorithm == true
            clusterData(ii).Analysis.ripley(Clusterparamstruct.ripley.radius,...
                Clusterparamstruct.ripley.maxDistance, 0, Clusterparamstruct.showPlots);
            if Clusterparamstruct.compareRandomData == true
                clusterData(ii).Analysis.ripley(Clusterparamstruct.ripley.radius,...
                Clusterparamstruct.ripley.maxDistance, 1, Clusterparamstruct.showPlots);
            end
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compaction Parameter
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if Clusterparamstruct.ccp.algorithm == true
            clusterData(ii).Analysis.compactionParameter(Clusterparamstruct.ccp.parameters);
            if Clusterparamstruct.compareRandomData == true
                warning('Random data are not supported for compaction parameter analysis.');
            end
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Voronoi
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if Clusterparamstruct.voronoi.algorithm == true
			Clusterparamstruct.voronoi.parameters.isRandom = false;
            clusterData(ii).Analysis.voronoiCluster(Clusterparamstruct.voronoi.parameters);
            if Clusterparamstruct.compareRandomData == true
				Clusterparamstruct.voronoi.parameters.isRandom = true;
                clusterData(ii).Analysis.voronoiCluster(Clusterparamstruct.voronoi.parameters);
            end
		end
    catch ME
        disp(['An error occured during evaluation of file: ' sampleID'.']);
        disp(getReport(ME,'extended'));
        fprintf('\n');
        disp('Continue execution...');
        failedFiles = [failedFiles  sampleID ', ' ];
        Clusterparamstruct.displayError(['Could not evaluate files: ' failedFiles]);
    end
end
% remove entries which are empty
keep = []; 
for ii = 1:FILESmax
    if ~isfield(clusterData(ii),'Analysis') || isempty(clusterData(ii).Analysis)
        continue
    else
        keep = [keep ii];
    end
end
if isempty(keep)
    disp('Data is empty. Can not proceed.')
    return
end
clusterData = clusterData(keep);
%% create results report and plots
multiWaitbar( 'Finishing computation! Please wait a few mintues.', 'Busy');
% Nearest Neighbor analysis
visKNNDistances(Clusterparamstruct, clusterData)
% DBSCAN algorithm
visDBSCAN(Clusterparamstruct, clusterData)
% distance analysis
visDistances(Clusterparamstruct, clusterData)
% radial density function
visRDF(Clusterparamstruct, clusterData)
% Ripley's function
visRipley(Clusterparamstruct, clusterData)
%% save data
Clusterparamstruct.ClusterPathname = [Clusterparamstruct.DIR_output filesep];
Clusterparamstruct.ClusterFilename = [Clusterparamstruct.treatmentName, '.mat'];
treatmentName.(Clusterparamstruct.treatmentName) = clusterData;
save([Clusterparamstruct.ClusterPathname Clusterparamstruct.ClusterFilename], '-struct', 'treatmentName');
savefig(Clusterparamstruct.handles.figure1, [Clusterparamstruct.ClusterPathname 'settings'])
multiWaitbar('CLOSEALL');
end