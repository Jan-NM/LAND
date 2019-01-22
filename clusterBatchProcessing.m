function [] = clusterBatchProcessing(Clusterparamstruct)
%clusterBatchProcessing batch processes multiple localization files for
%cluster analysis
%   Detailed explanation goes here
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
        % grid analysis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Clusterparamstruct.gridAna.algorithm == true
            clusterData(ii).Analysis.gridAnalysis(10, 0);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nearest Neighbor analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Clusterparamstruct.kNNDistance.algorithm == true
    totalkNNDistance = [];
    totalkNNDistanceRandom = [];
    totalInnerkNNDistance = [];
    totalOuterkNNDistance = [];
    totalInnerkNNDistanceRandom = [];
    totalOuterkNNDistanceRandom = [];
    for kk=1:size(clusterData, 2)
        totalkNNDistance = [clusterData(kk).Analysis.clusterStruct(2).kNNDistance; totalkNNDistance];
    end
    % if DBSCAN was done before
    if isfield(clusterData(kk).Analysis.clusterStruct, 'clusterDBSCAN')
        for kk=1:size(clusterData, 2)
            totalInnerkNNDistance = [clusterData(kk).Analysis.clusterStruct(3).kNNDistance; totalInnerkNNDistance];
            totalOuterkNNDistance = [clusterData(kk).Analysis.clusterStruct(4).kNNDistance; totalOuterkNNDistance];
        end
    end
    % for random data
    if Clusterparamstruct.compareRandomData == true
        for kk=1:size(clusterData, 2)
            totalkNNDistanceRandom = [clusterData(kk).Analysis.randomClusterStruct(2).kNNDistance; totalkNNDistanceRandom];
        end
    end
    % if DBSCAN was done before
    if isfield(clusterData(kk).Analysis.randomClusterStruct, 'clusterDBSCAN')
        for kk=1:size(clusterData, 2)
            totalInnerkNNDistanceRandom = [clusterData(kk).Analysis.randomClusterStruct(3).kNNDistance; totalInnerkNNDistanceRandom];
            totalOuterkNNDistanceRandom = [clusterData(kk).Analysis.randomClusterStruct(4).kNNDistance; totalOuterkNNDistanceRandom];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% create NN-Plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h = figure( 'Name', 'kNN-Distance' );
    histogram(totalkNNDistance, ceil(max(totalkNNDistance(:, 1))), 'Normalization', 'probability');
    hold on
    if Clusterparamstruct.compareRandomData == true  
        histogram(totalkNNDistanceRandom, ceil(max(totalkNNDistanceRandom(:, 1))), 'Normalization', 'probability');
    end
    grid on;
    title([num2str(Clusterparamstruct.kNNDistance.k) '-Nearest Neighbor Distance']);
    xlabel('distance [nm]');
    ylabel('normalized frequency');
    if Clusterparamstruct.compareRandomData == true
        legend(['mean k-NN-Distance = ' num2str(mean(totalkNNDistance)) ' nm \pm ' num2str(std(totalkNNDistance)) ' nm'],...
            ['mean k-NN-Distance = ' num2str(mean(totalkNNDistanceRandom)) ' nm \pm ' num2str(std(totalkNNDistanceRandom)) ' nm']);
    else
        legend(['mean k-NN-Distance = ' num2str(mean(totalkNNDistance)) ' nm \pm ' num2str(std(totalkNNDistance)) ' nm']);
    end
    ax = gca;
    ax.XLim = [0 Clusterparamstruct.kNNDistance.maxDistance];
    hold off
    savefig(h, [Clusterparamstruct.DIR_output filesep 'kNNDistance.fig']);
    
    % inner kNN distances
    if isfield(clusterData(kk).Analysis.clusterStruct, 'clusterDBSCAN')
        h = figure( 'Name', 'kNN-Distance of points within a cluster' );
        histogram(totalInnerkNNDistance, ceil(max(totalInnerkNNDistance(:, 1))), 'Normalization', 'probability');
        hold on
        if Clusterparamstruct.compareRandomData == true && isfield(clusterData(kk).Analysis.randomClusterStruct, 'clusterDBSCAN')
            histogram(totalInnerkNNDistanceRandom, ceil(max(totalInnerkNNDistanceRandom(:, 1))), 'Normalization', 'probability');
        end
        grid on;
        title([num2str(Clusterparamstruct.kNNDistance.k) '-Nearest Neighbor Distance of points within a cluster']);
        xlabel('distance [nm]');
        ylabel('normalized frequency');
        if Clusterparamstruct.compareRandomData == true
            legend(['mean k-NN-Distance = ' num2str(mean(totalInnerkNNDistance)) ' nm \pm ' num2str(std(totalInnerkNNDistance)) ' nm'],...
                ['mean k-NN-Distance = ' num2str(mean(totalInnerkNNDistanceRandom)) ' nm \pm ' num2str(std(totalInnerkNNDistanceRandom)) ' nm']);
        else
            legend(['mean k-NN-Distance = ' num2str(mean(totalInnerkNNDistance)) ' nm \pm ' num2str(std(totalInnerkNNDistance)) ' nm']);
        end    
        ax = gca;
        ax.XLim = [0 Clusterparamstruct.kNNDistance.maxDistance];
        hold off
        savefig(h, [Clusterparamstruct.DIR_output filesep 'innerkNNDistance.fig']);
    end
    % outer kNNdistances
    if isfield(clusterData(kk).Analysis.clusterStruct, 'clusterDBSCAN')
        h = figure( 'Name', 'kNN-Distance of points outside a cluster' );
        histogram(totalOuterkNNDistance, ceil(max(totalOuterkNNDistance(:, 1))), 'Normalization', 'probability');
        hold on
        if Clusterparamstruct.compareRandomData == true && isfield(clusterData(kk).Analysis.randomClusterStruct, 'clusterDBSCAN')
            histogram(totalOuterkNNDistanceRandom, ceil(max(totalOuterkNNDistanceRandom(:, 1))), 'Normalization', 'probability');
        end
        grid on;
        title([num2str(Clusterparamstruct.kNNDistance.k) '-Nearest Neighbor Distance of points outside a cluster']);
        xlabel('distance [nm]');
        ylabel('normalized frequency');
        if Clusterparamstruct.compareRandomData == true
            legend(['mean k-NN-Distance = ' num2str(mean(totalOuterkNNDistance)) ' nm \pm ' num2str(std(totalOuterkNNDistance)) ' nm'],...
                ['mean k-NN-Distance = ' num2str(mean(totalOuterkNNDistanceRandom)) ' nm \pm ' num2str(std(totalOuterkNNDistanceRandom)) ' nm']);
        else
            legend(['mean k-NN-Distance = ' num2str(mean(totalOuterkNNDistance)) ' nm \pm ' num2str(std(totalOuterkNNDistance)) ' nm']);
        end    
        ax = gca;
        ax.XLim = [0 Clusterparamstruct.kNNDistance.maxDistance];
        hold off
        savefig(h, [Clusterparamstruct.DIR_output filesep 'outerkNNDistance.fig']);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DBSCAN algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Clusterparamstruct.DBSCAN.algorithm == true
    % init algorithm specific values
    noOfCluster = 0;
    pointsInCluster = 0;
    totalPoints = 0;
    totalSignalsPerCluster = [];
    totalDiameter = [];
    totalDensity = [];
    totalNNDistance = [];
    noOfClusterRandom = 0;
    pointsInClusterRandom = 0;
    totalPointsRandom = 0;
    totalSignalsPerClusterRandom = [];
    totalDiameterRandom = [];
    totalDensityRandom = [];
    totalNNDistanceRandom = [];
    % measure number of clusters
    for kk=1:size(clusterData, 2)
        noOfCluster = noOfCluster + size(clusterData(kk).Analysis.clusterStruct(2).clusterDBSCAN, 2);
        pointsInCluster = pointsInCluster + clusterData(kk).Analysis.clusterStruct(1).clusterDBSCAN;
        totalPoints = totalPoints + clusterData(kk).Analysis.NoOfPoints;
        if size(clusterData(kk).Analysis.clusterStruct(2).clusterDBSCAN, 2) ~= 0
            totalSignalsPerCluster = [clusterData(kk).Analysis.clusterStruct(2).clusterDBSCAN.Molecules totalSignalsPerCluster];
            totalDiameter = [clusterData(kk).Analysis.clusterStruct(2).clusterDBSCAN.Diameter totalDiameter];
            totalDensity = [clusterData(kk).Analysis.clusterStruct(2).clusterDBSCAN.ClusterDensity totalDensity];
            if isfield(clusterData(kk).Analysis.clusterStruct(2).clusterDBSCAN, 'NNCluster')% image consists of more then one cluster
                totalNNDistance = [clusterData(kk).Analysis.clusterStruct(2).clusterDBSCAN.NNCluster totalNNDistance];
            end
        end
    end
    if Clusterparamstruct.compareRandomData == true
        for kk=1:size(clusterData, 2)
            noOfClusterRandom = noOfClusterRandom + size(clusterData(kk).Analysis.randomClusterStruct(2).clusterDBSCAN, 2);
            pointsInClusterRandom = pointsInClusterRandom + clusterData(kk).Analysis.randomClusterStruct(1).clusterDBSCAN;
            totalPointsRandom = totalPointsRandom + clusterData(kk).Analysis.NoOfPoints;
            if size(clusterData(kk).Analysis.randomClusterStruct(2).clusterDBSCAN, 2) ~= 0
                totalSignalsPerClusterRandom = [clusterData(kk).Analysis.randomClusterStruct(2).clusterDBSCAN.Molecules totalSignalsPerClusterRandom];
                totalDiameterRandom = [clusterData(kk).Analysis.randomClusterStruct(2).clusterDBSCAN.Diameter totalDiameterRandom];
                totalDensityRandom = [clusterData(kk).Analysis.randomClusterStruct(2).clusterDBSCAN.ClusterDensity totalDensityRandom];
                if isfield(clusterData(kk).Analysis.randomClusterStruct(2).clusterDBSCAN, 'NNCluster')% image consists of more then one cluster
                    totalNNDistanceRandom = [clusterData(kk).Analysis.randomClusterStruct(2).clusterDBSCAN.NNCluster totalNNDistanceRandom];
                end
            end
        end
    end
    totalSignalsPerCluster = single(totalSignalsPerCluster.');
    totalDiameter = single(totalDiameter.');
    totalDensity = single(totalDensity.');
    totalNNDistance = single(totalNNDistance.');
    totalSignalsPerClusterRandom = single(totalSignalsPerClusterRandom.');
    totalDiameterRandom = single(totalDiameterRandom.');
    totalDensityRandom = single(totalDensityRandom.');
    totalNNDistanceRandom = single(totalNNDistanceRandom.');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% create cluster plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h = figure( 'Name', 'Results of Cluster Analysis by DBSCAN Algorithm' );
    subplot(2, 2, 1);
    h1 = histogram(totalSignalsPerCluster);
    hold on
    if Clusterparamstruct.compareRandomData == true  
        histogram(totalSignalsPerClusterRandom, 'BinWidth', h1.BinWidth);
    end
    title('Signals per Cluster');
    % Label axes
    xlabel('# of detected signals per cluster');
    ylabel ('frequency');
    if Clusterparamstruct.compareRandomData == true
        legend(['mean number of signals per cluster = ' num2str(mean(totalSignalsPerCluster)) ' \pm ' num2str(std(totalSignalsPerCluster))],...
            ['mean number of signals per cluster = ' num2str(mean(totalSignalsPerClusterRandom)) ' \pm ' num2str(std(totalSignalsPerClusterRandom))]);
    else
        legend(['mean number of signals per cluster = ' num2str(mean(totalSignalsPerCluster)) ' \pm ' num2str(std(totalSignalsPerCluster))]);
    end    
    grid on
    hold off
    
    subplot(2, 2, 2);
    h2 = histogram(totalDiameter);
    hold on
    if Clusterparamstruct.compareRandomData == true  
        histogram(totalDiameterRandom, 'BinWidth', h2.BinWidth);
    end
    title('Cluster Diameter');  %( h, 'NeNa distances', ['Fit \sigma = ' num2str(fitresult.a) ' nm'], 'Location', 'NorthEast' );
    % Label axes
    xlabel('diameter [nm]');
    ylabel ('frequency');
    if Clusterparamstruct.compareRandomData == true
        legend(['mean cluster diameter = ' num2str(mean(totalDiameter)) ' nm \pm ' num2str(std(totalDiameter)) ' nm'],...
            ['mean cluster diameter = ' num2str(mean(totalDiameterRandom)) ' nm \pm ' num2str(std(totalDiameterRandom)) ' nm']);
    else
        legend(['mean cluster diameter = ' num2str(mean(totalDiameter)) ' nm \pm ' num2str(std(totalDiameter)) ' nm']);
    end    
    grid on
    hold off
    
    subplot(2, 2, 3);
    h3 = histogram(totalDensity);
    hold on
    if Clusterparamstruct.compareRandomData == true  
        histogram(totalDensityRandom, 'BinWidth', h3.BinWidth);
    end
    title('Cluster Density');  %( h, 'NeNa distances', ['Fit \sigma = ' num2str(fitresult.a) ' nm'], 'Location', 'NorthEast' );
    % Label axes
    xlabel('density [1/µm^2]');
    ylabel ('frequency');
    if Clusterparamstruct.compareRandomData == true
        legend(['mean cluster density = ' num2str(mean(totalDensity)) ' 1/µm^2 \pm ' num2str(std(totalDensity)) ' 1/µm^2'],...
            ['mean cluster density = ' num2str(mean(totalDensityRandom)) ' 1/µm^2 \pm ' num2str(std(totalDensityRandom)) ' 1/µm^2']);
    else
        legend(['mean cluster density = ' num2str(mean(totalDensity)) ' 1/µm^2 \pm ' num2str(std(totalDensity)) ' 1/µm^2']);
    end    
    grid on
    hold off
    
    subplot(2, 2, 4);
    h4 = histogram(totalNNDistance);
    hold on
    if Clusterparamstruct.compareRandomData == true  
        histogram(totalNNDistanceRandom, 'BinWidth', h4.BinWidth);
    end
    grid on
    title('Distance to next neighboring cluster');  %( h, 'NeNa distances', ['Fit \sigma = ' num2str(fitresult.a) ' nm'], 'Location', 'NorthEast' );
    % Label axes
    xlabel('distance [nm]');
    ylabel ('frequency');
    if Clusterparamstruct.compareRandomData == true
        legend(['mean distance = ' num2str(mean(totalNNDistance)) ' nm \pm ' num2str(std(totalNNDistance)) ' nm'],...
            ['mean distance = ' num2str(mean(totalNNDistanceRandom)) ' nm \pm ' num2str(std(totalNNDistanceRandom)) ' nm']);
    else
        legend(['mean distance = ' num2str(mean(totalNNDistance)) ' nm \pm ' num2str(std(totalNNDistance)) ' nm']);
    end
    hold off
    savefig(h, [Clusterparamstruct.DIR_output filesep 'DBSCANAlgo.fig']);
    
    % disp more statistics
    disp('DBSCAN Cluster algorithm statistics:');
    disp('The following values will not be saved!');
    disp(['Points: ' num2str(totalPoints)]);
    disp(['Points within Cluster: ' num2str(pointsInCluster)]);
    disp(['Number of Cluster: ' num2str(noOfCluster)]);
    disp(['Points within Cluster [random data]: ' num2str(pointsInClusterRandom)]);
    disp(['Number of Cluster [random data]: ' num2str(noOfClusterRandom)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distance analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Clusterparamstruct.Distance.algorithm == true
    totalDistance = [];
    totalDistanceRandom = [];
    totalInnerDistance = [];
    totalOuterDistance = [];
    totalInnerDistanceRandom = [];
    totalOuterDistanceRandom = []; 
    for kk=1:size(clusterData, 2)
        totalDistance = [clusterData(kk).Analysis.clusterStruct(1).allDistances; totalDistance];
    end
    % if DBSCAN was done before
    if isfield(clusterData(kk).Analysis.clusterStruct, 'clusterDBSCAN')
        for kk=1:size(clusterData, 2)
            totalInnerDistance = [clusterData(kk).Analysis.clusterStruct(2).allDistances; totalInnerDistance];
            totalOuterDistance = [clusterData(kk).Analysis.clusterStruct(3).allDistances; totalOuterDistance];
        end
    end
    if Clusterparamstruct.compareRandomData == true
        for kk=1:size(clusterData, 2)
            totalDistanceRandom = [clusterData(kk).Analysis.randomClusterStruct(1).allDistances; totalDistanceRandom];
        end
    end
    % if DBSCAN was done before
    if isfield(clusterData(kk).Analysis.randomClusterStruct, 'clusterDBSCAN')
        for kk=1:size(clusterData, 2)
            totalInnerDistanceRandom = [clusterData(kk).Analysis.randomClusterStruct(2).allDistances; totalInnerDistanceRandom];
            totalOuterDistanceRandom = [clusterData(kk).Analysis.randomClusterStruct(3).allDistances; totalOuterDistanceRandom];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% create distance plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h = figure( 'Name', 'Distance' );
    histogram(totalDistance, ceil(max(totalDistance(:, 1))), 'Normalization', 'probability');
    hold on
    if Clusterparamstruct.compareRandomData == true  
        histogram(totalDistanceRandom, ceil(max(totalDistanceRandom(:, 1))), 'Normalization', 'probability');
    end
    grid on;
    title('Distance Analysis');
    xlabel('distance [nm]');
    ylabel('normalized frequency');
    if Clusterparamstruct.compareRandomData == true
        legend(['mean distance = ' num2str(mean(totalDistance)) ' nm \pm ' num2str(std(totalDistance)) ' nm'],...
            ['mean distance = ' num2str(mean(totalDistanceRandom)) ' nm \pm ' num2str(std(totalDistanceRandom)) ' nm']);
    else
        legend(['mean distance = ' num2str(mean(totalDistance)) ' nm \pm ' num2str(std(totalDistance)) ' nm']);
    end    
    ax = gca;
    ax.XLim = [0 Clusterparamstruct.Distance.maxDistance];
    hold off
    savefig(h, [Clusterparamstruct.DIR_output filesep 'distance.fig']);
    
     % inner distances
    if isfield(clusterData(kk).Analysis.clusterStruct, 'clusterDBSCAN')
        h = figure( 'Name', 'Distances of points within a cluster' );
        histogram(totalInnerDistance, ceil(max(totalInnerDistance(:, 1))), 'Normalization', 'probability');
        hold on
        if Clusterparamstruct.compareRandomData == true && isfield(clusterData(kk).Analysis.randomClusterStruct, 'clusterDBSCAN')
            histogram(totalInnerDistanceRandom, ceil(max(totalInnerDistanceRandom(:, 1))), 'Normalization', 'probability');
        end
        grid on;
        title('Distance Analysis of points within a cluster');
        xlabel('distance [nm]');
        ylabel('normalized frequency');
        if Clusterparamstruct.compareRandomData == true
            legend(['mean distance = ' num2str(mean(totalInnerDistance)) ' nm \pm ' num2str(std(totalInnerDistance)) ' nm'],...
                ['mean distance = ' num2str(mean(totalInnerDistanceRandom)) ' nm \pm ' num2str(std(totalInnerDistanceRandom)) ' nm']);
        else
            legend(['mean distance = ' num2str(mean(totalInnerDistance)) ' nm \pm ' num2str(std(totalInnerDistance)) ' nm']);
        end    
        ax = gca;
        ax.XLim = [0 Clusterparamstruct.Distance.maxDistance];
        hold off
        savefig(h, [Clusterparamstruct.DIR_output filesep 'innerDistance.fig']);
    end
    % outer distances
    if isfield(clusterData(kk).Analysis.clusterStruct, 'clusterDBSCAN')
        h = figure( 'Name', 'Distances of points outside a cluster' );
        histogram(totalOuterDistance, ceil(max(totalOuterDistance(:, 1))), 'Normalization', 'probability');
        hold on
        if Clusterparamstruct.compareRandomData == true && isfield(clusterData(kk).Analysis.randomClusterStruct, 'clusterDBSCAN')
            histogram(totalOuterDistanceRandom, ceil(max(totalOuterDistanceRandom(:, 1))), 'Normalization', 'probability');
        end
        grid on;
        title('Distance Analysis of points outside a cluster');
        xlabel('distance [nm]');
        ylabel('normalized frequency');
        if Clusterparamstruct.compareRandomData == true
            legend(['mean distance = ' num2str(mean(totalOuterDistance)) ' nm \pm ' num2str(std(totalOuterDistance)) ' nm'],...
                ['mean distance = ' num2str(mean(totalOuterDistanceRandom)) ' nm \pm ' num2str(std(totalOuterDistanceRandom)) ' nm']);
        else
            legend(['mean distance = ' num2str(mean(totalOuterDistance)) ' nm \pm ' num2str(std(totalOuterDistance)) ' nm']);
        end    
        ax = gca;
        ax.XLim = [0 Clusterparamstruct.Distance.maxDistance];
        hold off
        savefig(h, [Clusterparamstruct.DIR_output filesep 'outerDistance.fig']);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% radial density function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Clusterparamstruct.RDF.algorithm == true
    % get bin size - should be the same for all samples
    binSize = clusterData(1).Analysis.clusterStruct(2).RDF;
    binValues = [];
    for kk=1:size(clusterData, 2)
        binValues = cat(2, clusterData(kk).Analysis.clusterStruct(1).RDF, binValues);
    end
    binValues = mean(binValues, 2);
    % random data
    binValuesRandom = [];
    if Clusterparamstruct.compareRandomData == true
        binSizeRandom = clusterData(1).Analysis.randomClusterStruct(2).RDF;
        for kk=1:size(clusterData, 2)
            binValuesRandom = cat(2, clusterData(kk).Analysis.randomClusterStruct(1).RDF, binValuesRandom);
        end
         binValuesRandom = mean(binValuesRandom, 2);
    end
    h = figure( 'Name', 'Radial Density Function' );
    plot( binSize, binValues);
    hold on
    if Clusterparamstruct.compareRandomData == true
        plot( binSizeRandom, binValuesRandom);
    end
    grid on;
    title('radial density function');
    xlabel('distance [nm]');
    ylabel('g(r)');
    if Clusterparamstruct.compareRandomData == true
        legend('experimental data', 'random data');
    else 
        legend('experimental data');
    end
    hold off
    savefig(h, [Clusterparamstruct.DIR_output filesep 'RDF.fig']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ripley's function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Clusterparamstruct.ripley.algorithm == true
    % get bin size - should be the same for all samples
    binSize = clusterData(1).Analysis.clusterStruct(2).ripley;
    binValues = [];
    for kk=1:size(clusterData, 2)
        binValues = cat(2, clusterData(kk).Analysis.clusterStruct(1).ripley, binValues);
    end
    binValues = mean(binValues, 2);
    % random data
    binValuesRandom = [];
    if Clusterparamstruct.compareRandomData == true
        binSizeRandom = clusterData(1).Analysis.randomClusterStruct(2).ripley;
        for kk=1:size(clusterData, 2)
            binValuesRandom = cat(2, clusterData(kk).Analysis.randomClusterStruct(1).ripley, binValuesRandom);
        end
         binValuesRandom = mean(binValuesRandom, 2);
    end
    h = figure( 'Name', 'Ripley Function' );
    plot( binSize, binValues);
    hold on
    if Clusterparamstruct.compareRandomData == true
        plot( binSizeRandom, binValuesRandom);
    end
    grid on;
    title('ripley function');
    xlabel('distance [nm]');
    ylabel('L(r)-r');
    if Clusterparamstruct.compareRandomData == true
        legend('experimental data', 'random data');
    else
        legend('experimental data');
    end    
    hold off
    savefig(h, [Clusterparamstruct.DIR_output filesep 'ripley.fig']);
end
%% save data
Clusterparamstruct.ClusterPathname = [Clusterparamstruct.DIR_output filesep];
Clusterparamstruct.ClusterFilename = [Clusterparamstruct.treatmentName, '.mat'];
treatmentName.(Clusterparamstruct.treatmentName) = clusterData;
save([Clusterparamstruct.ClusterPathname Clusterparamstruct.ClusterFilename], '-struct', 'treatmentName');
savefig(Clusterparamstruct.handles.figure1, [Clusterparamstruct.ClusterPathname 'settings'])
multiWaitbar('CLOSEALL');
end