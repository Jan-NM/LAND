function visKNN_Distances(Clusterparamstruct, clusterData)
%visKNN_Distances Summary of this function goes here
%   Detailed explanation goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Clusterparamstruct.kNNDistance.algorithm == true
    nSamples = size(clusterData, 2);
    totalkNNDistance = cell(1, nSamples);
    totalkNNDistanceRandom = cell(1, nSamples);
    totalInnerkNNDistance = cell(1, nSamples);
    totalOuterkNNDistance = cell(1, nSamples);
    totalInnerkNNDistanceRandom = cell(1, nSamples);
    totalOuterkNNDistanceRandom = cell(1, nSamples);
    for kk=1:nSamples
        if isfield(clusterData(kk).Analysis.clusterStruct, 'kNNDistance')
            totalkNNDistance{kk} = clusterData(kk).Analysis.clusterStruct(2).kNNDistance;
        else
            totalkNNDistance{kk} = 0;
        end
    end
    totalkNNDistance = cell2mat(totalkNNDistance.');
    % if DBSCAN was done before
    if isfield(clusterData(kk).Analysis.clusterStruct, 'clusterDBSCAN')
        for kk=1:nSamples
            if isfield(clusterData(kk).Analysis.clusterStruct, 'kNNDistance')
                totalInnerkNNDistance{kk} = clusterData(kk).Analysis.clusterStruct(3).kNNDistance;
                totalOuterkNNDistance{kk} = clusterData(kk).Analysis.clusterStruct(4).kNNDistance;
            else
                totalInnerkNNDistance{kk} = 0;
                totalOuterkNNDistance{kk} = 0;
            end
        end
    end
    totalInnerkNNDistance = cell2mat(totalInnerkNNDistance.');
    totalOuterkNNDistance = cell2mat(totalOuterkNNDistance.');
    %% for random data
    if Clusterparamstruct.compareRandomData == true
        for kk=1:nSamples
            if isfield(clusterData(kk).Analysis.randomClusterStruct, 'kNNDistance')
                totalkNNDistanceRandom{kk} = clusterData(kk).Analysis.randomClusterStruct(2).kNNDistance;
            end
        end
        totalkNNDistanceRandom = cell2mat(totalkNNDistanceRandom.');
    end
    % if DBSCAN was done before
    if isfield(clusterData(kk).Analysis.randomClusterStruct, 'clusterDBSCAN')
        if ~isempty(clusterData(kk).Analysis.randomClusterStruct(2).clusterDBSCAN)
            for kk=1:nSamples
                if isfield(clusterData(kk).Analysis.randomClusterStruct, 'kNNDistance')
                    totalInnerkNNDistanceRandom{kk} = clusterData(kk).Analysis.randomClusterStruct(3).kNNDistance;
                    totalOuterkNNDistanceRandom{kk} = clusterData(kk).Analysis.randomClusterStruct(4).kNNDistance;
                end
            end
            totalInnerkNNDistanceRandom = cell2mat(totalInnerkNNDistanceRandom.');
            totalOuterkNNDistanceRandom = cell2mat(totalOuterkNNDistanceRandom.');
        end
    end
    %% plotting
    h = figure( 'Name', 'kNN-Distance' );
    histogram(totalkNNDistance, ceil(max(totalkNNDistance(:, 1))), 'Normalization', 'probability');
    hold on
    if Clusterparamstruct.compareRandomData == true && isfield(clusterData(kk).Analysis.randomClusterStruct, 'kNNDistance')
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
        if ~isempty(clusterData(kk).Analysis.clusterStruct(2).clusterDBSCAN)
            h = figure( 'Name', 'kNN-Distance of points within a cluster' );
            histogram(totalInnerkNNDistance, ceil(max(totalInnerkNNDistance(:, 1))), 'Normalization', 'probability');
            hold on
            if Clusterparamstruct.compareRandomData == true && isfield(clusterData(kk).Analysis.randomClusterStruct, 'kNNDistance')
                histogram(totalInnerkNNDistanceRandom, ceil(max(totalInnerkNNDistanceRandom(:, 1))), 'Normalization', 'probability');
            end
            grid on;
            title([num2str(Clusterparamstruct.kNNDistance.k) '-Nearest Neighbor Distance of points within a cluster']);
            xlabel('distance [nm]');
            ylabel('normalized frequency');
            if Clusterparamstruct.compareRandomData == true && isfield(clusterData(kk).Analysis.randomClusterStruct, 'kNNDistance')
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
    end
    % outer kNNdistances
    if isfield(clusterData(kk).Analysis.clusterStruct, 'clusterDBSCAN')
        if ~isempty(clusterData(kk).Analysis.clusterStruct(2).clusterDBSCAN)
            h = figure( 'Name', 'kNN-Distance of points outside a cluster' );
            histogram(totalOuterkNNDistance, ceil(max(totalOuterkNNDistance(:, 1))), 'Normalization', 'probability');
            hold on
            if Clusterparamstruct.compareRandomData == true && isfield(clusterData(kk).Analysis.randomClusterStruct, 'kNNDistance')
                histogram(totalOuterkNNDistanceRandom, ceil(max(totalOuterkNNDistanceRandom(:, 1))), 'Normalization', 'probability');
            end
            grid on;
            title([num2str(Clusterparamstruct.kNNDistance.k) '-Nearest Neighbor Distance of points outside a cluster']);
            xlabel('distance [nm]');
            ylabel('normalized frequency');
            if Clusterparamstruct.compareRandomData == true && isfield(clusterData(kk).Analysis.randomClusterStruct, 'kNNDistance')
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
end