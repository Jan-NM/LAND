function vis_Distances(Clusterparamstruct, clusterData)
%vis_Distances Summary of this function goes here
%   Detailed explanation goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Clusterparamstruct.Distance.algorithm == true
    nSamples = size(clusterData, 2);
    totalDistance = cell(1, nSamples);
    totalDistanceRandom = cell(1, nSamples);
    totalInnerDistance = cell(1, nSamples);
    totalOuterDistance = cell(1, nSamples);
    totalInnerDistanceRandom = cell(1, nSamples);
    totalOuterDistanceRandom = cell(1, nSamples);
    for kk=1:nSamples
        if isfield(clusterData(kk).Analysis.clusterStruct, 'allDistances')
            totalDistance{kk} = clusterData(kk).Analysis.clusterStruct(1).allDistances;
        else
            totalDistance{kk} = 0;
        end
    end
    totalDistance = cell2mat(totalDistance.');
    % if DBSCAN was done before
    if isfield(clusterData(kk).Analysis.clusterStruct, 'clusterDBSCAN')
        for kk=1:nSamples
            if isfield(clusterData(kk).Analysis.clusterStruct, 'allDistances')
                totalInnerDistance{kk} = clusterData(kk).Analysis.clusterStruct(2).allDistances;
                totalOuterDistance{kk} = clusterData(kk).Analysis.clusterStruct(3).allDistances;
            else
                totalInnerDistance{kk} = 0;
                totalOuterDistance{kk} = 0;
            end
        end
    end
    totalInnerDistance = cell2mat(totalInnerDistance.');
    totalOuterDistance = cell2mat(totalOuterDistance.');
    %% for random data
    if Clusterparamstruct.compareRandomData == true
        for kk=1:nSamples
            if isfield(clusterData(kk).Analysis.randomClusterStruct, 'allDistances')
                totalDistanceRandom{kk} = clusterData(kk).Analysis.randomClusterStruct(1).allDistances;
            end
        end
        totalDistanceRandom = cell2mat(totalDistanceRandom.');
    end
    % if DBSCAN was done before
    if isfield(clusterData(kk).Analysis.randomClusterStruct, 'clusterDBSCAN')
        if ~isempty(clusterData(kk).Analysis.randomClusterStruct(2).clusterDBSCAN)
            for kk=1:nSamples
                if isfield(clusterData(kk).Analysis.randomClusterStruct, 'allDistances')
                    totalInnerDistanceRandom{kk} = clusterData(kk).Analysis.randomClusterStruct(2).allDistances;
                    totalOuterDistanceRandom{kk} = clusterData(kk).Analysis.randomClusterStruct(3).allDistances;
                end
            end
            totalInnerDistanceRandom = cell2mat(totalInnerDistanceRandom.');
            totalOuterDistanceRandom = cell2mat(totalOuterDistanceRandom.');
        end
    end
    %% plotting
    h = figure( 'Name', 'Distance' );
    histogram(totalDistance, ceil(max(totalDistance(:, 1))), 'Normalization', 'probability');
    hold on
    if Clusterparamstruct.compareRandomData == true && isfield(clusterData(kk).Analysis.randomClusterStruct, 'allDistances')
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
        if ~isempty(clusterData(kk).Analysis.clusterStruct(2).clusterDBSCAN)
            h = figure( 'Name', 'Distances of points within a cluster' );
            histogram(totalInnerDistance, ceil(max(totalInnerDistance(:, 1))), 'Normalization', 'probability');
            hold on
            if Clusterparamstruct.compareRandomData == true && isfield(clusterData(kk).Analysis.randomClusterStruct, 'allDistances')
                histogram(totalInnerDistanceRandom, ceil(max(totalInnerDistanceRandom(:, 1))), 'Normalization', 'probability');
            end
            grid on;
            title('Distance Analysis of points within a cluster');
            xlabel('distance [nm]');
            ylabel('normalized frequency');
            if Clusterparamstruct.compareRandomData == true && isfield(clusterData(kk).Analysis.randomClusterStruct, 'allDistances')
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
    end
    % outer distances
    if isfield(clusterData(kk).Analysis.clusterStruct, 'clusterDBSCAN')
        if ~isempty(clusterData(kk).Analysis.clusterStruct(2).clusterDBSCAN)
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
            if Clusterparamstruct.compareRandomData == true && isfield(clusterData(kk).Analysis.randomClusterStruct, 'clusterDBSCAN')
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
end