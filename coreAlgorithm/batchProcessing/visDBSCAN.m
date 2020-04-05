function visDBSCAN(Clusterparamstruct, clusterData)
%visDBSCAN Summary of this function goes here
%   Detailed explanation goes here
% integrate structfun
% integrate plotting func
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Clusterparamstruct.DBSCAN.algorithm == true
    nSamples = size(clusterData, 2);
    totalSignalsPerCluster = cell(1, nSamples);
    totalDiameter = cell(1, nSamples);
    totalDensity = cell(1, nSamples);
    totalNNDistance = cell(1, nSamples);
    totalSignalsPerClusterRandom = cell(1, nSamples);
    totalDiameterRandom = cell(1, nSamples);
    totalDensityRandom = cell(1, nSamples);
    totalNNDistanceRandom = cell(1, nSamples);
    % measure number of clusters
    for kk=1:nSamples
        if size(clusterData(kk).Analysis.clusterStruct(2).clusterDBSCAN, 2) ~= 0
            for ii=1:size(clusterData(kk).Analysis.clusterStruct(2).clusterDBSCAN, 2)
                totalSignalsPerCluster{kk} = cat(2, clusterData(kk).Analysis.clusterStruct(2).clusterDBSCAN(ii).Molecules, totalSignalsPerCluster{kk});
                totalDiameter{kk} = cat(2, clusterData(kk).Analysis.clusterStruct(2).clusterDBSCAN(ii).Diameter, totalDiameter{kk});
                totalDensity{kk} = cat(2, clusterData(kk).Analysis.clusterStruct(2).clusterDBSCAN(ii).ClusterDensity, totalDensity{kk});
                if isfield(clusterData(kk).Analysis.clusterStruct(2).clusterDBSCAN, 'NNCluster')% image consists of more then one cluster
                    totalNNDistance{kk} = cat(2, clusterData(kk).Analysis.clusterStruct(2).clusterDBSCAN(ii).NNCluster, totalNNDistance{kk});
                end
            end
        end
    end
    if Clusterparamstruct.compareRandomData == true
        for kk=1:nSamples
            if size(clusterData(kk).Analysis.randomClusterStruct(2).clusterDBSCAN, 2) ~= 0
                for ii=1:size(clusterData(kk).Analysis.randomClusterStruct(2).clusterDBSCAN, 2)
                    totalSignalsPerClusterRandom{kk} = cat(2, clusterData(kk).Analysis.randomClusterStruct(2).clusterDBSCAN(ii).Molecules, totalSignalsPerClusterRandom{kk});
                    totalDiameterRandom{kk} = cat(2, clusterData(kk).Analysis.randomClusterStruct(2).clusterDBSCAN(ii).Diameter, totalDiameterRandom{kk});
                    totalDensityRandom{kk} = cat(2, clusterData(kk).Analysis.randomClusterStruct(2).clusterDBSCAN(ii).ClusterDensity, totalDensityRandom{kk});
                    if isfield(clusterData(kk).Analysis.randomClusterStruct(2).clusterDBSCAN, 'NNCluster')% image consists of more then one cluster
                        totalNNDistanceRandom{kk} = cat(2, clusterData(kk).Analysis.randomClusterStruct(2).clusterDBSCAN(ii).NNCluster, totalNNDistanceRandom{kk});
                    end
                end
            end
        end
    end
    totalSignalsPerCluster = cell2mat(totalSignalsPerCluster);
    totalDiameter = cell2mat(totalDiameter);
    totalDensity = cell2mat(totalDensity);
    totalNNDistance = cell2mat(totalNNDistance);
    totalSignalsPerClusterRandom = cell2mat(totalSignalsPerClusterRandom);
    totalDiameterRandom = cell2mat(totalDiameterRandom);
    totalDensityRandom = cell2mat(totalDensityRandom);
    totalNNDistanceRandom = cell2mat(totalNNDistanceRandom);
    %% plotting
    h = figure( 'Name', 'Results of Cluster Analysis by DBSCAN Algorithm' );
    % signals per cluster
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
    % diameter
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
    % density of signals within clusters
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
    % distance to next neighboring cluster
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
end
end