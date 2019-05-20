function visRDF(Clusterparamstruct, clusterData)
%visRDF Summary of this function goes here
%   Detailed explanation goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    %% plotting
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
end