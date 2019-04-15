function [] = parameterEstimation( obj, method)
%parameterEstimation helps to estimate appropriate parameters for DBSCAN 
%   Jan Neumann, 27.06.17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch method
    case 'DBSCAN'
        figure('Name', 'DBSCAN parameter estimation');
        set(gcf,'units', 'normalized', 'outerposition', [0 0 1 1])
        subplot(2,2,1)
        %% calculate nearest neighbor distance (k = 2)
        positions = obj.positionTable;
        treeData = obj.queryData;
        switch obj.flagSearchTree
            case 'cTree' % for knn search, use Matlab's internal function
                if obj.dimension == 3
                    [~, nnDistance] = knnsearch(positions(:, 1:3), positions(:, 1:3), 'K', 3);
                else
                    [~, nnDistance] = knnsearch(positions(:, 1:2), positions(:, 1:2), 'K', 3);
                end
            case 'matlabTree'
                [~, nnDistance] = knnsearch(treeData.X, treeData.X, 'K', 3);
        end
        nnDistances = mean(nnDistance(:, 2:end), 2);
        % plot histogram
        histogram(nnDistances, 'Normalization', 'probability');
        grid on
        xlabel('distance [nm]');
        ylabel('normalized frequency');
        title('Nearest neighbor distances');
        %% enter radius and do rangesearch with corresponding radius
        prompt = {'Enter radius value [nm]:'};
        dlg_title = 'Input';
        num_lines = 1;
        defaultans = {'2'};
        radius = inputdlg(prompt,dlg_title,num_lines,defaultans);
        radius = str2double(radius{1});
        switch obj.flagSearchTree  % should be done on a cropped dataset for large position tables
            case 'cTree'
                if obj.dimension == 3
                    [~, distance] = rangesearch(positions(:, 1:3), positions(:, 1:3), radius);
                else
                    [~, distance] = rangesearch(positions(:, 1:2), positions(:, 1:2), radius);
                end
            case 'matlabTree'
                [~, distance] = rangesearch(positions(:, 1:2), treeData.X, radius);
        end
        minPoints = cellfun('length', distance);
        subplot(2,2,2)
        histogram(minPoints)
        title(['Frequency of number of neighbors for radius ' num2str(radius)]);
        xlabel('Number of neighbors');
        ylabel('frequency');
        grid on
        %% second method according to DBSCAN paper
        % enter k value
        prompt = {'Enter k value:'};
        dlg_title = 'Input';
        num_lines = 1;
        defaultans = {'4'};
        k = inputdlg(prompt,dlg_title,num_lines,defaultans);
        k = str2double(k{1});
        % calculate k-NN distance
        switch obj.flagSearchTree
            case 'cTree' % for knn search, use Matlab's internal function
                if obj.dimension == 3
                    [~, nnDistance] = knnsearch(positions(:, 1:3), positions(:, 1:3), 'K', k);
                else
                    [~, nnDistance] = knnsearch(positions(:, 1:2), positions(:, 1:2), 'K', k);
                end
            case 'matlabTree'
                [~, nnDistance] = knnsearch(treeData.X, treeData.X, 'K', k);
        end
        % mean over NN distance
        nnDistances = mean(nnDistance(:, 2:end), 2);
        %% create k-distance graph
        % sort k-distances in ascending order
        nnDistances = sortrows(nnDistances, 'ascend');
        % plot figures for estimation of radius parameter
        subplot(2,2,3)
        histogram(nnDistances, 'Normalization', 'probability');
        title([num2str(k), '-nearest neighbor distance']);
        grid on
        xlabel('distance [nm]');
        ylabel('normalized frequency');
        subplot(2,2,4)
        plot(nnDistances);
        title([num2str(k), '-distance graph']);
        grid on;
        xlabel('points');
        ylabel('radius [nm]');
    otherwise
        disp('Choose DBSCAN as method!')
end
end
