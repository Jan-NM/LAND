function [] = scatterPlot(obj)
%scatterPlot creates scatterplot of experimental and random data
%
% Jan Neumann, 18.04.18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% random data
figure( 'Name', num2str(obj.sampleID) );
if obj.dimension == 3
    scatter3(obj.randomTable(:, 2), obj.randomTable(:, 1), obj.randomTable(:, 3), '.'); 
    axis equal vis3d
    zlabel('z [nm]')
else
    scatter(obj.randomTable(:, 2), obj.randomTable(:, 1), '.');
    axis equal
end
title('Random data')
set(gca,'Ydir','reverse')
xlabel('x [nm]')
ylabel('y [nm]')
%% experimental data
figure( 'Name', num2str(obj.sampleID) );
if obj.dimension == 3
    scatter3(obj.positionTable(:, 2), obj.positionTable(:, 1), obj.positionTable(:, 3), '.'); 
    axis equal vis3d
    zlabel('z [nm]')
else
    scatter(obj.positionTable(:, 2), obj.positionTable(:, 1), '.');
    axis equal
end
title('Experimental data')
set(gca,'Ydir','reverse')
xlabel('x [nm]')
ylabel('y [nm]')
end