function [] = scatterPlot(obj)
%scatterPlot creates scatterplot of experimental and random data
%
% Jan Neumann, 18.04.18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure( 'Name', num2str(obj.sampleID) );
plot(obj.randomTable(:, 2), obj.randomTable(:, 1), '.');
title('Random data')
set(gca,'Ydir','reverse')
axis equal
axis off
figure( 'Name', num2str(obj.sampleID) );
title('Experimental data')
plot(obj.positionTable(:, 2), obj.positionTable(:, 1), '.');
title('Experimental data')
set(gca,'Ydir','reverse')
axis equal
axis off
end