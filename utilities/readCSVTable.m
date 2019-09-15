[FileName, PathName] = uigetfile('*.csv', 'Select data');
% read .csv file
tempLocTable = csvread([PathName FileName],1,0);
Orte = zeros(size(tempLocTable, 1), 10);

% x 
Orte(:, 2) = tempLocTable(:, 3) - min(tempLocTable(:, 3));

% y
Orte(:, 3) = tempLocTable(:, 4) - min(tempLocTable(:, 4));

% y
Orte(:, 4) = tempLocTable(:, 5) - min(tempLocTable(:, 5));

%Orte(:, 11) = tempLocTable(:, 5);
%Orte(:, 12) = 0;
Orte(:, 9) = tempLocTable(:, 2);
% localization precision
%Orte(:, 4) = tempLocTable(:, 10);
%Orte(:, 5) = tempLocTable(:, 10);
% sort signals according to frame number in ascending order
Orte = sortrows(Orte, 9);
[savefile, savepath] = uiputfile('*.mat','Save file as...');

fullFileName = fullfile(savepath, savefile);
save(fullFileName, 'Orte');