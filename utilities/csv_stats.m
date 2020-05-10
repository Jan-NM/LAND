%% run script for every csv file in a folder
%% write the result to a csv file

% d = uigetdir(pwd, 'Select a folder');
title_class = "DC";
d = "/extra/home/ole/Hiwi/Data/Mohan_raw/"+title_class+"/filtered";
out_dir = "/extra/home/ole/Hiwi/matlab_plots/Mohan/"+title_class;

%d = "/home/ole/Hiwi/Data";
%out_dir = "/home/ole/Hiwi/matlab_plots/Pres/small";

disp(d)
files = dir(fullfile(d, 'GSD*.csv'));

filenames = {files(:).name};

x_col = 3;
y_col = 4;
z_col = 5;
dims = 2;

max_area = 20*10^(-3);
max_cum_area= 50;
for k = 1 : length(filenames)
    [filename] = convertCharsToStrings(filenames{k});
    folder_names = split(d,"/");
    folder_name = folder_names(end);
    filename_no_ext = extractBetween(filename,1,strlength(filename)-4);
  
    disp(filename)
    file = csvread(d+"/"+filename, 1, 0);
    
    Orte = zeros(size(file, 1), 12);
    Orte(:, 2) = file(:, x_col) - min(file(:, x_col));
    Orte(:, 3) = file(:, y_col) - min(file(:, y_col));
    
    %Orte(:, 11) = file(:, z_col) - min(file(:, z_col));
    disp("size of orte");
    disp(size(Orte));
    Orte = datasample(Orte, 180000, 'Replace', false);
    disp("size of orte after");
    disp(size(Orte));
    
    cell = ClusterAnalysis(Orte,'001',dims);
    [smallest_vol, biggest_vol, lowest_plot_density, highest_plot_density] = ...
        cell.voronoiCluster();
        %cell.voronoiCluster(0,0,out_dir, folder_name+"_pres1_"+filename_no_ext, title_class, max_area, max_cum_area);        
        
    fprintf(fid_out,'%s,%s,%.5f,%.5f,%.15f,%.8f\n',d, filename, ...
                    smallest_vol, biggest_vol, lowest_plot_density, highest_plot_density);
        
end
