%% run script for every csv file in a folder
%% write the result to a csv file

% d = uigetdir(pwd, 'Select a folder');
d = "/home/ole/Hiwi/Data/1A-filtered-zpm300";
out_dir = "/home/ole/Hiwi/matlab_plots";
disp(d)
files = dir(fullfile(d, '*.csv'));

filenames = {files(:).name};

x_col = 3;
y_col = 4;
z_col = 5;
dims = 3;

fid_out = fopen('file_out_2.csv','w'); 
fprintf(fid_out,'%s,%s,%s,%s,%s,%s\n', 'path', 'filename','smallest_vol','biggest_vol','0.5% quant density', '99.5% quant density');

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
    Orte(:, 11) = file(:, z_col) - min(file(:, z_col));
    
    cell = ClusterAnalysis(Orte,'001',dims);
    [smallest_vol, biggest_vol, lowest_plot_density, highest_plot_density] = cell.voronoiCluster(0,dims,0,out_dir, folder_name+"_"+filename_no_ext);
    fprintf(fid_out,'%s,%s,%.5f,%.5f,%.15f,%.8f\n',d, filename, ...
                    smallest_vol, biggest_vol, lowest_plot_density, highest_plot_density);
        
end

fclose(fid_out);
