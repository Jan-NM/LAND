%% run script for every csv file in a folder
%% write the result to a csv file

% d = uigetdir(pwd, 'Select a folder');
d = "/home/ole/Hiwi/Data/Mohan_filtered/Sperm/PDC/filtered/cropped";
out_dir = "/home/ole/Hiwi/matlab_plots/Mohan/PDC";
disp(d)
files = dir(fullfile(d, '*cropped.mat'));

filenames = {files(:).name};

dims = 2;

fid_out = fopen('stats_mohan_PDC.csv','w'); 
fprintf(fid_out,'%s,%s,%s,%s,%s,%s\n', 'path', 'filename','smallest_area','biggest_area','0.5% quant density', '99.5% quant density');

for k = 1 : length(filenames)
    [filename] = convertCharsToStrings(filenames{k});
    folder_names = split(d,"/");
    folder_name = folder_names(end);
    filename_no_ext = extractBetween(filename,1,strlength(filename)-4);
  
    disp(filename)
    load(d+"/"+filename);
    
    cell = ClusterAnalysis(Orte,'001',dims);
    [smallest_vol, biggest_vol, lowest_plot_density, highest_plot_density] = cell.voronoiCluster(0,dims,0,out_dir, folder_name+"_"+filename_no_ext);
    fprintf(fid_out,'%s,%s,%.5f,%.5f,%.15f,%.8f\n',d, filename, ...
                    smallest_vol, biggest_vol, lowest_plot_density, highest_plot_density);
        
end

fclose(fid_out);
