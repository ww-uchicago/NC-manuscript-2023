% This serves as the main fu
close all; clear all;clearvars; close all; clc;

folder = 'J:\Wei Lab HD\recordings\New Rig\20190304\Cell 1';


moviefiles = dir([folder  '\Time Trace(s)*.xls']);


for i = 1:size(moviefiles,1)
    close all
    clearvars -except paras folder moviefiles i retinaInfo cell_area num_stimuli compare pixel_xcorr smoothing unmixing 
    movieNum = str2num(moviefiles(i,1).name(14:end-4));
    moviestring = num2str(movieNum);
    txtfile_path = [folder  '\Time Trace(s)' moviestring '.xls'];
    txtfile = csvread(txtfile_path,1,0);
    dlmwrite([folder '\CellSignals\CellSignals' moviestring '.txt'],txtfile(:,1:end-2), 'delimiter','\t');
    
end
    
%dlmwrite([folder date '\MouseDatasetArea.txt'],dataset, 'delimiter','\t');
