%% get prominence, or peak amplitdue of Vm events during baseline period or during moving bar stim for all cells
% first same abf file names in a excel file and import into matlab to save
% as a workspace file, e.g. ctrl.mat

clear all
close all

savedir = 'C:\Users\WW\Downloads\test\';
save_bl_binCount_File = [savedir 'ctrl_baseline_Events_per_bin'];
save_bar_binCount_File = [savedir 'ctrl_bar_Events_per_bin'];

filena = readtable([savedir 'ctrl.xlsx']);
fileNames = table2array(filena(:,1));
time_bin = 1000; % 100 ms each bin


for m = 1:size(fileNames,1)
  xx = fileNames{m};
    xx = char(xx);
    sl = strlength(xx);
    slc = sl-4;
  abfFile = xx(1:slc);
savedir = 'C:\Users\WW\Downloads\test\';
  filename = [savedir abfFile];
  load(filename);
% pull out event counts in new time axis  
baseline_index = find(time_bin_axis<-2 | time_bin_axis>=20);
bar_index = find(time_bin_axis >= 0 & time_bin_axis <16);

baseline_counts = mean(meanEventCounts(baseline_index));
bar_counts = mean(meanEventCounts(bar_index));
bar_events_per_bin{m} = bar_counts;
baseline_events_per_bin{m} = baseline_counts;
end

save(save_bl_binCount_File ,'baseline_events_per_bin');

save(save_bar_binCount_File ,'bar_events_per_bin');

