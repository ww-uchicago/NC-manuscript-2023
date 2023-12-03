%% Align start times for histogram of SAC Vm events during moving bar - WW 10/30/2023 - workspace of
%% each abf file is now appended with an adjusted time axis "time_bin_axis'
clear all
close all


savedir = 'E:\Dropbox\work_new\manuscripts\Hector mGluR2\data files\Vm_matlab_analysis\Vm Files\Control\';
filena = readtable([savedir 'ctrl.xlsx']);
fileNames = table2array(filena(:,1));
time_bin = 1000; % 100 ms each bin

for i = 1:size(fileNames,1)
   x = fileNames{i};
    x = char(x);
    sl = strlength(x);
    slc = sl-4;
  abfFile = x(1:slc);
  filename = [savedir abfFile];
load(filename,'eventCounts');
num_bin = load(filename,'num_bin');
num_bin =num_bin.num_bin;
%% find the start time of the abf file
start_t = filena{i,2};
start_bin = ceil(start_t*10/time_bin);
%re-align time bin axis based on start bin
time_bin_axis = [1-start_bin:1:num_bin-start_bin];
% save mean event counts of this file with new time bin axis into a cell
% array
save(filename,'time_bin_axis',"-append");
end


