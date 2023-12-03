%% get prominence, or peak amplitdue of Vm events during baseline period or during moving bar stim for all cells
% first same abf file names in a excel file and import into matlab to save
% as a workspace file, e.g. ctrl.mat

clear all
close all

savedir = 'E:\Dropbox\work_new\manuscripts\Hector mGluR2\data files\Vm_matlab_analysis\Vm Files\Cs Internal\';
saveblFile = [savedir 'cs_baselineEvents'];
savebarFile = [savedir 'cs_barEvents'];

filena = readtable([savedir 'cs.xlsx']);
fileNames = table2array(filena(:,1));
time_bin = 1000; % 100 ms each bin


for m = 1:size(fileNames,1)
  xx = fileNames{m};
    xx = char(xx);
    sl = strlength(xx);
    slc = sl-4;
  abfFile = xx(1:slc);
  filename = [savedir abfFile];
  load(filename);
% pull out event counts in new time axis  
  EvCounts{m} = meanEventCounts;
Time_axis{m} = time_bin_axis;

% pool baseline events, and bar evoked events for amplitude, freq and
% distribution analysis

%  baseline bin is -2 before or +20 or after for time bin axis
all_bl_p = []; % baseline prominence
all_bl_a = []; % baseline absolute Vm
all_b_p = []; % bar prominence
all_b_a = []; %bar amplitude Vm absolute

baseline_index = find(time_bin_axis<-2 | time_bin_axis>=20);
bar_index = find(time_bin_axis >= 0 & time_bin_axis <16);
for i = 1:size(bin_event_ampl,2)
    for j = 1:size(baseline_index,2)
        p = bin_event_ampl{1,i}{baseline_index(j)};
        all_bl_p = [all_bl_p, p'];
        a= bin_event_absolute_amp{1,i}{baseline_index(j)};
        all_bl_a = [all_bl_a, a'];
        
    end
      for l = 1:size(bar_index,2)
        bp = bin_event_ampl{1,i}{bar_index(l)};
        all_b_p = [all_b_p, bp'];
        ba= bin_event_absolute_amp{1,i}{bar_index(l)};
        all_b_a = [all_b_a, ba'];
      end
end
% %% remove aberrent outliers with ammplitude > 10 mV
% remov_index = find(all_b_p > 15);
% all_b_p_final = all_b_p;
% all_b_p_final(remov_index)=[];
% 
% remov_bl_index = find(all_bl_p > 15);
% all_bl_p_final = all_bl_p;
% all_bl_p_final(remov_bl_index)=[];

figure(m)
plot(all_b_p,'r')
hold on
% plot(all_b_p_final,'k')
% hold on
plot(all_bl_p,'b')
% hold on
% plot(all_bl_p_final,'g')
hold off
shg
bar_prom{m} = all_b_p;
baseline_prom{m} = all_bl_p;
end
tilefigs

% put all baseline events into a vector
baseline_amp = cell2mat(baseline_prom);
bar_amp= cell2mat(bar_prom);

[bl,bl_stats] = cdfplot(baseline_amp)
hold on
[bar,bar_stats] = cdfplot(bar_amp)
legend([bl bar],'baseline','bar');
hold off


save(saveblFile,'baseline_amp');
save(savebarFile,'bar_amp');
save([savedir 'cs_ev_workspace']);
