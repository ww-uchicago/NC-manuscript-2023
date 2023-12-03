%plot traces and histograms for spike,epsc and ipsc, and do cross correlation among
%them
%input dir:
close all
clear all
pref_direction = 180;
channel = 1;
save_output_plot = 1;
spike_data =1;%change if have spike data

prompt = 'input e baseline start in data points'
e_bl_start = input(prompt);
prompt = 'input e baseline end in data points'
e_bl_end = input(prompt);

prompt = 'input i baseline start in data points'
i_bl_start = input(prompt);
prompt = 'input i baseline end in data points'
i_bl_end = input(prompt);

% specify paths: change for every file
common_path = 'F:\Wei Lab HD\recordings\New Rig\20170330\'; %change
date_info = '17330';
e_index_num = '004';%changebv       
e_data_fn = [date_info e_index_num '.abf'];%change

i_index_num = '003';%changebv       
i_data_fn = [date_info i_index_num '.abf'];%change

if spike_data
s_index_num = '001';%changebv   

s_data_fn = [date_info s_index_num '.abf'];%change
spike_threshold = -4;
data_start = 0; %seconds
data_middle = 2;%seconds for dividing on and off responses-Wei 20141107
data_stop = inf; %seconds
start_rep = 1; 
stop_rep = 3;
end

if pref_direction < 180
        null_direction = 90+pref_direction;
    else
        null_direction = pref_direction-90;
  end
% null_direction = 0;% optional: define arbitory null and pref dirs


if spike_data
s_stimulus_fn= ['stim' s_index_num '.txt'];
s_stimulus_path = [common_path s_stimulus_fn];

% load stimulus file and find pref and null dirs
s_stimulus_data=load('-ASCII',s_stimulus_path);
s_directions = s_stimulus_data(:,1);
        
s_pref_inds = find(s_directions == pref_direction)';
 s_null_inds = find(s_directions == null_direction)';
 
% % load spike traces
 [s_Data,SampInt]=abfload([common_path s_data_fn]);
s_Data = s_Data(:,channel,:);
s_Data = squeeze(s_Data);
end
    
%load stimulus data

e_stimulus_fn= ['stim' e_index_num '.txt'];
e_stimulus_path = [common_path e_stimulus_fn];

i_stimulus_fn= ['stim' i_index_num '.txt'];
i_stimulus_path = [common_path i_stimulus_fn];


e_stimulus_data=load('-ASCII',e_stimulus_path);
        e_directions = e_stimulus_data(:,1);
        
e_pref_inds = find(e_directions == pref_direction)'
e_pref_inds = e_pref_inds(start_rep:stop_rep);
  
i_stimulus_data=load('-ASCII',i_stimulus_path);
        i_directions = i_stimulus_data(:,1);
        
i_pref_inds = find(i_directions == pref_direction)'
  
i_pref_inds = i_pref_inds(start_rep:stop_rep);
    
e_null_inds = find(e_directions == null_direction)'
e_null_inds = e_null_inds(start_rep:stop_rep);

i_null_inds = find(i_directions == null_direction)'
i_null_inds = i_null_inds(start_rep:stop_rep);

% % load whole cell traces and stim directions
[e_Data,SampInt]=abfload([common_path e_data_fn]);
e_Data = e_Data(:,channel,:);
e_Data = squeeze(e_Data);

[i_Data,SampInt]=abfload([common_path i_data_fn]);
i_Data = i_Data(:,channel,:);
i_Data = squeeze(i_Data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot raw data
%now set the right scale in ms in the x axis, e.g. sampling interval 100 us
data_length = size(e_Data,1);
%data_length = 15900;
time_axis = 0:0.1:(data_length-1)/10;%in milliseconds

% % optional: manually input pref and/or null indices
% e_pref_inds = e_pref_inds(:,[2 3]);
% e_null_inds = e_null_inds(:,[2 3]);
% i_pref_inds = i_pref_inds(:,[2 3]);
% i_null_inds = i_null_inds(:,[2 3]);

figure(1)
subplot(2,1,1)
plot(time_axis, e_Data(:,e_pref_inds),'r');
xlabel('Time', 'FontSize',10);
title([common_path e_index_num 'e_individual_traces']);
hold on
plot(time_axis, e_Data(:,e_null_inds),'k');
xlabel('Time', 'FontSize',10);
hold off

subplot(2,1,2)
plot(time_axis, i_Data(:,i_pref_inds),'r');
xlabel('Time', 'FontSize',10);
title([common_path i_index_num 'i_individual_traces']);
hold on
plot(time_axis, i_Data(:,i_null_inds),'k');
xlabel('Time', 'FontSize',10);
legend('red-pref','black-null');
hold off


figure(2)
subplot(2,1,1)
% substract baseline from each sweep using the first 100 ms (1000 data point)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e_bs = e_Data([e_bl_start:e_bl_end],:); % take the first 1000 data points of all sweeps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e_baseline = mean(e_bs,1); %take the mean of each sweep
e_Data_bs_adjusted = bsxfun(@minus,e_Data,e_baseline);%substract the baselien value from the raw data

e_mean_null_trace = mean(e_Data_bs_adjusted(:,e_null_inds)')';
e_mean_pref_trace = mean(e_Data_bs_adjusted(:,e_pref_inds)')';

%Optional: if not using baseline adjustment, just average the raw traces:
% e_mean_null_trace = mean(e_Data(:,e_null_inds)')';
% e_mean_pref_trace = mean(e_Data(:,e_pref_inds)')';


plot(time_axis, e_mean_pref_trace,'r');
xlabel('Time', 'FontSize',10);
title([common_path e_index_num 'e_mean_traces']);
legend('red-pref','black-null');
hold on
plot(time_axis, e_mean_null_trace,'k');
xlabel('Time', 'FontSize',10);
hold off

subplot(2,1,2)
% substract baseline from each sweep using the first 100 ms (1000 data point)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_bs = i_Data([i_bl_start:i_bl_end],:); % take the first 1000 data points of all sweeps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_baseline = mean(i_bs,1); %take the mean of each sweep
i_Data_bs_adjusted = bsxfun(@minus,i_Data,i_baseline);%substract the baselien value from the raw data

i_mean_null_trace = mean(i_Data_bs_adjusted(:,i_null_inds)')';
i_mean_pref_trace = mean(i_Data_bs_adjusted(:,i_pref_inds)')';
plot(time_axis, i_mean_pref_trace,'r');
xlabel('Time', 'FontSize',10);
title([common_path i_index_num 'i_mean_traces']);
hold on
plot(time_axis, i_mean_null_trace,'k');
xlabel('Time', 'FontSize',10);
hold off


% fig 3 plot individual sweeps after baseline correction
figure(3)
subplot(2,1,1)
plot(time_axis, e_Data_bs_adjusted(:,e_pref_inds),'r');
xlabel('Time', 'FontSize',10);
title([common_path e_index_num 'e_individual_traces_bs_adjusted']);
hold on
plot(time_axis, e_Data_bs_adjusted(:,e_null_inds),'k');
xlabel('Time', 'FontSize',10);
hold off

subplot(2,1,2)
plot(time_axis, i_Data_bs_adjusted(:,i_pref_inds),'r');
xlabel('Time', 'FontSize',10);
title([common_path i_index_num 'i_individual_traces_bs_adjusted']);
hold on
plot(time_axis, i_Data_bs_adjusted(:,i_null_inds),'k');
xlabel('Time', 'FontSize',10);
legend('red-pref','black-null');
hold off

prompt = 'use current baseline?';
answer = input(prompt); %use 1 for yes, 0 no.
if answer == 1

%fig4, plot spike, Ie,Ii for pref vs null, in the form of traces
figure (4)
subplot(2,1,1)
plot(time_axis, e_mean_pref_trace,'r');
xlabel('Time (ms)', 'FontSize',10);
title([common_path e_index_num i_index_num 'pref s Ie Ii']);
hold on
plot(time_axis, i_mean_pref_trace,'b');
if spike_data
    hold on
    plot(time_axis, s_Data(:,s_pref_inds)-200);
    hold off
end

subplot (2,1,2)
plot(time_axis, e_mean_null_trace,'r');
xlabel('Time (ms)', 'FontSize',10);
title([common_path e_index_num i_index_num 'null s Ie Ii']);
hold on
plot(time_axis, i_mean_null_trace,'b');
if spike_data
    hold on
    plot(time_axis, s_Data(:,s_null_inds)-200);
    hold off
end

%fig5, plot ge-gi for pref vs null, in the form of traces
pref_relative_difference = -1*e_mean_pref_trace - i_mean_pref_trace;
null_relative_difference = -1*e_mean_null_trace - i_mean_null_trace;
figure (5)
subplot(2,1,1)
plot(time_axis, pref_relative_difference,'r');
xlabel('Time (ms)', 'FontSize',10);
title([common_path e_index_num i_index_num 'pref ge-gi']);
if spike_data
    hold on
    plot(time_axis, s_Data(:,s_pref_inds)-200);
    hold off
end

subplot (2,1,2)
plot(time_axis, null_relative_difference,'k');
xlabel('Time (ms)', 'FontSize',10);
title([common_path e_index_num i_index_num 'null ge-gi']);
if spike_data
    hold on
    plot(time_axis, s_Data(:,s_null_inds)-200);
    hold off
end

% Next make histograms: first spike data
if spike_data
filtData = zeros(size(s_Data));

%filter data, and load spiketimes into a structure
SpikeTimes = cell(size(s_Data,2),1); %creat a cell array of number of rep x 1
for c = 1:size(s_Data,2) %for every repetition
    filtData(:,c) = filterData(s_Data(:,c),SampInt);
    SpikeTimes{c} = findSpikeTimes(filtData(:,c),SampInt,...
        'thresh',spike_threshold,'show_plot',0);
    
end

        %remove first 0.12 secs form each trial if necessary
        clippedSpikeTimes = subseqSpikeTimes(SpikeTimes,'start',data_start,'stop',data_stop);
        %sort data according to directions and plot if DS data
%         s_stimulus_data=load('-ASCII',stimulus_path);
%         s_directions = stimulus_data(:,1);
        SortedDSData.total = sortDSdatabeta(clippedSpikeTimes,s_directions,'first_repetition',start_rep,'last_repetition',stop_rep);
        %On spikes
        clippedSpikeTimesOn = subseqSpikeTimes(SpikeTimes,'start',data_start,'stop',data_middle);
        SortedDSData.On = sortDSdatabeta(clippedSpikeTimesOn,s_directions,'first_repetition',start_rep,'last_repetition',stop_rep);
         %Off spikes
         clippedSpikeTimesOff = subseqSpikeTimes(SpikeTimes,'start',data_middle,'stop',data_stop);
        SortedDSData.Off = sortDSdatabeta(clippedSpikeTimesOff,s_directions,'first_repetition',start_rep,'last_repetition',stop_rep)
     
    duration = (size(s_Data,1)-1)/10000;%in seconds
   

% %find index for pref and null dirs
null_ind = find(abs((s_directions-SortedDSData.total.null_dir)) < eps('single')); %hack to find match
pref_ind = find(abs((s_directions-SortedDSData.total.pref_dir)) < eps('single')); %hack to find match

%select spikes times:
nullSpikeTimes = {clippedSpikeTimes{null_ind'}};
prefSpikeTimes = {clippedSpikeTimes{pref_ind'}};


% make PSTH plots for null and pref dir
%this is figure 6
figure(6) %PSTH spike pref and null
subplot(6,1,1)%PSTH spike pref
all_pref_spike_times = cat(2,prefSpikeTimes{:});
s_edges = 0:0.1:duration;
[pref_PSTH pref_indices] = histc(all_pref_spike_times,s_edges);
%choose total number of spikes or firing rate
 %pref_PSTH = pref_PSTH./length(all_pref_spike_times);%normalized firing rate
 bar(s_edges,pref_PSTH,1,'r')
% makePSTH(prefSpikeTimes,0,duration,0.1);%in seconds
legend('spike_pref')

subplot(6,1,2)% PSTH spike null
%this is figure 6
all_null_spike_times = cat(2,nullSpikeTimes{:});
[null_PSTH null_indices] = histc(all_null_spike_times,s_edges);
 %null_PSTH = null_PSTH./length(all_null_spike_times);
 bar(s_edges,null_PSTH,1,'k')%convert to normalized firing rate
% makePSTH(prefSpikeTimes,0,duration,0.1);%in seconds
legend('spike_null')
         
end

%now make bar graph for whole cell currents
%first take integral for each bin
e_edges = 1:1000:(data_length-1); %in data points
e_pref_bin_data = [];
for i = 1:size(e_edges,2) %for every repetition
    if i < size(e_edges,2)
    e_pref_bin_data(i)= sum(e_mean_pref_trace(e_edges(i):e_edges(i+1)));
    else
      e_pref_bin_data(i)= sum(e_mean_pref_trace(e_edges(i):data_length));
    end
    
end
%figure 7 binned pref e and i
subplot(6,1,3)
bar(e_edges,e_pref_bin_data,1,'r')
hold on

i_pref_bin_data = [];
for i = 1:size(e_edges,2) %for every repetition
    if i < size(e_edges,2)
    i_pref_bin_data(i)= sum(i_mean_pref_trace(e_edges(i):e_edges(i+1)));
    else
      i_pref_bin_data(i)= sum(i_mean_pref_trace(e_edges(i):data_length));
    end
    
end

bar(e_edges,i_pref_bin_data,1,'b')
title('current pref')
legend('red-e; blue I')
hold off

%figure 8 binned null e and i
e_null_bin_data = [];
for i = 1:size(e_edges,2) %for every repetition
    if i < size(e_edges,2)
    e_null_bin_data(i)= sum(e_mean_null_trace(e_edges(i):e_edges(i+1)));
    else
      e_null_bin_data(i)= sum(e_mean_null_trace(e_edges(i):data_length));
    end
    
end
%figure 8 binned null e and i
subplot(6,1,4)
bar(e_edges,e_null_bin_data,1,'r')
hold on

i_null_bin_data = [];
for i = 1:size(e_edges,2) %for every repetition
    if i < size(e_edges,2)
    i_null_bin_data(i)= sum(i_mean_null_trace(e_edges(i):e_edges(i+1)));
    else
      i_null_bin_data(i)= sum(i_mean_null_trace(e_edges(i):data_length));
    end
    
end

bar(e_edges,i_null_bin_data,1,'b')
title('current null')
legend('e-red, i-blue')
hold off

%figure 9 plot relative binned difference between e and i
relative_pref_bin_data = [];
for i = 1:size(e_edges,2) %for every repetition
    if i < size(e_edges,2)
    relative_pref_bin_data(i)= sum(pref_relative_difference(e_edges(i):e_edges(i+1)));
    else
      relative_pref_bin_data(i)= sum(pref_relative_difference(e_edges(i):data_length));
    end
    
end

subplot(6,1,5)
bar(e_edges,relative_pref_bin_data,1,'r');
title('pref e-i')

subplot(6,1,6)
relative_null_bin_data = [];
for i = 1:size(e_edges,2) %for every repetition
    if i < size(e_edges,2)
    relative_null_bin_data(i)= sum(null_relative_difference(e_edges(i):e_edges(i+1)));
    else
      relative_null_bin_data(i)= sum(null_relative_difference(e_edges(i):data_length));
    end
    
end
bar(e_edges,relative_null_bin_data,1,'r');
title('null e-i')

%examine cross correlation - pref, x-axis is seconds
[p_s_e_r,p_s_e_lags]=xcorr(pref_PSTH,-1*e_pref_bin_data);
figure(7)
norm_p_s_e_r=p_s_e_r/(sqrt(sum(pref_PSTH.^2))*sqrt(sum(e_pref_bin_data.^2)));
plot(0.1*p_s_e_lags,norm_p_s_e_r,'k')
hold on

[p_s_i_r,p_s_i_lags]=xcorr(pref_PSTH,i_pref_bin_data);
norm_p_s_i_r=p_s_i_r/(sqrt(sum(pref_PSTH.^2))*sqrt(sum(i_pref_bin_data.^2)));
plot(0.1*p_s_i_lags,norm_p_s_i_r,'b')
hold on

[p_s_s_r,p_s_s_lags]=xcorr(pref_PSTH,relative_pref_bin_data);
norm_p_s_s_r=p_s_s_r/(sqrt(sum(pref_PSTH.^2))*sqrt(sum(relative_pref_bin_data.^2)));

plot(0.1*p_s_s_lags,norm_p_s_s_r,'c')
hold off

%plot xcorr for null
figure(8)
[n_s_e_r,n_s_e_lags]=xcorr(null_PSTH,-1*e_null_bin_data);
norm_n_s_e_r = n_s_e_r/(sqrt(sum(null_PSTH.^2))*sqrt(sum(e_null_bin_data.^2)));
plot(0.1*n_s_e_lags,norm_n_s_e_r,'k')
hold on
[n_s_i_r,n_s_i_lags]=xcorr(null_PSTH,i_null_bin_data);
norm_n_s_i_r = n_s_i_r/(sqrt(sum(null_PSTH.^2))*sqrt(sum(i_null_bin_data.^2)));
plot(0.1*n_s_i_lags,norm_n_s_i_r,'b')
hold on

[n_s_s_r,n_s_s_lags]=xcorr(null_PSTH,relative_null_bin_data);
norm_n_s_s_r = n_s_s_r/(sqrt(sum(null_PSTH.^2))*sqrt(sum(relative_null_bin_data.^2)));
plot(0.1*n_s_s_lags,norm_n_s_s_r,'c')
hold off


tilefigs

