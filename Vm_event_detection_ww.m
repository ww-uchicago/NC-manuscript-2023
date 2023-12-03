%% event detection for SAC Vm traces in control, 1 mM TEA or internal Cesium WW 10/30/2023

clear all
close all

%% load data
channel = 1; %channel to analyze on digidata
data_start = 0; %seconds
data_stop = inf; %seconds
common_path = 'E:\Dropbox\work_new\manuscripts\Hector mGluR2\data files\Vm_matlab_analysis\Vm Files\Control\';
abfFile = '20722011';
data_fn = [abfFile '.abf'];
savedir = 'E:\Dropbox\work_new\manuscripts\Hector mGluR2\data files\Vm_matlab_analysis\Vm Files\Control\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_threshold = 2; %  mV peak amplitude/prominence
% histogram of events counts
time_bin = 1000; % 100 ms each bin


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(data_stop,'inf');
    data_stop = inf;
end
data_path = [common_path data_fn];

%added this line to choose channel: Justin 2009-08-26
[Data,SampInt]=abfload(data_path);
Data = Data(:,channel,:);
Data = squeeze(Data);



%% filter data and find peaks and their amplitudes after removing slow depol
x0=zeros(size(Data)); %filtered data
slow_x = zeros(size(Data)); % slow component of vm
x=zeros(size(Data));% high freq component after subtraction
pks=cell(1,size(Data,2)); % absolute amplitude in trace, Vm
locs=cell(1,size(Data,2));% time for each peak
width=cell(1,size(Data,2));% half max width of each peak
prom=cell(1,size(Data,2)); % amplutude of peaks delata Vm

time = [1:size(Data,1)];% unit is 0.1 ms at 10kHz sampling frequency of pClamp
num_bin = floor(size(time,2)/time_bin);% how many bins in each sweep

for n = 1: size(Data, 2)
    x0(:,n) =smoothdata(Data(:,n),'gaussian',10);% check if filter is good
    slow_x(:,n) = smoothdata(Data(:,n),'movmean',1000);% check
    x(:,n)= x0(:,n) - slow_x(:,n);
  % plot(time,Data(:,n),'k',time,x(:,n),'r',time,slow_x(:,n),'b')
% plot(time,x(:,n))

  % find peaks
    [pks{n},locs{n},width{n},prom{n}]= findpeaks(x(:,n),'MinPeakProminence',v_threshold);

%  peakTime{n}= pktime(n,:);
%  peakAmplitude{n}=pkamp(n,:);
figure(n)
plot(time,Data(:,n),time,x(:,n),locs{n},pks{n},'o',locs{n},prom{n},'diamond')% 
hold on
xline(locs{n})
hold off
shg
% savefig([savedir abfFile '_' num2str(n)]); %optional save figs

%%% histogram of events counts
for k = 1:num_bin
    index{n}{k} = find(locs{n}>=time_bin*(k-1)+1 & locs{n}<k*time_bin);
bin_event_loc{n}{k} = locs{n}(index{n}{k});
bin_event_ampl{n}{k} = prom{n}(index{n}{k});
bin_event_absolute_amp{n}{k} = pks{n}(index{n}{k});
bin_event_width{n}{k} = width{n}(index{n}{k});
ev_count{n}(k) = size(index{n}{k},1);
end

end

% plot histogram of event counts
eventCounts = reshape (cell2mat(ev_count),[size(Data,2),num_bin]);% 10 rows, each row is histogram counts
meanEventCounts=mean(eventCounts);
stdEventCounts=std(eventCounts);

figure(n+1)
errorbar(meanEventCounts,stdEventCounts)
% savefig([savedir abfFile 'events']);
save([savedir abfFile]);
tilefigs