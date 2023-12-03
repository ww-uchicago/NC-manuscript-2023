clear all
wt_filenames = ['J:\WeiLabHD\recordings\NewRig\20180726\18726011.abf';					
'J:\WeiLabHD\recordings\NewRig\20180729\18729002.abf';					
'J:\WeiLabHD\recordings\NewRig\20180729\18729014.abf';
'J:\WeiLabHD\recordings\NewRig\20180828\18828004.abf';					
'J:\WeiLabHD\recordings\NewRig\20180828\18828026.abf';					
'J:\WeiLabHD\recordings\NewRig\20180828\18828049.abf'];

ant_filenames = ['J:\WeiLabHD\recordings\NewRig\20180729\18729007.abf';					
'J:\WeiLabHD\recordings\NewRig\20180729\18729017.abf';					
'J:\WeiLabHD\recordings\NewRig\20180731\18731019.abf';					
'J:\WeiLabHD\recordings\NewRig\20180826\18826016.abf';					
'J:\WeiLabHD\recordings\NewRig\20180826\18826017.abf';					
'J:\WeiLabHD\recordings\NewRig\20180826\18826018.abf';										
'J:\WeiLabHD\recordings\NewRig\20180827\18827035.abf'];

tea_filenames = ['J:\WeiLabHD\recordings\NewRig\20180808\18808011.abf';
'J:\WeiLabHD\recordings\NewRig\20180828\18828008.abf';
'J:\WeiLabHD\recordings\NewRig\20180828\18828009.abf';
'J:\WeiLabHD\recordings\NewRig\20180828\18828030.abf';
'J:\WeiLabHD\recordings\NewRig\20180828\18828031.abf';
'J:\WeiLabHD\recordings\NewRig\20180828\18828054.abf';
'J:\WeiLabHD\recordings\NewRig\20180828\18828055.abf'];

drugs_filenames = ['J:\WeiLabHD\recordings\NewRig\20180808\18808029.abf';								
'J:\WeiLabHD\recordings\NewRig\20180826\18826023.abf';					
'J:\WeiLabHD\recordings\NewRig\20180826\18826024.abf';					
'J:\WeiLabHD\recordings\NewRig\20180826\18826026.abf';					
'J:\WeiLabHD\recordings\NewRig\20180827\18827017.abf';					
'J:\WeiLabHD\recordings\NewRig\20180827\18827015.abf';					
'J:\WeiLabHD\recordings\NewRig\20180827\18827016.abf';										
'J:\WeiLabHD\recordings\NewRig\20180827\18827022.abf';					
'J:\WeiLabHD\recordings\NewRig\20180827\18827023.abf';					
'J:\WeiLabHD\recordings\NewRig\20180827\18827026.abf';																						
'J:\WeiLabHD\recordings\NewRig\20180828\18828020.abf';					
'J:\WeiLabHD\recordings\NewRig\20180828\18828021.abf';								
'J:\WeiLabHD\recordings\NewRig\20180828\18828034.abf';					
'J:\WeiLabHD\recordings\NewRig\20180828\18828035.abf';					
'J:\WeiLabHD\recordings\NewRig\20180828\18828036.abf';					
'J:\WeiLabHD\recordings\NewRig\20180828\18828037.abf';								
'J:\WeiLabHD\recordings\NewRig\20180828\18828043.abf';									
'J:\WeiLabHD\recordings\NewRig\20180828\18828045.abf';															
'J:\WeiLabHD\recordings\NewRig\20180828\18828060.abf'; 
'J:\WeiLabHD\recordings\NewRig\20191217\19d17005.abf';															
'J:\WeiLabHD\recordings\NewRig\20191217\19d17009.abf';															
'J:\WeiLabHD\recordings\NewRig\20191217\19d17014.abf';															
'J:\WeiLabHD\recordings\NewRig\20191217\19d17015.abf';															
'J:\WeiLabHD\recordings\NewRig\20191217\19d17016.abf';															
'J:\WeiLabHD\recordings\NewRig\20191217\19d17020.abf';															
'J:\WeiLabHD\recordings\NewRig\20191217\19d17022.abf';															
'J:\WeiLabHD\recordings\NewRig\20191217\19d17024.abf';															
'J:\WeiLabHD\recordings\NewRig\20191217\19d17025.abf';															
'J:\WeiLabHD\recordings\NewRig\20191217\19d17026.abf'];			

wt_baseline_diff = [];
wt_area_activeBL = [];
wt_area_minBL = [];
wt_amp_activeBL = [];
wt_amp_minBL = [];
for i = 1:size(wt_filenames,1)
    clearvars folder index stimulus_fn stim_data Data Data_mean Data_sem minPks minLocs avg_baseline std_baseline minPks minLocs time_diff min_baseline
    folder = wt_filenames(i,1:(end-13));
    index = wt_filenames(i,(end-6):(end-4));
    stimulus_fn=strcat(folder,'\stim',index,'.txt');
    stim_data = load('-ASCII',stimulus_fn);
    
    if (stim_data(1,2) == 400.00) && (stim_data(1,3) == 400.00) && (stim_data(1,end) >= 100.00)
        wt_filenames(i,:)
        [Data, Samp] = abfload(wt_filenames(i,:));
        Data = Data(1:40000,1,:);
        Data = squeeze(Data);
        Data_mean = mean(Data,2);
        Data_sem = std(Data,0,2)/sqrt(size(Data,2));
        
        avg_baseline = mean(Data_mean(35000:40000,1));
        std_baseline = std(Data_mean(35000:40000,1));
        [minPks, minLocs] = findpeaks(-Data_mean,'MinPeakHeight',(-avg_baseline - 0.2*std_baseline),  'MinPeakDistance', 1500, 'MinPeakWidth',500);
        minPks = -minPks;
        
        time_diff = minLocs - 5000;
        min_baseline = Data_mean(minLocs(find(time_diff == min(time_diff))),1);
        
        wt_baseline_diff = [wt_baseline_diff ; abs(avg_baseline - min_baseline)];
        
        time_range = minLocs(find(time_diff == min(time_diff))): 23000;
        
        wt_area_activeBL = [wt_area_activeBL; trapz(time_range, (Data_mean(time_range,1)-avg_baseline))];
        wt_area_minBL = [wt_area_minBL; trapz(time_range, (Data_mean(time_range,1)-min_baseline))];
        
        wt_amp_activeBL = [wt_amp_activeBL; (max(Data_mean(time_range,1))-avg_baseline)];
        wt_amp_minBL = [wt_amp_minBL; (max(Data_mean(time_range,1))-min_baseline)];
        
        figure
        shadedErrorBar([1:size(Data,1)],Data_mean ,Data_sem,'lineprops',  '-k');
        title(wt_filenames(i,:))
        hold on
        plot(minLocs, Data_mean(minLocs), '*r')
        
        fig = gcf;
        saveas(fig, [folder '\MeanTrace' index], 'fig')
        close all
    end
end



ant_baseline_diff = [];
ant_area_activeBL = [];
ant_area_minBL = [];
ant_amp_activeBL = [];
ant_amp_minBL = [];
for i = 1:size(ant_filenames,1)
    clearvars folder index stimulus_fn stim_data Data Data_mean Data_sem avg_baseline std_baseline minPks minLocs time_diff min_baseline
    folder = ant_filenames(i,1:(end-13));
    index = ant_filenames(i,(end-6):(end-4));
    stimulus_fn=strcat(folder,'\stim',index,'.txt');
    stim_data = load('-ASCII',stimulus_fn);
    
    if (stim_data(1,2) == 400.00) && (stim_data(1,3) == 400.00)
        [Data, Samp] = abfload(ant_filenames(i,:));
        Data = Data(1:40000,1,:);
        Data = squeeze(Data);
        Data_mean = mean(Data,2);
        Data_sem = std(Data,0,2)/sqrt(size(Data,2));
        
        avg_baseline = mean(Data_mean(35000:40000,1));
        std_baseline = std(Data_mean(35000:40000,1));
        [minPks, minLocs] = findpeaks(-Data_mean,'MinPeakHeight',(-avg_baseline - 0.2*std_baseline),  'MinPeakDistance', 1500, 'MinPeakWidth',500);
        minPks = -minPks;
        
        time_diff = minLocs - 5000;
        min_baseline = Data_mean(minLocs(find(time_diff == min(time_diff))),1);
        
        ant_baseline_diff = [ant_baseline_diff ; abs(avg_baseline - min_baseline)];
        
        time_range = minLocs(find(time_diff == min(time_diff))): 23000;
        
        ant_area_activeBL = [ant_area_activeBL; trapz(time_range, (Data_mean(time_range,1)-avg_baseline))];
        ant_area_minBL = [ant_area_minBL; trapz(time_range, (Data_mean(time_range,1)-min_baseline))];
        
        ant_amp_activeBL = [ant_amp_activeBL; (max(Data_mean(time_range,1))-avg_baseline)];
        ant_amp_minBL = [ant_amp_minBL; (max(Data_mean(time_range,1))-min_baseline)];
        
        figure
        shadedErrorBar([1:size(Data,1)],Data_mean ,Data_sem,'lineprops',  '-k');
        title(ant_filenames(i,:))
        hold on
        plot(minLocs, Data_mean(minLocs), '*r')
        
     
        
        fig = gcf;
        saveas(fig, [folder '\MeanTrace' index], 'fig')
        close all
    end
end

tea_baseline_diff = [];
tea_area_activeBL = [];
tea_area_minBL = [];
tea_amp_activeBL = [];
tea_amp_minBL = [];
for i = 1:size(tea_filenames,1)
    clearvars folder index stimulus_fn stim_data Data Data_mean Data_sem avg_baseline std_baseline minPks minLocs time_diff min_baseline
    folder = tea_filenames(i,1:(end-13));
    index = tea_filenames(i,(end-6):(end-4));
    stimulus_fn=strcat(folder,'\stim',index,'.txt');
    stim_data = load('-ASCII',stimulus_fn);
    
    if (stim_data(1,2) == 400.00) && (stim_data(1,3) == 400.00)
        [Data, Samp] = abfload(tea_filenames(i,:));
        Data = Data(1:40000,1,:);
        Data = squeeze(Data);
        Data_mean = mean(Data,2);
        Data_sem = std(Data,0,2)/sqrt(size(Data,2));
       
        avg_baseline = mean(Data_mean(35000:40000,1));
        std_baseline = std(Data_mean(35000:40000,1));
        [minPks, minLocs] = findpeaks(-Data_mean,'MinPeakHeight',(-avg_baseline - 0.2*std_baseline),  'MinPeakDistance', 1500, 'MinPeakWidth',500);
        minPks = -minPks;
        
        time_diff = minLocs - 5000;
        min_baseline = Data_mean(minLocs(find(time_diff == min(time_diff))),1);
        
        tea_baseline_diff = [tea_baseline_diff ; abs(avg_baseline - min_baseline)];
        
        time_range = minLocs(find(time_diff == min(time_diff))): 23000;
        
        tea_area_activeBL = [tea_area_activeBL; trapz(time_range, (Data_mean(time_range,1)-avg_baseline))];
        tea_area_minBL = [tea_area_minBL; trapz(time_range, (Data_mean(time_range,1)-min_baseline))];
        
        tea_amp_activeBL = [tea_amp_activeBL; (max(Data_mean(time_range,1))-avg_baseline)];
        tea_amp_minBL = [tea_amp_minBL; (max(Data_mean(time_range,1))-min_baseline)];
        
        figure
        shadedErrorBar([1:size(Data,1)],Data_mean ,Data_sem,'lineprops',  '-k');
        title(tea_filenames(i,:))
        hold on
        plot(minLocs, Data_mean(minLocs), '*r')
        
        fig = gcf;
        saveas(fig, [folder '\MeanTrace' index], 'fig')
        close all
    end
end


drugs_baseline_diff = [];
drugs_area_activeBL = [];
drugs_area_minBL = [];
drugs_amp_activeBL = [];
drugs_amp_minBL = [];
for i = 1:size(drugs_filenames,1)
    clearvars folder index stimulus_fn stim_data Data Data_mean Data_sem avg_baseline std_baseline minPks minLocs time_diff min_baseline
    folder = drugs_filenames(i,1:(end-13));
    index = drugs_filenames(i,(end-6):(end-4));
    stimulus_fn=strcat(folder,'\stim',index,'.txt');
    stim_data = load('-ASCII',stimulus_fn);
    
    if (stim_data(1,2) == 400.00) && (stim_data(1,3) == 400.00)
        [Data, Samp] = abfload(drugs_filenames(i,:));
        Data = Data(1:40000,1,:);
        Data = squeeze(Data);
        Data_mean = mean(Data,2);
        Data_sem = std(Data,0,2)/sqrt(size(Data,2));
        
        avg_baseline = mean(Data_mean(35000:40000,1));
        std_baseline = std(Data_mean(35000:40000,1));
        [minPks, minLocs] = findpeaks(-Data_mean,'MinPeakHeight',(-avg_baseline - 0.2*std_baseline),  'MinPeakDistance', 1500, 'MinPeakWidth',500);
        minPks = -minPks;
        
        time_diff = minLocs - 5000;
        min_baseline = Data_mean(minLocs(find(time_diff == min(time_diff))),1);
        
        drugs_baseline_diff = [drugs_baseline_diff ; abs(avg_baseline - min_baseline)];
        
        time_range = minLocs(find(time_diff == min(time_diff))): 23000;
        
        drugs_area_activeBL = [drugs_area_activeBL; trapz(time_range, (Data_mean(time_range,1)-avg_baseline))];
        drugs_area_minBL = [drugs_area_minBL; trapz(time_range, (Data_mean(time_range,1)-min_baseline))];
        
        drugs_amp_activeBL = [drugs_amp_activeBL; (max(Data_mean(time_range,1))-avg_baseline)];
        drugs_amp_minBL = [drugs_amp_minBL; (max(Data_mean(time_range,1))-min_baseline)];
        
        figure
        shadedErrorBar([1:size(Data,1)],Data_mean ,Data_sem,'lineprops',  '-k');
        title(drugs_filenames(i,:))
        hold on
        plot(minLocs, Data_mean(minLocs), '*r')
        
        fig = gcf;
        saveas(fig, [folder '\MeanTrace' index], 'fig')
        close all
    end
end

averages = [mean(wt_area_minBL), mean(wt_area_activeBL); mean(ant_area_minBL), mean(ant_area_activeBL); mean(tea_area_minBL), mean(tea_area_activeBL); mean(drugs_area_minBL), mean(drugs_area_activeBL)];
errors = [std(wt_area_minBL)/sqrt(length(wt_area_minBL)),std(wt_area_activeBL)/sqrt(length(wt_area_activeBL)); std(ant_area_minBL)/sqrt(length(ant_area_minBL)), std(ant_area_activeBL)/sqrt(length(ant_area_activeBL)); std(tea_area_minBL)/sqrt(length(tea_area_minBL)), std(tea_area_activeBL)/sqrt(length(tea_area_activeBL)); std(drugs_area_minBL)/sqrt(length(drugs_area_minBL)), std(drugs_area_activeBL)/sqrt(length(drugs_area_activeBL))];

bar_xtick = errorbar_groups(averages, errors,'bar_width',0.75,'errorbar_width',2, 'bar_colors', [0 0 0; 0 1 0; 1 0 0; 0 0 1]);

averages = [mean(wt_amp_minBL), mean(wt_amp_activeBL); mean(ant_amp_minBL), mean(ant_amp_activeBL); mean(tea_amp_minBL), mean(tea_amp_activeBL); mean(drugs_amp_minBL), mean(drugs_amp_activeBL)];
errors = [std(wt_amp_minBL)/sqrt(length(wt_amp_minBL)),std(wt_amp_activeBL)/sqrt(length(wt_amp_activeBL)); std(ant_amp_minBL)/sqrt(length(ant_amp_minBL)), std(ant_amp_activeBL)/sqrt(length(ant_amp_activeBL)); std(tea_amp_minBL)/sqrt(length(tea_amp_minBL)), std(tea_amp_activeBL)/sqrt(length(tea_amp_activeBL)); std(drugs_amp_minBL)/sqrt(length(drugs_amp_minBL)), std(drugs_amp_activeBL)/sqrt(length(drugs_amp_activeBL))];

bar_xtick = errorbar_groups(averages, errors,'bar_width',0.75,'errorbar_width',2, 'bar_colors', [0 0 0; 0 1 0; 1 0 0; 0 0 1]);

averages = [mean(wt_baseline_diff), mean(ant_baseline_diff), mean(tea_baseline_diff), mean(drugs_baseline_diff)];
errors = [std(wt_baseline_diff)/sqrt(length(wt_baseline_diff)), std(ant_baseline_diff)/sqrt(length(ant_baseline_diff)), std(tea_baseline_diff)/sqrt(length(tea_baseline_diff)), std(drugs_baseline_diff)/sqrt(length(drugs_baseline_diff))];

plot([1 2 3 4],averages, 'Color',[0 0 0], 'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0 0 0]);
hold on
errorbar([1 2 3 4],averages,errors, 'Color',[0 0 0], 'linewidth', 0.5)

wt_baseline_diff = [];
wt_area_activeBL = [];
wt_area_minBL = [];
wt_amp_activeBL = [];
wt_amp_minBL = [];
wt_avg_dataset = [];
for i = 1:size(wt_filenames,1)
    clearvars folder index stimulus_fn stim_data Data Data_mean Data_sem minPks minLocs avg_baseline std_baseline minPks minLocs time_diff min_baseline
    folder = wt_filenames(i,1:(end-13));
    index = wt_filenames(i,(end-6):(end-4));
    stimulus_fn=strcat(folder,'\stim',index,'.txt');
    stim_data = load('-ASCII',stimulus_fn);
    
    if (stim_data(1,2) == 400.00) && (stim_data(1,3) == 400.00) && (stim_data(1,end) >= 100.00)
        [Data, Samp] = abfload(wt_filenames(i,:));
        Data = Data(:,1,:);
        Data = squeeze(Data);
        Data = Data(1:40000,:);
        Data = mean(Data,2);
        Data = lowpass(Data,20, 10000);
        wt_avg_dataset(:,i) = Data;
        
        avg_baseline = mean(Data(32000:38000,1));
        std_baseline = std(Data(32000:38000,1));
        [minPks, minLocs] = findpeaks(-Data,'MinPeakHeight',(-avg_baseline - 0.2*std_baseline),  'MinPeakDistance', 1500, 'MinPeakWidth',500);
        minPks = -minPks;
        
        time_diff = minLocs - 5000;
        min_baseline = Data(minLocs(find(time_diff == min(time_diff))),1);
        
        wt_baseline_diff = [wt_baseline_diff ; avg_baseline];
        
        time_range = minLocs(find(time_diff == min(time_diff))): 23000;
        
        wt_area_activeBL = [wt_area_activeBL; trapz((time_range/10000), (Data(time_range,1)-avg_baseline))];
        wt_area_minBL = [wt_area_minBL; trapz((time_range/10000), (Data(time_range,1)-min_baseline))];
        
        wt_amp_activeBL = [wt_amp_activeBL; (max(Data(time_range,1))-avg_baseline)];
        wt_amp_minBL = [wt_amp_minBL; (max(Data(time_range,1))-min_baseline)];
        
        figure
        plot([1:size(Data,1)],Data,'-k','LineWidth',2);
        title(wt_filenames(i,:))
        hold on
        plot(minLocs, Data(minLocs), '*r')
        
        fig = gcf;
        saveas(fig, [folder '\FilteredMeanTrace' index], 'fig')
        close all
    end
end
writematrix(wt_avg_dataset,'C:\Users\Hector Acaron\Dropbox\2019 Dendrites Meeting\SAC Recordings\Wt_avgTraces.csv', 'delimiter', '\,','newline','pc')

% 
% ant_baseline_diff = [];
% ant_area_activeBL = [];
% ant_area_minBL = [];
% ant_amp_activeBL = [];
% ant_amp_minBL = [];
% for i = 1:size(ant_filenames,1)
%     clearvars folder index stimulus_fn stim_data Data Data_mean Data_sem avg_baseline std_baseline minPks minLocs time_diff min_baseline
%     folder = ant_filenames(i,1:(end-13));
%     index = ant_filenames(i,(end-6):(end-4));
%     stimulus_fn=strcat(folder,'\stim',index,'.txt');
%     stim_data = load('-ASCII',stimulus_fn);
%     
%     if (stim_data(1,2) == 400.00) && (stim_data(1,3) == 400.00)
%         [Data, Samp] = abfload(ant_filenames(i,:));
%         Data = Data(:,1,:);
%         Data = squeeze(Data);
%         Data = lowpass(Data,12, 10000);
%         Data = Data(1:40000,:);
%         Data_mean = mean(Data,2);
%         Data_sem = std(Data,0,2)/sqrt(size(Data,2));
%         
%         avg_baseline = mean(Data_mean(32000:38000,1));
%         std_baseline = std(Data_mean(32000:38000,1));
%         [minPks, minLocs] = findpeaks(-Data_mean,'MinPeakHeight',(-avg_baseline - 0.2*std_baseline),  'MinPeakDistance', 1500, 'MinPeakWidth',500);
%         minPks = -minPks;
%         
%         time_diff = minLocs - 5000;
%         min_baseline = Data_mean(minLocs(find(time_diff == min(time_diff))),1);
%         
%         ant_baseline_diff = [wt_baseline_diff ; avg_baseline];
%         
%         time_range = minLocs(find(time_diff == min(time_diff))): 23000;
%         
%         ant_area_activeBL = [ant_area_activeBL; trapz((time_range/10000), (Data_mean(time_range,1)-avg_baseline))];
%         ant_area_minBL = [ant_area_minBL; trapz((time_range/10000), (Data_mean(time_range,1)-min_baseline))];
%         
%         ant_amp_activeBL = [ant_amp_activeBL; (max(Data_mean(time_range,1))-avg_baseline)];
%         ant_amp_minBL = [ant_amp_minBL; (max(Data_mean(time_range,1))-min_baseline)];
%         
%         figure
%         shadedErrorBar([1:size(Data,1)],Data_mean ,Data_sem,'lineprops',  '-k');
%         title(ant_filenames(i,:))
%         hold on
%         plot(minLocs, Data_mean(minLocs), '*r')
%         
%      
%         
%         fig = gcf;
%         saveas(fig, [folder '\FilteredMeanTrace' index], 'fig')
%         close all
%     end
% end
% 
% tea_baseline_diff = [];
% tea_area_activeBL = [];
% tea_area_minBL = [];
% tea_amp_activeBL = [];
% tea_amp_minBL = [];
% for i = 1:size(tea_filenames,1)
%     clearvars folder index stimulus_fn stim_data Data Data_mean Data_sem avg_baseline std_baseline minPks minLocs time_diff min_baseline
%     folder = tea_filenames(i,1:(end-13));
%     index = tea_filenames(i,(end-6):(end-4));
%     stimulus_fn=strcat(folder,'\stim',index,'.txt');
%     stim_data = load('-ASCII',stimulus_fn);
%     
%     if (stim_data(1,2) == 400.00) && (stim_data(1,3) == 400.00)
%         [Data, Samp] = abfload(tea_filenames(i,:));
%         Data = Data(:,1,:);
%         Data = squeeze(Data);
%         Data = lowpass(Data,12, 10000);
%         Data = Data(1:40000,:);
%         Data_mean = mean(Data,2);
%         Data_sem = std(Data,0,2)/sqrt(size(Data,2));
%        
%         avg_baseline = mean(Data_mean(32000:38000,1));
%         std_baseline = std(Data_mean(32000:38000,1));
%         [minPks, minLocs] = findpeaks(-Data_mean,'MinPeakHeight',(-avg_baseline - 0.2*std_baseline),  'MinPeakDistance', 1500, 'MinPeakWidth',500);
%         minPks = -minPks;
%         
%         time_diff = minLocs - 5000;
%         min_baseline = Data_mean(minLocs(find(time_diff == min(time_diff))),1);
%         
%         tea_baseline_diff = [wt_baseline_diff ; avg_baseline];
%         
%         time_range = minLocs(find(time_diff == min(time_diff))): 23000;
%         
%         tea_area_activeBL = [tea_area_activeBL; trapz((time_range/10000), (Data_mean(time_range,1)-avg_baseline))];
%         tea_area_minBL = [tea_area_minBL; trapz((time_range/10000), (Data_mean(time_range,1)-min_baseline))];
%         
%         tea_amp_activeBL = [tea_amp_activeBL; (max(Data_mean(time_range,1))-avg_baseline)];
%         tea_amp_minBL = [tea_amp_minBL; (max(Data_mean(time_range,1))-min_baseline)];
%         
%         figure
%         shadedErrorBar([1:size(Data,1)],Data_mean ,Data_sem,'lineprops',  '-k');
%         title(tea_filenames(i,:))
%         hold on
%         plot(minLocs, Data_mean(minLocs), '*r')
%         
%         fig = gcf;
%         saveas(fig, [folder '\FilteredMeanTrace' index], 'fig')
%         close all
%     end
% end
% 
% 
% drugs_baseline_diff = [];
% drugs_area_activeBL = [];
% drugs_area_minBL = [];
% drugs_amp_activeBL = [];
% drugs_amp_minBL = [];
% for i = 1:size(drugs_filenames,1)
%     clearvars folder index stimulus_fn stim_data Data Data_mean Data_sem avg_baseline std_baseline minPks minLocs time_diff min_baseline
%     folder = drugs_filenames(i,1:(end-13));
%     index = drugs_filenames(i,(end-6):(end-4));
%     stimulus_fn=strcat(folder,'\stim',index,'.txt');
%     stim_data = load('-ASCII',stimulus_fn);
%     
%     if (stim_data(1,2) == 400.00) && (stim_data(1,3) == 400.00)
%         [Data, Samp] = abfload(drugs_filenames(i,:));
%         Data = Data(:,1,:);
%         Data = squeeze(Data);
%         Data = lowpass(Data,12, 10000);
%         Data = Data(1:40000,:);
%         Data_mean = mean(Data,2);
%         Data_sem = std(Data,0,2)/sqrt(size(Data,2));
%         
%         avg_baseline = mean(Data_mean(32000:38000,1));
%         std_baseline = std(Data_mean(32000:38000,1));
%         [minPks, minLocs] = findpeaks(-Data_mean,'MinPeakHeight',(-avg_baseline - 0.2*std_baseline),  'MinPeakDistance', 1500, 'MinPeakWidth',500);
%         minPks = -minPks;
%         
%         time_diff = minLocs - 5000;
%         min_baseline = Data_mean(minLocs(find(time_diff == min(time_diff))),1);
%         
%         drugs_baseline_diff = [wt_baseline_diff ; avg_baseline];
%         
%         time_range = minLocs(find(time_diff == min(time_diff))): 23000;
%         
%         drugs_area_activeBL = [drugs_area_activeBL; trapz((time_range/10000), (Data_mean(time_range,1)-avg_baseline))];
%         drugs_area_minBL = [drugs_area_minBL; trapz((time_range/10000), (Data_mean(time_range,1)-min_baseline))];
%         
%         drugs_amp_activeBL = [drugs_amp_activeBL; (max(Data_mean(time_range,1))-avg_baseline)];
%         drugs_amp_minBL = [drugs_amp_minBL; (max(Data_mean(time_range,1))-min_baseline)];
%         
%         figure
%         shadedErrorBar([1:size(Data,1)],Data_mean ,Data_sem,'lineprops',  '-k');
%         title(drugs_filenames(i,:))
%         hold on
%         plot(minLocs, Data_mean(minLocs), '*r')
%         
%         fig = gcf;
%         saveas(fig, [folder '\FilteredMeanTrace' index], 'fig')
%         close all
%     end
% end
% 
% averages = [mean(wt_area_minBL), mean(wt_area_activeBL); mean(ant_area_minBL), mean(ant_area_activeBL); mean(tea_area_minBL), mean(tea_area_activeBL); mean(drugs_area_minBL), mean(drugs_area_activeBL)];
% errors = [std(wt_area_minBL)/sqrt(length(wt_area_minBL)),std(wt_area_activeBL)/sqrt(length(wt_area_activeBL)); std(ant_area_minBL)/sqrt(length(ant_area_minBL)), std(ant_area_activeBL)/sqrt(length(ant_area_activeBL)); std(tea_area_minBL)/sqrt(length(tea_area_minBL)), std(tea_area_activeBL)/sqrt(length(tea_area_activeBL)); std(drugs_area_minBL)/sqrt(length(drugs_area_minBL)), std(drugs_area_activeBL)/sqrt(length(drugs_area_activeBL))];
% 
% bar_xtick = errorbar_groups([mean(wt_area_minBL), mean(wt_area_activeBL); mean(ant_area_minBL), mean(ant_area_activeBL)] ,[std(wt_area_minBL)/sqrt(length(wt_area_minBL)),std(wt_area_activeBL)/sqrt(length(wt_area_activeBL)); std(ant_area_minBL)/sqrt(length(ant_area_minBL)), std(ant_area_activeBL)/sqrt(length(ant_area_activeBL))],'bar_width',0.75,'errorbar_width',2, 'bar_colors', [1 0 0; 0 0 1]);
% hold on
% set(gca,'XTickLabel',{'Minimum Baseline','Active Baseline'},'fontweight','bold')
% set(gcf,'color','white')
% ylabel(['Area'],'FontSize',12,'FontWeight','bold','Color','k')
% plot ([bar_xtick(1) - 0.375]*ones(length(wt_area_minBL),1), wt_area_minBL,'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0 0 0]);
% plot ([bar_xtick(1) + 0.375]*ones(length(ant_area_minBL),1), ant_area_minBL,'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0 0 0]);
% plot ([bar_xtick(2) - 0.375]*ones(length(wt_area_activeBL),1), wt_area_activeBL,'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0 0 0]);
% plot ([bar_xtick(2) + 0.375]*ones(length(ant_area_activeBL),1), ant_area_activeBL,'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0 0 0]);
% [h_area1, p_area1] = ttest2(wt_area_minBL, ant_area_minBL);
% [h_area2, p_area2] = ttest2(wt_area_activeBL, ant_area_activeBL);
% 
% 
% bar_xtick = errorbar_groups([mean(wt_amp_minBL), mean(wt_amp_activeBL); mean(ant_amp_minBL), mean(ant_amp_activeBL)] ,[std(wt_amp_minBL)/sqrt(length(wt_amp_minBL)),std(wt_amp_activeBL)/sqrt(length(wt_amp_activeBL)); std(ant_amp_minBL)/sqrt(length(ant_amp_minBL)), std(ant_amp_activeBL)/sqrt(length(ant_amp_activeBL))],'bar_width',0.75,'errorbar_width',2, 'bar_colors', [1 0 0; 0 0 1]);
% hold on
% set(gca,'XTickLabel',{'Minimum Baseline','Active Baseline'},'fontweight','bold')
% set(gcf,'color','white')
% ylabel(['amp'],'FontSize',12,'FontWeight','bold','Color','k')
% plot ([bar_xtick(1) - 0.375]*ones(length(wt_amp_minBL),1), wt_amp_minBL,'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0 0 0]);
% plot ([bar_xtick(1) + 0.375]*ones(length(ant_amp_minBL),1), ant_amp_minBL,'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0 0 0]);
% plot ([bar_xtick(2) - 0.375]*ones(length(wt_amp_activeBL),1), wt_amp_activeBL,'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0 0 0]);
% plot ([bar_xtick(2) + 0.375]*ones(length(ant_amp_activeBL),1), ant_amp_activeBL,'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0 0 0]);
% [h_amp1, p_amp1] = ttest2(wt_amp_minBL, ant_amp_minBL);
% [h_amp2, p_amp2] = ttest2(wt_amp_activeBL, ant_amp_activeBL);
% 
% averages = [mean(wt_area_minBL); mean(ant_area_minBL); mean(tea_area_minBL); mean(drugs_area_minBL)];
% errors = [std(wt_area_minBL)/sqrt(length(wt_area_minBL)); std(ant_area_minBL)/sqrt(length(ant_area_minBL)); std(tea_area_minBL)/sqrt(length(tea_area_minBL)); std(drugs_area_minBL)/sqrt(length(drugs_area_minBL))];
% 
% bar_xtick = errorbar_groups(averages, errors,'bar_width',1,'errorbar_width',0.75, 'bar_colors', [0 0 0; 0 1 0; 1 0 0; 0 0 1]);
% hold on
% plot ([1]*ones(length(wt_area_minBL),1), wt_area_minBL,'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0 0 0]);
% plot ([2]*ones(length(ant_area_minBL),1), ant_area_minBL,'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0 0 0]);
% plot ([3]*ones(length(tea_area_minBL),1), tea_area_minBL,'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0 0 0]);
% plot ([4]*ones(length(drugs_area_minBL),1), drugs_area_minBL,'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0 0 0]);
% hold off
% 
% 
% 
% 
% averages = [mean(wt_amp_minBL); mean(ant_amp_minBL); mean(tea_amp_minBL); mean(drugs_amp_minBL)];
% errors = [std(wt_amp_minBL)/sqrt(length(wt_amp_minBL)); std(ant_amp_minBL)/sqrt(length(ant_amp_minBL)); std(tea_amp_minBL)/sqrt(length(tea_amp_minBL)); std(drugs_amp_minBL)/sqrt(length(drugs_amp_minBL))];
% 
% bar_xtick = errorbar_groups(averages, errors,'bar_width',1,'errorbar_width',2, 'bar_colors', [0 0 0; 0 1 0; 1 0 0; 0 0 1]);
% hold on
% plot ([1]*ones(length(wt_amp_minBL),1), wt_amp_minBL,'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0 0 0]);
% plot ([2]*ones(length(ant_amp_minBL),1), ant_amp_minBL,'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0 0 0]);
% plot ([3]*ones(length(tea_amp_minBL),1), tea_amp_minBL,'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0 0 0]);
% plot ([4]*ones(length(drugs_amp_minBL),1), drugs_amp_minBL,'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0 0 0]);
% hold off
% 
% averages = [mean(wt_baseline_diff), mean(ant_baseline_diff), mean(tea_baseline_diff), mean(drugs_baseline_diff)];
% errors = [std(wt_baseline_diff)/sqrt(length(wt_baseline_diff)), std(ant_baseline_diff)/sqrt(length(ant_baseline_diff)), std(tea_baseline_diff)/sqrt(length(tea_baseline_diff)), std(drugs_baseline_diff)/sqrt(length(drugs_baseline_diff))];
% 
% plot([1 2 3 4],averages, 'Color',[0 0 0], 'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0 0 0]);
% hold on
% errorbar([1 2 3 4],averages,errors, 'Color',[0 0 0], 'linewidth', 0.5)
% 
