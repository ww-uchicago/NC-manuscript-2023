function [DSData1 DSData2 DSData3 DSData4 DSData5 DSData6 filtData] = dendloadClampexDataScript(varargin)
addpath('/Users/hectoracaron/Desktop/Research/Wei Lab HD/recordings/wei_codes/loadingTraces');
common_path = '/Users/hectoracaron/Desktop/Research/Wei Lab HD/recordings/New Rig/20150121/';
close all
index_250 = 012; 
index_50 = 012 ;
index_loc1 = 014 ;
index_loc2 = 016 ;
index_loc3 = 018 ;
index_loc4 = 020 ;
fidxs = [index_250, index_50 index_loc1, index_loc2, index_loc3, index_loc4];
for k = 1:6
    fidx = num2str(fidxs(k));
    stype=2;
    switch stype
        case 1
            stimulus_type = 'triggered_steps';
        case 2
            stimulus_type = 'moving_bars';
        case 3
            stimulus_type = 'drifting_gratings';
    end

    channel = 2; %channel to analyze on digidata
    show_repetitions = 1;
    data_start = 0; %seconds
    data_stop = inf; %seconds
    start_rep = 1;
    stop_rep = 3; %change to match # of repetions
    spike_threshold = -4.0; %standard deviation, -4 is good
    show_spike_finding_plots = 1;
    close_spike_finding_plots = 1;
    show_output_plots = 1;
    save_output_plot = 1;
    data_fn = strcat('151210',fidx,'.abf');
    stimulus_fn=strcat('stim0',fidx,'.txt');
    save_plot_fn=strcat('fig0',fidx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if strcmpi(data_stop,'inf');
        data_stop = inf;
    end

    data_path = [common_path data_fn];
    save_plot_path = [common_path save_plot_fn];
    stimulus_path = [common_path stimulus_fn];
    [Data,SampInt]=abfload(data_path);
    %added this line to choose channel: Justin 2009-08-26
    Data = Data(:,channel,:);
    Data = squeeze(Data);
    filtData = zeros(size(Data));

    %filter data, and load spiketimes into a structure
    SpikeTimes = cell(size(Data,2),1);
    for c = 1:size(Data,2)
        filtData(:,c) = filterData(Data(:,c),SampInt);
        SpikeTimes{c} = findSpikeTimes(filtData(:,c),SampInt,'order',2,...
            'thresh',spike_threshold,'show_plot',show_spike_finding_plots,...
            'show_plot_thresh',1);
        if show_spike_finding_plots
            legend(sprintf('trial %d',c));
        end
        %close after every 12 plots
        if close_spike_finding_plots && mod(c,13)==0
            close all
        end
    end
    %close all

    switch stimulus_type
        case 'drifting_gratings'
            %remove time if requested
            clippedSpikeTimes = subseqSpikeTimes(SpikeTimes,'start',data_start,'stop',data_stop);
            %sort data according to directions and plot if DS data
            stimulus_data=load('-ASCII',stimulus_path);
            directions = stimulus_data(:,1);
            SortedDSData = sortDSdata(clippedSpikeTimes,directions,'first_repetition',start_rep,'last_repetition',stop_rep);
            if show_output_plots
                makePSTH(clippedSpikeTimes,0,3,0.01);
                figure
                makePolarPlot(SortedDSData,'show_repetitions',1,'show_repetitions',show_repetitions);
            end

        case 'moving_bars'
            %remove first 0.3 secs form each trial
            clippedSpikeTimes = subseqSpikeTimes(SpikeTimes,'start',data_start,'stop',data_stop);
            %sort data according to directions and plot if DS data
            stimulus_data=load('-ASCII',stimulus_path);
            directions = stimulus_data(:,1);
            SortedDSData = sortDSdata(clippedSpikeTimes,directions,'first_repetition',start_rep,'last_repetition',stop_rep);
            if k == 1 
                DSData1 = SortedDSData;
            elseif k == 2
                DSData2 = SortedDSData;
            elseif k == 3
                DSData3 = SortedDSData;
            elseif k == 4
                DSData4 = SortedDSData;
            elseif k == 5
                DSData5 = SortedDSData;
            elseif k == 6
                DSData6 = SortedDSData;
            end
            if show_output_plots
                makePSTH(clippedSpikeTimes,0,2.4,0.01);
                figure
                makePolarPlot(SortedDSData,'show_repetitions',1);
            end

        case 'triggered_steps'
            %if step data, generate PSTH
            clippedSpikeTimes = subseqSpikeTimes(SpikeTimes,'start',data_start,'stop',data_stop);
            %[PSTH bins h] = makePSTH(clippedSpikeTimes,0,6,0.05,'y_range',[0 120]);
            [PSTH bins h] = makePSTH(clippedSpikeTimes,0,6,0.05);

        case 'continuous_trace'
            [h ax] = makeTracePlot(Data,SampInt,'font_size',18);

        otherwise
            error('stimulus_type does not match any cases')
    end
    save(strcat('exps_',fidx,'_data','SortedDsData'));
    if save_output_plot
        set(gcf,'PaperPositionMode','auto')
        print(gcf,save_plot_path,'-depsc2')
        saveas(gcf,save_plot_path,'fig')
    end
end

return 

end


