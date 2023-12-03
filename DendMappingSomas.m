% This serves as the main fu
close all; clear all;clearvars;  clc;
%% setup parameters
paras.vfreq = 1000;   % TTL sampling rate, 1000, or 10000
paras.rmrep=[2]; % remove certain reps of recording
%% main part
folder = 'C:\Users\Hector Acaron\Dropbox\2019 Dendrites Meeting';
%pname=uigetdir(folder);
%read the xml, and tif file and calculate the mean intensity value over
%time

directory = 'J';
%%%%%%%%%%%%%%%%%%%% dend %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen([folder '\WT_VoltageImaging_20200629.txt'], 'rt');
dend_Info = textscan(fid, '%f %f %f %f %f %f %f %s %f', 'HeaderLines', 0, 'CollectOutput', true);
fclose(fid);
dend_Info{1,2} = dend_Info{1,2}(find(dend_Info{1,1}(:,6) == 0),:);
dend_Info{1,3} = dend_Info{1,3}(find(dend_Info{1,1}(:,6) == 0),:);
dend_Info{1,1} = dend_Info{1,1}(find(dend_Info{1,1}(:,6) == 0),:);

dend_area = unique(dend_Info{1,1}(:,7));
xmlFactor = 0.928505; %microns per XML units

avgShadedPlots = input('Average Plots?\n');

PPT = input('Generate ppt?\n');

dend_events = input('Event Detection?\n');




for i = 1:length(dend_area)
    close all
    clearvars -except paras xmlFactor folder moviefiles i dend_area num_stimuli compare analysis avgShadedPlots  dend_events PPT date mouse age slidesFile slides directory var_pairedInfo var_Info dend_pairedInfo dend_Info dend_data
    movieInf = find(dend_Info{1,1}(:,7) == dend_area(i));
    paths = dend_Info{1,2}((dend_Info{1,1}(:,7) == dend_area(i)),:);
    for t = 1:length(movieInf)
        clearvars -except paras xmlFactor folder moviefiles i dend_area num_stimuli compare analysis avgShadedPlots dend_events PPT t movieInf cellSignals cp_traces cf_traces dend_PSTH edges window_size time_window SortedPeakData time date mouse age slidesFile slides directory var_pairedInfo var_Info dend_pairedInfo dend_Info paths dend_data
        movieNum = dend_Info{1,1}(movieInf(t),1);
        moviestring = num2str(movieNum);
        pathstr = paths(t,:);
        pname = dir([char(pathstr) '\TSeries*' moviestring]);
        tseriesxml = [pname.name '.xml'];
        pname = [pname.folder '\' pname.name];
        cellSignalString = strcat('CellSignals',moviestring,'.txt');

        %% IMPORT TIFF FILES%%%%
       
        roi_select = dend_Info{1,1}(movieInf(t),2);
        cp_dir = dend_Info{1,1}(movieInf(t),3);
        cf_dir = dend_Info{1,1}(movieInf(t),4);
        exp_type = dend_Info{1,1}(movieInf(t),5);
        dend_orient = dend_Info{1,3}(movieInf(t),1);
        eye = 1;
        piece = 1;
        barType = 2;
        interTrialWait = 0.5;
        %%IMPORT CELL SIGNALS%%%%%

        cellSignals{1,t} = dlmread([char(pathstr) '\CellSignals\' cellSignalString]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        if (barType == 1 || barType == 2 || barType == 3)
            %%%%%IMPORT VOLTAGE RECORDING & STIMULUS FILE%%%%%%%%%%%%%%%%
            fp.vol=dir([pname '\*.csv']);   fp.vol=([pname '\' fp.vol.name]);
            %fp.trace=dir([folder '\*.xls']); fp.trace=([folder '\' fp.trace.name]);
            fp.stim=dir([pname '\*stim*.txt']);
            vFreq=10000; 
            save_plots = 1;
            locs=[1 200 1400 800];
            % load the raw data
            numbg=1;   % number of background ROIs
            try
                data.stim=dlmread([pname '\' fp.stim(1).name]);
                fp.stim=([pname '\' fp.stim(1).name]);
            catch
                    data.stim=dlmread([pname '\' fp.stim(2).name]);
                    fp.stim=([pname '\' fp.stim(2).name]);
            end
            stimInf.speed=data.stim(1,2); stimInf.bw=data.stim(1,3);
            stimInf.bl=data.stim(1,4); stimInf.rad=data.stim(1,5); data.stim=data.stim(:,1); stimInf.dirs=data.stim;
            data.vol=csvread(fp.vol,1); data.vol_time = data.vol(:,1); data.vol=data.vol(:,2); 
            data.trace= cellSignals{1,t};
            fp.xml = dir([pname '\*' moviestring '.xml']); fp.xml =([pname '\' fp.xml.name]);

            data.bg=data.trace(:,end);
            %data.trace_smooth = DataSmooth(data.trace, 6, 'average');
            %data.trace_sub = data.trace_smooth(:,1)- data.trace_smooth(:,2);

            %%
            rois=1:(size(data.trace,2) -1 );
            %rois=[6];
            rm_rep= [ ];
            %process the voltage data
            [data.vlocs,data.pks]=peakdetect(data.vol,1,4000);


            if length(data.vlocs) < length(stimInf.dirs)
                loc_diff = diff(data.vlocs);
                location = find(loc_diff > (mean(loc_diff) + 2*std(loc_diff)));
                vlocs_2 = data.vlocs;
                data.vlocs = zeros(length(stimInf.dirs),1);
                data.vlocs(1:location,1) = vlocs_2(1:location,1);
                interval = floor(mean(diff(vlocs_2)));
                data.vlocs((location + 1),1) = vlocs_2(location,1) + interval;
                data.vlocs((location + 2):length(stimInf.dirs),1) = vlocs_2((location + 1):(length(stimInf.dirs) -1),1);
            end


            dlist=unique(data.stim);

            % determine the truncation window size

            fRate=length(data.vol)/size(data.trace,1);

            data.fFrame=round((data.vlocs/fRate),0); %data.fFrame(1)=[];

            %f = fit([data.fFrame(1):numel(data.trace_sub(:,1))]',data.trace_sub(data.fFrame(1):end,1), 'exp2');
            data.timestamps = data.vol_time(data.vlocs)*0.001;
            wsize=round(mean(diff(data.fFrame)));
            time_bin = mean(diff(data.timestamps));
            
            stim_time = time_bin;

            fid = fopen(fp.xml);
            tline=fgets(fid);
            while ischar(tline)
                 if ~isempty(strfind(tline, 'framePeriod')) 
                     expression = '\"';
                     tok=regexp(tline,expression, 'split');
                     FramePeriod = str2num(tok{1,4});
                     tline = 0;
                 else
                     tline=fgets(fid);
                 end
            end
            clearvars fid tline tok expression

            framerate = 1/FramePeriod;

            wsize_stim = wsize - round(framerate*(0.8)*interTrialWait);

            %%%%%%   FIT   %%%%%%%
            rng(0,'twister');
            
            if exp_type == 1
                smooth_window = 3;
            elseif exp_type == 2
                smooth_window = 12;
            end


            data.raw = data.trace(1:(data.fFrame(end)+wsize),:);
            data.raw(:,end) = DataSmooth(data.raw(:,end),(smooth_window*2), 'average');

            data.trace = [data.raw(:,(1:(end-1))) - data.raw(:,end)];
            y_fit = data.trace;
            y_fit2 = y_fit;

            y_fit2 = DataSmooth(y_fit2,smooth_window, 'average');

            for j = 1: length(data.fFrame)
                y_fit2(data.fFrame(j):(data.fFrame(j)+wsize_stim),:)= NaN;
            end


            if exp_type == 1
                y_fit2(1:(round(framerate*4)),:)= NaN;
            end

            x_fit = 1:length(y_fit2);
            x_fit = x_fit(find(isnan(y_fit2(:,1)) == 0));
           
            y_fit3 = y_fit2(find(isnan(y_fit2(:,1)) == 0),:);

            for i = 1:size(data.trace,2)

                 try
                    [f2, g2, o2] = fit([x_fit]', [y_fit3(:,i)], 'exp2');
                    [f2_e1, g2_e1, o2_e1] = fit([x_fit]', [y_fit3(:,i)], 'exp1');
                    if (g2_e1.rsquare > g2.rsquare) || (g2_e1.adjrsquare > g2.adjrsquare)
                        clearvars f2 g2 o2
                        f2 = f2_e1;
                        g2 = g2_e1;
                        o2 = o2_e1;
                    end
                catch

                    try 
                        [f2, g2, o2] = fit([x_fit]', [y_fit3(:,i)], 'exp2');
                    catch
                        [f2, g2, o2] = fit([x_fit]', [y_fit3(:,i)], 'exp1');
                    end
                end

                figure
                x_trace = [1:length(data.trace(:,i))];
                plot(x_trace, data.trace(:,i), '-b', 'Linewidth', 1)
                hold on
                plot(x_trace, f2(x_trace), '-r', 'Linewidth', 2)
                hold off



                if exp_type == 2
                    baseline = f2(x_trace);
                    clipframe = find(baseline <= 0.1);
                    if isempty(clipframe) == 1
                        clip_trace(i) = 0;
                    elseif isempty(clipframe) == 0
                        clip_trace(i) = clipframe(1);
                    end
                else
                    clip_trace(i) = 0;
                end


                close all

                for k = 1:length(x_trace)

                    data.trace(k,i) = (data.trace(k,i) - f2(x_trace(k)))/f2(x_trace(k));
                end


                clearvars f2 g2 o2 baseline clipframe
            end


            data.signal=data.trace(:,roi_select);
            

            if exp_type == 1
               new_FR = 75;
            elseif exp_type == 2
               new_FR = 1000;
            end

            time_tx = linspace(FramePeriod, (FramePeriod*length(data.signal)), length(data.signal));
            for u = 1: size(data.signal,2)

                data.signalSmooth(:,u) = resample(data.signal(:,u),time_tx,  new_FR, 1, 1, 'pchip');
                data.signalSmooth(:,u) = DataSmooth(data.signalSmooth(:,u),(2*smooth_window), 'average');
                data.signalSmooth(:,u) = DataSmooth(data.signalSmooth(:,u),(smooth_window), 'average');
                data.signalSmooth(:,u) = DataSmooth(data.signalSmooth(:,u),(smooth_window), 'average');
            end

            loc_diff = round(mean(diff(data.vlocs)));

            data.vol = data.vol(1:(data.vlocs(end)+loc_diff));

            fRate_resamp=length(data.vol)/size(data.signalSmooth,1);

            data.fFrameResamp=round((data.vlocs/fRate_resamp),0); %data.fFrame(1)=[];

            %f = fit([data.fFrame(1):numel(data.trace_sub(:,1))]',data.trace_sub(data.fFrame(1):end,1), 'exp2');

            wsize_resamp=round(mean(diff(data.fFrameResamp)));

            new_FR = wsize_resamp/(time_bin);


            % CLIP DATA SET FOR EACH SWEEP%%%%
            
            data.signal=clip(data.signal,data.fFrame,wsize);
            data.bg=clip(data.bg,data.fFrame,wsize);

            data.signalSmooth=clip(data.signalSmooth,data.fFrameResamp,wsize_resamp);
            
            
            if clip_trace(roi_select) > 0
                frame_diffs = data.fFrame - clip_trace(roi_select);
                sweep_number = find(frame_diffs > 0);
                sweep_number = sweep_number(1);
                if mod(sweep_number,2) == 0
                    sweep_number = sweep_number - 1;
                    data.signal{1,1}(:,sweep_number:end) = [];
                    reps = [1:(sweep_number-1)]';
                elseif mod(sweep_number,2) == 1
                    data.signal{1,1}(:,sweep_number:end) = [];
                    reps = [1:(sweep_number-1)]';
                end
            elseif clip_trace(roi_select) == 0
                reps=reshape(1:length(stimInf.dirs),length(dlist),(length(stimInf.dirs))/length(dlist)); 
                reps=reshape(reps,size(reps,1)*size(reps,2),1); 
            end

            



            %data.fFrame=data.fFrame(reps);
            data.stim=data.stim(reps);

            % ORGANIZE BY DIRECTION
            
            dlist=unique(data.stim);
            sortd=zeros(length(data.stim)/length(dlist),length(dlist));
            for i=1:length(dlist)
                sortd(:,i)= find(dlist(i)==data.stim);
            end
            sortd=reshape(sortd,1,size(sortd,1)*size(sortd,2));
            data.signal=cellfun(@(x) x(:,sortd),data.signal,'uniformoutput',false);
            data.bg=cellfun(@(x) x(:,sortd),data.bg,'uniformoutput',false);

            data.signalSmooth=cellfun(@(x) x(:,sortd),data.signalSmooth,'uniformoutput',false);


            %%%%%%%%%%%%%%%%%%% CALCULATE QI (QUALITY INDEX) %%%%%%%%%%%%%%%%%%%%%%%
            dend_data.window_size{1,t} = wsize_resamp;
            dend_data.time_window{1,t} = (2*stimInf.rad +stimInf.bw)/stimInf.speed + interTrialWait;
            
            cp_inf = find(dlist == cp_dir);
            cf_inf = find(dlist == cf_dir);
            
            cp_end = cp_inf*(length(data.stim)/2);
            cp_start = cp_end - (length(data.stim)/2 -1);
            
            cf_end = cf_inf*(length(data.stim)/2);
            cf_start = cf_end - (length(data.stim)/2 -1);
            if exp_type == 1
                mult_deltaF = 1;
            elseif exp_type == 2
                mult_deltaF = -1;
            end
            %%%%%% CP and CF traces %%%%%%%%%%%%%%
            dend_data.cp_traces{1,t} = mult_deltaF*data.signalSmooth{1,1}(:,[cp_start:cp_end]);
            dend_data.cf_traces{1,t} = mult_deltaF*data.signalSmooth{1,1}(:,[cf_start:cf_end]);

            dend_data.cp_mean{1,t} = mean(dend_data.cp_traces{1,t},2);
            dend_data.cp_err{1,t} = std(dend_data.cp_traces{1,t},0,2)/sqrt(size(dend_data.cp_traces{1,t},2));
                
            dend_data.cf_mean{1,t} = mean(dend_data.cf_traces{1,t},2);
            dend_data.cf_err{1,t} = std(dend_data.cf_traces{1,t},0,2)/sqrt(size(dend_data.cf_traces{1,t},2));
            
            [r, c] = xcorr(dend_data.cp_mean{1,t}, dend_data.cf_mean{1,t});
            dend_data.shift(movieInf(t),1) = c(find(r == max(r)));
            spatial_shift = ((dend_data.shift(movieInf(t),1)/new_FR)*(stimInf.speed*1.07))/2;
            fid = fopen(fp.xml);
            tline=fgets(fid);
            while ischar(tline)
                 if ~isempty(strfind(tline, 'positionCurrent')) 
                     tline=fgets(fid);
                     tline=fgets(fid);
                     expression = '\"';
                     tok=regexp(tline,expression, 'split');
                     x_pos = str2num(tok{1,4});
                     tline=fgets(fid);
                     tline=fgets(fid);
                     tline=fgets(fid);
                     tok=regexp(tline,expression, 'split');
                     y_pos = str2num(tok{1,4});
                     tline = 0;
                 else
                     tline=fgets(fid);
                 end
            end
            dend_data.position(movieInf(t),1:2) = [x_pos y_pos];
            
            if spatial_shift > 0
                if dend_orient < 180
                    angle_shift = dend_orient + 180;
                elseif dend_orient > 180
                    angle_shift = dend_orient - 180;
                end
            elseif spatial_shift <= 0
                angle_shift = dend_orient;
            end
            x_shift = abs(spatial_shift)*cosd(angle_shift);
            y_shift = abs(spatial_shift)*sind(angle_shift);
            
            y_shift = y_shift*(1/xmlFactor);
            x_shift = x_shift*(1/xmlFactor);
            dend_data.center(movieInf(t),1:2) = [(x_pos + x_shift) (y_pos + y_shift)];
            dend_data.orient(movieInf(t),1) = dend_orient;
            clearvars r c tline tok 
            
            dend_data.QI(movieInf(t),1) = (var(dend_data.cf_mean{1,t}))/(mean(var(dend_data.cf_traces{1,t},0,1)));
            
            dend_data.time_axis{1,t} = linspace(0, dend_data.time_window{1,t}, dend_data.window_size{1,t})';
            dend_data.search_indexes{1,t} = intersect(find(dend_data.time_axis{1,t} >= 0.35), find(dend_data.time_axis{1,t} <= 2));
            
            [dend_data.cp_pk(movieInf(t),1), dend_data.cp_loc(movieInf(t),1)] = findpeaks(dend_data.cp_mean{1,t}(dend_data.search_indexes{1,t},1),'SortStr', 'descend', 'NPeaks',1);
            [dend_data.cf_pk(movieInf(t),1), dend_data.cf_loc(movieInf(t),1)] = findpeaks(dend_data.cf_mean{1,t}(dend_data.search_indexes{1,t},1),'SortStr', 'descend', 'NPeaks',1);
            dend_data.cf_pkindex(movieInf(t),1) = dend_data.search_indexes{1,t}(1) + dend_data.cf_loc(movieInf(t),1) - 1;
            dend_data.cp_pkindex(movieInf(t),1) = dend_data.search_indexes{1,t}(1) + dend_data.cp_loc(movieInf(t),1) - 1;
            
            [dend_data.cp_min(movieInf(t),1), dend_data.cp_minloc(movieInf(t),1)] = findpeaks((-1)*dend_data.cp_mean{1,t}((dend_data.search_indexes{1,t}(1):dend_data.cp_pkindex(movieInf(t),1)),1),'SortStr', 'descend', 'NPeaks',1);
            [dend_data.cf_min(movieInf(t),1), dend_data.cf_minloc(movieInf(t),1)] = findpeaks((-1)*dend_data.cf_mean{1,t}((dend_data.search_indexes{1,t}(1):dend_data.cf_pkindex(movieInf(t),1)),1),'SortStr', 'descend', 'NPeaks',1);
            
            dend_data.cp_min(movieInf(t),1) = (-1)*dend_data.cp_min(movieInf(t),1);
            dend_data.cf_min(movieInf(t),1) = (-1)*dend_data.cf_min(movieInf(t),1);
            
            dend_data.cf_minindex(movieInf(t),1) = dend_data.search_indexes{1,t}(1) + dend_data.cf_minloc(movieInf(t),1) - 1;
            dend_data.cp_minindex(movieInf(t),1) = dend_data.search_indexes{1,t}(1) + dend_data.cp_minloc(movieInf(t),1) - 1;
            
            dend_data.cp_deltaF(movieInf(t),1) = abs(dend_data.cp_pk(movieInf(t),1) - (-1)*dend_data.cp_min(movieInf(t),1));
            dend_data.cf_deltaF(movieInf(t),1) = abs(dend_data.cf_pk(movieInf(t),1) - (-1)*dend_data.cf_min(movieInf(t),1));
            
            dend_data.cp_areaTrace{1,t} = dend_data.cp_mean{1,t}(dend_data.search_indexes{1,t},:) - dend_data.cp_min(movieInf(t),1);
            dend_data.cf_areaTrace{1,t} = dend_data.cf_mean{1,t}(dend_data.search_indexes{1,t},:) - dend_data.cf_min(movieInf(t),1);
            
            dend_data.cp_Area(movieInf(t),1) = trapz(dend_data.cp_areaTrace{1,t})*FramePeriod;
            dend_data.cf_Area(movieInf(t),1) = trapz(dend_data.cf_areaTrace{1,t})*FramePeriod;
            
            
            dend_data.dsi(movieInf(t),1) = (dend_data.cf_deltaF(movieInf(t),1) - dend_data.cp_deltaF(movieInf(t),1))/(dend_data.cf_deltaF(movieInf(t),1) + dend_data.cp_deltaF(movieInf(t),1));
            
            dend_data.dsiArea(movieInf(t),1) = (dend_data.cf_Area(movieInf(t),1) - dend_data.cp_Area(movieInf(t),1))/(dend_data.cf_Area(movieInf(t),1) + dend_data.cp_Area(movieInf(t),1));
            
            
            dend_data.cf_rise{1,t} = intersect(find(dend_data.cf_mean{1,t}(dend_data.cf_minindex(movieInf(t),1):dend_data.cf_pkindex(movieInf(t),1),1) >= (0.10)*dend_data.cf_pk(movieInf(t),1)), find(dend_data.cf_mean{1,t}(dend_data.cf_minindex(movieInf(t),1):dend_data.cf_pkindex(movieInf(t),1),1) <= (0.90)*dend_data.cf_pk(movieInf(t),1)));
            if isempty(dend_data.cf_rise{1,t})
                dend_data.cf_rise{1,t} = dend_data.cf_minindex(movieInf(t),1):dend_data.cf_pkindex(movieInf(t),1);
            end
            
            dend_data.cp_rise{1,t} = intersect(find(dend_data.cp_mean{1,t}(dend_data.cp_minindex(movieInf(t),1):dend_data.cp_pkindex(movieInf(t),1),1) >= (0.10)*dend_data.cp_pk(movieInf(t),1)), find(dend_data.cp_mean{1,t}(dend_data.cp_minindex(movieInf(t),1):dend_data.cp_pkindex(movieInf(t),1),1) <= (0.90)*dend_data.cp_pk(movieInf(t),1)));
            if isempty(dend_data.cp_rise{1,t})
                dend_data.cp_rise{1,t} = dend_data.cp_minindex(movieInf(t),1):dend_data.cp_pkindex(movieInf(t),1);
            end
            
            dend_data.cf_risestart(movieInf(t),1) = dend_data.cf_minindex(movieInf(t),1) + dend_data.cf_rise{1,t}(1) - 1;
            dend_data.cf_risestop(movieInf(t),1) = dend_data.cf_minindex(movieInf(t),1) + dend_data.cf_rise{1,t}(end) - 1;
            dend_data.cf_risetime(movieInf(t),1) = dend_data.time_axis{1,t}(dend_data.cf_risestop(movieInf(t),1)) - dend_data.time_axis{1,t}(dend_data.cf_risestart(movieInf(t),1));
            dend_data.cf_riseslope(movieInf(t),1) = (0.90*dend_data.cf_pk(movieInf(t),1) - 0.10*dend_data.cf_pk(movieInf(t),1))/dend_data.cf_risetime(movieInf(t),1);
            dend_data.cf_risestart(movieInf(t),1) = dend_data.time_axis{1,t}(dend_data.cf_risestart(movieInf(t),1));
            
            dend_data.cp_risestart(movieInf(t),1) = dend_data.cp_minindex(movieInf(t),1) + dend_data.cp_rise{1,t}(1) - 1;
            dend_data.cp_risestop(movieInf(t),1) = dend_data.cp_minindex(movieInf(t),1) + dend_data.cp_rise{1,t}(end) - 1;
            dend_data.cp_risetime(movieInf(t),1) = dend_data.time_axis{1,t}(dend_data.cp_risestop(movieInf(t),1)) - dend_data.time_axis{1,t}(dend_data.cp_risestart(movieInf(t),1));
            dend_data.cp_riseslope(movieInf(t),1) = (0.90*dend_data.cp_pk(movieInf(t),1) - 0.10*dend_data.cp_pk(movieInf(t),1))/dend_data.cp_risetime(movieInf(t),1);
            dend_data.cp_risestart(movieInf(t),1) = dend_data.time_axis{1,t}(dend_data.cp_risestart(movieInf(t),1));
            
            if dend_Info{1,1}(movieInf(t),6) == 0
                dend_data.bin_num(movieInf(t),1) = 1;
            elseif dend_Info{1,1}(movieInf(t),6) < 0.2
                dend_data.bin_num(movieInf(t),1) = 2;
            elseif (dend_Info{1,1}(movieInf(t),6) >= 0.2) && (dend_Info{1,1}(movieInf(t),6) < 0.4)
                dend_data.bin_num(movieInf(t),1) = 3;
            elseif (dend_Info{1,1}(movieInf(t),6) >= 0.4) && (dend_Info{1,1}(movieInf(t),6) < 0.6)
                dend_data.bin_num(movieInf(t),1) = 4;
            elseif (dend_Info{1,1}(movieInf(t),6) >= 0.6) && (dend_Info{1,1}(movieInf(t),6) < 0.8)
                dend_data.bin_num(movieInf(t),1) = 5;
            elseif (dend_Info{1,1}(movieInf(t),6) >= 0.8) && (dend_Info{1,1}(movieInf(t),6) <= 1)
                dend_data.bin_num(movieInf(t),1) = 6;
            end
        end
    end

        if avgShadedPlots == 1

            figure('color',[1 1 1],'position',locs.*[1 2 1 .5]);
            for b = 1:length(movieInf)

                if rem(length(movieInf),2) == 1
                    row_plots = ceil(length(movieInf)/2);
                elseif rem(length(movieInf),2) == 0
                    row_plots = floor(length(movieInf)/2) + 1;
                end
            
                sbplt(1,b) = subplot(2,row_plots,[b]);

                cp_plot = shadedErrorBar(dend_data.time_axis{1,b},dend_data.cp_mean{1,b},dend_data.cp_err{1,b},'lineprops', {'-k','Linewidth',1});
                                
                hold on
                
                cf_plot = shadedErrorBar(dend_data.time_axis{1,b},dend_data.cf_mean{1,b},dend_data.cf_err{1,b},'lineprops', {'-r','Linewidth',1});
                                
                title(['Position ' num2str(b)])
                legend([cp_plot.mainLine, cf_plot.mainLine],'CP', 'CF');
                ylabel('DeltaF/F0')
                xlabel('Time (seconds)')
                
                hold off
                
                sbplt(1,b) = subplot(2,row_plots,[b]);

                cp_plot = shadedErrorBar(dend_data.time_axis{1,b},dend_data.cp_mean{1,b},dend_data.cp_err{1,b},'lineprops', {'-k','Linewidth',1});
                                
                hold on
                
                cf_plot = shadedErrorBar(dend_data.time_axis{1,b},dend_data.cf_mean{1,b},dend_data.cf_err{1,b},'lineprops', {'-r','Linewidth',1});
                                
                title(['Position ' num2str(b)])
                legend([cp_plot.mainLine, cf_plot.mainLine],'CP', 'CF');
                ylabel('DeltaF/F0')
                xlabel('Time (seconds)')
                
                hold off

            end

            %fullfig(gcf) 

            fig = gcf;
            
            linkaxes(sbplt,'xy')
            saveas(fig,strcat(folder,'\','CompiledWTSoma',num2str(dend_Info{1,1}(movieInf(t),7))),'fig') 
            saveas(fig,strcat(folder,'\','CompiledWTSoma',num2str(dend_Info{1,1}(movieInf(t),7))),'eps') 

            saveas(fig,strcat(folder,'\','CompiledWTSoma',num2str(dend_Info{1,1}(movieInf(t),7))),'png') 

            close all

            
        end
            
            
end
dend_data.info = dend_Info{1,1}(:,:);
dend_data.paths = [dend_Info{1,2}(:,:)];
fid = fopen([folder '\WT_CalciumMappingResults_Soma20200525.txt'], 'w');
for i = 1:length(dend_data.QI)
    fprintf(fid,'\n %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f' ,dend_data.info(i,:), dend_data.QI(i), dend_data.dsi(i), dend_data.cp_deltaF(i), dend_data.cf_deltaF(i),dend_data.dsiArea(i), dend_data.cp_Area(i), dend_data.cf_Area(i), dend_data.cf_riseslope(i), dend_data.bin_num(i), dend_data.cf_risestart(i), dend_data.cp_risestart(i), dend_data.shift(i), dend_data.position(i,1), dend_data.position(i,2), dend_data.orient(i), dend_data.center(i,1), dend_data.center(i,2));
end
fclose(fid);

%dend_data.results = [dend_data.info, dend_data.QI, dend_data.dsi, dend_data.cp_deltaF, dend_data.cf_deltaF, dend_data.riseslope];

recipients = {'787-244-4805'};
subject    = 'Function has finished successfully';
message    = 'Now go out and have fun!';
carrier    = 'att';
send_msg(recipients, subject, message, carrier);