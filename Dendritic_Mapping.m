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

directory = 'J:\';
%%%%%%%%%%%%%%%%%%%% dend %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen([folder '\WT_VoltageImaging_20200629.txt'], 'rt');
dend_Info = textscan(fid, '%f %f %f %f %f %f %f %s %f %f %f ', 'HeaderLines', 0, 'CollectOutput', true);
fclose(fid);

dend_area = unique(dend_Info{1,1}(:,7));


avgShadedPlots = input('Average Plots?\n');

PPT = input('Generate ppt?\n');

dend_events = input('Event Detection?\n');




for i = 1:length(dend_area)
    close all
    clearvars -except paras folder moviefiles i dend_area num_stimuli compare analysis avgShadedPlots  dend_events PPT date mouse age slidesFile slides directory var_pairedInfo var_Info dend_pairedInfo dend_Info dend_data
    movieInf = find(dend_Info{1,1}(:,7) == dend_area(i));
    paths = dend_Info{1,2}((dend_Info{1,1}(:,7) == dend_area(i)),:);
    for t = 1:length(movieInf)
        clearvars -except paras folder moviefiles i dend_area num_stimuli compare analysis avgShadedPlots dend_events PPT t movieInf cellSignals cp_traces cf_traces dend_PSTH edges window_size time_window SortedPeakData time date mouse age slidesFile slides directory var_pairedInfo var_Info dend_pairedInfo dend_Info paths dend_data mean_baseline
        movieNum = dend_Info{1,1}(movieInf(t),1);
        moviestring = num2str(movieNum);
        pathstr = paths(t,:);
        pname = dir([directory char(pathstr) '\TSeries*' moviestring]);
        tseriesxml = [pname.name '.xml'];
        pname = [pname.folder '\' pname.name];
        cellSignalString = strcat('CellSignals',moviestring,'.txt');

        %% IMPORT TIFF FILES%%%%
       
        roi_select = dend_Info{1,1}(movieInf(t),2);
        cp_dir = dend_Info{1,1}(movieInf(t),3);
        cf_dir = dend_Info{1,1}(movieInf(t),4);
        exp_type = dend_Info{1,1}(movieInf(t),5);
        
        soma_center = [dend_Info{1,3}(movieInf(t),2) dend_Info{1,3}(movieInf(t),3)];
        dend_orient = [dend_Info{1,3}(movieInf(t),1)];
        
        eye = 1;
        piece = 1;
        barType = 2;
        interTrialWait = 0.5;
        %%IMPORT CELL SIGNALS%%%%%

        cellSignals{1,t} = dlmread([directory char(pathstr) '\CellSignals\' cellSignalString]);
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

            y_fit2 = DataSmooth(y_fit2,(smooth_window), 'average');

            for j = 1: length(data.fFrame)
                y_fit2(data.fFrame(j):(data.fFrame(j)+wsize_stim),:)= NaN;
            end


            if exp_type == 1
                y_fit2(1:(round(framerate*4)),:)= NaN;
            end

            x_fit = 1:length(y_fit2);
            x_fit = x_fit(find(isnan(y_fit2(:,1)) == 0));
           
            y_fit3 = y_fit2(find(isnan(y_fit2(:,1)) == 0),:);
            
            dend_data.mean_baseline(movieInf(t)) = mean(y_fit3((round(framerate*5)):end, roi_select));

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
            time_1000 = [0:(1/new_FR):time_tx(end)];
            for u = 1: size(data.signal,2)
                data.signalSmooth(:,u) = interp1(time_tx, data.signal(:,u),time_1000, 'pchip');
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
            
            dend_data.new_FR(movieInf(t),1) = new_FR;


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
            for q=1:length(dlist)
                sortd(:,q)= find(dlist(q)==data.stim);
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
            clearvars r c tline tok 
            
            theta_PosVector = atand(((soma_center(1,2) - y_pos)/(soma_center(1,1) - x_pos)));
            vector_distance = sqrt((soma_center(1,1) - x_pos)^2 + (soma_center(1,2) - y_pos)^2);
            if theta_PosVector < 0
                theta_PosVector = 360 + theta_PosVector;
            end
            delta_theta = theta_PosVector - dend_orient;
            dend_data.spatial_shift(movieInf(t),1) = vector_distance*cosd(delta_theta);
            if (x_pos == soma_center(1,1)) && (y_pos == soma_center(1,2))
                dend_data.time_shift(movieInf(t),1) = 0;
            else
                dend_data.time_shift(movieInf(t),1) = (dend_data.spatial_shift(movieInf(t),1)/(stimInf.speed*1.07));
            end
            
                              
            dend_data.QI(movieInf(t),1) = (var(dend_data.cf_mean{1,t}))/(mean(var(dend_data.cf_traces{1,t},0,1)));
            
            dend_data.time_axis{1,t} = linspace(0, dend_data.time_window{1,t}, dend_data.window_size{1,t})';
            dend_data.search_indexes{1,t} = intersect(find(dend_data.time_axis{1,t} >= 0.35), find(dend_data.time_axis{1,t} <= 2));
            
            [dend_data.cp_pk(movieInf(t),1), dend_data.cp_loc(movieInf(t),1)] = findpeaks(dend_data.cp_mean{1,t}(dend_data.search_indexes{1,t},1),'SortStr', 'descend', 'NPeaks',1);
            [dend_data.cf_pk(movieInf(t),1), dend_data.cf_loc(movieInf(t),1)] = findpeaks(dend_data.cf_mean{1,t}(dend_data.search_indexes{1,t},1),'SortStr', 'descend', 'NPeaks',1);
            dend_data.cf_pkindex(movieInf(t),1) = dend_data.search_indexes{1,t}(1) + dend_data.cf_loc(movieInf(t),1) - 1;
            dend_data.cp_pkindex(movieInf(t),1) = dend_data.search_indexes{1,t}(1) + dend_data.cp_loc(movieInf(t),1) - 1;
            
            try
                [dend_data.cp_min(movieInf(t),1), dend_data.cp_minloc(movieInf(t),1)] = findpeaks((-1)*dend_data.cp_mean{1,t}((dend_data.search_indexes{1,t}(1):dend_data.cp_pkindex(movieInf(t),1)),1),'SortStr', 'descend', 'NPeaks',1);
            catch
                [dend_data.cp_min(movieInf(t),1), dend_data.cp_minloc(movieInf(t),1)] = findpeaks((-1)*dend_data.cp_mean{1,t}((1:dend_data.cp_pkindex(movieInf(t),1)),1),'SortStr', 'descend', 'NPeaks',1);
            end
            
            try
                [dend_data.cf_min(movieInf(t),1), dend_data.cf_minloc(movieInf(t),1)] = findpeaks((-1)*dend_data.cf_mean{1,t}((dend_data.search_indexes{1,t}(1):dend_data.cf_pkindex(movieInf(t),1)),1),'SortStr', 'descend', 'NPeaks',1);
            catch
                [dend_data.cf_min(movieInf(t),1), dend_data.cf_minloc(movieInf(t),1)] = findpeaks((-1)*dend_data.cf_mean{1,t}((1:dend_data.cf_pkindex(movieInf(t),1)),1),'SortStr', 'descend', 'NPeaks',1);
            end
            
            
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
            
            
            
            dend_data.dtpt_shift(movieInf(t),1) = round(dend_data.time_shift(movieInf(t),1)*new_FR);
                
            if dend_data.dtpt_shift(movieInf(t),1) > 0
                  cf_mean_val = dend_data.cf_mean{1,t}(1)*ones(abs(dend_data.dtpt_shift(movieInf(t))),1);
                  cf_err_val = dend_data.cf_err{1,t}(1)*ones(abs(dend_data.dtpt_shift(movieInf(t))),1);
                  cp_mean_val = dend_data.cp_mean{1,t}(end)*ones(abs(dend_data.dtpt_shift(movieInf(t))),1);
                  cp_err_val = dend_data.cp_err{1,t}(end)*ones(abs(dend_data.dtpt_shift(movieInf(t))),1);
                  dend_data.cf_mean_corr{1,t} = [cf_mean_val; dend_data.cf_mean{1,t}(1:(end-abs(dend_data.dtpt_shift(movieInf(t)))))];
                  dend_data.cf_err_corr{1,t} = [cf_err_val; dend_data.cf_err{1,t}(1:(end-abs(dend_data.dtpt_shift(movieInf(t)))))];
                  dend_data.cp_mean_corr{1,t} = [dend_data.cp_mean{1,t}((abs(dend_data.dtpt_shift(movieInf(t)))+1):end); cp_mean_val];
                  dend_data.cp_err_corr{1,t} = [dend_data.cp_err{1,t}((abs(dend_data.dtpt_shift(movieInf(t)))+1):end); cp_err_val];
            elseif dend_data.dtpt_shift(movieInf(t),1) < 0
                cf_mean_val = dend_data.cf_mean{1,t}(end)*ones(abs(dend_data.dtpt_shift(movieInf(t))),1);
                  cf_err_val = dend_data.cf_err{1,t}(end)*ones(abs(dend_data.dtpt_shift(movieInf(t))),1);
                  cp_mean_val = dend_data.cp_mean{1,t}(1)*ones(abs(dend_data.dtpt_shift(movieInf(t))),1);
                  cp_err_val = dend_data.cp_err{1,t}(1)*ones(abs(dend_data.dtpt_shift(movieInf(t))),1);
                  dend_data.cp_mean_corr{1,t} = [cp_mean_val; dend_data.cp_mean{1,t}(1:(end-abs(dend_data.dtpt_shift(movieInf(t)))))];
                  dend_data.cp_err_corr{1,t} = [cp_err_val; dend_data.cp_err{1,t}(1:(end-abs(dend_data.dtpt_shift(movieInf(t)))))];
                  dend_data.cf_mean_corr{1,t} = [dend_data.cf_mean{1,t}((abs(dend_data.dtpt_shift(movieInf(t)))+1):end); cf_mean_val];
                  dend_data.cf_err_corr{1,t} = [dend_data.cf_err{1,t}((abs(dend_data.dtpt_shift(movieInf(t)))+1):end); cf_err_val];
            elseif dend_data.dtpt_shift(movieInf(t),1) == 0
                dend_data.cp_mean_corr{1,t} = dend_data.cp_mean{1,t};
                dend_data.cp_err_corr{1,t} = dend_data.cp_err{1,t};
                dend_data.cf_mean_corr{1,t} = dend_data.cf_mean{1,t};
                dend_data.cf_err_corr{1,t} = dend_data.cf_err{1,t};
            end
            
            if (exp_type == 2) && (dend_data.QI(movieInf(t)) >=0.1350)
                out_fn = strcat(folder, '\CellTraces\Bin', num2str(dend_data.bin_num(movieInf(t),1)),'\Speed',num2str(stimInf.speed),'\Cell',num2str(dend_area(i)),'Movie',moviestring,'.txt');
                dlmwrite(out_fn, [dend_data.cp_mean_corr{1,t},dend_data.cf_mean_corr{1,t}], 'delimiter', '\t','newline','pc')
            end
                
            dend_data.meanTrace_length(movieInf(t),1) = length(dend_data.cf_mean_corr{1,t});
            
            dend_data.cp_min(movieInf(t),1) = min(dend_data.cp_mean_corr{1,t}(find(dend_data.time_axis{1,t} <= 0.35),1));
            dend_data.cf_min(movieInf(t),1) = min(dend_data.cf_mean_corr{1,t}(find(dend_data.time_axis{1,t} <= 0.35),1));
            
            dend_data.cp_mean_corr{1,t} = dend_data.cp_mean_corr{1,t} - dend_data.cp_min(movieInf(t),1);
            dend_data.cp_mean{1,t} = dend_data.cp_mean{1,t} - dend_data.cp_min(movieInf(t),1);
            
            dend_data.cf_mean_corr{1,t} = dend_data.cf_mean_corr{1,t} - dend_data.cf_min(movieInf(t),1);
            dend_data.cf_mean{1,t} = dend_data.cf_mean{1,t} - dend_data.cf_min(movieInf(t),1);
            
            dend_data.cp_baseline(movieInf(t),1) = mean(dend_data.cp_mean_corr{1,t}(find(dend_data.time_axis{1,t} <= 0.35),1));
            dend_data.cp_baselineErr(movieInf(t),1) = std(dend_data.cp_mean_corr{1,t}(find(dend_data.time_axis{1,t} <= 0.35),1));

            dend_data.cf_baseline(movieInf(t),1) = mean(dend_data.cf_mean_corr{1,t}(find(dend_data.time_axis{1,t} <= 0.35),1));
            dend_data.cf_baselineErr(movieInf(t),1) = std(dend_data.cf_mean_corr{1,t}(find(dend_data.time_axis{1,t} <= 0.35),1));

            dend_data.cp_ThrA(movieInf(t),1) = dend_data.cp_baseline(movieInf(t),1) + 3*dend_data.cp_baselineErr(movieInf(t),1);
            dend_data.cp_ThrB(movieInf(t),1) = dend_data.cp_baseline(movieInf(t),1) + 1.5*dend_data.cp_baselineErr(movieInf(t),1);

            dend_data.cf_ThrA(movieInf(t),1) = dend_data.cf_baseline(movieInf(t),1) + 3*dend_data.cf_baselineErr(movieInf(t),1);
            dend_data.cf_ThrB(movieInf(t),1) = dend_data.cf_baseline(movieInf(t),1) + 1.5*dend_data.cf_baselineErr(movieInf(t),1);

            dend_data.lowpass1_cp{1,t} = lowpass(dend_data.cp_mean_corr{1,t}, 2, round(new_FR));
            dend_data.lowpass1_cf{1,t} = lowpass(dend_data.cf_mean_corr{1,t}, 2, round(new_FR));

            [cf_pks, cf_locs] = findpeaks(dend_data.lowpass1_cf{1,t}(intersect(find(dend_data.time_axis{1,t} >= 0.35), find(dend_data.time_axis{1,t} <= 2))),'MinPeakHeight', dend_data.cf_ThrA(movieInf(t),1),'SortStr', 'descend', 'NPeaks',1);
            if isempty(cf_pks) == 1
                dend_data.cf_loc(movieInf(t),1) = 0;
                dend_data.cf_pk(movieInf(t),1) = 0;
                dend_data.cf_risestart(movieInf(t),1) = 0;
                dend_data.cf_riseslope(movieInf(t),1) = 0;
                %dend_data.cf_peaktime(movieInf(t),1) = 0;
            elseif isempty(cf_pks) == 0
                indices = find(dend_data.time_axis{1,t} <= 0.35);
                dend_data.cf_peaktime(movieInf(t),1) = (indices(end) + cf_locs)/new_FR;
                dend_data.cf_risestart(movieInf(t),1) = dsearchn(dend_data.lowpass1_cf{1,t}(intersect(find(dend_data.time_axis{1,t} >= 0.35), find(dend_data.time_axis{1,t} <= dend_data.cf_peaktime(movieInf(t),1)))),dend_data.cf_ThrB(movieInf(t),1));
                dend_data.cf_risestart(movieInf(t),1) = (indices(end) + dend_data.cf_risestart(movieInf(t),1))/new_FR;
                dend_data.cf_riseslope(movieInf(t),1) = (cf_pks - dend_data.cf_ThrB(movieInf(t),1))/(dend_data.cf_peaktime(movieInf(t),1) - dend_data.cf_risestart(movieInf(t),1));
                indices = [];
            end
            cf_pks = [];
            cf_locs = [];

            [cp_pks, cp_locs] = findpeaks(dend_data.lowpass1_cp{1,t}(intersect(find(dend_data.time_axis{1,t} >= 0.35), find(dend_data.time_axis{1,t} <= 2))),'MinPeakHeight', dend_data.cp_ThrA(movieInf(t),1),'SortStr', 'descend', 'NPeaks',1);
            if isempty(cp_pks) == 1
                dend_data.cp_loc(movieInf(t),1) = 0;
                dend_data.cp_pk(movieInf(t),1) = 0;
                dend_data.cp_risestart(movieInf(t),1) = 0;
                dend_data.cp_riseslope(movieInf(t),1) = 0;
                %dend_data.cp_peaktime(movieInf(t),1) = 0;
            elseif isempty(cp_pks) == 0
                indices = find(dend_data.time_axis{1,t} <= 0.35);
                dend_data.cp_peaktime(movieInf(t),1) = (indices(end) + cp_locs)/new_FR;
                dend_data.cp_risestart(movieInf(t),1) = dsearchn(dend_data.lowpass1_cp{1,t}(intersect(find(dend_data.time_axis{1,t} >= 0.35), find(dend_data.time_axis{1,t} <= dend_data.cp_peaktime(movieInf(t),1)))),dend_data.cp_ThrB(movieInf(t),1));
                dend_data.cp_risestart(movieInf(t),1) = (indices(end) + dend_data.cp_risestart(movieInf(t),1))/new_FR;
                
                indices = [];
            end
            cp_pks = [];
            cp_locs = [];
            
            
            
            
        end
    end

        if avgShadedPlots == 1
            locs=[1 200 1400 800];
            figure('color',[1 1 1],'position',locs.*[1 2 1 .5]);
            for d = 1:length(movieInf)

                if rem(length(movieInf),2) == 1
                    row_plots = ceil(length(movieInf)/2);
                elseif rem(length(movieInf),2) == 0
                    row_plots = floor(length(movieInf)/2) + 1;
                end
            
                sbplt(1,d) = subplot(2,row_plots,[d]);

                cp_plot = shadedErrorBar(dend_data.time_axis{1,d},dend_data.cp_mean{1,d},dend_data.cp_err{1,d},'lineprops', {'-k','Linewidth',1});
                                
                hold on
                
                cf_plot = shadedErrorBar(dend_data.time_axis{1,d},dend_data.cf_mean{1,d},dend_data.cf_err{1,d},'lineprops', {'-r','Linewidth',1});
                                
                title(['Position ' num2str(d)])
                legend([cp_plot.mainLine, cf_plot.mainLine],'CP', 'CF');
                ylabel('DeltaF/F0')
                xlabel('Time (seconds)')
                
                hold off

            end

            %fullfig(gcf) 

            fig = gcf;
            
            linkaxes(sbplt,'xy')
            saveas(fig,strcat(folder,'\','CompiledWTVoltage',num2str(dend_Info{1,1}(movieInf(t),7))),'fig') 
            saveas(fig,strcat(folder,'\','CompiledWTVoltage',num2str(dend_Info{1,1}(movieInf(t),7))),'eps') 

            saveas(fig,strcat(folder,'\','CompiledWTVoltage',num2str(dend_Info{1,1}(movieInf(t),7))),'png') 

            close all
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            locs=[1 200 1400 800];
            figure('color',[1 1 1],'position',locs.*[1 2 1 .5]);
            for d = 1:length(movieInf)

                if rem(length(movieInf),2) == 1
                    row_plots = ceil(length(movieInf)/2);
                elseif rem(length(movieInf),2) == 0
                    row_plots = floor(length(movieInf)/2) + 1;
                end

                sbplt(1,d) = subplot(2,row_plots,[d]);

                cp_plot = shadedErrorBar(dend_data.time_axis{1,d},dend_data.cp_mean_corr{1,d},dend_data.cp_err_corr{1,d},'lineprops', {'-k','Linewidth',1});
                                
                hold on
                
                cf_plot = shadedErrorBar(dend_data.time_axis{1,d},dend_data.cf_mean_corr{1,d},dend_data.cf_err_corr{1,d},'lineprops', {'-r','Linewidth',1});
                                
                title(['Position ' num2str(d)])
                legend([cp_plot.mainLine, cf_plot.mainLine],'CP', 'CF');
                ylabel('DeltaF/F0')
                xlabel('Time (seconds)')
                
                hold off

            end

            %fullfig(gcf) 

            fig = gcf;
            
            linkaxes(sbplt,'xy')
            saveas(fig,strcat(folder,'\','CompiledWTVoltageCORR',num2str(dend_Info{1,1}(movieInf(t),7))),'fig') 
            saveas(fig,strcat(folder,'\','CompiledWTVoltageCORR',num2str(dend_Info{1,1}(movieInf(t),7))),'eps') 
            saveas(fig,strcat(folder,'\','CompiledWTVoltageCORR',num2str(dend_Info{1,1}(movieInf(t),7))),'png') 

            close all

            
        end
            
            
end
dend_data.info = dend_Info{1,1}(:,:);
dend_data.paths = [dend_Info{1,2}(:,:)];
fid = fopen([folder '\WT_VoltageMappingResults_20200629_corr.txt'], 'w');
for i = 1:length(dend_data.QI)
    fprintf(fid,'\n %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f' ,dend_data.info(i,:), dend_data.new_FR(i), dend_data.meanTrace_length(i), dend_data.QI(i), dend_data.mean_baseline(i), dend_data.dsi(i), dend_data.cp_deltaF(i), dend_data.cf_deltaF(i),dend_data.dsiArea(i), dend_data.cp_Area(i), dend_data.cf_Area(i), dend_data.cf_riseslope(i), dend_data.bin_num(i), dend_data.cf_risestart(i), dend_data.cp_risestart(i), dend_data.cf_peaktime(i), dend_data.cp_peaktime(i), dend_data.position(i,1), dend_data.position(i,2));
end
fclose(fid);
% 
% %dend_data.results = [dend_data.info, dend_data.QI, dend_data.dsi, dend_data.cp_deltaF, dend_data.cf_deltaF, dend_data.riseslope];
% 
% recipients = {'787-244-4805'};
% subject    = 'Function has finished successfully';
% message    = 'Now go out and have fun!';
% carrier    = 'att';
% send_msg(recipients, subject, message, carrier);