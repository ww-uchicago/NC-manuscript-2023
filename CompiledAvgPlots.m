% This serves as the main fu
close all; clear all;clearvars;  clc;
%% setup parameters
paras.vfreq = 1000;   % TTL sampling rate, 1000, or 10000
paras.rmrep=[2]; % remove certain reps of recording
%% main part
folder = 'J:\WeiLabHD\recordings\NewRig\20170508';
%pname=uigetdir(folder);
%read the xml, and tif file and calculate the mean intensity value over
%time
retinaInfo = dlmread([folder '\RetinaInfo.txt']);
moviefiles = dir([folder '\*.zip']);
cell_area = unique(retinaInfo(:,6));

avgShadedPlots = input('Average Plots?\n');

PPT = input('Generate ppt?\n');

if PPT
    date = input('YR-MNTH-DAY\n', 's');
    mouse = input('Mouse strain?\n', 's');
    age = input('Age?\n', 's');

    import mlreportgen.ppt.*;

    slidesFile = [folder '\' date '.pptx'];
    slides = Presentation(slidesFile,'myTemplate');
end


for i = 1:length(cell_area)
    close all
    clearvars -except paras folder moviefiles i retinaInfo cell_area num_stimuli compare analysis avgShadedPlots  PPT date mouse age slidesFile slides
    movieNums = retinaInfo((retinaInfo(:,6) == cell_area(i)),1);
    for t = 1:length(movieNums)
        clearvars -except paras folder moviefiles i retinaInfo cell_area num_stimuli compare analysis avgShadedPlots PPT t movieNums cellSignals ptrace SortedPeakData time date mouse age slidesFile slides
        moviestring = num2str(movieNums(t));
        pname = dir([folder '\TSeries*' moviestring]);
        tseriesxml = [pname.name '.xml'];
        pname = [folder '\' pname.name];
        cellSignalString = strcat('CellSignals',num2str(movieNums(t)),'.txt');

        %% IMPORT TIFF FILES%%%%

        movieInf = find(retinaInfo(:,1) == movieNums(t));
        eye = retinaInfo(movieInf, 2);
        piece = retinaInfo(movieInf, 3);
        barType = retinaInfo(movieInf, 4);
        interTrialWait = retinaInfo(movieInf, 5);

        %%IMPORT CELL SIGNALS%%%%%

        cellSignals{1,t} = dlmread([folder '\CellSignals\' cellSignalString]);
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
            data.vol=csvread(fp.vol,1); data.vol=data.vol(:,2); 
            data.trace= cellSignals{1,t};
            data.signal=data.trace(:,1:end-1);
            data.bg=data.trace(:,end);
            %%
            rois=1:size(data.signal,2);
            %rois=[6];
            rm_rep= [ ];
            %process the voltage data
            [data.vlocs,data.pks]=peakdetect(data.vol,2,7000);



            if length(data.vlocs) < length(stimInf.dirs)
                loc_diff = diff(data.vlocs);
                location = find(loc_diff > 100000);
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

            data.fFrame=ceil(data.vlocs/fRate); %data.fFrame(1)=[];
            wsize=min(diff(data.fFrame));

            reps=reshape(1:length(stimInf.dirs),length(dlist),(length(stimInf.dirs))/length(dlist)); 
            reps=reshape(reps,size(reps,1)*size(reps,2),1); nreps=length(reps)/length(dlist);

            % CLIP DATA SET FOR EACH SWEEP%%%%
            data.signal=clip(data.signal,data.fFrame,wsize);
            data.bg=clip(data.bg,data.fFrame,wsize);

            % remove certain reps
            data.signal=cellfun(@(x) x(:,reps),data.signal,'uniformoutput',false);
            data.bg=cellfun(@(x) x(:,reps),data.bg,'uniformoutput',false);
            data.fFrame=data.fFrame(reps);
            data.stim=data.stim(reps);

            % ORGANIZE BY DIRECTION

            

            %%%%SMOOTHING OF RGC CALCIUM SIGNALS



            %%% CREATE A SMOOTH AVERAGED BACKGROUNG TRACE & BACKGROUND SUBSTRACTION
            meanBgd = mean(data.bg{1,1},2);
            meanBgd = DataSmooth(meanBgd, 6, 'average');


    %         for i = 1:size(data.signal,2);
    %             for k = 1:length(data.stim)
    %                 data.signalNew{1,i}(:,k) = data.signal{1,i}(:,k) - meanBgd;
    %             end
    %         end
    %         
            data.signalNew = data.signal;

            data.signalNew = cellfun(@(x) DataSmooth(x,6,'average'),data.signalNew,'uniformoutput',false);
            data.signalNew = cellfun(@(x) DataSmooth(x,4,'average'),data.signalNew,'uniformoutput',false);
            data.signalNew = cellfun(@(x) DataSmooth(x,3,'average'),data.signalNew,'uniformoutput',false);
            %data.signalNew = cellfun(@(x) DataSmooth(x,3,'average'),data.signalNew,'uniformoutput',false);

            %%%PLOT THE SORTED RAW DATA FOR EACH CELL
            ptrace{1,t} =cellfun(@(x) reconstruct(x, dlist),data.signalNew,'uniformoutput',false);
            for i=1:nreps
                rlegends{i}=['Trial',' ',num2str(i)];
            end 

            % frames = [1;34;107;141;214;248;321;355;428;462;535;569;642;676;749;783];
            % 
            % ptrace{5}(:,6) = 0;
            % for w = 1:length(frames)
            %     ptrace{5}(frames(w),6) = 800;
            % end

            %%%%%%%%%%%%%%%%%%% CALCULATE QI (QUALITY INDEX) %%%%%%%%%%%%%%%%%%%%%%%

            QI = zeros(size(ptrace{1,t},2),1);
            for i = 1:size(ptrace{1,t},2)
                var1 = var(mean(ptrace{1,t}{1,i},2));
                var2 = mean(var(ptrace{1,t}{1,i}));
                QI(i) = var1/var2;
            end



            %calculate the average traces across sweeps
            if nreps>1
                data.signalAve=cellfun(@(x) sweepMean(x,nreps),data.signalNew,'uniformoutput',false);
                data.bgAve=cellfun(@(x) sweepMean(x,nreps),data.bg,'uniformoutput',false);
                data.signalErr=cellfun(@(x) sweepErr(x,nreps),data.signalNew,'uniformoutput',false);
                data.bgErr=cellfun(@(x) sweepErr(x,nreps),data.bg,'uniformoutput',false);
            else
                data.signalAve=data.signalNew;
                data.bgAve=data.bg;

            end


            [nframes, num_reps] = size(data.signalNew{1,1});
            if stimInf.speed == 4
                time = stimInf.speed + interTrialWait;
            end
            if stimInf.speed > 50
                time(t) = (2*stimInf.rad +stimInf.bw)/stimInf.speed + interTrialWait;
            end
            framerate = nframes/time(t);
            for i = 1: size(ptrace{1,t},2)
                data.f0(i) = abs(mean(mean(data.signalNew{1,i}((nframes-10):nframes,:),2)));
            end

            twoP_response = data.trace(1:(ceil(10*framerate)),:);
            twoP_response = DataSmooth(twoP_response, 10, 'average');

            % CALCULATE DELTA_F/F0
            for i=1:size(data.signalNew,2)
                [data.f{i},~]=max(data.signalNew{i});
                [r, c] = size(data.signalNew{i});
                time_axis = [0:(1/framerate):((r-1)/framerate)];
                %[data.signalSub{i}] = data.signal{i} - data.bg{1,1};
                for k = 1:length(stimInf.dirs)
                    areaTrace = data.signalNew{i}(:,k) - data.f0(i);
                    [data.area{i}(:,k)] = trapz(time_axis, areaTrace);
                    if data.area{i}(:,k) < 1
                        data.area{i}(:,k) = 1;
                    end
                end
                data.df{i}=(data.f{i}-data.f0(i))/data.f0(i); 
                data.df{i}=reshape(data.df{i},nreps,length(dlist));
                data.area{i}=reshape(data.area{i},nreps,length(dlist));
                if nreps>1
                    data.mamp{i}=mean(data.df{i});
                else
                    data.mamp{i}=data.df{i};
                end
            end

            close all

            polarityIndex = zeros(size(ptrace{1,t},2),1);
            polarityIndex2 = zeros(size(ptrace{1,t},2),1);
            dir = [0 45 90 135 180 225 270 315];

            OOi = zeros(size(ptrace{1,t},2),1);
            POi = zeros(size(ptrace{1,t},2),1);
            [nframes, num_reps] = size(data.signalNew{1,1});
            pFrate = nframes/time(t);
            trial_startFrames = 1;
            if (barType == 1 || barType == 2)
                on_startFrame = trial_startFrames;
                on_endFrame = ceil(pFrate*(2*stimInf.rad/stimInf.speed));
                on_endtime = 2*stimInf.rad/stimInf.speed;
                off_startFrame = ceil(pFrate*(stimInf.bw/stimInf.speed));
                off_endFrame = ceil(pFrate*((2*stimInf.rad +stimInf.bw)/stimInf.speed));
                off_endtime = (2*stimInf.rad +stimInf.bw)/stimInf.speed;
            elseif barType == 4
                off_startFrame = trial_startFrames;
                off_endFrame = ceil(pFrate*1.3);
                on_startFrame = floor(pFrate*1.49);
                on_endFrame = ceil(pFrate*(1.5+1.3));
            end

            for w = 1:size(ptrace{1,t},2)
                [nframes, num_reps] = size(data.signalNew{1,w});
                if stimInf.speed == 4
                    time = stimInf.speed + interTrialWait;
                end
                if stimInf.speed > 50
                    time(t) = (2*stimInf.rad +stimInf.bw)/stimInf.speed + interTrialWait;
                end

                pFrate = nframes/time(t);
                trial_startFrames = 1;


                if (barType == 1 || barType == 2)
                    on_resp{1,w} = data.signalNew{1,w}(on_startFrame:on_endFrame,:);
                    off_resp{1,w} = data.signalNew{1,w}(off_startFrame:off_endFrame,:);
                    on_peaks(w,:) = (max(on_resp{1,w}) - on_resp{1,w}(1,:))/data.f0(w);
                    off_peaks(w,:) = (max(off_resp{1,w}) - off_resp{1,w}(1,:))/data.f0(w);

                    on_rThr = on_peaks(w,:) > 0.4;
                    off_rThr = off_peaks(w,:) > 0.4;

                    sum_on = sum(sum(on_resp{1,w}));
                    sum_off = sum(sum(off_resp{1,w}));
                    if barType ~= 1
                        POi(w, 1) = (sum_on - sum_off)/(sum_on + sum_off);
                    end

                    OOi(w,1) = (mean(on_peaks(w,find(on_rThr))-mean(off_peaks(w, find(off_rThr)))))/(mean(on_peaks(w,find(on_rThr))) + mean(off_peaks(w, find(off_rThr)))); 
                end
                for a = 2:nframes
                    dFdf{1,w}((a-1),:) = (data.signalNew{1,w}(a,:) - data.signalNew{1,w}((a-1),:))/(1/pFrate);
                end

                if (barType == 1 || barType == 2)
                    on_dFdt{1,w} = dFdf{1,w}(1:(on_endFrame-1),:);
                    off_dFdt{1,w} = dFdf{1,w}(off_startFrame:(off_endFrame-1),:);

                    on_onset = sum(on_dFdt{1,w}(on_dFdt{1,w} > 0));
                    off_onset = sum(off_dFdt{1,w}(off_dFdt{1,w} > 0));

                    POi(w,1) = (on_onset - off_onset)/(on_onset + off_onset);
                end

                sbins=4;
                dFdf=cellfun(@(x) DataSmooth(x,sbins,'average'),dFdf,'uniformoutput',false);

                Ave_dFdf{1,w} = mean(dFdf{1,w},2);


                dF.Baseline{1,w} = mean(Ave_dFdf{1,w}(floor(pFrate*3.2):end));
                dF.std{1,w} = std(Ave_dFdf{1,w}(floor(pFrate*3.2):end));


                if OOi(w,1) >= 0.5 && OOi(w,1) <=1.0 
                    polarityIndex(w,1) = 2;
                elseif OOi(w,1) <= -0.5 && OOi(w,1) >= -1.0
                    polarityIndex(w,1) = 3;
                elseif OOi(w,1) > -0.3 && OOi(w,1) < 0.3
                    polarityIndex(w,1) = 1;
                elseif OOi(w,1) == 30
                    polarityIndex(w,1) = 0;
                else
                    polarityIndex(w,1) = 4;
                end  
            end

            %%%%PLOT TUNING CURVES, CALCULATE DSI OSI VECTOR SUM

            [SortedPeakData{1,t}, SortedAreaData, outputMatrixPeaks, outputMatrixArea, max_valPeaks, max_valArea] = loadImagingDataScript2(movieNums(t), data.df, data.area, length(rois), folder, dlist);
        end
    end

        if avgShadedPlots == 1
            for b = 1:length(movieNums)
                    stim_numCells(b) = size(ptrace{1,b},2);
            end

            if PPT
                h = retinaInfo((retinaInfo(:,1) == movieNums(1)),6);

                if h == 1
                    % Add a title slide
                    presentationTitleSlide = add(slides,'Title Slide');
                    replace(presentationTitleSlide,'Title',mouse);

                    subtitleText = Paragraph(age);
                    replace(presentationTitleSlide,'Subtitle',subtitleText);
                end

                
                % Add a text slide
                textSlide = add(slides,'Title and Content');
                area_title = ['Area ' num2str(h) ': Movies ' num2str(movieNums(1)) ' - ' num2str(movieNums(length(movieNums)))];    
                titleText2 = Paragraph(area_title);
                contents = find(textSlide,'Title');
                replace(contents(1),titleText2);
            end

                        
            for c=1:min(stim_numCells)
                figure('color',[1 1 1],'position',locs.*[1 2 1 .5]);
                for b = 1:length(movieNums)
                    y = mean(ptrace{1,b}{1,c},2);
                    err = std(ptrace{1,b}{1,c},0,2)/sqrt(size(ptrace{1,b}{1,c},2));
                    
                    if b == 1
                        subplot(length(movieNums),4,[1 2 3])

                        shadedErrorBar(linspace(0, length(dlist)*time(b), size(ptrace{1,b}{1,1},1)),y,err, 'lineprops', {'-k','Linewidth',2})
                        ylabel('a.u.')
                        xlabel('Time (seconds)')

                        subplot(length(movieNums),4,[4])

                        fig=polar(degToRad(SortedPeakData{1,b}.directions),SortedPeakData{1,b}.meanResponses{1,c},'k');
                        set(fig,'linewidth',4);
                    elseif b == 2 
                        subplot(length(movieNums),4,[5 6 7])

                        shadedErrorBar(linspace(0, length(dlist)*time(b), size(ptrace{1,b}{1,1},1)),y,err, 'lineprops', {'-k','Linewidth',2})
                        ylabel('a.u.')
                        xlabel('Time (seconds)')

                        subplot(length(movieNums),4,[8])

                        fig=polar(degToRad(SortedPeakData{1,b}.directions),SortedPeakData{1,b}.meanResponses{1,c},'k');
                        set(fig,'linewidth',4);
                    elseif b == 3
                        subplot(length(movieNums),4,[9 10 11])

                        shadedErrorBar(linspace(0, length(dlist)*time(b), size(ptrace{1,b}{1,1},1)),y,err, 'lineprops', {'-k','Linewidth',2})
                        ylabel('a.u.')
                        xlabel('Time (seconds)')

                        subplot(length(movieNums),4,[12])

                        fig=polar(degToRad(SortedPeakData{1,b}.directions),SortedPeakData{1,b}.meanResponses{1,c},'k');
                        set(fig,'linewidth',4);
                    end
                    
                end
                
                fullfig(gcf) 
                
                fig = gcf;
                

                saveas(fig,strcat(folder,'\AvgPlots\','CompiledArea',num2str(retinaInfo((retinaInfo(:,1) == movieNums(b)),6)),'Cell',num2str(c)),'fig') 
                saveas(fig,strcat(folder,'\AvgPlots\','CompiledArea',num2str(retinaInfo((retinaInfo(:,1) == movieNums(b)),6)),'Cell',num2str(c)),'eps') 
                
                saveas(fig,strcat(folder,'\AvgPlots\','CompiledArea',num2str(retinaInfo((retinaInfo(:,1) == movieNums(b)),6)),'Cell',num2str(c)),'png') 
                
                
                if PPT
                    picture_name = strcat(folder,'\AvgPlots\','CompiledArea',num2str(retinaInfo((retinaInfo(:,1) == movieNums(b)),6)),'Cell',num2str(c),'.png');

                    plot1 = Picture(picture_name);
                    cell_title = ['Cell ' num2str(c) ': Movies ' num2str(movieNums(1)) ' - ' num2str(movieNums(length(movieNums)))];  
                    pictureSlide = add(slides,'Title and Content');
                    replace(pictureSlide,'Title',cell_title);
                    contents = find(pictureSlide,'Content');
                    replace(contents(1),plot1);
                end
                
                
                close all
            end
            
        end
            
            
end
    
% Generate and open the presentation
if PPT
    close(slides);

    if ispc
        winopen(slidesFile);
    end
end
    
 
