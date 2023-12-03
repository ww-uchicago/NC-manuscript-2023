% This serves as the main fu
close all; clear all;clearvars;  clc;
%% setup parameters
paras.vfreq = 1000;   % TTL sampling rate, 1000, or 10000
paras.rmrep=[2]; % remove certain reps of recording
%% main part
folder = 'D:\Wei Lab HD\recordings\New Rig\20180309';
%pname=uigetdir(folder);
%read the xml, and tif file and calculate the mean intensity value over
%time
retinaInfo = dlmread([folder '\RetinaInfo.txt']);
moviefiles = dir([folder '\*.zip']);
cell_area = unique(retinaInfo(:,6));

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
 
            save_plots = 1;
            locs=[1 200 1400 800];
            % load the raw data
            numbg=1;   % number of background ROIs

            data.trace= cellSignals{1,t};
            data.signal=data.trace(:,1:end-1);
            data.bg=data.trace(:,end);
            %%
            rois=1:size(data.signal,2);
            %rois=[6];
            
            %process the voltage data
            



            %%% CREATE A SMOOTH AVERAGED BACKGROUNG TRACE & BACKGROUND SUBSTRACTION
            data.signal = DataSmooth(data.signal, 3, 'average');


    %         for i = 1:size(data.signal,2);
    %             for k = 1:length(data.stim)
    %                 data.signalNew{1,i}(:,k) = data.signal{1,i}(:,k) - meanBgd;
    %             end
    %         end
    %         
            data.signalNew = data.signal;
            
            rho_corr{1,t} = zeros(size(data.signal,2));
            
            for b = 1:size(data.signal,2)
                for a = 1:size(data.signal,2)
                    [crosscorr] = xcorr(data.signal(:,b), data.signal(:,a), 0,'coeff');
                    rho_corr{1,t}(b,a) = crosscorr;
                end
            end
            
            colormap('jet')
            imagesc(rho_corr{1,t}, [0.2 1])
            colorbar
            
            save_plot_fn=strcat('Movie-',num2str(movieNums(t)));
            save_plot_path = [folder '\' save_plot_fn];
            saveas(gcf,save_plot_path,'fig') 
            saveas(gcf,save_plot_path,'eps') 

            %%%PLOT THE SORTED RAW DATA FOR EACH CELL
            [rois] = ReadImageJROI([folder '\RoiSet' num2str(movieNums(t)) '.zip']);

            coords = zeros(length(rois),4);

            for q = 1:length(rois)
                coords(q,:) = rois{1,q}.vnRectBounds;
            end

            coords(coords < 1) = 1;
            for q = 1:length(rois)
                pos(q,1:2) = [mean(coords(q,1),coords(q,3)), mean(coords(q,2),coords(q,4))];
            end
            
            for b = 1:size(data.signal,2)
                for a = 1:size(data.signal,2)
                    distance_matrix(b,a) = pdist([pos(b,:);pos(a,:)],'euclidean');
                end
            end
            close all
            figure
            scatter(distance_matrix(:), rho_corr{1,t}(:))
            
            save_plot_fn=strcat('Movie-',num2str(movieNums(t)),'DC');
            save_plot_path = [folder '\' save_plot_fn];
            saveas(gcf,save_plot_path,'fig') 
            saveas(gcf,save_plot_path,'eps') 
           

        end
    end

    
            
            
end
    

 