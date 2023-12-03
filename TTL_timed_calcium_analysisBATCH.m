
% Read the files, make sure only stimulus text, time series excel, and voltage .csv file present in the folder
folder1= 'C:\Users\admin\Desktop\Wei Lab HD\recordings\New Rig\';
fp_strings = ['20151105\TSeries-11052015-1204-';
    '20151105\TSeries-11052015-1204-';
    '20151105\TSeries-11052015-1204-';
    '20151105\TSeries-11052015-1204-';
    '20151105\TSeries-11052015-1204-';
    '20151105\TSeries-11052015-1204-';
    '20151105\TSeries-11052015-1204-';
    '20151106\TSeries-11062015-1332-';
    '20151106\TSeries-11062015-1332-';
    '20151106\TSeries-11062015-1332-';
    '20151109\TSeries-11092015-1229-';
    '20151109\TSeries-11092015-1229-';
    '20151109\TSeries-11092015-1229-';
    '20151109\TSeries-11092015-1229-';
    '20151109\TSeries-11092015-1229-';
    '20151111\TSeries-11112015-1145-';
    '20151111\TSeries-11112015-1145-';
    '20151111\TSeries-11112015-1145-';
    '20151111\TSeries-11112015-1145-';
    '20151111\TSeries-11112015-1145-';
    '20151111\TSeries-11112015-1145-';
    '20151111\TSeries-11112015-1145-';
    '20151111\TSeries-11112015-1145-';
    '20151201\TSeries-12012015-1113-';
    '20151201\TSeries-12012015-1113-';
    '20151201\TSeries-12012015-1113-';
    '20151201\TSeries-12012015-1113-';
    '20151201\TSeries-12012015-1113-';
    '20151201\TSeries-12012015-1113-';
    '20151201\TSeries-12012015-1113-';
    '20151204\TSeries-12042015-1502-';
    '20151204\TSeries-12042015-1502-';
    '20151204\TSeries-12042015-1502-';
    '20151204\TSeries-12042015-1502-';
    '20151204\TSeries-12042015-1502-';
    '20151207\TSeries-12072015-1212-';
    '20151207\TSeries-12072015-1212-';
    '20151207\TSeries-12072015-1212-';
    '20151207\TSeries-12072015-1212-';
    '20151207\TSeries-12072015-1212-';
    '20151207\TSeries-12072015-1212-';
    '20151207\TSeries-12072015-1212-';
    '20151207\TSeries-12072015-1212-';
    '20151209\TSeries-12092015-1259-';
    '20151209\TSeries-12092015-1259-';
    '20151209\TSeries-12092015-1450-';
    '20151216\TSeries-12162015-1310-';
    '20151216\TSeries-12162015-1310-';
    '20151216\TSeries-12162015-1310-';
    '20151216\TSeries-12162015-1310-';
    '20151217\TSeries-12172015-1324-';
    '20151217\TSeries-12172015-1324-';
    '20151217\TSeries-12172015-1324-';
    '20151217\TSeries-12172015-1324-';
    '20151217\TSeries-12172015-1324-';
    '20160120\TSeries-01202016-1222-';
    '20160120\TSeries-01202016-1222-';
    '20160120\TSeries-01202016-1222-';
    '20160120\TSeries-01202016-1222-';
    '20160120\TSeries-01202016-1222-';
    '20160120\TSeries-01202016-1222-';
    '20160120\TSeries-01202016-1222-';
    '20160123\TSeries-01232016-1550-';
    '20160123\TSeries-01232016-1550-';
    '20160123\TSeries-01232016-1550-';
    '20160123\TSeries-01232016-1550-';
    '20160123\TSeries-01232016-1550-';
    '20160123\TSeries-01232016-1550-';
    '20160123\TSeries-01232016-1550-';
    '20160123\TSeries-01232016-1550-';
    '20160123\TSeries-01232016-1550-';
    '20160125\TSeries-01252016-1335-';
    '20160125\TSeries-01252016-1335-';
    '20160125\TSeries-01252016-1335-';
    '20160125\TSeries-01252016-1335-';
    '20160125\TSeries-01252016-1335-';
    '20160125\TSeries-01252016-1335-';
    '20160127\TSeries-01272016-1537-';
    '20160127\TSeries-01272016-1537-';
    '20160127\TSeries-01272016-1537-';
    '20160127\TSeries-01272016-1537-';
    '20160127\TSeries-01272016-1537-';
    '20160127\TSeries-01272016-1537-';
    '20160127\TSeries-01272016-1537-';
    '20160203\TSeries-02032016-1257-';
    '20160203\TSeries-02032016-1257-';
    '20160203\TSeries-02032016-1257-';
    '20160203\TSeries-02032016-1257-';
    '20160207\TSeries-02072016-1505-';
    '20160207\TSeries-02072016-1505-';
    '20160207\TSeries-02072016-1505-';
    '20160207\TSeries-02072016-1505-';
    '20160207\TSeries-02072016-1505-';
    '20160207\TSeries-02072016-1505-';
    '20160207\TSeries-02072016-1505-';
    '20160210\TSeries-02102016-1409-';
    '20160210\TSeries-02102016-1409-';
    '20160210\TSeries-02102016-1409-';
    '20160210\TSeries-02102016-1409-';
    '20160210\TSeries-02102016-1409-';
    '20160210\TSeries-02102016-1409-';
    '20160210\TSeries-02102016-1409-';
    '20160210\TSeries-02102016-1409-';
    '20160210\TSeries-02102016-1658-';
    '20160210\TSeries-02102016-1658-';
    '20160210\TSeries-02102016-1658-';
    '20160210\TSeries-02102016-1658-';
    '20160210\TSeries-02102016-1658-';
    '20160210\TSeries-02102016-1658-';
    '20160210\TSeries-02102016-1658-';
    '20160210\TSeries-02102016-1658-';
    '20160210\TSeries-02102016-1658-';
    '20160210\TSeries-02102016-1658-';
    '20160401\TSeries-04012016-1447-';
    '20160401\TSeries-04012016-1447-';
    '20160401\TSeries-04012016-1447-';
    '20160405\TSeries-04052016-1426-';
    '20160405\TSeries-04052016-1426-';
    '20160405\TSeries-04052016-1426-';
    '20160405\TSeries-04052016-1426-';
    '20160405\TSeries-04052016-1426-';
    '20160405\TSeries-04052016-1426-';
    '20160405\TSeries-04052016-1426-';
    '20160405\TSeries-04052016-1426-';
    '20160405\TSeries-04052016-1426-';
    '20160503\TSeries-05032016-0604-';
    '20160503\TSeries-05032016-0604-';
    '20160503\TSeries-05032016-0604-';
    '20160503\TSeries-05032016-0604-';
    '20160503\TSeries-05032016-0604-';
    '20160503\TSeries-05032016-0604-';
    '20160504\TSeries-05042016-0639-';
    '20160504\TSeries-05042016-0639-';
    '20160504\TSeries-05042016-0639-';
    '20160504\TSeries-05042016-0639-';
    '20160504\TSeries-05042016-0639-';
    '20160504\TSeries-05042016-0639-';
    '20160504\TSeries-05042016-0639-';
    '20160504\TSeries-05042016-0639-';
    '20160504\TSeries-05042016-0639-';
    '20160504\TSeries-05042016-0639-';
    '20160504\TSeries-05042016-0639-';
    '20160504\TSeries-05042016-0639-'];


fp_strings = cellstr(fp_strings);

movieNum = [226;229;232;235;238;241;245;251;255;258;265;269;272;275;278;295;298;301;304;308;311;314;317;376;379;382;385;388;391;394;404;407;410;413;415;419;426;429;432;435;444;448;451;454;457;463;509;513;524;527;538;541;544;547;550;301;304;307;311;316;319;322;350;353;354;356;357;359;360;362;363;367;369;370;371;373;377;380;383;384;385;387;388;389;391;392;394;395;397;398;401;403;404;406;410;412;413;415;416;418;419;420;421;423;424;426;427;428;430;432;433;435;436;595;596;598;601;602;604;606;609;611;612;614;615;625;627;629;631;634;636;647;648;649;650;651;652;654;655;656;657;659;660];

for p = 104:length(movieNum)
    close all
    clearvars -except folder1 fp_strings movieNum p
    folder = char(strcat(folder1, fp_strings(p),num2str(movieNum(p))));
    fp.vol=dir([folder '\*.csv']);   fp.vol=([folder '\' fp.vol.name]);
    fp.trace=dir([folder '\*.xls']); fp.trace=([folder '\' fp.trace.name]);
    fp.stim=dir([folder '\*stim*.txt']);
    vFreq=10000; 
    save_plots = 1;
    locs=[1 200 1400 800];
    % load the raw data
    numbg=1;   % number of background ROIs
    try
        data.stim=dlmread([folder '\' fp.stim(1).name]);
        fp.stim=([folder '\' fp.stim(1).name]);
    catch
            data.stim=dlmread([folder '\' fp.stim(2).name]);
            fp.stim=([folder '\' fp.stim(2).name]);
    end
    stimInf.speed=data.stim(1,2); stimInf.bw=data.stim(1,3);
    stimInf.bl=data.stim(1,4); stimInf.rad=data.stim(1,5); data.stim=data.stim(:,1); stimInf.dirs=data.stim;
    data.vol=csvread(fp.vol,1); data.vol=data.vol(:,2); 
    data.trace=xlsread(fp.trace,1);  data.trace(:,1)=[];% data.trace(:,end-1:end)=[]; % first column is index, last two columns are average and error
    data.signal=data.trace(:,1:end-1);
    data.bg=data.trace(:,end);
    %%
    rois=1:size(data.signal,2);
    cell_nums = [1;2;4;5;12];
    %rois=[6];
    rm_rep= [ ];
    %process the voltage data
    [data.vlocs,data.pks]=peakdetect(data.vol,2,25000);
    % remove some rois and reps
    tr=1:size(data.signal,2);rm_roi=tr(~ismember(tr,rois));
    data.signal(:,rm_roi)=[];

    % visualize the whole raw data
    scnsize=get(0,'screensize');
    figure('color',[1 1 1],'position',locs);  hold all;
    for i=1:size(data.signal,2)
        plot(data.signal(:,i)+i*15*size(data.signal,2),'linewidth',1.5);
        legends{i}=['Roi' ' ' num2str(rois(i))];
    end
    legend(legends);
    for i=1:numbg
        plot((data.bg(:,i)-70)*2,'k','linewidth',2); box off; axis off;
    end
    if save_plots == 1
        saveas(gcf,[folder 'RawTraces'],'fig')
    end

    hold off;

    % determine the truncation window size
    fRate=length(data.vol)/size(data.trace,1);
    framerate=vFreq/fRate;
    data.fFrame=ceil(data.vlocs/fRate); data.fFrame(1)=[];
    wsize=min(diff(data.fFrame));
    % data.bg=data.bg-mean(data.bg(data.fFrame(1)-100:data.fFrame(1)-1));
    % data.signal=data.signal-repmat(data.bg,1,size(data.signal,2));
    % remove certain reps from the data
    reps=reshape(1:40,8,5);% reps(:,rm_rep)=[]; 
    reps=reshape(reps,size(reps,1)*size(reps,2),1); nreps=length(reps)/8;
    % clip the dataset
    data.signal=clip(data.signal,data.fFrame,wsize);
    data.bg=clip(data.bg,data.fFrame,wsize);
    data.f0=mean(data.trace(data.fFrame(1)-100:data.fFrame(1)-1,1:end-numbg));
    bgd.f0=mean(data.trace(data.fFrame(1)-100:data.fFrame(1)-1,end-numbg+1:end));
    data.f0=data.f0-bgd.f0;
    % remove certain reps
    data.signal=cellfun(@(x) x(:,reps),data.signal,'uniformoutput',false);
    data.bg=cellfun(@(x) x(:,reps),data.bg,'uniformoutput',false);
    data.fFrame=data.fFrame(reps);
    data.stim=data.stim(reps);

    % analysis the tuning of the data
    dlist=unique(data.stim);
    sortd=zeros(length(data.stim)/8,8);
    for i=1:length(dlist);
        sortd(:,i)= find(dlist(i)==data.stim);
    end
    sortd=reshape(sortd,1,size(sortd,1)*size(sortd,2));
    data.signal=cellfun(@(x) x(:,sortd),data.signal,'uniformoutput',false);
    data.bg=cellfun(@(x) x(:,sortd),data.bg,'uniformoutput',false);

    ptrace=cellfun(@(x) reconstruct(x),data.signal,'uniformoutput',false);
    btrace=cellfun(@(x) reconstruct(x),data.bg,'uniformoutput',false);
    colors=varycolor(nreps);
    for i=1:nreps
        rlegends{i}=['Trial',' ',num2str(i)];
    end 
    for i=1:size(ptrace,2)
        figure('color',[1 1 1],'position',locs.*[1 2 1 .5]);
        for j=1:size(ptrace{i},2)
            hold all;
            plot(ptrace{i}(:,j),'color',colors(j,:),'linewidth',2);
        end
        legend(rlegends); legend boxoff;
        hold off;
        title(['ROI', ' ', num2str(rois(i))]);
        box off;
        save_plot_fn=strcat('Cell',num2str(i));
        save_plot_path = [folder save_plot_fn];
        if save_plots == 1
            saveas(gcf,save_plot_path,'fig') 
        end
    end

    sbins=4;
    data.signal=cellfun(@(x) DataSmooth(x,sbins,'average'),data.signal,'uniformoutput',false);
    data.bg=cellfun(@(x) DataSmooth(x,sbins*2,'average'),data.bg,'uniformoutput',false);

    % subtract background from signal
    data.signal=cellfun(@(x) x-data.bg{1}, data.signal,'uniformoutput',false);

    %calculate the average traces across sweeps
    if nreps>1
        data.signalAve=cellfun(@(x) sweepMean(x,nreps),data.signal,'uniformoutput',false);
        data.bgAve=cellfun(@(x) sweepMean(x,nreps),data.bg,'uniformoutput',false);
        data.signalErr=cellfun(@(x) sweepErr(x,nreps),data.signal,'uniformoutput',false);
        data.bgErr=cellfun(@(x) sweepErr(x,nreps),data.bg,'uniformoutput',false);
    else
        data.signalAve=data.signal;
        data.bgAve=data.bg;

    end
    %analyze  background recording
    for i=1:size(data.bg,2)
        [bgd.f{i},~]=max(data.bg{i});
        bgd.df{i}=(bgd.f{i}-bgd.f0(i))/bgd.f0(i);
        bgd.df{i}=reshape(bgd.df{i},nreps,8);
        [~,bgd.mt{i}]=max(data.bgAve{i});
        if nreps>1
            bgd.mamp{i}=mean(bgd.df{i});
        else
            bgd.mamp{i}=bgd.df{i};
        end
    end

    %framerate = input('What is the framerate for this movie?\n');
    %framerate = framerate/1000;

    % do the same for signal
    for i=1:size(data.signal,2)
        [data.f{i},~]=max(data.signal{i});
        [r, c] = size(data.signal{i});
        time_axis = [0:(1/framerate):((r-1)/framerate)];
        [data.signalSub{i}] = data.signal{i} - data.bg{1,1};
        for k = 1:40
            [data.area{i}(:,k)] = trapz(time_axis, data.signalSub{i}(:,k));
        end
        data.df{i}=(data.f{i}-data.f0(i))/data.f0(i); 
        data.df{i}=reshape(data.df{i},nreps,8);
        data.area{i}=reshape(data.area{i},nreps,8);
        [~,data.mt{i}]=max(data.signalAve{i});
        data.rmt{i}=data.mt{i}-bgd.mt{1};
        if nreps>1
            data.mamp{i}=mean(data.df{i});
        else
            data.mamp{i}=data.df{i};
        end
    end
    [SortedPeakData, SortedAreaData] = loadImagingDataScript2(movieNum(p), data.df, data.area, length(rois), folder, dlist);

    % plot the tuning curves
    result=zeros(size(data.signal,2),11);
    % the column of results are ROI INDEX ;DSI	VS 	PREF ANGLE	CF dF/F0	CP
    % dF/F0	CF Peak Time	CP Peak Time ,Relative CF peak tiem, relative cp
    % time, F0 ... Bg CP time bg CF time
    theta = linspace(0,2*pi,9);
    figuret=figure('color',[1 1 1],'position',locs );
    idxs=zeros(1,size(data.mamp,2));
    nidxs=zeros(1,size(data.mamp,2));

    for i=1:size(data.mamp,2)
        subplot(4,ceil(size(data.mamp,2)/4),i);
        fig=polar(theta,[data.mamp{i},data.mamp{i}(1)],'--k');set(fig,'linewidth',3);
        hold all
        for j=1:size(data.df{1,i},1)
            fig=polar(theta,[data.df{1,i}(j,:) data.df{1,i}(j,1)]);
            set(fig,'linewidth',1.5);
        end
        [ ll,dsi, osi, angle,idx ,nidx, o1dx, o2dx, os_angle, global_osi]=vectorsum(data.mamp{i});
        title(['ROI' ' '  num2str(rois(i))]);
        idxs(i)=idx; nidxs(i)=nidx;
        result(i,:)=[rois(i), dsi,ll, angle, data.mamp{i}(idx),  data.mamp{i}(nidx) , ...,
            data.mt{i}(idx)*fRate/10000, data.mt{i}(nidx)*fRate/10000, data.f0(i), bgd.mt{1}(idx)*fRate/10000,bgd.mt{1}(nidx)*fRate/10000];
    end



    % plot the CF versus CP traces
    figurec=figure('color',[1 1 1],'position',locs);
    t=(1:wsize)*fRate/10000;
    overlay.cf=[];
    overlay.cp=[];
    for i=1:size(result,1)
        idx=idxs(i); nidx=nidxs(i);
        subplot(4,ceil(size(result,1)/4),i);
        plot(t,data.signal{1,i}(:, (idx-1)*nreps+1:idx*nreps),'k','linewidth',2);
        hold on;
        plot(t,data.signal{1,i}(:, (nidx-1)*nreps+1:nidx*nreps),'r','linewidth',2);
        % add the lines for the timing
        l1=line([data.mt{i}(idx)*fRate/10000 data.mt{i}(idx)*fRate/10000],[ data.f0(i)  max(max((data.signal{1,i}(:, (idx-1)*nreps+1:idx*nreps))))] ...,
            ,'color','k','linewidth',2);
        l2=line([data.mt{i}(nidx)*fRate/10000 data.mt{i}(nidx)*fRate/10000],[ data.f0(i)  max(max(data.signal{1,i}(:, (idx-1)*nreps+1:idx*nreps)))] ...,
            ,'color'  ,'r','linewidth',2);
        hold off;
        xlabel('Time/s'); ylabel('Fluorescence');
        title(['ROI', ' ', num2str(rois(i))]); xlim([0 t(end)]); box off;
        overlay.cf=[overlay.cf data.signalAve{1,i}(:,idx)];
        overlay.cp=[overlay.cp data.signalAve{1,i}(:,nidx)];
    end
    %normalize and  plot the overlay
    overlay.cf= overlay.cf./repmat(max(overlay.cf),size(overlay.cf,1),1);
    overlay.cp= overlay.cp./repmat(max(overlay.cp),size(overlay.cp,1),1);
    figureo=figure('color',[1 1 1],'position',locs);
    subplot(1,2,1);
    plot(t,overlay.cf,'linewidth',2);
    xlim([0,t(end)]); xlabel('Centrifugal'); legend(legends); legend boxoff;
    ylabel('Scaled Fluo'); box off;
    subplot(1,2,2);
    plot(t,overlay.cp,'linewidth',2)
    xlim([0,t(end)]); xlabel('Centrapetal');
    ylabel('Scaled Fluo');box off;legend(legends); legend boxoff;


    % % plot the tuning curve of response time
    % figure('color',[1 1 1],'position',locs);
    % for i=1:size(data.mt,2)
    %     idx=idxs(i); nidx=nidxs(i);
    %     subplot(4,ceil(size(data.mamp,2)/4),i);
    %     fig=polar(theta,[data.rmt{i}-data.rmt{i}(idx),data.rmt{i}(1)-data.rmt{i}(idx)],'--k');set(fig,'linewidth',3);
    %     title(['ROI', ' ', num2str(rois(i))]);
    % end


    % plot the background trace
    figureb=figure('color',[1 1 1],'position',locs.*[1 2 1 .5]);
    idx=mode(idxs); nidx=mode(nidxs);
    for i=1:size(data.bg,2)
        plot(t,data.bg{1}(:, (idx-1)*nreps+1:idx*nreps),'k','linewidth',1);
        hold all;
        plot(t,data.bg{1}(:, (nidx-1)*nreps+1:nidx*nreps),'r','linewidth',1);
        hold off;
        xlabel('Time/s'); ylabel('Fluorescence');
        title(['Background', ' ', num2str(i)]); xlim([0 t(end)]); box off;
        % add the lines for the timing
        l1=line([bgd.mt{i}(idx)*fRate/10000 bgd.mt{i}(idx)*fRate/10000],[ bgd.f0(i)  max(max((data.bg{1,i}(:, (idx-1)*nreps+1:idx*nreps))))] ...,
            ,'color','g','linewidth',2);
        l2=line([bgd.mt{i}(nidx)*fRate/10000 bgd.mt{i}(nidx)*fRate/10000],[ bgd.f0(i)  max(max(data.bg{1,i}(:, (idx-1)*nreps+1:idx*nreps)))] ...,
            ,'color'  ,'b','linewidth',2);
    end
    hold on;
    plot(t,data.bgAve{1}(:,idx) ,'k','linewidth',3);
    plot(t,data.bgAve{1}(:,nidx) ,'r','linewidth',3);

    if save_plots == 1
        saveas(figuret,strrep(fp.stim,'.txt','_Tuning.fig'));
        saveas(figurec,strrep(fp.stim,'.txt','_CF_CP_Traces.fig'));
        saveas(figureb,strrep(fp.stim,'.txt','_Background.fig'));
        saveas(figureo,strrep(fp.stim,'.txt','_CF_CP_Scaled_Overlay.fig'));
    end
    figures=figure('color',[1 1 1],'position',locs.*[1 2 1 .5]);
    subplot(1,2,1);
    plot(t,overlay.cf(:,[1 floor(size(overlay.cf,2)/2) size(overlay.cf,2)]),'linewidth',2.5);legend(['ROI',' ',num2str(result(1,1))], ...,
        ['ROI',' ',num2str(result(floor(size(overlay.cp,2)/2),1))],['ROI',' ',num2str(result(size(overlay.cp,2),1))]);
    legend boxoff;
    xlabel('Centrifugal'); xlim([0 t(end)]); box off;
    subplot(1,2,2);
    plot(t,overlay.cp(:,[1 floor(size(overlay.cp,2)/2) size(overlay.cp,2)]),'linewidth',2.5);legend(['ROI',' ',num2str(result(1,1))], ...,
       ['ROI',' ',num2str(result(floor(size(overlay.cp,2)/2),1))],['ROI',' ',num2str(result(size(overlay.cp,2),1))]);
    xlabel('Centrapetal'); xlim([0 t(end)]); box off; legend  boxoff;
    if save_plots == 1
        saveas(figures,strrep(fp.stim,'.txt','_Region_Specific_Kinetics.fig'));
    end

end
    
    
