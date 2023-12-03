%Import Data
clear all
clc
close all

controlD = uiimport;
vgatKO_D = uiimport;

control = controlD.WT_Dataset_10132016;
vgat_KO = vgatKO_D.KO_Dataset06212016;

[m, n] = size(control);
WT_RGCs = m;
[u, v] = size(vgat_KO);
vgat_RGCs = u;

% DATA PROCESSING
%control
for i = 1:m
    %correct index values
    if control(i,14) == 1
        control(i,3) = 0.150;
    elseif control(i, 14) == 2
        control(i,2) = 0.100;
    elseif control(i, 14) == 0
        control(i, 2) = 0.100;
        control(i, 3) = 0.150;
    end
    
    %correct angles
    if control(i,13) == 2
        if (control(i, 5) > 210 && control(i, 5) < 280) || (control(i, 5) > 30 && control(i, 5) < 100)
            if control(i, 5) >= 180
                control(i, 5) = control(i,5) - 180;
            elseif control(i, 5) < 180
                control(i, 5) = control(i,5) + 180;
            end
        elseif (control(i, 10) > 210 && control(i, 10) < 280) || (control(i, 10) > 30 && control(i, 10) < 100)
            if control(i, 10) >= 180
                control(i, 10) = control(i,10) - 180;
            elseif control(i, 10) < 180
                control(i, 10) = control(i,10) + 180;
            end
        end
   end
end

%VgatKO
for i = 1:u
    %correct index values
    if vgat_KO(i,14) == 1
        vgat_KO(i,3) = 0.150;
    elseif vgat_KO(i, 14) == 2
        vgat_KO(i,2) = 0.100;
    elseif vgat_KO(i, 14) == 0
        vgat_KO(i, 2) = 0.100;
        vgat_KO(i, 3) = 0.150;
    end
    
    %correct angles
    if vgat_KO(i,13) == 2
        if (vgat_KO(i, 5) > 210 && vgat_KO(i, 5) < 280) || (vgat_KO(i, 5) > 30 && vgat_KO(i, 5) < 100)
            if vgat_KO(i, 5) >= 180
                vgat_KO(i, 5) = vgat_KO(i,5) - 180;
            elseif vgat_KO(i, 5) < 180
                vgat_KO(i, 5) = vgat_KO(i,5) + 180;
            end
        elseif (vgat_KO(i, 10) > 210 && vgat_KO(i, 10) < 280) || (vgat_KO(i, 10) > 30 && vgat_KO(i, 10) < 100)
            if vgat_KO(i, 10) >= 180
                vgat_KO(i, 10) = vgat_KO(i,10) - 180;
            elseif vgat_KO(i, 10) < 180
                vgat_KO(i, 10) = vgat_KO(i,10) + 180;
            end
        end
    end
end

%Separation of data
for i = 1:3
    WT_data{1,i} = control((control(:,12) == i), :);
    KO_data{1,i} = vgat_KO((vgat_KO(:,12) == i),:);
end

for i = 1:3
    WT_DSdata{1,i} =  WT_data{1,i}(( WT_data{1,i}(:,14) == 1), :);
    KO_DSdata{1,i} =  KO_data{1,i}(( KO_data{1,i}(:,14) == 1), :);
end

for i = 1:3
    WT_OSdata{1,i} =  WT_data{1,i}(( WT_data{1,i}(:,14) == 2), :);
    KO_OSdata{1,i} =  KO_data{1,i}(( KO_data{1,i}(:,14) == 2), :);
end



%statistics: 1)avg dsi 2)std dsi 3)avg osi 4)std osi 5)avg VS 6)std VS
for i = 1:3
    statsWT_DS{1,i} =[mean( WT_DSdata{1,i}(:,2)), std( WT_DSdata{1,i}(:,2)),mean( WT_DSdata{1,i}(:,3)), std( WT_DSdata{1,i}(:,3)), mean( WT_DSdata{1,i}(:,4)), std( WT_DSdata{1,i}(:,4))];
    statsKO_DS{1,i} =[mean( KO_DSdata{1,i}(:,2)), std( KO_DSdata{1,i}(:,2)),mean( KO_DSdata{1,i}(:,3)), std( KO_DSdata{1,i}(:,3)), mean( KO_DSdata{1,i}(:,4)), std( KO_DSdata{1,i}(:,4))];
    
    statsWT_OS{1,i} =[mean( WT_OSdata{1,i}(:,2)), std( WT_OSdata{1,i}(:,2)),mean( WT_OSdata{1,i}(:,3)), std( WT_OSdata{1,i}(:,3)), mean( WT_OSdata{1,i}(:,4)), std( WT_OSdata{1,i}(:,4))];
    statsKO_OS{1,i} =[mean( KO_OSdata{1,i}(:,2)), std( KO_OSdata{1,i}(:,2)),mean( KO_OSdata{1,i}(:,3)), std( KO_OSdata{1,i}(:,3)), mean( KO_OSdata{1,i}(:,4)), std( KO_OSdata{1,i}(:,4))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DS_averages = [ statsWT_DS{1,1}(1),  statsWT_DS{1,2}(1),  statsWT_DS{1,3}(1)];
DS_errors = [ statsWT_DS{1,1}(2),  statsWT_DS{1,2}(2),  statsWT_DS{1,3}(2)];

OS_averages = [ statsWT_OS{1,1}(3),  statsWT_OS{1,2}(3), statsWT_OS{1,3}(3)];
OS_errors = [ statsWT_OS{1,1}(4),  statsWT_OS{1,2}(4),  statsWT_OS{1,3}(4)];

DSpercentages = [size(WT_DSdata{1,1},1)/m*100, size(WT_DSdata{1,2},1)/m*100, size(WT_DSdata{1,3},1)/m*100; size(KO_DSdata{1,1},1)/u*100, size(KO_DSdata{1,2},1)/u*100, size(KO_DSdata{1,3},1)/u*100];
OSpercentages = [size(WT_OSdata{1,1},1)/m*100, size(WT_OSdata{1,2},1)/m*100, size(WT_OSdata{1,3},1)/m*100; size(KO_OSdata{1,1},1)/u*100, size(KO_OSdata{1,2},1)/u*100, size(KO_OSdata{1,3},1)/u*100];

vgatDS_averages = [ statsKO_DS{1,1}(1),  statsKO_DS{1,2}(1),  statsKO_DS{1,3}(1)];
vgatDS_errors = [ statsKO_DS{1,1}(2),  statsKO_DS{1,2}(2),  statsKO_DS{1,3}(2)];

vgatOS_averages = [ statsKO_OS{1,1}(3),  statsKO_OS{1,2}(3), statsKO_OS{1,3}(3)];
vgatOS_errors = [ statsKO_OS{1,1}(4),  statsKO_OS{1,2}(4),  statsKO_OS{1,3}(4)];


versus_DSaverages = [DS_averages; vgatDS_averages];
versus_DSerrors = [DS_errors; vgatDS_errors];

versus_OSaverages = [OS_averages; vgatOS_averages];
versus_OSerrors = [OS_errors; vgatOS_errors];

h_vgat = errorbar_groups(versus_DSaverages,versus_DSerrors,'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'ON-OFF DSGCs','ON DSGCs', 'OFF DSGCs'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Average DSI value','FontSize',12,'FontWeight','bold','Color','k')
legend('WT', 'KO')


h_vgat = errorbar_groups(versus_OSaverages,versus_OSerrors,'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'ON-OFF OSGCs','ON OSGCs','OFF OSGCs'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Average OSI value','FontSize',12,'FontWeight','bold','Color','k')
legend('WT', 'KO')

WT_mice = max(control(:,15));
KO_mice = max(vgat_KO(:,15));

for i = 1:WT_mice
    WT_miceDSData{1,i} = WT_DSdata{1,1}((WT_DSdata{1,1}(:,15) == i),:);
    WT_miceDSData{2,i} = WT_DSdata{1,2}((WT_DSdata{1,2}(:,15) == i),:);
    WT_miceDSData{3,i} = WT_DSdata{1,3}((WT_DSdata{1,3}(:,15) == i),:);
    
    WT_miceOSData{1,i} = WT_OSdata{1,1}((WT_OSdata{1,1}(:,15) == i),:);
    WT_miceOSData{2,i} = WT_OSdata{1,2}((WT_OSdata{1,2}(:,15) == i),:);
    WT_miceOSData{3,i} = WT_OSdata{1,3}((WT_OSdata{1,3}(:,15) == i),:);
end

for i = 1:KO_mice
    KO_miceDSData{1,i} = KO_DSdata{1,1}((KO_DSdata{1,1}(:,15) == i),:);
    KO_miceDSData{2,i} = KO_DSdata{1,2}((KO_DSdata{1,2}(:,15) == i),:);
    KO_miceDSData{3,i} = KO_DSdata{1,3}((KO_DSdata{1,3}(:,15) == i),:);
    
    KO_miceOSData{1,i} = KO_OSdata{1,1}((KO_OSdata{1,1}(:,15) == i),:);
    KO_miceOSData{2,i} = KO_OSdata{1,2}((KO_OSdata{1,2}(:,15) == i),:);
    KO_miceOSData{3,i} = KO_OSdata{1,3}((KO_OSdata{1,3}(:,15) == i),:);
end

%percentages
perc_errors = zeros(2, 3);
h_perc= errorbar_groups(DSpercentages, perc_errors, 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'ON-OFF DSGCs','ON DSGCs','OFF DSGCs'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of RGCs (%)','FontSize',12,'FontWeight','bold','Color','k')
legend('WT', 'KO')

h_perc= errorbar_groups(OSpercentages, perc_errors,'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'ON-OFF OSGCs','ON OSGCs','OFF OSGCs'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of RGCs (%)','FontSize',12,'FontWeight','bold','Color','k')
legend('WT', 'KO')

%%%%%%%%%%%%Percentages by Mouse%%%%%%%%%%%%%%%%%%%%

WT_DSperc_byMouse = zeros(WT_mice, 3);
WT_OSperc_byMouse = zeros(WT_mice, 3);

KO_DSperc_byMouse = zeros(KO_mice, 3);
KO_OSperc_byMouse = zeros(KO_mice, 3);

for i = 1:WT_mice
    WT_DSperc_byMouse(i, 1) = size(WT_miceDSData{1,i},1)/length(find(control(:,15)== i))*100;
    WT_DSperc_byMouse(i, 2) = size(WT_miceDSData{2,i},1)/length(find(control(:,15)== i))*100;
    WT_DSperc_byMouse(i, 3) = size(WT_miceDSData{3,i},1)/length(find(control(:,15)== i))*100;

    WT_OSperc_byMouse(i, 1) = size(WT_miceOSData{1,i},1)/length(find(control(:,15)== i))*100;
    WT_OSperc_byMouse(i, 2) = size(WT_miceOSData{2,i},1)/length(find(control(:,15)== i))*100;
    WT_OSperc_byMouse(i, 3) = size(WT_miceOSData{3,i},1)/length(find(control(:,15)== i))*100;
end

for i = 1:KO_mice
    KO_DSperc_byMouse(i, 1) = size(KO_miceDSData{1,i},1)/length(find(vgat_KO(:,15)== i))*100;
    KO_DSperc_byMouse(i, 2) = size(KO_miceDSData{2,i},1)/length(find(vgat_KO(:,15)== i))*100;
    KO_DSperc_byMouse(i, 3) = size(KO_miceDSData{3,i},1)/length(find(vgat_KO(:,15)== i))*100;

    KO_OSperc_byMouse(i, 1) = size(KO_miceOSData{1,i},1)/length(find(vgat_KO(:,15)== i))*100;
    KO_OSperc_byMouse(i, 2) = size(KO_miceOSData{2,i},1)/length(find(vgat_KO(:,15)== i))*100;
    KO_OSperc_byMouse(i, 3) = size(KO_miceOSData{3,i},1)/length(find(vgat_KO(:,15)== i))*100;
end

stats_DSbyM = [mean(WT_DSperc_byMouse,1), mean(KO_DSperc_byMouse,1); std(WT_DSperc_byMouse,1), std(KO_DSperc_byMouse,1)];

stats_OSbyM = [mean(WT_OSperc_byMouse,1), mean(KO_OSperc_byMouse,1); std(WT_OSperc_byMouse,1), std(KO_OSperc_byMouse,1)];


h_perc= errorbar_groups([stats_DSbyM(1, 1:3); stats_DSbyM(1, 4:6)] , [stats_DSbyM(2, 1:3); stats_DSbyM(2, 4:6)], 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'ON-OFF DSGCs','ON DSGCs','OFF DSGCs'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of RGCs (%, mouse averages)','FontSize',12,'FontWeight','bold','Color','k')
legend('WT', 'KO')

h_perc= errorbar_groups([stats_OSbyM(1, 1:3); stats_OSbyM(1, 4:6)] , [stats_OSbyM(2, 1:3); stats_OSbyM(2, 4:6)], 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'ON-OFF OSGCs','ON OSGCs','OFF OSGCs'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of RGCs (%, mouse averages)','FontSize',12,'FontWeight','bold','Color','k')
legend('WT', 'KO')

%%%%%%%%%%%%%%%%%%%%%%Cummulative
%%%%%%%%%%%%%%%%%%%%%%Distributions%%%%%%%%%%%%

for i = 1:2
     [f_dsgcs{1,i}, x_dsgcs{1,i}] = ecdf( WT_DSdata{1,i}(:,2));
     [f_KOdsgcs{1,i}, x_KOdsgcs{1,i}] = ecdf( KO_DSdata{1,i}(:,2));
end

for i = 1:3
     [f_osgcs{1,i}, x_osgcs{1,i}] = ecdf( WT_OSdata{1,i}(:,3));
     [f_KOosgcs{1,i}, x_KOosgcs{1,i}] = ecdf( KO_OSdata{1,i}(:,3));
end


for i = 1:2
    figure
    hold on
    scatter(x_dsgcs{1,i},f_dsgcs{1,i},'o', 'filled', 'b')
    scatter(x_KOdsgcs{1,i},f_KOdsgcs{1,i},'o', 'filled', 'r')
    set(gcf,'color','white')
    xlabel('DSI value','FontSize',12,'FontWeight','bold','Color','k')
    gname = 'Cummulative frequency f(x)';
    if i == 1
        gname = [gname, ' ON-OFF DSGCs'];
    elseif i == 2
        gname = [gname, ' ON DSGCs'];
    elseif i == 3
        gname = [gname, ' OFF DSGCs'];
    end
    ylabel(gname,'FontSize',12,'FontWeight','bold','Color','k')
    legend('WT', 'KO')
end

for i = 1:3
    figure
    hold on
    scatter(x_osgcs{1,i},f_osgcs{1,i},'o', 'filled', 'b')
    scatter(x_KOosgcs{1,i},f_KOosgcs{1,i},'o', 'filled', 'r')
    set(gcf,'color','white')
    xlabel('OSI value','FontSize',12,'FontWeight','bold','Color','k')
    gname = 'Cummulative frequency f(x)';
    if i == 1
        gname = [gname, ' ON-OFF OSGCs'];
    elseif i == 2
        gname = [gname, ' ON OSGCs'];
    elseif i == 3
        gname = [gname, ' OFF OSGCs'];
    end
    ylabel(gname,'FontSize',12,'FontWeight','bold','Color','k')
    legend('WT', 'KO')
end

%%%%%%%%%%%%%%%%Vector Sum Alignments%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:3
    figure
    set(gcf,'color','white')
    hold on
    gname = 'gDSI Alignments WT';
    for k = 1:size(WT_DSdata{1,i},1)
        WT_DS_polarVS{1,i}(k,:) = [WT_DSdata{1,i}(k,4)*cosd(WT_DSdata{1,i}(k,5)) WT_DSdata{1,i}(k,4)*sind(WT_DSdata{1,i}(k,5))];
        PlotAxisAtOrigin([0 WT_DS_polarVS{1,i}(k,1)],[0 WT_DS_polarVS{1,i}(k,2)],'m-o')
        hold on
    end
    if i == 1
        gname = [gname, ' ON-OFF DSGCs'];
    elseif i == 2
        gname = [gname, ' ON DSGCs'];
    elseif i == 3
        gname = [gname, ' OFF DSGCs'];
    end
    title(gname,'FontSize',12,'FontWeight','bold','Color','k');
end

for i = 1:3
    figure
    set(gcf,'color','white')
    hold on
    gname = 'gOSI Alignments WT';
    for k = 1:size(WT_OSdata{1,i},1)
        WT_OS_polarVS{1,i}(k,:) = [WT_OSdata{1,i}(k,11)*cosd(WT_OSdata{1,i}(k,10)) WT_OSdata{1,i}(k,11)*sind(WT_OSdata{1,i}(k,10))];
        PlotAxisAtOrigin([0 WT_OS_polarVS{1,i}(k,1)],[0 WT_OS_polarVS{1,i}(k,2)],'m-o')
        hold on
    end
    if i == 1
        gname = [gname, ' ON-OFF OSGCs'];
    elseif i == 2
        gname = [gname, ' ON OSGCs'];
    elseif i == 3
        gname = [gname, ' OFF OSGCs'];
    end
    title(gname,'FontSize',12,'FontWeight','bold','Color','k');
end

for i = 1:3
    figure
    set(gcf,'color','white')
    hold on
    gname = 'gDSI Alignments KO';
    for k = 1:size(KO_DSdata{1,i},1)
        KO_DS_polarVS{1,i}(k,:) = [KO_DSdata{1,i}(k,4)*cosd(KO_DSdata{1,i}(k,5)) KO_DSdata{1,i}(k,4)*sind(KO_DSdata{1,i}(k,5))];
        PlotAxisAtOrigin([0 KO_DS_polarVS{1,i}(k,1)],[0 KO_DS_polarVS{1,i}(k,2)],'m-o')
        hold on
    end
    if i == 1
        gname = [gname, ' ON-OFF DSGCs'];
    elseif i == 2
        gname = [gname, ' ON DSGCs'];
    elseif i == 3
        gname = [gname, ' OFF DSGCs'];
    end
    title(gname,'FontSize',12,'FontWeight','bold','Color','k');
end

for i = 1:3
    figure
    set(gcf,'color','white')
    hold on
    gname = 'gOSI Alignments KO';
    for k = 1:size(KO_OSdata{1,i},1)
        KO_OS_polarVS{1,i}(k,:) = [KO_OSdata{1,i}(k,11)*cosd(KO_OSdata{1,i}(k,10)) KO_OSdata{1,i}(k,11)*sind(KO_OSdata{1,i}(k,10))];
        PlotAxisAtOrigin([0 KO_OS_polarVS{1,i}(k,1)],[0 KO_OS_polarVS{1,i}(k,2)],'m-o')
        hold on
    end
    if i == 1
        gname = [gname, ' ON-OFF OSGCs'];
    elseif i == 2
        gname = [gname, ' ON OSGCs'];
    elseif i == 3
        gname = [gname, ' OFF OSGCs'];
    end
    title(gname,'FontSize',12,'FontWeight','bold','Color','k');
end


