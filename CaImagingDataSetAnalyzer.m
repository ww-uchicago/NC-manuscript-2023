%Import Data

controlD = uiimport;

control = controlD.data;

%ko = uiimport;

%ko = ko.data;
[m, n] = size(control);

% DATA PROCESSING
for i = 1:m
    %correct index values
    if control(i,12) == 1
        control(i,3) = 0.150;
    elseif control(i, 12) == 2
        control(i,2) = 0.100;
    elseif control(i, 12) == 0
        control(i, 2) = 0.100;
        control(i, 3) = 0.150;
    end
    
    
    %correct angles
    if control(i,11) == 2
        if (control(i, 5) > 200 && control(i, 5) < 290) || (control(i, 5) > 20 && control(i, 5) < 110)
            if control(i, 5) >= 180
                control(i, 5) = control(i,5) - 180;
            elseif control(i, 5) < 180
                control(i, 5) = control(i,5) + 180;
            end
        end
    end
end

%Separation of data
ON_indexes = find(control(:,10) == 2);
OFF_indexes = find(control(:,10) == 3);
OO_indexes = find(control(:,10) == 1);
DS_indexes = find(control(:,12) == 1);
OS_indexes = find(control(:,12) == 2);
ON_DS_indexes = find((control(:,10) == 2) & (control(:,12) == 1));
OFF_DS_indexes = find(control(:,10) == 3 & control(:,12) == 1);
OO_DS_indexes = find(control(:,10) == 1 & control(:,12) == 1);
ON_OS_indexes = find(control(:,10) == 2 & control(:,12) == 2);
OFF_OS_indexes = find(control(:,10) == 3 & control(:,12) == 2);
OO_OS_indexes = find(control(:,10) == 1 & control(:,12) == 2);

ON_cells = length(ON_indexes);
OFF_cells = length(OFF_indexes);
OO_cells = length(OO_indexes);
DS_cells = length(DS_indexes);
OS_cells = length(OS_indexes);
ON_DS_cells = length(ON_DS_indexes);
OFF_DS_cells = length(OFF_DS_indexes);
OO_DS_cells = length(OO_DS_indexes);
ON_OS_cells = length(ON_OS_indexes);
OFF_OS_cells = length(OFF_OS_indexes);
OO_OS_cells = length(OO_OS_indexes);

ON_Data = control(ON_indexes,:);
OFF_Data = control(OFF_indexes,:);
OO_Data = control(OO_indexes,:);
DS_Data = control(DS_indexes,:);
OS_Data = control(OS_indexes,:);
ON_DS_Data = control(ON_DS_indexes,:);
OFF_DS_Data = control(OFF_DS_indexes,:);
OO_DS_Data = control(OO_DS_indexes,:);
ON_OS_Data = control(ON_OS_indexes,:);
OFF_OS_Data = control(OFF_OS_indexes,:);
OO_OS_Data = control(OO_OS_indexes,:);

%statistics: 1)avg dsi 2)std dsi 3)avg osi 4)std osi 5)avg VS 6)std VS
%7)percentage
stats_allcells = [mean(control(:,2)), std(control(:,2)),mean(control(:,3)), std(control(:,3)), mean(control(:,4)), std(control(:,4))];
stats_ONcells = [mean(ON_Data(:,2)), std(ON_Data(:,2)),mean(ON_Data(:,3)), std(ON_Data(:,3)), mean(ON_Data(:,4)), std(ON_Data(:,4)), ON_cells/m];
stats_OFFcells = [mean(OFF_Data(:,2)), std(OFF_Data(:,2)),mean(OFF_Data(:,3)), std(OFF_Data(:,3)), mean(OFF_Data(:,4)), std(OFF_Data(:,4)),  OFF_cells/m];
stats_OOcells = [mean(OO_Data(:,2)), std(OO_Data(:,2)),mean(OO_Data(:,3)), std(OO_Data(:,3)), mean(OO_Data(:,4)), std(OO_Data(:,4)),  OO_cells/m];
stats_DScells = [mean(DS_Data(:,2)), std(DS_Data(:,2)),mean(DS_Data(:,3)), std(DS_Data(:,3)), mean(DS_Data(:,4)), std(DS_Data(:,4)),  DS_cells/m];
stats_OScells = [mean(OS_Data(:,2)), std(OS_Data(:,2)),mean(OS_Data(:,3)), std(OS_Data(:,3)), mean(OS_Data(:,4)), std(OS_Data(:,4)),  OS_cells/m];
stats_ONDScells = [mean(ON_DS_Data(:,2)), std(ON_DS_Data(:,2)),mean(ON_DS_Data(:,3)), std(ON_DS_Data(:,3)), mean(ON_DS_Data(:,4)), std(ON_DS_Data(:,4)),  ON_DS_cells/m];
stats_OFFDScells = [mean(OFF_DS_Data(:,2)), std(OFF_DS_Data(:,2)),mean(OFF_DS_Data(:,3)), std(OFF_DS_Data(:,3)), mean(OFF_DS_Data(:,4)), std(OFF_DS_Data(:,4)),  OFF_DS_cells/m];
stats_OODScells = [mean(OO_DS_Data(:,2)), std(OO_DS_Data(:,2)),mean(OO_DS_Data(:,3)), std(OO_DS_Data(:,3)), mean(OO_DS_Data(:,4)), std(OO_DS_Data(:,4)),  OO_DS_cells/m];
stats_ONOScells = [mean(ON_OS_Data(:,2)), std(ON_OS_Data(:,2)),mean(ON_OS_Data(:,3)), std(ON_OS_Data(:,3)), mean(ON_OS_Data(:,4)), std(ON_OS_Data(:,4)),  ON_OS_cells/m];
stats_OFFOScells = [mean(OFF_OS_Data(:,2)), std(OFF_OS_Data(:,2)),mean(OFF_OS_Data(:,3)), std(OFF_OS_Data(:,3)), mean(OFF_OS_Data(:,4)), std(OFF_OS_Data(:,4)),  OFF_OS_cells/m];
stats_OOOScells = [mean(OO_OS_Data(:,2)), std(OO_OS_Data(:,2)),mean(OO_OS_Data(:,3)), std(OO_OS_Data(:,3)), mean(OO_OS_Data(:,4)), std(OO_OS_Data(:,4)),  OO_OS_cells/m];

%Plot Stats

DS_averages = [stats_DScells(1), stats_ONDScells(1), stats_OFFDScells(1), stats_OODScells(1)];
DS_errors = [stats_DScells(2), stats_ONDScells(2), stats_OFFDScells(2), stats_OODScells(2)];
h = errorbar_groups(DS_averages, DS_errors);
hold on
set(gca,'XTickLabel',{'All DSGCs','ON DSGCs','OFF DSGCs', 'ON-OFF DSGCs'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Average DSI value','FontSize',12,'FontWeight','bold','Color','k')

OS_averages = [stats_OScells(3), stats_ONOScells(3), stats_OFFOScells(3), stats_OOOScells(3)];
OS_errors = [stats_OScells(4), stats_ONOScells(4), stats_OFFOScells(4), stats_OOOScells(4)];
h = errorbar_groups(OS_averages, OS_errors);
hold on
set(gca,'XTickLabel',{'All OSGCs','ON OSGCs','OFF OSGCs', 'ON-OFF OSGCs'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Average OSI value','FontSize',12,'FontWeight','bold','Color','k')

%Cummulative Distributions
[f_DSIrgcs, x_DSIrgcs] = ecdf(control(:,2));
[f_dsgcs, x_dsgcs] = ecdf(DS_Data(:,2));
[f_onds, x_onds] = ecdf(ON_DS_Data(:,2));
[f_offds, x_offds] = ecdf(OFF_DS_Data(:,2));
[f_oods, x_oods] = ecdf(OO_DS_Data(:,2));

[f_OSIrgcs, x_OSIrgcs] = ecdf(control(:,3));
[f_osgcs, x_osgcs] = ecdf(OS_Data(:,3));
[f_onos, x_onos] = ecdf(ON_OS_Data(:,3));
[f_offos, x_offos] = ecdf(OFF_OS_Data(:,3));
[f_ooos, x_ooos] = ecdf(OO_OS_Data(:,3));

figure
hold on
scatter(x_dsgcs,f_dsgcs,'o', 'filled', 'k')
scatter(x_onds,f_onds,'o', 'filled', 'b')
scatter(x_offds,f_offds,'o', 'filled', 'r')
scatter(x_oods,f_oods,'o', 'filled', 'g')
set(gcf,'color','white')
xlabel('DSI value','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Cummulative frequency f(x)','FontSize',12,'FontWeight','bold','Color','k')
legend('All DSGCs', 'ON DSGCs', 'OFF DSGCs', 'ON-OFF DSGCs')

figure
hold on
scatter(x_osgcs,f_osgcs,'o', 'filled', 'k')
scatter(x_onos,f_onos,'o', 'filled', 'b')
scatter(x_offos,f_offos,'o', 'filled', 'r')
scatter(x_ooos,f_ooos,'o', 'filled', 'g')
set(gcf,'color','white')
xlabel('OSI value','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Cummulative frequency f(x)','FontSize',12,'FontWeight','bold','Color','k')
legend('All OSGCs', 'ON OSGCs', 'OFF OSGCs', 'ON-OFF OSGCs')

%Histogram Analysis
figure
edges = (0:0.2:1);
set(gcf,'color','white')
h_rgcs = histc(control(:,2), edges);

bar(edges, h_rgcs, 'histc')
set(gcf,'color','w');
xlabel('DSI','FontSize',12,'FontWeight','bold')
ylabel('Number of Cells (WT)','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'fontWeight','bold')
legend('All RGCs')
axis([0 1 0 inf])

figure
set(gcf,'color','white')
h_dsgcs = histc(DS_Data(:,2), edges);

bar(edges, h_dsgcs, 'histc')
set(gcf,'color','w');
xlabel('DSI','FontSize',12,'FontWeight','bold')
ylabel('Number of Cells (WT)','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'fontWeight','bold')
legend('DSGCs')
axis([0 1 0 inf])

figure
set(gcf,'color','white')
h_rgcs = histc(control(:,3), edges);

bar(edges, h_rgcs, 'histc')
set(gcf,'color','w');
xlabel('OSI','FontSize',12,'FontWeight','bold')
ylabel('Number of Cells (WT)','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'fontWeight','bold')
legend('All RGCs')
axis([0 1 0 inf])

figure
set(gcf,'color','white')
h_osgcs = histc(OS_Data(:,3), edges);

bar(edges, h_osgcs, 'histc')
set(gcf,'color','w');
xlabel('OSI','FontSize',12,'FontWeight','bold')
ylabel('Number of Cells (WT)','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'fontWeight','bold')
legend('OSGCs')
axis([0 1 0 inf])

%Plot VS
OO_vsAngles = zeros(OO_DS_cells,1);
for k = 1:OO_DS_cells
    OO_vsAngles(k,1) = angleCorrection(OO_DS_Data(k,5),OO_DS_Data(k,6), OO_DS_Data(k,7), OO_DS_Data(k,8), OO_DS_Data(k,9));
    if OO_vsAngles(k,1)< 0
        OO_vsAngles(k,1) = OO_vsAngles(k,1) + 360;
    end
end

ON_vsAngles = zeros(ON_DS_cells,1);
for k = 1:ON_DS_cells
    ON_vsAngles(k,1) = angleCorrection(ON_DS_Data(k,5),ON_DS_Data(k,6), ON_DS_Data(k,7), ON_DS_Data(k,8), ON_DS_Data(k,9));
    if ON_vsAngles(k,1)< 0
        ON_vsAngles(k,1) = ON_vsAngles(k,1) + 360;
    end
end


ONDS_polarVS = zeros(ON_DS_cells,2);
for i = 1:ON_DS_cells
    ONDS_polarVS(i,:) = [ON_DS_Data(i,4)*cosd(ON_vsAngles(i)) ON_DS_Data(i,4)*sind(ON_vsAngles(i))];
end

OFFDS_polarVS = zeros(OFF_DS_cells,2);
for i = 1:OFF_DS_cells
    OFFDS_polarVS(i,:) = [OFF_DS_Data(i,4)*cosd(OFF_DS_Data(i,5)) OFF_DS_Data(i,4)*sind(OFF_DS_Data(i,5))];
end

OODS_polarVS = zeros(OO_DS_cells,2);
for i = 1:OO_DS_cells
    OODS_polarVS(i,:) = [OO_DS_Data(i,4)*cosd(OO_vsAngles(i)) OO_DS_Data(i,4)*sind(OO_vsAngles(i))];
end

figure
set(gcf,'color','white')
hold on
for i=1:ON_DS_cells
    PlotAxisAtOrigin([0 ONDS_polarVS(i,1)],[0 ONDS_polarVS(i,2)],'b-o')
end

figure
set(gcf,'color','white')
hold on
for i=1:OO_DS_cells
    PlotAxisAtOrigin([0 OODS_polarVS(i,1)],[0 OODS_polarVS(i,2)],'g-o')
end

