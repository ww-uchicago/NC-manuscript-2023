%Import Data
clear all
clc
close all

controlD = uiimport;
vgatKO_D = uiimport;

control = controlD.WT_Dataset06212016;
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
ON_indexes = find(control(:,12) == 2);
OFF_indexes = find(control(:,12) == 3);
OO_indexes = find(control(:,12) == 1);
DS_indexes = find(control(:,14) == 1);
OS_indexes = find(control(:,14) == 2);
ON_DS_indexes = find((control(:,12) == 2) & (control(:,14) == 1));
OFF_DS_indexes = find(control(:,12) == 3 & control(:,14) == 1);
OO_DS_indexes = find(control(:,12) == 1 & control(:,14) == 1);
ON_OS_indexes = find(control(:,12) == 2 & control(:,14) == 2);
OFF_OS_indexes = find(control(:,12) == 3 & control(:,14) == 2);
OO_OS_indexes = find(control(:,12) == 1 & control(:,14) == 2);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vgatON_indexes = find(vgat_KO(:,12) == 2);
vgatOFF_indexes = find(vgat_KO(:,12) == 3);
vgatOO_indexes = find(vgat_KO(:,12) == 1);
vgatDS_indexes = find(vgat_KO(:,14) == 1);
vgatOS_indexes = find(vgat_KO(:,14) == 2);
vgatON_DS_indexes = find((vgat_KO(:,12) == 2) & (vgat_KO(:,14) == 1));
vgatOFF_DS_indexes = find(vgat_KO(:,12) == 3 & vgat_KO(:,14) == 1);
vgatOO_DS_indexes = find(vgat_KO(:,12) == 1 & vgat_KO(:,14) == 1);
vgatON_OS_indexes = find(vgat_KO(:,12) == 2 & vgat_KO(:,14) == 2);
vgatOFF_OS_indexes = find(vgat_KO(:,12) == 3 & vgat_KO(:,14) == 2);
vgatOO_OS_indexes = find(vgat_KO(:,12) == 1 & vgat_KO(:,14) == 2);

vgatON_cells = length(vgatON_indexes);
vgatOFF_cells = length(vgatOFF_indexes);
vgatOO_cells = length(vgatOO_indexes);
vgatDS_cells = length(vgatDS_indexes);
vgatOS_cells = length(vgatOS_indexes);
vgatON_DS_cells = length(vgatON_DS_indexes);
vgatOFF_DS_cells = length(vgatOFF_DS_indexes);
vgatOO_DS_cells = length(vgatOO_DS_indexes);
vgatON_OS_cells = length(vgatON_OS_indexes);
vgatOFF_OS_cells = length(vgatOFF_OS_indexes);
vgatOO_OS_cells = length(vgatOO_OS_indexes);

vgatON_Data = vgat_KO(vgatON_indexes,:);
vgatOFF_Data = vgat_KO(vgatOFF_indexes,:);
vgatOO_Data = vgat_KO(vgatOO_indexes,:);
vgatDS_Data = vgat_KO(vgatDS_indexes,:);
vgatOS_Data = vgat_KO(vgatOS_indexes,:);
vgatON_DS_Data = vgat_KO(vgatON_DS_indexes,:);
vgatOFF_DS_Data = vgat_KO(vgatOFF_DS_indexes,:);
vgatOO_DS_Data = vgat_KO(vgatOO_DS_indexes,:);
vgatON_OS_Data = vgat_KO(vgatON_OS_indexes,:);
vgatOFF_OS_Data = vgat_KO(vgatOFF_OS_indexes,:);
vgatOO_OS_Data = vgat_KO(vgatOO_OS_indexes,:);


%statistics: 1)avg dsi 2)std dsi 3)avg osi 4)std osi 5)avg VS 6)std VS
stats_allcells = [mean(control(:,2)), std(control(:,2)),mean(control(:,3)), std(control(:,3)), mean(control(:,4)), std(control(:,4))];
stats_ONcells = [mean(ON_Data(:,2)), std(ON_Data(:,2)),mean(ON_Data(:,3)), std(ON_Data(:,3)), mean(ON_Data(:,4)), std(ON_Data(:,4))];
stats_OFFcells = [mean(OFF_Data(:,2)), std(OFF_Data(:,2)),mean(OFF_Data(:,3)), std(OFF_Data(:,3)), mean(OFF_Data(:,4)), std(OFF_Data(:,4))];
stats_OOcells = [mean(OO_Data(:,2)), std(OO_Data(:,2)),mean(OO_Data(:,3)), std(OO_Data(:,3)), mean(OO_Data(:,4)), std(OO_Data(:,4))];
stats_DScells = [mean(DS_Data(:,2)), std(DS_Data(:,2)),mean(DS_Data(:,3)), std(DS_Data(:,3)), mean(DS_Data(:,4)), std(DS_Data(:,4))];
stats_OScells = [mean(OS_Data(:,2)), std(OS_Data(:,2)),mean(OS_Data(:,3)), std(OS_Data(:,3)), mean(OS_Data(:,4)), std(OS_Data(:,4))];
stats_ONDScells = [mean(ON_DS_Data(:,2)), std(ON_DS_Data(:,2)),mean(ON_DS_Data(:,3)), std(ON_DS_Data(:,3)), mean(ON_DS_Data(:,4)), std(ON_DS_Data(:,4))];
stats_OFFDScells = [mean(OFF_DS_Data(:,2)), std(OFF_DS_Data(:,2)),mean(OFF_DS_Data(:,3)), std(OFF_DS_Data(:,3)), mean(OFF_DS_Data(:,4)), std(OFF_DS_Data(:,4))];
stats_OODScells = [mean(OO_DS_Data(:,2)), std(OO_DS_Data(:,2)),mean(OO_DS_Data(:,3)), std(OO_DS_Data(:,3)), mean(OO_DS_Data(:,4)), std(OO_DS_Data(:,4))];
stats_ONOScells = [mean(ON_OS_Data(:,2)), std(ON_OS_Data(:,2)),mean(ON_OS_Data(:,3)), std(ON_OS_Data(:,3)), mean(ON_OS_Data(:,4)), std(ON_OS_Data(:,4))];
stats_OFFOScells = [mean(OFF_OS_Data(:,2)), std(OFF_OS_Data(:,2)),mean(OFF_OS_Data(:,3)), std(OFF_OS_Data(:,3)), mean(OFF_OS_Data(:,4)), std(OFF_OS_Data(:,4))];
stats_OOOScells = [mean(OO_OS_Data(:,2)), std(OO_OS_Data(:,2)),mean(OO_OS_Data(:,3)), std(OO_OS_Data(:,3)), mean(OO_OS_Data(:,4)), std(OO_OS_Data(:,4))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vgatstats_allcells = [mean(vgat_KO(:,2)), std(vgat_KO(:,2)),mean(vgat_KO(:,3)), std(vgat_KO(:,3)), mean(vgat_KO(:,4)), std(vgat_KO(:,4))];
vgatstats_ONcells = [mean(vgatON_Data(:,2)), std(vgatON_Data(:,2)),mean(vgatON_Data(:,3)), std(vgatON_Data(:,3)), mean(vgatON_Data(:,4)), std(vgatON_Data(:,4))];
vgatstats_OFFcells = [mean(vgatOFF_Data(:,2)), std(vgatOFF_Data(:,2)),mean(vgatOFF_Data(:,3)), std(vgatOFF_Data(:,3)), mean(vgatOFF_Data(:,4)), std(vgatOFF_Data(:,4))];
vgatstats_OOcells = [mean(vgatOO_Data(:,2)), std(vgatOO_Data(:,2)),mean(vgatOO_Data(:,3)), std(vgatOO_Data(:,3)), mean(vgatOO_Data(:,4)), std(vgatOO_Data(:,4))];
vgatstats_DScells = [mean(vgatDS_Data(:,2)), std(vgatDS_Data(:,2)),mean(vgatDS_Data(:,3)), std(vgatDS_Data(:,3)), mean(vgatDS_Data(:,4)), std(vgatDS_Data(:,4))];
vgatstats_OScells = [mean(vgatOS_Data(:,2)), std(vgatOS_Data(:,2)),mean(vgatOS_Data(:,3)), std(vgatOS_Data(:,3)), mean(vgatOS_Data(:,4)), std(vgatOS_Data(:,4))];
vgatstats_ONDScells = [mean(vgatON_DS_Data(:,2)), std(vgatON_DS_Data(:,2)),mean(vgatON_DS_Data(:,3)), std(vgatON_DS_Data(:,3)), mean(vgatON_DS_Data(:,4)), std(vgatON_DS_Data(:,4))];
vgatstats_OFFDScells = [mean(vgatOFF_DS_Data(:,2)), std(vgatOFF_DS_Data(:,2)),mean(vgatOFF_DS_Data(:,3)), std(vgatOFF_DS_Data(:,3)), mean(vgatOFF_DS_Data(:,4)), std(vgatOFF_DS_Data(:,4))];
vgatstats_OODScells = [mean(vgatOO_DS_Data(:,2)), std(vgatOO_DS_Data(:,2)),mean(vgatOO_DS_Data(:,3)), std(vgatOO_DS_Data(:,3)), mean(vgatOO_DS_Data(:,4)), std(vgatOO_DS_Data(:,4))];
vgatstats_ONOScells = [mean(vgatON_OS_Data(:,2)), std(vgatON_OS_Data(:,2)),mean(vgatON_OS_Data(:,3)), std(vgatON_OS_Data(:,3)), mean(vgatON_OS_Data(:,4)), std(vgatON_OS_Data(:,4))];
vgatstats_OFFOScells = [mean(vgatOFF_OS_Data(:,2)), std(vgatOFF_OS_Data(:,2)),mean(vgatOFF_OS_Data(:,3)), std(vgatOFF_OS_Data(:,3)), mean(vgatOFF_OS_Data(:,4)), std(vgatOFF_OS_Data(:,4))];
vgatstats_OOOScells = [mean(vgatOO_OS_Data(:,2)), std(vgatOO_OS_Data(:,2)),mean(vgatOO_OS_Data(:,3)), std(vgatOO_OS_Data(:,3)), mean(vgatOO_OS_Data(:,4)), std(vgatOO_OS_Data(:,4))];

wtAve_gOSI = [mean(OS_Data(:,11)), mean(ON_OS_Data(:,11)), mean(OFF_OS_Data(:,11)), mean(OO_OS_Data(:,11)); std(OS_Data(:,11)), std(ON_OS_Data(:,11)), std(OFF_OS_Data(:,11)), std(OO_OS_Data(:,11))];
koAve_gOSI = [mean(vgatOS_Data(:,11)), mean(vgatON_OS_Data(:,11)), mean(vgatOFF_OS_Data(:,11)), mean(vgatOO_OS_Data(:,11)); std(vgatOS_Data(:,11)), std(vgatON_OS_Data(:,11)), std(vgatOFF_OS_Data(:,11)), std(vgatOO_OS_Data(:,11))];
h = errorbar_groups([wtAve_gOSI(1,:); koAve_gOSI(1,:)], [wtAve_gOSI(2,:); koAve_gOSI(2,:)],'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'All OSGCs','ON OSGCs','OFF OSGCs', 'ON-OFF OSGCs'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Average gOSI value','FontSize',12,'FontWeight','bold','Color','k')
legend('WT', 'KO')


%Plot Stats

DS_averages = [stats_DScells(1), stats_ONDScells(1), stats_OFFDScells(1), stats_OODScells(1)];
DS_errors = [stats_DScells(2), stats_ONDScells(2), stats_OFFDScells(2), stats_OODScells(2)];
h = errorbar_groups(DS_averages, DS_errors,'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'All DSGCs','ON DSGCs','OFF DSGCs', 'ON-OFF DSGCs'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Average DSI value (WT)','FontSize',12,'FontWeight','bold','Color','k')

OS_averages = [stats_OScells(3), stats_ONOScells(3), stats_OFFOScells(3), stats_OOOScells(3)];
OS_errors = [stats_OScells(4), stats_ONOScells(4), stats_OFFOScells(4), stats_OOOScells(4)];
h = errorbar_groups(OS_averages, OS_errors,'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'All OSGCs','ON OSGCs','OFF OSGCs', 'ON-OFF OSGCs'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Average OSI value (WT)','FontSize',12,'FontWeight','bold','Color','k')

%percentages
DSpercentages = [DS_cells/m*100, ON_DS_cells/m*100, OFF_DS_cells/m*100, OO_DS_cells/m*100; vgatDS_cells/u*100, vgatON_DS_cells/u*100, vgatOFF_DS_cells/u*100, vgatOO_DS_cells/u*100];
OSpercentages = [OS_cells/m*100, ON_OS_cells/m*100, OFF_OS_cells/m*100, OO_OS_cells/m*100; vgatOS_cells/u*100, vgatON_OS_cells/u*100, vgatOFF_OS_cells/u*100, vgatOO_OS_cells/u*100];
perc_errors = zeros(2, 4);
h_perc= errorbar_groups(DSpercentages, perc_errors, 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'DSGCs','ON DSGCs','OFF DSGCs', 'ON-OFF DSGCs'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of RGCs (%)','FontSize',12,'FontWeight','bold','Color','k')
legend('WT', 'KO')

h_perc= errorbar_groups([OO_DS_cells/m*100;vgatOO_DS_cells/u*100], [0;0], 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'ON-OFF DSGCs'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of RGCs (%)','FontSize',12,'FontWeight','bold','Color','k')
legend('WT', 'KO')

h_perc= errorbar_groups(OSpercentages, perc_errors,'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'OSGCs','ON OSGCs','OFF OSGCs', 'ON-OFF OSGCs'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of RGCs (%)','FontSize',12,'FontWeight','bold','Color','k')
legend('WT', 'KO')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vgatDS_averages = [vgatstats_DScells(1), vgatstats_ONDScells(1), vgatstats_OFFDScells(1), vgatstats_OODScells(1)];
vgatDS_errors = [vgatstats_DScells(2), vgatstats_ONDScells(2), vgatstats_OFFDScells(2), vgatstats_OODScells(2)];
h_vgat = errorbar_groups(vgatDS_averages,vgatDS_errors,'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'All DSGCs','ON DSGCs','OFF DSGCs', 'ON-OFF DSGCs'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Average DSI value (KO)','FontSize',12,'FontWeight','bold','Color','k')

vgatOS_averages = [vgatstats_OScells(3), vgatstats_ONOScells(3), vgatstats_OFFOScells(3), vgatstats_OOOScells(3)];
vgatOS_errors = [vgatstats_OScells(4), vgatstats_ONOScells(4), vgatstats_OFFOScells(4), vgatstats_OOOScells(4)];
h_vgat = errorbar_groups(vgatOS_averages, vgatOS_errors,'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'All OSGCs','ON OSGCs','OFF OSGCs', 'ON-OFF OSGCs'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Average OSI value (KO)','FontSize',12,'FontWeight','bold','Color','k')

%%%%%%%COMPARISON%%%%%%%%%%%%%%%%

versus_DSaverages = [DS_averages; vgatDS_averages];
versus_DSerrors = [DS_errors; vgatDS_errors];
h_vgat = errorbar_groups(versus_DSaverages,versus_DSerrors,'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'All DSGCs','ON DSGCs','OFF DSGCs', 'ON-OFF DSGCs'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Average DSI value','FontSize',12,'FontWeight','bold','Color','k')
legend('WT', 'KO')

versus_OSaverages = [OS_averages; vgatOS_averages];
versus_OSerrors = [OS_errors; vgatOS_errors];
h_vgat = errorbar_groups(versus_OSaverages,versus_OSerrors,'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'All OSGCs','ON OSGCs','OFF OSGCs', 'ON-OFF OSGCs'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Average OSI value','FontSize',12,'FontWeight','bold','Color','k')
legend('WT', 'KO')


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

[f_ongOSI, x_ongOSI] = ecdf(ON_OS_Data(:,11));
[f_vgatongOSI, x_vgatongOSI] = ecdf(vgatON_OS_Data(:,11));

figure
hold on
scatter(x_ongOSI,f_ongOSI,'o', 'filled', 'b')
scatter(x_vgatongOSI,f_vgatongOSI,'o', 'filled', 'r')
set(gcf,'color','white')
xlabel('gOSI value','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Cummulative frequency f(x)','FontSize',12,'FontWeight','bold','Color','k')
legend('WT ON OSGCs', 'KO ON OSGCs')


figure
hold on
scatter(x_dsgcs,f_dsgcs,'o', 'filled', 'k')
scatter(x_onds,f_onds,'o', 'filled', 'b')
scatter(x_offds,f_offds,'o', 'filled', 'r')
scatter(x_oods,f_oods,'o', 'filled', 'g')
set(gcf,'color','white')
xlabel('DSI value','FontSize',12,'FontWeight','bold','Color','k')
ylabel('WT Cummulative frequency f(x)','FontSize',12,'FontWeight','bold','Color','k')
legend('All DSGCs', 'ON DSGCs', 'OFF DSGCs', 'ON-OFF DSGCs')

figure
hold on
scatter(x_osgcs,f_osgcs,'o', 'filled', 'k')
scatter(x_onos,f_onos,'o', 'filled', 'b')
scatter(x_offos,f_offos,'o', 'filled', 'r')
scatter(x_ooos,f_ooos,'o', 'filled', 'g')
set(gcf,'color','white')
xlabel('OSI value','FontSize',12,'FontWeight','bold','Color','k')
ylabel('WT Cummulative frequency f(x)','FontSize',12,'FontWeight','bold','Color','k')
legend('All OSGCs', 'ON OSGCs', 'OFF OSGCs', 'ON-OFF OSGCs')

%%%%%%%%%%%%%%%%%%VGAT%%%%%%%%%%%%%%%%%%%%%%%%
[f_vgatDSIrgcs, x_vgatDSIrgcs] = ecdf(vgat_KO(:,2));
[f_vgatdsgcs, x_vgatdsgcs] = ecdf(vgatDS_Data(:,2));
[f_vgatonds, x_vgatonds] = ecdf(vgatON_DS_Data(:,2));
%[f_vgatoffds, x_vgatoffds] = ecdf(vgatOFF_DS_Data(:,2));
[f_vgatoods, x_vgatoods] = ecdf(vgatOO_DS_Data(:,2));

[f_vgatOSIrgcs, x_vgatOSIrgcs] = ecdf(vgat_KO(:,3));
[f_vgatosgcs, x_vgatosgcs] = ecdf(vgatOS_Data(:,3));
[f_vgatonos, x_vgatonos] = ecdf(vgatON_OS_Data(:,3));
%[f_vgatoffos, x_vgatoffos] = ecdf(vgatOFF_OS_Data(:,3));
[f_vgatooos, x_vgatooos] = ecdf(vgatOO_OS_Data(:,3));

figure
hold on
scatter(x_vgatdsgcs,f_vgatdsgcs,'o', 'filled', 'k')
scatter(x_vgatonds,f_vgatonds,'o', 'filled', 'b')
%scatter(x_vgatoffds,f_vgatoffds,'o', 'filled', 'r')
scatter(x_vgatoods,f_vgatoods,'o', 'filled', 'g')
set(gcf,'color','white')
xlabel('DSI value','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Vgat KO Cummulative frequency f(x)','FontSize',12,'FontWeight','bold','Color','k')
legend('All DSGCs', 'ON DSGCs', 'ON-OFF DSGCs')

figure
hold on
scatter(x_vgatosgcs,f_vgatosgcs,'o', 'filled', 'k')
scatter(x_vgatonos,f_vgatonos,'o', 'filled', 'b')
%scatter(x_vgatoffos,f_vgatoffos,'o', 'filled', 'r')
scatter(x_vgatooos,f_vgatooos,'o', 'filled', 'g')
set(gcf,'color','white')
xlabel('OSI value','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Vgat KO Cummulative frequency f(x)','FontSize',12,'FontWeight','bold','Color','k')
legend('All OSGCs', 'ON OSGCs', 'ON-OFF OSGCs')

%%%%%%%%%%%%%%%%%%%%COMPARISON%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
scatter(x_dsgcs,f_dsgcs,'o', 'filled', 'b')
scatter(x_vgatdsgcs,f_vgatdsgcs,'o', 'filled', 'r')
set(gcf,'color','white')
xlabel('DSI value','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Cummulative frequency f(x)','FontSize',12,'FontWeight','bold','Color','k')
legend('WT DSGCs', 'KO DSGCs')


figure
hold on
scatter(x_oods,f_oods,'o', 'filled', 'b')
scatter(x_vgatoods,f_vgatoods,'o', 'filled', 'r')
set(gcf,'color','white')
xlabel('DSI value','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Cummulative frequency f(x)','FontSize',12,'FontWeight','bold','Color','k')
legend('WT ON-OFF DSGCs', 'KO ON-OFF DSGCs')


figure
hold on
scatter(x_onds,f_onds,'o', 'filled', 'b')
scatter(x_vgatonds,f_vgatonds,'o', 'filled', 'r')
set(gcf,'color','white')
xlabel('DSI value','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Cummulative frequency f(x)','FontSize',12,'FontWeight','bold','Color','k')
legend('WT ON DSGCs', 'KO ON DSGCs')

% figure
% hold on
% scatter(x_offds,f_offds,'o', 'filled', 'b')
% scatter(x_vgatoffds,f_vgatoffds,'o', 'filled', 'r')
% set(gcf,'color','white')
% xlabel('DSI value','FontSize',12,'FontWeight','bold','Color','k')
% ylabel('Cummulative frequency f(x)','FontSize',12,'FontWeight','bold','Color','k')
% legend('WT OFF DSGCs', 'KO OFF DSGCs')


figure
hold on
scatter(x_osgcs,f_osgcs,'o', 'filled', 'b')
scatter(x_vgatosgcs,f_vgatosgcs,'o', 'filled', 'r')
set(gcf,'color','white')
xlabel('OSI value','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Cummulative frequency f(x)','FontSize',12,'FontWeight','bold','Color','k')
legend('WT OSGCs', 'KO OSGCs')


figure
hold on
scatter(x_ooos,f_ooos,'o', 'filled', 'b')
scatter(x_vgatooos,f_vgatooos,'o', 'filled', 'r')
set(gcf,'color','white')
xlabel('OSI value','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Cummulative frequency f(x)','FontSize',12,'FontWeight','bold','Color','k')
legend('WT ON-OFF OSGCs', 'KO ON-OFF OSGCs')


figure
hold on
scatter(x_onos,f_onos,'o', 'filled', 'b')
scatter(x_vgatonos,f_vgatonos,'o', 'filled', 'r')
set(gcf,'color','white')
xlabel('OSI value','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Cummulative frequency f(x)','FontSize',12,'FontWeight','bold','Color','k')
legend('WT ON OSGCs', 'KO ON OSGCs')
% 
% figure
% hold on
% scatter(x_offos,f_offos,'o', 'filled', 'b')
% %scatter(x_vgatoffos,f_vgatoffos,'o', 'filled', 'r')
% set(gcf,'color','white')
% xlabel('OSI value','FontSize',12,'FontWeight','bold','Color','k')
% ylabel('Cummulative frequency f(x)','FontSize',12,'FontWeight','bold','Color','k')
% legend('WT OFF OSGCs')





%Histogram Analysis
figure
edges = (0:0.05:1);
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
axis([0 1 0 max(h_osgcs)])

figure
set(gcf,'color','white')
h_gOSI = histc(OS_Data(:,11), edges);

bar(edges, h_gOSI, 'histc')
set(gcf,'color','w');
xlabel('gOSI','FontSize',12,'FontWeight','bold')
ylabel('Number of Cells (WT)','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'fontWeight','bold')
axis([0 1 0 inf])

%%%%%%%%%%%%%%%%%%%%%%%%%VGAT%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
edges = (0:0.05:1);
set(gcf,'color','white')
h_vgatrgcs = histc(vgat_KO(:,2), edges);

bar(edges, h_vgatrgcs, 'histc')
set(gcf,'color','w');
xlabel('DSI','FontSize',12,'FontWeight','bold')
ylabel('Number of Cells (KO)','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'fontWeight','bold')
legend('All RGCs')
axis([0 1 0 inf])

figure
set(gcf,'color','white')
h_vgatdsgcs = histc(vgatDS_Data(:,2), edges);

bar(edges, h_vgatdsgcs, 'histc')
set(gcf,'color','w');
xlabel('DSI','FontSize',12,'FontWeight','bold')
ylabel('Number of Cells (KO)','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'fontWeight','bold')
legend('DSGCs')
axis([0 1 0 inf])

figure
set(gcf,'color','white')
h_vgatrgcs = histc(vgat_KO(:,3), edges);

bar(edges, h_vgatrgcs, 'histc')
set(gcf,'color','w');
xlabel('OSI','FontSize',12,'FontWeight','bold')
ylabel('Number of Cells (KO)','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'fontWeight','bold')
legend('All RGCs')
axis([0 1 0 inf])

figure
set(gcf,'color','white')
h_vgatosgcs = histc(vgatOS_Data(:,3), edges);

bar(edges, h_vgatosgcs, 'histc')
set(gcf,'color','w');
xlabel('OSI','FontSize',12,'FontWeight','bold')
ylabel('Number of Cells (KO)','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'fontWeight','bold')
legend('OSGCs')
axis([0 1 0 max(h_osgcs)])

figure
set(gcf,'color','white')
h_vgatgOSI = histc(vgatOS_Data(:,11), edges);

bar(edges, h_vgatgOSI, 'histc')
set(gcf,'color','w');
xlabel('gOSI','FontSize',12,'FontWeight','bold')
ylabel('Number of Cells (KO)','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'fontWeight','bold')
axis([0 1 0 inf])

%Plot VS


ONDS_polarVS = zeros(ON_DS_cells,2);
for i = 1:ON_DS_cells
    ONDS_polarVS(i,:) = [ON_DS_Data(i,4)*cosd(ON_DS_Data(i,5)) ON_DS_Data(i,4)*sind(ON_DS_Data(i,5))];
end

OFFDS_polarVS = zeros(OFF_DS_cells,2);
for i = 1:OFF_DS_cells
    OFFDS_polarVS(i,:) = [OFF_DS_Data(i,4)*cosd(OFF_DS_Data(i,5)) OFF_DS_Data(i,4)*sind(OFF_DS_Data(i,5))];
end

OODS_polarVS = zeros(OO_DS_cells,2);
for i = 1:OO_DS_cells
    OODS_polarVS(i,:) = [OO_DS_Data(i,4)*cosd(OO_DS_Data(i,5)) OO_DS_Data(i,4)*sind(OO_DS_Data(i,5))];
end

figure
set(gcf,'color','white')
hold on
for i=1:ON_DS_cells
    PlotAxisAtOrigin([0 ONDS_polarVS(i,1)],[0 ONDS_polarVS(i,2)],'b-o')
end
title('Vector Sum Alignments WT ON DSGCs')

figure
set(gcf,'color','white')
hold on
for i=1:OO_DS_cells
    PlotAxisAtOrigin([0 OODS_polarVS(i,1)],[0 OODS_polarVS(i,2)],'g-o')
end
title('Vector Sum Alignments WT ON-OFF DSGCs')

figure
set(gcf,'color','white')
hold on
for i=1:OO_DS_cells
    PlotAxisAtOrigin([0 cosd(OO_DS_Data(i,5))],[0 sind(OO_DS_Data(i,5))],'g-o')
end
title('Vector Sum Alignments WT ON-OFF DSGCs')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%VGAT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


vgatONDS_polarVS = zeros(vgatON_DS_cells,2);
for i = 1:vgatON_DS_cells
    vgatONDS_polarVS(i,:) = [vgatON_DS_Data(i,4)*cosd(vgatON_DS_Data(i,5)) vgatON_DS_Data(i,4)*sind(vgatON_DS_Data(i,5))];
end


vgatOODS_polarVS = zeros(vgatOO_DS_cells,2);
for i = 1:vgatOO_DS_cells
    vgatOODS_polarVS(i,:) = [vgatOO_DS_Data(i,4)*cosd(vgatOO_DS_Data(i,5)) vgatOO_DS_Data(i,4)*sind(vgatOO_DS_Data(i,5))];
end

figure
set(gcf,'color','white')
hold on
for i=1:vgatON_DS_cells
    PlotAxisAtOrigin([0 vgatONDS_polarVS(i,1)],[0 vgatONDS_polarVS(i,2)],'m-o')
end
title('Vector Sum Alignments KO ON DSGCs','FontSize',12,'FontWeight','bold','Color','k')

figure
set(gcf,'color','white')
hold on
for i=1:vgatOO_DS_cells
    PlotAxisAtOrigin([0 vgatOODS_polarVS(i,1)],[0 vgatOODS_polarVS(i,2)],'c-o')
end
title('Vector Sum Alignments KO ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')


%%%%%%%Directional
%%%%%%%Preferences%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WTPreference = zeros(OO_DS_cells,1);
for i=1:OO_DS_cells
    if (OO_DS_Data(i,5) >= 110 && OO_DS_Data(i,5) <= 200) && OO_DS_Data(i,14)==1
        WTPreference(i) = 1;
    elseif (OO_DS_Data(i,5) >= 200 && OO_DS_Data(i,5) <= 290) && OO_DS_Data(i,14)==1
        WTPreference(i) = 2;
    elseif (OO_DS_Data(i,5) >= 20 && OO_DS_Data(i,5) <= 110) && OO_DS_Data(i,14)==1
        WTPreference(i) = 3;
    elseif (OO_DS_Data(i,5) >= 290 && OO_DS_Data(i,5) <= 360) && OO_DS_Data(i,14)==1
        WTPreference(i) = 4;
    elseif (OO_DS_Data(i,5) >= 0 && OO_DS_Data(i,5) <= 20) && OO_DS_Data(i,14)==1
        WTPreference(i) = 4;
    end
end




WT_ONPreference = zeros(ON_DS_cells,1);
for i=1:ON_DS_cells
    if (ON_DS_Data(i,5) >= 65 && ON_DS_Data(i,5) <= 140) && ON_DS_Data(i,14)==1
        WT_ONPreference(i) = 1;
    elseif (ON_DS_Data(i,5) >= 205 && ON_DS_Data(i,5) <= 287) && ON_DS_Data(i,14)==1
        WT_ONPreference(i) = 2;
    elseif (ON_DS_Data(i,5) >= 295 && ON_DS_Data(i,5) <= 360) && ON_DS_Data(i,14)==1
        WT_ONPreference(i) = 3;
    elseif (ON_DS_Data(i,5) >= 0 && ON_DS_Data(i,5) <= 15) && ON_DS_Data(i,14)==1
        WT_ONPreference(i) = 3;
    end
end

OO_DS_Data = [OO_DS_Data, WTPreference];
ON_DS_Data = [ON_DS_Data, WT_ONPreference];

WT_DR_OODS = find(OO_DS_Data(:,16)==1);
WT_VR_OODS = find(OO_DS_Data(:,16)==2);

WT_DR_OODS_Data = OO_DS_Data(WT_DR_OODS,:);
WT_VR_OODS_Data = OO_DS_Data(WT_VR_OODS,:);

figure
set(gcf,'color','white')
hold on
for i=1:length(WT_DR_OODS)
    PlotAxisAtOrigin([0 cosd(WT_DR_OODS_Data(i,5))],[0 sind(WT_DR_OODS_Data(i,5))],'b-o')
end
title('Vector Sum Alignments WT ON-OFF DSGCs from dorsal retina')

figure
set(gcf,'color','white')
hold on
for i=1:length(WT_VR_OODS)
    PlotAxisAtOrigin([0 cosd(WT_VR_OODS_Data(i,5))],[0 sind(WT_VR_OODS_Data(i,5))],'b-o')
end
title('Vector Sum Alignments WT ON-OFF DSGCs from ventral retina')

tempIndexes = find(OO_DS_Data(:,17) == 1);
nasalIndexes = find(OO_DS_Data(:,17) == 4);
ventralIndexes = find(OO_DS_Data(:,17) == 2);
dorsalIndexes = find(OO_DS_Data(:,17) == 3);

supIndexes = find(ON_DS_Data(:,17) == 1);
antIndexes = find(ON_DS_Data(:,17) == 3);
infIndexes = find(ON_DS_Data(:,17) == 2);

tempOODScells = length(tempIndexes);
nasalOODScells = length(nasalIndexes);
ventralOODScells = length(ventralIndexes);
dorsalOODScells = length(dorsalIndexes);


tempOO_DS_Data = OO_DS_Data(tempIndexes,:);
nasalOO_DS_Data = OO_DS_Data(nasalIndexes,:);
ventralOO_DS_Data = OO_DS_Data(ventralIndexes,:);
dorsalOO_DS_Data = OO_DS_Data(dorsalIndexes,:);

antON_DS_Data = ON_DS_Data(antIndexes,:);
supON_DS_Data = ON_DS_Data(supIndexes,:);
infON_DS_Data = ON_DS_Data(infIndexes,:);

stats_tempOODScells = [mean(tempOO_DS_Data(:,2)), std(tempOO_DS_Data(:,2)),mean(tempOO_DS_Data(:,3)), std(tempOO_DS_Data(:,3)), mean(tempOO_DS_Data(:,4)), std(tempOO_DS_Data(:,4))];
stats_nasalOODScells = [mean(nasalOO_DS_Data(:,2)), std(nasalOO_DS_Data(:,2)),mean(nasalOO_DS_Data(:,3)), std(nasalOO_DS_Data(:,3)), mean(nasalOO_DS_Data(:,4)), std(nasalOO_DS_Data(:,4))];
stats_ventralOODScells = [mean(ventralOO_DS_Data(:,2)), std(ventralOO_DS_Data(:,2)),mean(ventralOO_DS_Data(:,3)), std(ventralOO_DS_Data(:,3)), mean(ventralOO_DS_Data(:,4)), std(ventralOO_DS_Data(:,4))];
stats_dorsalOODScells = [mean(dorsalOO_DS_Data(:,2)), std(dorsalOO_DS_Data(:,2)),mean(dorsalOO_DS_Data(:,3)), std(dorsalOO_DS_Data(:,3)), mean(dorsalOO_DS_Data(:,4)), std(dorsalOO_DS_Data(:,4))];



vgatPreference = zeros(vgatOO_DS_cells,1);
for i=1:vgatOO_DS_cells
    if (vgatOO_DS_Data(i,5) >= 110 && vgatOO_DS_Data(i,5) <= 200) && vgatOO_DS_Data(i,14)==1
        vgatPreference(i) = 1;
    elseif (vgatOO_DS_Data(i,5) >= 200 && vgatOO_DS_Data(i,5) <= 290) && vgatOO_DS_Data(i,14)==1
        vgatPreference(i) = 2;
    elseif (vgatOO_DS_Data(i,5) >= 20 && vgatOO_DS_Data(i,5) <= 110) && vgatOO_DS_Data(i,14)==1
        vgatPreference(i) = 3;
    elseif (vgatOO_DS_Data(i,5) >= 290 && vgatOO_DS_Data(i,5) <= 360) && vgatOO_DS_Data(i,14)==1
        vgatPreference(i) = 4;
    elseif (vgatOO_DS_Data(i,5) >= 0 && vgatOO_DS_Data(i,5) <= 20) && vgatOO_DS_Data(i,14)==1
        vgatPreference(i) = 4;
    end
end

vgatOO_DS_Data = [vgatOO_DS_Data, vgatPreference];

vgat_DR_OODS = find(vgatOO_DS_Data(:,16)==1);
vgat_VR_OODS = find(vgatOO_DS_Data(:,16)==2);

vgat_DR_OODS_Data = vgatOO_DS_Data(vgat_DR_OODS,:);
vgat_VR_OODS_Data = vgatOO_DS_Data(vgat_VR_OODS,:);

figure
set(gcf,'color','white')
hold on
for i=1:length(vgat_DR_OODS)
    PlotAxisAtOrigin([0 cosd(vgat_DR_OODS_Data(i,5))],[0 sind(vgat_DR_OODS_Data(i,5))],'r-o')
end
title('Vector Sum Alignments KO ON-OFF DSGCs from dorsal retina')

figure
set(gcf,'color','white')
hold on
for i=1:length(vgat_VR_OODS)
    PlotAxisAtOrigin([0 cosd(vgat_VR_OODS_Data(i,5))],[0 sind(vgat_VR_OODS_Data(i,5))],'r-o')
end
title('Vector Sum Alignments WT ON-OFF DSGCs from ventral retina')


vgattempIndexes = find(vgatOO_DS_Data(:,17) == 1);
vgatnasalIndexes = find(vgatOO_DS_Data(:,17) == 4);
vgatventralIndexes = find(vgatOO_DS_Data(:,17) == 2);
vgatdorsalIndexes = find(vgatOO_DS_Data(:,17) == 3);

vgattempOO_DS_Data = vgatOO_DS_Data(vgattempIndexes,:);
vgatnasalOO_DS_Data = vgatOO_DS_Data(vgatnasalIndexes,:);
vgatventralOO_DS_Data = vgatOO_DS_Data(vgatventralIndexes,:);
vgatdorsalOO_DS_Data = vgatOO_DS_Data(vgatdorsalIndexes,:);

stats_vgattempOODScells = [mean(vgattempOO_DS_Data(:,2)), std(vgattempOO_DS_Data(:,2)),mean(vgattempOO_DS_Data(:,3)), std(vgattempOO_DS_Data(:,3)), mean(vgattempOO_DS_Data(:,4)), std(vgattempOO_DS_Data(:,4))];
stats_vgatnasalOODScells = [mean(vgatnasalOO_DS_Data(:,2)), std(vgatnasalOO_DS_Data(:,2)),mean(vgatnasalOO_DS_Data(:,3)), std(vgatnasalOO_DS_Data(:,3)), mean(vgatnasalOO_DS_Data(:,4)), std(vgatnasalOO_DS_Data(:,4))];
stats_vgatventralOODScells = [mean(vgatventralOO_DS_Data(:,2)), std(vgatventralOO_DS_Data(:,2)),mean(vgatventralOO_DS_Data(:,3)), std(vgatventralOO_DS_Data(:,3)), mean(vgatventralOO_DS_Data(:,4)), std(vgatventralOO_DS_Data(:,4))];
stats_vgatdorsalOODScells = [mean(vgatdorsalOO_DS_Data(:,2)), std(vgatdorsalOO_DS_Data(:,2)),mean(vgatdorsalOO_DS_Data(:,3)), std(vgatdorsalOO_DS_Data(:,3)), mean(vgatdorsalOO_DS_Data(:,4)), std(vgatdorsalOO_DS_Data(:,4))];

vgattempOODScells = length(vgattempIndexes);
vgatnasalOODScells = length(vgatnasalIndexes);
vgatventralOODScells = length(vgatventralIndexes);
vgatdorsalOODScells = length(vgatdorsalIndexes);

KO_ONPreference = zeros(vgatON_DS_cells,1);
for i=1:vgatON_DS_cells
    if (vgatON_DS_Data(i,5) >= 65 && vgatON_DS_Data(i,5) <= 140) && vgatON_DS_Data(i,14)==1
        KO_ONPreference(i) = 1;
    elseif (vgatON_DS_Data(i,5) >= 205 && vgatON_DS_Data(i,5) <= 287) && vgatON_DS_Data(i,14)==1
        KO_ONPreference(i) = 2;
    elseif (vgatON_DS_Data(i,5) >= 295 && vgatON_DS_Data(i,5) <= 360) && vgatON_DS_Data(i,14)==1
        KO_ONPreference(i) = 3;
    elseif (vgatON_DS_Data(i,5) >= 0 && vgatON_DS_Data(i,5) <= 15) && vgatON_DS_Data(i,14)==1
        KO_ONPreference(i) = 3;
    end
end

vgatON_DS_Data = [vgatON_DS_Data, KO_ONPreference];

vgatsupIndexes = find(vgatON_DS_Data(:,17) == 1);
vgatantIndexes = find(vgatON_DS_Data(:,17) == 3);
vgatinfIndexes = find(vgatON_DS_Data(:,17) == 2);

vgatantON_DS_Data = vgatON_DS_Data(vgatantIndexes,:);
vgatsupON_DS_Data = vgatON_DS_Data(vgatsupIndexes,:);
vgatinfON_DS_Data = vgatON_DS_Data(vgatinfIndexes,:);

%%%%%%%Orientation
%%%%%%%Preferences%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WTON_OSPreference = zeros(ON_OS_cells,1);
WTOO_OSPreference = zeros(OO_OS_cells,1);
WTOFF_OSPreference = zeros(OFF_OS_cells,1);

for i=1:ON_OS_cells
    if (ON_OS_Data(i,10) >= 110 && ON_OS_Data(i,10) <= 200) && ON_OS_Data(i,14)==2
        WTON_OSPreference(i) = 1;
    elseif (ON_OS_Data(i,10) >= 200 && ON_OS_Data(i,10) <= 290) && ON_OS_Data(i,14)==2
        WTON_OSPreference(i) = 2;
    elseif (ON_OS_Data(i,10) >= 20 && ON_OS_Data(i,10) <= 110) && ON_OS_Data(i,14)==2
        WTON_OSPreference(i) = 2;
    elseif (ON_OS_Data(i,10) >= 290 && ON_OS_Data(i,10) <= 360) && ON_OS_Data(i,14)==2
        WTON_OSPreference(i) = 1;
    elseif (ON_OS_Data(i,10) >= 0 && ON_OS_Data(i,10) <= 20) && ON_OS_Data(i,14)==2
        WTON_OSPreference(i) = 1;
    end
end

for i=1:OO_OS_cells
    if (OO_OS_Data(i,10) >= 110 && OO_OS_Data(i,10) <= 200) && OO_OS_Data(i,14)==2
        WTOO_OSPreference(i) = 1;
    elseif (OO_OS_Data(i,10) >= 200 && OO_OS_Data(i,10) <= 290) && OO_OS_Data(i,14)==2
        WTOO_OSPreference(i) = 2;
    elseif (OO_OS_Data(i,10) >= 20 && OO_OS_Data(i,10) <= 110) && OO_OS_Data(i,14)==2
        WTOO_OSPreference(i) = 2;
    elseif (OO_OS_Data(i,10) >= 290 && OO_OS_Data(i,10) <= 360) && OO_OS_Data(i,14)==2
        WTOO_OSPreference(i) = 1;
    elseif (OO_OS_Data(i,10) >= 0 && OO_OS_Data(i,10) <= 20) && OO_OS_Data(i,14)==2
        WTOO_OSPreference(i) = 1;
    end
end

for i=1:OFF_OS_cells
    if (OFF_OS_Data(i,10) >= 110 && OFF_OS_Data(i,10) <= 200) && OFF_OS_Data(i,14)==2
        WTOFF_OSPreference(i) = 1;
    elseif (OFF_OS_Data(i,10) >= 200 && OFF_OS_Data(i,10) <= 290) && OFF_OS_Data(i,14)==2
        WTOFF_OSPreference(i) = 2;
    elseif (OFF_OS_Data(i,10) >= 20 && OFF_OS_Data(i,10) <= 110) && OFF_OS_Data(i,14)==2
        WTOFF_OSPreference(i) = 2;
    elseif (OFF_OS_Data(i,10) >= 290 && OFF_OS_Data(i,10) <= 360) && OFF_OS_Data(i,14)==2
        WTOFF_OSPreference(i) = 1;
    elseif (OFF_OS_Data(i,10) >= 0 && OFF_OS_Data(i,10) <= 20) && OFF_OS_Data(i,14)==2
        WTOFF_OSPreference(i) = 1;
    end
end

ON_OS_Data = [ON_OS_Data, WTON_OSPreference];
OO_OS_Data = [OO_OS_Data, WTOO_OSPreference];
OFF_OS_Data = [OFF_OS_Data, WTOFF_OSPreference];

WT_DR_OOOS = find(OO_OS_Data(:,16)==1);
WT_VR_OOOS = find(OO_OS_Data(:,16)==2);
WT_DR_OOOS_Data = OO_OS_Data(WT_DR_OOOS,:);
WT_VR_OOOS_Data = OO_OS_Data(WT_VR_OOOS,:);

WT_DR_ONOS = find(ON_OS_Data(:,16)==1);
WT_VR_ONOS = find(ON_OS_Data(:,16)==2);
WT_DR_ONOS_Data = ON_OS_Data(WT_DR_ONOS,:);
WT_VR_ONOS_Data = ON_OS_Data(WT_VR_ONOS,:);

WT_DR_OFFOS = find(OFF_OS_Data(:,16)==1);
WT_VR_OFFOS = find(OFF_OS_Data(:,16)==2);
WT_DR_OFFOS_Data = OFF_OS_Data(WT_DR_OFFOS,:);
WT_VR_OFFOS_Data = OFF_OS_Data(WT_VR_OFFOS,:);

figure
set(gcf,'color','white')
hold on
for i=1:length(WT_DR_OOOS)
    PlotAxisAtOrigin([cosd(WT_DR_OOOS_Data(i,10)+180) cosd(WT_DR_OOOS_Data(i,10))],[sind(WT_DR_OOOS_Data(i,10)+180) sind(WT_DR_OOOS_Data(i,10))],'b-o')
end
title('Vector Sum Alignments WT ON-OFF OSGCs from dorsal retina')

figure
set(gcf,'color','white')
hold on
for i=1:length(WT_VR_OOOS)
    PlotAxisAtOrigin([cosd(WT_VR_OOOS_Data(i,10)+180) cosd(WT_VR_OOOS_Data(i,10))],[sind(WT_VR_OOOS_Data(i,10)+180) sind(WT_VR_OOOS_Data(i,10))],'b-o')
end
title('Vector Sum Alignments WT ON-OFF OSGCs from ventral retina')

figure
set(gcf,'color','white')
hold on
for i=1:length(WT_DR_ONOS)
    PlotAxisAtOrigin([cosd(WT_DR_ONOS_Data(i,10)+180) cosd(WT_DR_ONOS_Data(i,10))],[sind(WT_DR_ONOS_Data(i,10)+180) sind(WT_DR_ONOS_Data(i,10))],'b-o')
end
title('Vector Sum Alignments WT ON OSGCs from dorsal retina')

figure
set(gcf,'color','white')
hold on
for i=1:length(WT_VR_ONOS)
    PlotAxisAtOrigin([cosd(WT_VR_ONOS_Data(i,10)+180) cosd(WT_VR_ONOS_Data(i,10))],[sind(WT_VR_ONOS_Data(i,10)+180) sind(WT_VR_ONOS_Data(i,10))],'b-o')
end
title('Vector Sum Alignments WT ON OSGCs from ventral retina')

figure
set(gcf,'color','white')
hold on
for i=1:length(WT_DR_OFFOS)
    PlotAxisAtOrigin([cosd(WT_DR_OFFOS_Data(i,10)+180) cosd(WT_DR_OFFOS_Data(i,10))],[sind(WT_DR_OFFOS_Data(i,10)+180) sind(WT_DR_OFFOS_Data(i,10))],'b-o')
end
title('Vector Sum Alignments WT OFF OSGCs from dorsal retina')

figure
set(gcf,'color','white')
hold on
for i=1:length(WT_VR_OFFOS)
    PlotAxisAtOrigin([cosd(WT_VR_OFFOS_Data(i,10)+180) cosd(WT_VR_OFFOS_Data(i,10))],[sind(WT_VR_OFFOS_Data(i,10)+180) sind(WT_VR_OFFOS_Data(i,10))],'b-o')
end
title('Vector Sum Alignments WT OFF OSGCs from ventral retina')


hOS_ONIndexes = find(ON_OS_Data(:,17) == 1);
vOS_ONIndexes = find(ON_OS_Data(:,17) == 2);
hOS_OOIndexes = find(OO_OS_Data(:,17) == 1);
vOS_OOIndexes = find(OO_OS_Data(:,17) == 2);
hOS_OFFIndexes = find(OFF_OS_Data(:,17) == 1);
vOS_OFFIndexes = find(OFF_OS_Data(:,17) == 2);

hON_OS_Data = ON_OS_Data(hOS_ONIndexes,:);
vON_OS_Data = ON_OS_Data(vOS_ONIndexes,:);
hOO_OS_Data = OO_OS_Data(hOS_OOIndexes,:);
vOO_OS_Data = OO_OS_Data(vOS_OOIndexes,:);
hOFF_OS_Data = OFF_OS_Data(hOS_OFFIndexes,:);
vOFF_OS_Data = OFF_OS_Data(vOS_OFFIndexes,:);

KOON_OSPreference = zeros(vgatON_OS_cells,1);
KOOO_OSPreference = zeros(vgatOO_OS_cells,1);
KOOFF_OSPreference = zeros(vgatOFF_OS_cells,1);

for i=1:vgatON_OS_cells
    if (vgatON_OS_Data(i,10) >= 110 && vgatON_OS_Data(i,10) <= 200) && vgatON_OS_Data(i,14)==2
        KOON_OSPreference(i) = 1;
    elseif (vgatON_OS_Data(i,10) >= 200 && vgatON_OS_Data(i,10) <= 290) && vgatON_OS_Data(i,14)==2
        KOON_OSPreference(i) = 2;
    elseif (vgatON_OS_Data(i,10) >= 20 && vgatON_OS_Data(i,10) <= 110) && vgatON_OS_Data(i,14)==2
        KOON_OSPreference(i) = 2;
    elseif (vgatON_OS_Data(i,10) >= 290 && vgatON_OS_Data(i,10) <= 360) && vgatON_OS_Data(i,14)==2
        KOON_OSPreference(i) = 1;
    elseif (vgatON_OS_Data(i,10) >= 0 && vgatON_OS_Data(i,10) <= 20) && vgatON_OS_Data(i,14)==2
        KOON_OSPreference(i) = 1;
    end
end

for i=1:vgatOO_OS_cells
    if (vgatOO_OS_Data(i,10) >= 110 && vgatOO_OS_Data(i,10) <= 200) && vgatOO_OS_Data(i,14)==2
        KOOO_OSPreference(i) = 1;
    elseif (vgatOO_OS_Data(i,10) >= 200 && vgatOO_OS_Data(i,10) <= 290) && vgatOO_OS_Data(i,14)==2
        KOOO_OSPreference(i) = 2;
    elseif (vgatOO_OS_Data(i,10) >= 20 && vgatOO_OS_Data(i,10) <= 110) && vgatOO_OS_Data(i,14)==2
        KOOO_OSPreference(i) = 2;
    elseif (vgatOO_OS_Data(i,10) >= 290 && vgatOO_OS_Data(i,10) <= 360) && vgatOO_OS_Data(i,14)==2
        KOOO_OSPreference(i) = 1;
    elseif (vgatOO_OS_Data(i,10) >= 0 && vgatOO_OS_Data(i,10) <= 20) && vgatOO_OS_Data(i,14)==2
        KOOO_OSPreference(i) = 1;
    end
end

for i=1:vgatOFF_OS_cells
    if (vgatOFF_OS_Data(i,10) >= 110 && vgatOFF_OS_Data(i,10) <= 200) && vgatOFF_OS_Data(i,14)==2
        KOOFF_OSPreference(i) = 1;
    elseif (vgatOFF_OS_Data(i,10) >= 200 && vgatOFF_OS_Data(i,10) <= 290) && vgatOFF_OS_Data(i,14)==2
        KOOFF_OSPreference(i) = 2;
    elseif (vgatOFF_OS_Data(i,10) >= 20 && vgatOFF_OS_Data(i,10) <= 110) && vgatOFF_OS_Data(i,14)==2
        KOOFF_OSPreference(i) = 2;
    elseif (vgatOFF_OS_Data(i,10) >= 290 && vgatOFF_OS_Data(i,10) <= 360) && vgatOFF_OS_Data(i,14)==2
        KOOFF_OSPreference(i) = 1;
    elseif (vgatOFF_OS_Data(i,10) >= 0 && vgatOFF_OS_Data(i,10) <= 20) && vgatOFF_OS_Data(i,14)==2
        KOOFF_OSPreference(i) = 1;
    end
end

vgatON_OS_Data = [vgatON_OS_Data, KOON_OSPreference];
vgatOO_OS_Data = [vgatOO_OS_Data, KOOO_OSPreference];
vgatOFF_OS_Data = [vgatOFF_OS_Data, KOOFF_OSPreference];

vgat_DR_OOOS = find(vgatOO_OS_Data(:,16)==1);
vgat_VR_OOOS = find(vgatOO_OS_Data(:,16)==2);
vgat_DR_OOOS_Data = vgatOO_OS_Data(vgat_DR_OOOS,:);
vgat_VR_OOOS_Data = vgatOO_OS_Data(vgat_VR_OOOS,:);

vgat_DR_ONOS = find(vgatON_OS_Data(:,16)==1);
vgat_VR_ONOS = find(vgatON_OS_Data(:,16)==2);
vgat_DR_ONOS_Data = vgatON_OS_Data(vgat_DR_ONOS,:);
vgat_VR_ONOS_Data = vgatON_OS_Data(vgat_VR_ONOS,:);

vgat_DR_OFFOS = find(vgatOFF_OS_Data(:,16)==1);
vgat_VR_OFFOS = find(vgatOFF_OS_Data(:,16)==2);
vgat_DR_OFFOS_Data = vgatOFF_OS_Data(vgat_DR_OFFOS,:);
vgat_VR_OFFOS_Data = vgatOFF_OS_Data(vgat_VR_OFFOS,:);


figure
set(gcf,'color','white')
hold on
for i=1:length(vgat_DR_OOOS)
    PlotAxisAtOrigin([cosd(vgat_DR_OOOS_Data(i,10)+180) cosd(vgat_DR_OOOS_Data(i,10))],[sind(vgat_DR_OOOS_Data(i,10)+180) sind(vgat_DR_OOOS_Data(i,10))],'r-o')
end
title('Vector Sum Alignments KO ON-OFF OSGCs from dorsal retina')

figure
set(gcf,'color','white')
hold on
for i=1:length(vgat_VR_OOOS)
    PlotAxisAtOrigin([cosd(vgat_VR_OOOS_Data(i,10)+180) cosd(vgat_VR_OOOS_Data(i,10))],[sind(vgat_VR_OOOS_Data(i,10)+180) sind(vgat_VR_OOOS_Data(i,10))],'r-o')
end
title('Vector Sum Alignments KO ON-OFF OSGCs from ventral retina')

figure
set(gcf,'color','white')
hold on
for i=1:length(vgat_DR_ONOS)
    PlotAxisAtOrigin([cosd(vgat_DR_ONOS_Data(i,10)+180) cosd(vgat_DR_ONOS_Data(i,10))],[sind(vgat_DR_ONOS_Data(i,10)+180) sind(vgat_DR_ONOS_Data(i,10))],'r-o')
end
title('Vector Sum Alignments KO ON OSGCs from dorsal retina')

figure
set(gcf,'color','white')
hold on
for i=1:length(vgat_VR_ONOS)
    PlotAxisAtOrigin([cosd(vgat_VR_ONOS_Data(i,10)+180) cosd(vgat_VR_ONOS_Data(i,10))],[sind(vgat_VR_ONOS_Data(i,10)+180) sind(vgat_VR_ONOS_Data(i,10))],'r-o')
end
title('Vector Sum Alignments KO ON OSGCs from ventral retina')

figure
set(gcf,'color','white')
hold on
for i=1:length(vgat_DR_OFFOS)
    PlotAxisAtOrigin([cosd(vgat_DR_OFFOS_Data(i,10)+180) cosd(vgat_DR_OFFOS_Data(i,10))],[sind(vgat_DR_OFFOS_Data(i,10)+180) sind(vgat_DR_OFFOS_Data(i,10))],'r-o')
end
title('Vector Sum Alignments KO OFF OSGCs from dorsal retina')

figure
set(gcf,'color','white')
hold on
for i=1:length(vgat_VR_OFFOS)
    PlotAxisAtOrigin([cosd(vgat_VR_OFFOS_Data(i,10)+180) cosd(vgat_VR_OFFOS_Data(i,10))],[sind(vgat_VR_OFFOS_Data(i,10)+180) sind(vgat_VR_OFFOS_Data(i,10))],'r-o')
end
title('Vector Sum Alignments KO OFF OSGCs from ventral retina')


vgathOS_ONIndexes = find(vgatON_OS_Data(:,17) == 1);
vgatvOS_ONIndexes = find(vgatON_OS_Data(:,17) == 2);
vgathOS_OOIndexes = find(vgatOO_OS_Data(:,17) == 1);
vgatvOS_OOIndexes = find(vgatOO_OS_Data(:,17) == 2);
vgathOS_OFFIndexes = find(vgatOFF_OS_Data(:,17) == 1);
vgatvOS_OFFIndexes = find(vgatOFF_OS_Data(:,17) == 2);

vgathON_OS_Data = vgatON_OS_Data(vgathOS_ONIndexes,:);
vgatvON_OS_Data = vgatON_OS_Data(vgatvOS_ONIndexes,:);
vgathOO_OS_Data = vgatOO_OS_Data(vgathOS_OOIndexes,:);
vgatvOO_OS_Data = vgatOO_OS_Data(vgatvOS_OOIndexes,:);
vgathOFF_OS_Data = vgatOFF_OS_Data(vgathOS_OFFIndexes,:);
vgatvOFF_OS_Data = vgatOFF_OS_Data(vgatvOS_OFFIndexes,:);

%%%%%%%%%%%%%%%%%%%% Average DSi by directional
%%%%%%%%%%%%%%%%%%%% Preference%%%%%%%%%%%%%%%%%%%

dir_DSaverages = [stats_tempOODScells(1),stats_nasalOODScells(1),stats_ventralOODScells(1), stats_dorsalOODScells(1) ;...
    stats_vgattempOODScells(1),stats_vgatnasalOODScells(1),stats_vgatventralOODScells(1), stats_vgatdorsalOODScells(1)];
dir_DSerrors = [stats_tempOODScells(2),stats_nasalOODScells(2),stats_ventralOODScells(2), stats_dorsalOODScells(2) ;...
    stats_vgattempOODScells(2),stats_vgatnasalOODScells(2),stats_vgatventralOODScells(2), stats_vgatdorsalOODScells(2)];
h_vgat = errorbar_groups(dir_DSaverages,dir_DSerrors,'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'Temporal','Nasal','Ventral', 'Dorsal'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Average DSI value by Directional Preferrence','FontSize',12,'FontWeight','bold','Color','k')
legend('WT', 'KO')

%%%%%%%%%%%%%%%%%%%%%%%%%Percentages by directional
%%%%%%%%%%%%%%%%%%%%%%%%%preferrence%%%%%%%%%%%%%%%%%%

dirDSpercentages = [tempOODScells/m*100, nasalOODScells/m*100, ventralOODScells/m*100, dorsalOODScells/m*100; vgattempOODScells/u*100, vgatnasalOODScells/u*100, vgatventralOODScells/u*100, vgatdorsalOODScells/u*100];
perc_errors = zeros(2, 4);
h_perc= errorbar_groups(dirDSpercentages, perc_errors, 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'Temporal','Nasal','Ventral', 'Dorsal'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of ON-OFF DSGCs(%) by Directional Preferrence','FontSize',12,'FontWeight','bold','Color','k')
legend('WT', 'KO')


typeDSpercentages = [OO_DS_cells/OO_cells*100, ON_DS_cells/ON_cells*100, OFF_DS_cells/OFF_cells*100; vgatOO_DS_cells/vgatOO_cells*100, vgatON_DS_cells/vgatON_cells*100, vgatOFF_DS_cells/vgatOFF_cells*100];
perc_errors = zeros(2, 3);
h_perc= errorbar_groups(typeDSpercentages, perc_errors, 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'ON-OFF','ON','OFF'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of RGC type(%) showing direction selectivity','FontSize',12,'FontWeight','bold','Color','k')
legend('WT', 'KO')

%%%%%%%%%%%%%%%%%%%% WT Percentage Plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirOODSpercentages = [tempOODScells/m*100, nasalOODScells/m*100, ventralOODScells/m*100, dorsalOODScells/m*100];
perc_errors = zeros(1, 4);
h_perc= errorbar_groups(dirOODSpercentages, perc_errors, 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'Temporal','Nasal','Ventral', 'Dorsal'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of ON-OFF DSGCs(%) by Directional Preferrence','FontSize',12,'FontWeight','bold','Color','k')
legend('WT')

dirONDSpercentages = [length(antIndexes)/m*100, length(supIndexes)/m*100, length(infIndexes)/m*100];
perc_errors = zeros(1, 3);
h_perc= errorbar_groups(dirONDSpercentages, perc_errors, 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'Anterior','Superior','Inferior'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of ON DSGCs(%) by Directional Preferrence','FontSize',12,'FontWeight','bold','Color','k')
legend('WT')


oriONOSpercentages = [length(hOS_ONIndexes)/m*100, length(vOS_ONIndexes)/m*100];
perc_errors = zeros(1, 2);
h_perc= errorbar_groups(oriONOSpercentages, perc_errors, 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'Horizontal Preferring','Vertical Preferring'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of ON OSGCs(%) by Orientation Preferrence','FontSize',12,'FontWeight','bold','Color','k')
legend('WT')

oriOOOSpercentages = [length(hOS_OOIndexes)/m*100, length(vOS_OOIndexes)/m*100];
perc_errors = zeros(1, 2);
h_perc= errorbar_groups(oriOOOSpercentages, perc_errors, 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'Horizontal Preferring','Vertical Preferring'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of ON-OFF OSGCs(%) by Orientation Preferrence','FontSize',12,'FontWeight','bold','Color','k')
legend('WT')

oriOFFOSpercentages = [length(hOS_OFFIndexes)/m*100, length(vOS_OFFIndexes)/m*100];
perc_errors = zeros(1, 2);
h_perc= errorbar_groups(oriOFFOSpercentages, perc_errors, 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'Horizontal Preferring','Vertical Preferring'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of OFF OSGCs(%) by Orientation Preferrence','FontSize',12,'FontWeight','bold','Color','k')
legend('WT')

%%%%%%%%%%%%%%%%%%%% Vgat Percentage Plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vgatdirOODSpercentages = [vgattempOODScells/u*100, vgatnasalOODScells/u*100, vgatventralOODScells/u*100, vgatdorsalOODScells/u*100];
perc_errors = zeros(1, 4);
h_perc= errorbar_groups(vgatdirOODSpercentages, perc_errors, 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'Temporal','Nasal','Ventral', 'Dorsal'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of ON-OFF DSGCs(%) by Directional Preferrence','FontSize',12,'FontWeight','bold','Color','k')
legend('KO')

vgatdirONDSpercentages = [length(vgatantIndexes)/u*100, length(vgatsupIndexes)/u*100, length(vgatinfIndexes)/u*100];
perc_errors = zeros(1, 3);
h_perc= errorbar_groups(vgatdirONDSpercentages, perc_errors, 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'Anterior','Superior','Inferior'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of ON DSGCs(%) by Directional Preferrence','FontSize',12,'FontWeight','bold','Color','k')
legend('KO')


vgatoriONOSpercentages = [length(vgathOS_ONIndexes)/u*100, length(vgatvOS_ONIndexes)/u*100];
perc_errors = zeros(1, 2);
h_perc= errorbar_groups(vgatoriONOSpercentages, perc_errors, 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'Horizontal Preferring','Vertical Preferring'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of ON OSGCs(%) by Orientation Preferrence','FontSize',12,'FontWeight','bold','Color','k')
legend('KO')

vgatoriOOOSpercentages = [length(vgathOS_OOIndexes)/u*100, length(vgatvOS_OOIndexes)/u*100];
perc_errors = zeros(1, 2);
h_perc= errorbar_groups(vgatoriOOOSpercentages, perc_errors, 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'Horizontal Preferring','Vertical Preferring'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of ON-OFF OSGCs(%) by Orientation Preferrence','FontSize',12,'FontWeight','bold','Color','k')
legend('KO')

vgatoriOFFOSpercentages = [length(vgathOS_OFFIndexes)/u*100, length(vgatvOS_OFFIndexes)/u*100];
perc_errors = zeros(1, 2);
h_perc= errorbar_groups(vgatoriOFFOSpercentages, perc_errors, 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'Horizontal Preferring','Vertical Preferring'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of OFF OSGCs(%) by Orientation Preferrence','FontSize',12,'FontWeight','bold','Color','k')
legend('KO')

%%%%%%%%%%%%%%COMPARISON%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

compdirOODSpercentages = [tempOODScells/m*100, nasalOODScells/m*100, ventralOODScells/m*100, dorsalOODScells/m*100; vgattempOODScells/u*100, vgatnasalOODScells/u*100, vgatventralOODScells/u*100, vgatdorsalOODScells/u*100];
perc_errors = zeros(2, 4);
h_perc= errorbar_groups(compdirOODSpercentages, perc_errors, 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'Temporal','Nasal','Ventral', 'Dorsal'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of ON-OFF DSGCs(%) by Directional Preferrence','FontSize',12,'FontWeight','bold','Color','k')
legend('WT','KO')

compdirONDSpercentages = [length(antIndexes)/m*100, length(supIndexes)/m*100, length(infIndexes)/m*100; length(vgatantIndexes)/u*100, length(vgatsupIndexes)/u*100, length(vgatinfIndexes)/u*100];
perc_errors = zeros(2, 3);
h_perc= errorbar_groups(compdirONDSpercentages, perc_errors, 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'Anterior','Superior','Inferior'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of ON DSGCs(%) by Directional Preferrence','FontSize',12,'FontWeight','bold','Color','k')
legend('WT','KO')

comporiONOSpercentages = [length(hOS_ONIndexes)/m*100, length(vOS_ONIndexes)/m*100; length(vgathOS_ONIndexes)/u*100, length(vgatvOS_ONIndexes)/u*100];
perc_errors = zeros(2, 2);
h_perc= errorbar_groups(comporiONOSpercentages, perc_errors, 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'Horizontal Preferring','Vertical Preferring'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of ON OSGCs(%) by Orientation Preferrence','FontSize',12,'FontWeight','bold','Color','k')
legend('WT', 'KO')


comporiOOOSpercentages = [length(hOS_OOIndexes)/m*100, length(vOS_OOIndexes)/m*100; length(vgathOS_OOIndexes)/u*100, length(vgatvOS_OOIndexes)/u*100];
perc_errors = zeros(2, 2);
h_perc= errorbar_groups(comporiOOOSpercentages, perc_errors, 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'Horizontal Preferring','Vertical Preferring'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of ON-OFF OSGCs(%) by Orientation Preferrence','FontSize',12,'FontWeight','bold','Color','k')
legend('WT', 'KO')

comporiOFFOSpercentages = [length(hOS_OFFIndexes)/m*100, length(vOS_OFFIndexes)/m*100; length(vgathOS_OFFIndexes)/u*100, length(vgatvOS_OFFIndexes)/u*100];
perc_errors = zeros(2, 2);
h_perc= errorbar_groups(comporiOFFOSpercentages, perc_errors, 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'Horizontal Preferring','Vertical Preferring'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percentage of OFF OSGCs(%) by Orientation Preferrence','FontSize',12,'FontWeight','bold','Color','k')
legend('WT','KO')

OS_perchange = [(((comporiONOSpercentages(2,1)-comporiONOSpercentages(1,1))/comporiONOSpercentages(1,1))*100),(((comporiOFFOSpercentages(2,1)-comporiOFFOSpercentages(1,1))/comporiOFFOSpercentages(1,1))*100), (((comporiOOOSpercentages(2,1)-comporiOOOSpercentages(1,1))/comporiOOOSpercentages(1,1))*100); (((comporiONOSpercentages(2,2)-comporiONOSpercentages(1,2))/comporiONOSpercentages(1,2))*100), (((comporiOFFOSpercentages(2,2)-comporiOFFOSpercentages(1,2))/comporiOFFOSpercentages(1,2))*100), (((comporiOOOSpercentages(2,2)-comporiOOOSpercentages(1,2))/comporiOOOSpercentages(1,2))*100)];
perc_errors = zeros(2, 3);
h_perc= errorbar_groups(OS_perchange, perc_errors, 'bar_width',0.75,'errorbar_width',0.5);
hold on
set(gca,'XTickLabel',{'ON','OFF','ON-OFF'},'fontweight','bold')
set(gcf,'color','white')
ylabel('Percent Change in Vgat KO by Orientation Preferrence ','FontSize',12,'FontWeight','bold','Color','k')
legend('Horizontal','Vertical')

%%%%%%%%%%%%%%%%%%%%%%%%Angle Histogram%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_OODSWT = OO_DS_Data(:,5)*(pi/180);
figure
rose(theta_OODSWT,15)
set(gcf,'color','white')
title('Angle Histogram of WT ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')

theta_OODSWT_DR = WT_DR_OODS_Data(:,5)*(pi/180);
figure
rose(theta_OODSWT_DR,15)
set(gcf,'color','white')
title('Angle Histogram of WT ON-OFF DSGCs from dorsal retina','FontSize',12,'FontWeight','bold','Color','k')

theta_OODSWT_VR = WT_VR_OODS_Data(:,5)*(pi/180);
figure
rose(theta_OODSWT_VR,15)
set(gcf,'color','white')
title('Angle Histogram of WT ON-OFF DSGCs from ventral retina','FontSize',12,'FontWeight','bold','Color','k')

theta_ONDSWT = ON_DS_Data(:,5)*(pi/180);
figure
rose(theta_ONDSWT,15)
set(gcf,'color','white')
title('Angle Histogram of WT ON DSGCs','FontSize',12,'FontWeight','bold','Color','k')

theta_ONOSWT = ON_OS_Data(:,10)*(pi/180);
theta_ONOSorth = theta_ONOSWT;
for i = 1:length(theta_ONOSWT)
    if theta_ONOSWT(i) <= pi
        theta_ONOSorth(i) = theta_ONOSWT(i) + pi;
    else
        theta_ONOSorth(i) = theta_ONOSWT(i) - pi;
    end
end
theta_ONOSWT = [theta_ONOSWT; theta_ONOSorth];
figure
rose(theta_ONOSWT,15)
set(gcf,'color','white')
title('Angle Histogram of WT ON OSGCs','FontSize',12,'FontWeight','bold','Color','k')

theta_OOOSWT = OO_OS_Data(:,10)*(pi/180);
figure
rose(theta_OOOSWT,15)
set(gcf,'color','white')
title('Angle Histogram of WT ON-OFF OSGCs','FontSize',12,'FontWeight','bold','Color','k')

theta_OFFOSWT = OFF_OS_Data(:,10)*(pi/180);
figure
rose(theta_OFFOSWT,15)
set(gcf,'color','white')
title('Angle Histogram of WT OFF OSGCs','FontSize',12,'FontWeight','bold','Color','k')


theta_OODSKO = vgatOO_DS_Data(:,5)*(pi/180);
figure
rose(theta_OODSKO,24)
set(gcf,'color','white')
title('Angle Histogram of KO ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')

theta_OODSKO_DR = vgat_DR_OODS_Data(:,5)*(pi/180);
figure
rose(theta_OODSKO_DR,15)
set(gcf,'color','white')
title('Angle Histogram of KO ON-OFF DSGCs from dorsal retina','FontSize',12,'FontWeight','bold','Color','k')

theta_OODSKO_VR = vgat_VR_OODS_Data(:,5)*(pi/180);
figure
rose(theta_OODSKO_VR,15)
set(gcf,'color','white')
title('Angle Histogram of KO ON-OFF DSGCs from ventral retina','FontSize',12,'FontWeight','bold','Color','k')

theta_ONDSKO = vgatON_DS_Data(:,5)*(pi/180);
figure
rose(theta_ONDSKO,15)
set(gcf,'color','white')
title('Angle Histogram of KO ON DSGCs','FontSize',12,'FontWeight','bold','Color','k')
%%%%%%%%%%%%%%%%%%%Spiking vs CaImaging%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%spiking = uiimport;
%spikeData = spiking.WTspikeDSI;

%[f_nasaloods, x_nasaloods] = ecdf(nasalOO_DS_Data(:,2));
%[f_spiking, x_spiking] = ecdf(spikeData(:,2));


% figure
% hold on
% scatter(x_nasaloods,f_nasaloods,'o', 'filled', 'b')
% scatter(x_spiking,f_spiking,'o', 'filled', 'r')
% set(gcf,'color','white')
% xlabel('DSI value posterior preferring ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')
% ylabel('Cummulative frequency f(x)','FontSize',12,'FontWeight','bold','Color','k')
% legend('Ca++ Imaging', 'Spiking')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Vector Sum Alignments by Mouse%%%%%%%%%%%%%%%
WTmouse1Indexes = find(OO_DS_Data(:,15) == 1);
WTmouse2Indexes = find(OO_DS_Data(:,15) == 2);
WTmouse3Indexes = find(OO_DS_Data(:,15) == 3);
WTmouse4Indexes = find(OO_DS_Data(:,15) == 4);
WTmouse5Indexes = find(OO_DS_Data(:,15) == 5);
WTmouse6Indexes = find(OO_DS_Data(:,15) == 6);
WTmouse7Indexes = find(OO_DS_Data(:,15) == 7);
WTmouse8Indexes = find(OO_DS_Data(:,15) == 8);

mouse1OO_DS_Data = OO_DS_Data(WTmouse1Indexes,:);
mouse2OO_DS_Data = OO_DS_Data(WTmouse2Indexes,:);
mouse3OO_DS_Data = OO_DS_Data(WTmouse3Indexes,:);
mouse4OO_DS_Data = OO_DS_Data(WTmouse4Indexes,:);
mouse5OO_DS_Data = OO_DS_Data(WTmouse5Indexes,:);
mouse6OO_DS_Data = OO_DS_Data(WTmouse6Indexes,:);
mouse7OO_DS_Data = OO_DS_Data(WTmouse7Indexes,:);
mouse8OO_DS_Data = OO_DS_Data(WTmouse8Indexes,:);

figure
set(gcf,'color','white')
hold on
for i=1:length(WTmouse1Indexes)
    PlotAxisAtOrigin([0 1*cosd(mouse1OO_DS_Data(i,5))],[0 1*sind(mouse1OO_DS_Data(i,5))],'b-o');
end
title('Angle Alignments WT mouse #1 ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')

figure
set(gcf,'color','white')
hold on
for i=1:length(WTmouse2Indexes)
    PlotAxisAtOrigin([0 1*cosd(mouse2OO_DS_Data(i,5))],[0 1*sind(mouse2OO_DS_Data(i,5))],'b-o');
end
title('Angle Alignments WT mouse #2 ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')

figure
set(gcf,'color','white')
hold on
for i=1:length(WTmouse3Indexes)
    PlotAxisAtOrigin([0 1*cosd(mouse3OO_DS_Data(i,5))],[0 1*sind(mouse3OO_DS_Data(i,5))],'b-o');
end
title('Angle Alignments WT mouse #3 ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')

figure
set(gcf,'color','white')
hold on
for i=1:length(WTmouse4Indexes)
    PlotAxisAtOrigin([0 1*cosd(mouse4OO_DS_Data(i,5))],[0 1*sind(mouse4OO_DS_Data(i,5))],'b-o');
end
title('Angle Alignments WT mouse #4 ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')

figure
set(gcf,'color','white')
hold on
for i=1:length(WTmouse5Indexes)
    PlotAxisAtOrigin([0 1*cosd(mouse5OO_DS_Data(i,5))],[0 1*sind(mouse5OO_DS_Data(i,5))],'b-o');
end
title('Angle Alignments WT mouse #5 ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')

figure
set(gcf,'color','white')
hold on
for i=1:length(WTmouse6Indexes)
    PlotAxisAtOrigin([0 1*cosd(mouse6OO_DS_Data(i,5))],[0 1*sind(mouse6OO_DS_Data(i,5))],'b-o');
end
title('Angle Alignments WT mouse #6 ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')

figure
set(gcf,'color','white')
hold on
for i=1:length(WTmouse7Indexes)
    PlotAxisAtOrigin([0 1*cosd(mouse7OO_DS_Data(i,5))],[0 1*sind(mouse7OO_DS_Data(i,5))],'b-o');
end
title('Angle Alignments WT mouse #7 ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')

figure
set(gcf,'color','white')
hold on
for i=1:length(WTmouse8Indexes)
    PlotAxisAtOrigin([0 1*cosd(mouse8OO_DS_Data(i,5))],[0 1*sind(mouse8OO_DS_Data(i,5))],'b-o');
end
title('Angle Alignments WT mouse #8 ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')


KOmouse1Indexes = find(vgatOO_DS_Data(:,15) == 1);
KOmouse2Indexes = find(vgatOO_DS_Data(:,15) == 2);
KOmouse3Indexes = find(vgatOO_DS_Data(:,15) == 3);
KOmouse4Indexes = find(vgatOO_DS_Data(:,15) == 4);
KOmouse5Indexes = find(vgatOO_DS_Data(:,15) == 5);
KOmouse6Indexes = find(vgatOO_DS_Data(:,15) == 6);
KOmouse7Indexes = find(vgatOO_DS_Data(:,15) == 7);
KOmouse8Indexes = find(vgatOO_DS_Data(:,15) == 8);
KOmouse9Indexes = find(vgatOO_DS_Data(:,15) == 9);
KOmouse10Indexes = find(vgatOO_DS_Data(:,15) == 10);
KOmouse11Indexes = find(vgatOO_DS_Data(:,15) == 11);
KOmouse12Indexes = find(vgatOO_DS_Data(:,15) == 12);

KOmouse1OO_DS_Data = vgatOO_DS_Data(KOmouse1Indexes,:);
KOmouse2OO_DS_Data = vgatOO_DS_Data(KOmouse2Indexes,:);
KOmouse3OO_DS_Data = vgatOO_DS_Data(KOmouse3Indexes,:);
KOmouse4OO_DS_Data = vgatOO_DS_Data(KOmouse4Indexes,:);
KOmouse5OO_DS_Data = vgatOO_DS_Data(KOmouse5Indexes,:);
KOmouse6OO_DS_Data = vgatOO_DS_Data(KOmouse6Indexes,:);
KOmouse7OO_DS_Data = vgatOO_DS_Data(KOmouse7Indexes,:);
KOmouse8OO_DS_Data = vgatOO_DS_Data(KOmouse8Indexes,:);
KOmouse9OO_DS_Data = vgatOO_DS_Data(KOmouse9Indexes,:);
KOmouse10OO_DS_Data = vgatOO_DS_Data(KOmouse10Indexes,:);
KOmouse11OO_DS_Data = vgatOO_DS_Data(KOmouse11Indexes,:);
KOmouse12OO_DS_Data = vgatOO_DS_Data(KOmouse12Indexes,:);

if isempty(KOmouse1Indexes) == 0
    figure
    set(gcf,'color','white')
    hold on
    for i=1:length(KOmouse1Indexes)
        PlotAxisAtOrigin([0 1*cosd(KOmouse1OO_DS_Data(i,5))],[0 1*sind(KOmouse1OO_DS_Data(i,5))],'r-o');
    end
    title('Angle Alignments KO mouse #1 ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')
end

if isempty(KOmouse2Indexes) == 0
    figure
    set(gcf,'color','white')
    hold on
    for i=1:length(KOmouse2Indexes)
        PlotAxisAtOrigin([0 1*cosd(KOmouse2OO_DS_Data(i,5))],[0 1*sind(KOmouse2OO_DS_Data(i,5))],'r-o');
    end
    title('Angle Alignments KO mouse #2 ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')
end

if isempty(KOmouse3Indexes) == 0
    figure
    set(gcf,'color','white')
    hold on
    for i=1:length(KOmouse3Indexes)
        PlotAxisAtOrigin([0 1*cosd(KOmouse3OO_DS_Data(i,5))],[0 1*sind(KOmouse3OO_DS_Data(i,5))],'r-o');
    end
    title('Angle Alignments KO mouse #3 ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')
end

if isempty(KOmouse4Indexes) == 0
    figure
    set(gcf,'color','white')
    hold on
    for i=1:length(KOmouse4Indexes)
        PlotAxisAtOrigin([0 1*cosd(KOmouse4OO_DS_Data(i,5))],[0 1*sind(KOmouse4OO_DS_Data(i,5))],'r-o');
    end
    title('Angle Alignments KO mouse #4 ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')
end

if isempty(KOmouse5Indexes) == 0
    figure
    set(gcf,'color','white')
    hold on
    for i=1:length(KOmouse5Indexes)
        PlotAxisAtOrigin([0 1*cosd(KOmouse5OO_DS_Data(i,5))],[0 1*sind(KOmouse5OO_DS_Data(i,5))],'r-o');
    end
    title('Angle Alignments KO mouse #5 ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')
end

if isempty(KOmouse6Indexes) == 0
    figure
    set(gcf,'color','white')
    hold on
    for i=1:length(KOmouse6Indexes)
        PlotAxisAtOrigin([0 1*cosd(KOmouse6OO_DS_Data(i,5))],[0 1*sind(KOmouse6OO_DS_Data(i,5))],'r-o');
    end
    title('Angle Alignments KO mouse #6 ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')
end

if isempty(KOmouse7Indexes) == 0
    figure
    set(gcf,'color','white')
    hold on
    for i=1:length(KOmouse7Indexes)
        PlotAxisAtOrigin([0 1*cosd(KOmouse7OO_DS_Data(i,5))],[0 1*sind(KOmouse7OO_DS_Data(i,5))],'r-o');
    end
    title('Angle Alignments KO mouse #7 ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')
end

if isempty(KOmouse8Indexes) == 0
    figure
    set(gcf,'color','white')
    hold on
    for i=1:length(KOmouse8Indexes)
        PlotAxisAtOrigin([0 1*cosd(KOmouse8OO_DS_Data(i,5))],[0 1*sind(KOmouse8OO_DS_Data(i,5))],'r-o');
    end
    title('Angle Alignments KO mouse #8 ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')
end

if isempty(KOmouse9Indexes) == 0
    figure
    set(gcf,'color','white')
    hold on
    for i=1:length(KOmouse9Indexes)
        PlotAxisAtOrigin([0 1*cosd(KOmouse9OO_DS_Data(i,5))],[0 1*sind(KOmouse9OO_DS_Data(i,5))],'r-o');
    end
    title('Angle Alignments KO mouse #9 ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')
end

if isempty(KOmouse10Indexes) == 0
    figure
    set(gcf,'color','white')
    hold on
    for i=1:length(KOmouse10Indexes)
        PlotAxisAtOrigin([0 1*cosd(KOmouse10OO_DS_Data(i,5))],[0 1*sind(KOmouse10OO_DS_Data(i,5))],'r-o');
    end
    title('Angle Alignments KO mouse #10 ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')
end

if isempty(KOmouse11Indexes) == 0
    figure
    set(gcf,'color','white')
    hold on
    for i=1:length(KOmouse11Indexes)
        PlotAxisAtOrigin([0 1*cosd(KOmouse11OO_DS_Data(i,5))],[0 1*sind(KOmouse11OO_DS_Data(i,5))],'r-o');
    end
    title('Angle Alignments KO mouse #11 ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')
end

if isempty(KOmouse12Indexes) == 0
    figure
    set(gcf,'color','white')
    hold on
    for i=1:length(KOmouse12Indexes)
        PlotAxisAtOrigin([0 1*cosd(KOmouse12OO_DS_Data(i,5))],[0 1*sind(KOmouse12OO_DS_Data(i,5))],'r-o');
    end
    title('Angle Alignments KO mouse #12 ON-OFF DSGCs','FontSize',12,'FontWeight','bold','Color','k')
end

theta_OODSWT_M8 = mouse8OO_DS_Data(:,5)*(pi/180);
figure
rose(theta_OODSWT_M8,24)
set(gcf,'color','white')
title('Angle Histogram of WT ON-OFF DSGCs Mouse #8','FontSize',12,'FontWeight','bold','Color','k')

theta_OODSKO_M6 = KOmouse6OO_DS_Data(:,5)*(pi/180);
figure
rose(theta_OODSKO_M6,24)
set(gcf,'color','white')
title('Angle Histogram of KO ON-OFF DSGCs Mouse #6','FontSize',12,'FontWeight','bold','Color','k')

%%%%%%%%%%%%%%%%%%OS axes%%%%%%%%%%%%%%%%%%

figure
set(gcf,'color','white')
hold on
for i=1:OS_cells
    PlotAxisAtOrigin([1*cosd(OS_Data(i,10)+180) 1*cosd(OS_Data(i,10))],[1*sind(OS_Data(i,10)+180) 1*sind(OS_Data(i,10))],'b-o');
end
title('OS Angle Alignments WT All OSGCs','FontSize',12,'FontWeight','bold','Color','k')

figure
set(gcf,'color','white')
hold on
for i=1:ON_OS_cells
    PlotAxisAtOrigin([1*cosd(ON_OS_Data(i,10)+180) 1*cosd(ON_OS_Data(i,10))],[1*sind(ON_OS_Data(i,10)+180) 1*sind(ON_OS_Data(i,10))],'b-o');
end
title('OS Angle Alignments WT ON OSGCs','FontSize',12,'FontWeight','bold','Color','k')

figure
set(gcf,'color','white')
hold on
for i=1:OFF_OS_cells
    PlotAxisAtOrigin([1*cosd(OFF_OS_Data(i,10)+180) 1*cosd(OFF_OS_Data(i,10))],[1*sind(OFF_OS_Data(i,10)+180) 1*sind(OFF_OS_Data(i,10))],'b-o');
end
title('OS Angle Alignments WT OFF OSGCs','FontSize',12,'FontWeight','bold','Color','k')

figure
set(gcf,'color','white')
hold on
for i=1:OO_OS_cells
    PlotAxisAtOrigin([1*cosd(OO_OS_Data(i,10)+180) 1*cosd(OO_OS_Data(i,10))],[1*sind(OO_OS_Data(i,10)+180) 1*sind(OO_OS_Data(i,10))],'b-o');
end
title('OS Angle Alignments WT ON-OFF OSGCs','FontSize',12,'FontWeight','bold','Color','k')

figure
set(gcf,'color','white')
hold on
for i=1:vgatOS_cells
    PlotAxisAtOrigin([1*cosd(vgatOS_Data(i,10)+180) 1*cosd(vgatOS_Data(i,10))],[1*sind(vgatOS_Data(i,10)+180) 1*sind(vgatOS_Data(i,10))],'r-o');
end
title('OS Angle Alignments KO All OSGCs','FontSize',12,'FontWeight','bold','Color','k')

figure
set(gcf,'color','white')
hold on
for i=1:vgatON_OS_cells
    PlotAxisAtOrigin([1*cosd(vgatON_OS_Data(i,10)+180) 1*cosd(vgatON_OS_Data(i,10))],[1*sind(vgatON_OS_Data(i,10)+180) 1*sind(vgatON_OS_Data(i,10))],'r-o');
end
title('OS Angle Alignments KO ON OSGCs','FontSize',12,'FontWeight','bold','Color','k')

figure
set(gcf,'color','white')
hold on
for i=1:vgatOO_OS_cells
    PlotAxisAtOrigin([1*cosd(vgatOO_OS_Data(i,10)+180) 1*cosd(vgatOO_OS_Data(i,10))],[1*sind(vgatOO_OS_Data(i,10)+180) 1*sind(vgatOO_OS_Data(i,10))],'r-o');
end
title('OS Angle Alignments KO ON-OFF OSGCs','FontSize',12,'FontWeight','bold','Color','k')

figure
set(gcf,'color','white')
hold on
for i=1:vgatOFF_OS_cells
    PlotAxisAtOrigin([1*cosd(vgatOFF_OS_Data(i,10)+180) 1*cosd(vgatOFF_OS_Data(i,10))],[1*sind(vgatOFF_OS_Data(i,10)+180) 1*sind(vgatOFF_OS_Data(i,10))],'r-o');
end
title('OS Angle Alignments KO OFF OSGCs','FontSize',12,'FontWeight','bold','Color','k')