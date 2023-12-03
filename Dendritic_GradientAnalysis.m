
close all
clearvars -except Alexa594 Alexa488 ASAP3b
fullfig
for i = 1:size(Alexa594,2)
    max_dist = [];
    for j = 1:size(Alexa594{1,i},2)
        Alexa594{1,i}{1,j} = round(Alexa594{1,i}{1,j},4);
        max_dist(1,j) = max(Alexa594{1,i}{1,j}(:,1));
    end
    max_dist = max(max_dist);
    step_size = Alexa594{1,i}{1,1}(2,1) - Alexa594{1,i}{1,1}(1,1);
    dend_length = [step_size:step_size:max_dist]';
    Alexa594_summary{1,i} = [];
    for k = 1:length(dend_length)
        datapts = [];
        for q = 1:size(Alexa594{1,i},2)
            index = intersect(find(Alexa594{1,i}{1,q}(:,1) >= (dend_length(k) - 0.08)), find(Alexa594{1,i}{1,q}(:,1) <= (dend_length(k) + 0.08)));
            datapts = [datapts, Alexa594{1,i}{1,q}(index,2)];
        end
        Alexa594_summary{1,i}(k,1:4) = [dend_length(k), mean(datapts), std(datapts), (std(datapts)/sqrt(length(datapts)))];

    end
    try
        [f2, g2, o2] = fit(dend_length, [Alexa594_summary{1,i}(:,2)], 'exp2');
        [f2_e1, g2_e1, o2_e1] = fit(dend_length, [Alexa594_summary{1,i}(:,2)], 'exp1');
        if (g2_e1.rsquare > g2.rsquare) || (g2_e1.adjrsquare > g2.adjrsquare)
            clearvars f2 g2 o2
            f2 = f2_e1;
            g2 = g2_e1;
            o2 = o2_e1;
        end
    catch

        try 
            [f2, g2, o2] = fit(dend_length, [Alexa594_summary{1,i}(:,2)], 'exp2');
        catch
            [f2, g2, o2] = fit(dend_length, [Alexa594_summary{1,i}(:,2)], 'exp1');
        end
    end         
    sbplt(1,i) = subplot(4, 5, i);
    for r = 1:size(Alexa594{1,i},2)
        plot(Alexa594{1,i}{1,r}(:,1), Alexa594{1,i}{1,r}(:,2), 'r')
        hold on
    end
    x_trace = [0:0.1:145]';
    Alexa594_fit(:,i) = f2(x_trace);
    plot(Alexa594_summary{1,i}(:,1), Alexa594_summary{1,i}(:,2), 'k', 'LineWidth', 2);
    plot(dend_length, f2(dend_length),  'b', 'LineWidth', 2);
    hold off
    title(['Cell ' num2str(i)])
    ylabel('a.u.')
    xlabel('dendritic distance (um)')
    sgtitle('Alexa 594')
    hold off   
    clearvars f2 g2 o2
end
sbplt(1, (i+1)) = subplot(4, 5, (i+1));
shadedErrorBar(x_trace, mean(Alexa594_fit,2), [std(Alexa594_fit,0,2)/sqrt(size(Alexa594_fit,2))],'lineprops', {'-b','Linewidth',2})
title(['Summary'])
ylabel('a.u.')
xlabel('dendritic distance (um)')


fullfig
for i = 1:size(Alexa594,2)
    for j = 1:size(Alexa594{1,i},2)
        max_dist(1,j) = max(Alexa594{1,i}{1,j}(:,1));
        Alexa594_relRad{1,i}{1,j} = [Alexa594{1,i}{1,j}(:,1) sqrt((Alexa594{1,i}{1,j}(:,2)/mean(Alexa594_summary{1,i}(1:3,2))))];
    end

    max_dist = max(max_dist);
    step_size = Alexa594{1,i}{1,1}(2,1) - Alexa594{1,i}{1,1}(1,1);
    dend_length = [step_size:step_size:max_dist]';
    Alexa594_radii{1,i} = [];
    for k = 1:length(dend_length)
        datapts = [];
        for q = 1:size(Alexa594_relRad{1,i},2)
            index = intersect(find(Alexa594_relRad{1,i}{1,q}(:,1) >= (dend_length(k) - 0.08)), find(Alexa594_relRad{1,i}{1,q}(:,1) <= (dend_length(k) + 0.08)));
            datapts = [datapts, Alexa594_relRad{1,i}{1,q}(index,2)];
        end
        Alexa594_radii{1,i}(k,1:4) = [dend_length(k), mean(datapts), std(datapts), (std(datapts)/sqrt(length(datapts)))];

    end
       
    try
        [f2, g2, o2] = fit(dend_length, [Alexa594_radii{1,i}(:,2)], 'exp2');
        [f2_e1, g2_e1, o2_e1] = fit(dend_length, [Alexa594_radii{1,i}(:,2)], 'exp1');
        if (g2_e1.rsquare > g2.rsquare) || (g2_e1.adjrsquare > g2.adjrsquare)
            clearvars f2 g2 o2
            f2 = f2_e1;
            g2 = g2_e1;
            o2 = o2_e1;
        end
    catch

        try 
            [f2, g2, o2] = fit(dend_length, [Alexa594_radii{1,i}(:,2)], 'exp2');
        catch
            [f2, g2, o2] = fit(dend_length, [Alexa594_radii{1,i}(:,2)], 'exp1');
        end
    end  
    subplot(4, 5, i)
    x_trace = [0:0.1:145]';
    Alexa594_radfit(:,i) = f2(x_trace);
    for r = 1:size(Alexa594_relRad{1,i},2)
        plot(Alexa594_relRad{1,i}{1,r}(:,1), Alexa594_relRad{1,i}{1,r}(:,2), 'r')
        hold on
    end
    plot(Alexa594_radii{1,i}(:,1), Alexa594_radii{1,i}(:,2), 'k', 'LineWidth', 2);
    plot(dend_length, f2(dend_length),  'b', 'LineWidth', 2);
    hold off
    title(['Cell ' num2str(i)])
    ylabel('Relative Radius')
    xlabel('dendritic distance (um)')
    sgtitle('Alexa 594')
    hold off        
end
sbplt(1, (i+1)) = subplot(4, 5, (i+1));
shadedErrorBar(x_trace, mean(Alexa594_radfit,2), [std(Alexa594_radfit,0,2)/sqrt(size(Alexa594_radfit,2))],'lineprops', {'-b','Linewidth',2})
title(['Summary'])
ylabel('Relative Radius')
xlabel('dendritic distance (um)')

clearvars -except Alexa594 Alexa488 ASAP3b Alexa594_summary Alexa594_radii Alexa594_relRad Alexa594_fit Alexa594_radfit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fullfig
for i = 1:size(Alexa488,2)
    max_dist = [];
    for j = 1:size(Alexa488{1,i},2)
        Alexa488{1,i}{1,j} = round(Alexa488{1,i}{1,j},4);
        max_dist(1,j) = max(Alexa488{1,i}{1,j}(:,1));
    end
    max_dist = max(max_dist);
    step_size = Alexa488{1,i}{1,1}(2,1) - Alexa488{1,i}{1,1}(1,1);
    dend_length = [step_size:step_size:max_dist]';
    Alexa488_summary{1,i} = [];
    for k = 1:length(dend_length)
        datapts = [];
        for q = 1:size(Alexa488{1,i},2)
            index = intersect(find(Alexa488{1,i}{1,q}(:,1) >= (dend_length(k) - 0.08)), find(Alexa488{1,i}{1,q}(:,1) <= (dend_length(k) + 0.08)));
            datapts = [datapts, Alexa488{1,i}{1,q}(index,2)];
        end
        Alexa488_summary{1,i}(k,1:4) = [dend_length(k), mean(datapts), std(datapts), (std(datapts)/sqrt(length(datapts)))];

    end
    
    try
        [f2, g2, o2] = fit(dend_length, [Alexa488_summary{1,i}(:,2)], 'exp2');
        [f2_e1, g2_e1, o2_e1] = fit(dend_length, [Alexa488_summary{1,i}(:,2)], 'exp1');
        if (g2_e1.rsquare > g2.rsquare) || (g2_e1.adjrsquare > g2.adjrsquare)
            clearvars f2 g2 o2
            f2 = f2_e1;
            g2 = g2_e1;
            o2 = o2_e1;
        end
    catch

        try 
            [f2, g2, o2] = fit(dend_length, [Alexa488_summary{1,i}(:,2)], 'exp2');
        catch
            [f2, g2, o2] = fit(dend_length, [Alexa488_summary{1,i}(:,2)], 'exp1');
        end
    end         
    
    subplot(4, 5, i)
    for r = 1:size(Alexa488{1,i},2)
        plot(Alexa488{1,i}{1,r}(:,1), Alexa488{1,i}{1,r}(:,2), 'g')
        hold on
    end
    x_trace = [0:0.1:145]';
    Alexa488_fit(:,i) = f2(x_trace);
    plot(Alexa488_summary{1,i}(:,1), Alexa488_summary{1,i}(:,2), 'k', 'LineWidth', 2);
    plot(dend_length, f2(dend_length),  'b', 'LineWidth', 2);
    hold off
    title(['Cell ' num2str(i)])
    ylabel('a.u.')
    xlabel('dendritic distance (um)')
    sgtitle('Alexa 488')
    hold off        
end
sbplt(1, (i+1)) = subplot(4, 5, (i+1));
shadedErrorBar(x_trace, mean(Alexa488_fit,2), [std(Alexa488_fit,0,2)/sqrt(size(Alexa488_fit,2))],'lineprops', {'-b','Linewidth',2})
title(['Summary'])
ylabel('a.u.')
xlabel('dendritic distance (um)')


fullfig
for i = 1:size(Alexa488,2)
    for j = 1:size(Alexa488{1,i},2)
        max_dist(1,j) = max(Alexa488{1,i}{1,j}(:,1));
        Alexa488_relRad{1,i}{1,j} = [Alexa488{1,i}{1,j}(:,1) sqrt((Alexa488{1,i}{1,j}(:,2)/mean(Alexa488_summary{1,i}(1:3,2))))];
    end

    max_dist = max(max_dist);
    step_size = Alexa488{1,i}{1,1}(2,1) - Alexa488{1,i}{1,1}(1,1);
    dend_length = [step_size:step_size:max_dist]';
    Alexa488_radii{1,i} = [];
    for k = 1:length(dend_length)
        datapts = [];
        for q = 1:size(Alexa488_relRad{1,i},2)
            index = intersect(find(Alexa488_relRad{1,i}{1,q}(:,1) >= (dend_length(k) - 0.08)), find(Alexa488_relRad{1,i}{1,q}(:,1) <= (dend_length(k) + 0.08)));
            datapts = [datapts, Alexa488_relRad{1,i}{1,q}(index,2)];
        end
        Alexa488_radii{1,i}(k,1:4) = [dend_length(k), mean(datapts), std(datapts), (std(datapts)/sqrt(length(datapts)))];

    end
    
     try
        [f2, g2, o2] = fit(dend_length, [Alexa488_radii{1,i}(:,2)], 'exp2');
        [f2_e1, g2_e1, o2_e1] = fit(dend_length, [Alexa488_radii{1,i}(:,2)], 'exp1');
        if (g2_e1.rsquare > g2.rsquare) || (g2_e1.adjrsquare > g2.adjrsquare)
            clearvars f2 g2 o2
            f2 = f2_e1;
            g2 = g2_e1;
            o2 = o2_e1;
        end
    catch

        try 
            [f2, g2, o2] = fit(dend_length, [Alexa488_radii{1,i}(:,2)], 'exp2');
        catch
            [f2, g2, o2] = fit(dend_length, [Alexa488_radii{1,i}(:,2)], 'exp1');
        end
    end  
             
    subplot(4, 5, i)
    x_trace = [0:0.1:145]';
    Alexa488_radfit(:,i) = f2(x_trace);
    for r = 1:size(Alexa488_relRad{1,i},2)
        plot(Alexa488_relRad{1,i}{1,r}(:,1), Alexa488_relRad{1,i}{1,r}(:,2), 'g')
        hold on
    end
    plot(Alexa488_radii{1,i}(:,1), Alexa488_radii{1,i}(:,2), 'k', 'LineWidth', 2);
    plot(dend_length, f2(dend_length),  'b', 'LineWidth', 2);
    hold off
    title(['Cell ' num2str(i)])
    ylabel('Relative Radius')
    xlabel('dendritic distance (um)')
    sgtitle('Alexa 488')
    hold off        
end

sbplt(1, (i+1)) = subplot(4, 5, (i+1));
shadedErrorBar(x_trace, mean(Alexa488_radfit,2), [std(Alexa488_radfit,0,2)/sqrt(size(Alexa488_radfit,2))],'lineprops', {'-b','Linewidth',2})
title(['Summary'])
ylabel('Relative Radius')
xlabel('dendritic distance (um)')

clearvars -except Alexa594 Alexa488 ASAP3b Alexa594_summary Alexa594_radii Alexa594_relRad Alexa488_summary Alexa488_radii Alexa488_relRad

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fullfig
for i = 1:size(ASAP3b,2)
    max_dist = [];
    for j = 1:size(ASAP3b{1,i},2)
        ASAP3b{1,i}{1,j} = round(ASAP3b{1,i}{1,j},4);
        max_dist(1,j) = max(ASAP3b{1,i}{1,j}(:,1));
    end
    max_dist = max(max_dist);
    step_size = ASAP3b{1,i}{1,1}(2,1) - ASAP3b{1,i}{1,1}(1,1);
    dend_length = [step_size:step_size:max_dist]';
    ASAP3b_summary{1,i} = [];
    for k = 1:length(dend_length)
        datapts = [];
        for q = 1:size(ASAP3b{1,i},2)
            index = intersect(find(ASAP3b{1,i}{1,q}(:,1) >= (dend_length(k) - 0.08)), find(ASAP3b{1,i}{1,q}(:,1) <= (dend_length(k) + 0.08)));
            datapts = [datapts, ASAP3b{1,i}{1,q}(index,2)];
        end
        ASAP3b_summary{1,i}(k,1:4) = [dend_length(k), mean(datapts), std(datapts), (std(datapts)/sqrt(length(datapts)))];

    end
    
    try
        [f2, g2, o2] = fit(dend_length, [ASAP3b_summary{1,i}(:,2)], 'exp2');
        [f2_e1, g2_e1, o2_e1] = fit(dend_length, [ASAP3b_summary{1,i}(:,2)], 'exp1');
        if (g2_e1.rsquare > g2.rsquare) || (g2_e1.adjrsquare > g2.adjrsquare)
            clearvars f2 g2 o2
            f2 = f2_e1;
            g2 = g2_e1;
            o2 = o2_e1;
        end
    catch

        try 
            [f2, g2, o2] = fit(dend_length, [ASAP3b_summary{1,i}(:,2)], 'exp2');
        catch
            [f2, g2, o2] = fit(dend_length, [ASAP3b_summary{1,i}(:,2)], 'exp1');
        end
    end  
             
    subplot(4, 5, i)
    for r = 1:size(ASAP3b{1,i},2)
        plot(ASAP3b{1,i}{1,r}(:,1), ASAP3b{1,i}{1,r}(:,2), 'g')
        hold on
    end
    x_trace = [0:0.1:80]';
    ASAP3b_fit(:,i) = f2(x_trace);
    plot(ASAP3b_summary{1,i}(:,1), ASAP3b_summary{1,i}(:,2), 'k', 'LineWidth', 2);
    plot(dend_length, f2(dend_length),  'b', 'LineWidth', 2);
    hold off
    title(['Cell ' num2str(i)])
    ylabel('a.u.')
    xlabel('dendritic distance (um)')
    sgtitle('ASAP3b')
    hold off        
end

sbplt(1, (i+1)) = subplot(4, 5, (i+1));
shadedErrorBar(x_trace, mean(ASAP3b_fit,2), [std(ASAP3b_fit,0,2)/sqrt(size(ASAP3b_fit,2))],'lineprops', {'-b','Linewidth',2})
title(['Summary'])
ylabel('a.u.')
xlabel('dendritic distance (um)')

fullfig
for i = 1:size(ASAP3b,2)
    for j = 1:size(ASAP3b{1,i},2)
        max_dist(1,j) = max(ASAP3b{1,i}{1,j}(:,1));
        ASAP3b_relRad{1,i}{1,j} = [ASAP3b{1,i}{1,j}(:,1) (ASAP3b{1,i}{1,j}(:,2)/mean(ASAP3b_summary{1,i}(1:3,2)))];
    end

    max_dist = max(max_dist);
    step_size = ASAP3b{1,i}{1,1}(2,1) - ASAP3b{1,i}{1,1}(1,1);
    dend_length = [step_size:step_size:max_dist]';
    ASAP3b_radii{1,i} = [];
    for k = 1:length(dend_length)
        datapts = [];
        for q = 1:size(ASAP3b_relRad{1,i},2)
            index = intersect(find(ASAP3b_relRad{1,i}{1,q}(:,1) >= (dend_length(k) - 0.08)), find(ASAP3b_relRad{1,i}{1,q}(:,1) <= (dend_length(k) + 0.08)));
            datapts = [datapts, ASAP3b_relRad{1,i}{1,q}(index,2)];
        end
        ASAP3b_radii{1,i}(k,1:4) = [dend_length(k), mean(datapts), std(datapts), (std(datapts)/sqrt(length(datapts)))];

    end
    
    try
        [f2, g2, o2] = fit(dend_length, [ASAP3b_radii{1,i}(:,2)], 'exp2');
        [f2_e1, g2_e1, o2_e1] = fit(dend_length, [ASAP3b_radii{1,i}(:,2)], 'exp1');
        if (g2_e1.rsquare > g2.rsquare) || (g2_e1.adjrsquare > g2.adjrsquare)
            clearvars f2 g2 o2
            f2 = f2_e1;
            g2 = g2_e1;
            o2 = o2_e1;
        end
    catch

        try 
            [f2, g2, o2] = fit(dend_length, [ASAP3b_radii{1,i}(:,2)], 'exp2');
        catch
            [f2, g2, o2] = fit(dend_length, [ASAP3b_radii{1,i}(:,2)], 'exp1');
        end
    end  
             
    subplot(4, 5, i)
    x_trace = [0:0.1:80]';
    ASAP3b_radfit(:,i) = f2(x_trace);
    for r = 1:size(ASAP3b_relRad{1,i},2)
        plot(ASAP3b_relRad{1,i}{1,r}(:,1), ASAP3b_relRad{1,i}{1,r}(:,2), 'g')
        hold on
    end
    plot(ASAP3b_radii{1,i}(:,1), ASAP3b_radii{1,i}(:,2), 'k', 'LineWidth', 2);
    plot(dend_length, f2(dend_length),  'b', 'LineWidth', 2);
    hold off
    title(['Cell ' num2str(i)])
    ylabel('Relative Radius')
    xlabel('dendritic distance (um)')
    sgtitle('ASAP3b')
    hold off        
end
sbplt(1, (i+1)) = subplot(4, 5, (i+1));
shadedErrorBar(x_trace, mean(ASAP3b_radfit,2), [std(ASAP3b_radfit,0,2)/sqrt(size(ASAP3b_radfit,2))],'lineprops', {'-b','Linewidth',2})
title(['Summary'])
ylabel('Relative Radius')
xlabel('dendritic distance (um)')
clearvars -except Alexa594 Alexa488 ASAP3b Alexa594_summary Alexa594_radii Alexa594_relRad Alexa488_summary Alexa488_radii Alexa488_relRad ASAP3b_summary ASAP3b_radii ASAP3b_relRad
