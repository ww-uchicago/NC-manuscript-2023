folder = 'C:\Users\Hector Acaron\Desktop\Wei Lab\Data\20170505';
%pname=uigetdir(folder);
%read the xml, and tif file and calculate the mean intensity value over
%time
moviefiles = dir([folder '\*.zip']);
cell_area = 0;


for q = 1:size(moviefiles,1)
    close all
    clearvars -except paras folder moviefiles q ratioCh1 ratioCh2
    movieNum = str2num(moviefiles(q,1).name(7:10));
    moviestring = moviefiles(q,1).name(7:10);
    pname = dir([folder '\ZSeries*' moviestring]);
    tseriesxml = [pname.name '.xml'];
    pname = [folder '\' pname.name];

    %% IMPORT TIFF FILES%%%%
  


    [~, tseries] = fileparts(pname);
    xmlfile = fileread([pname '\' tseriesxml]);
    chImgNames =  regexp(xmlfile, ...
        'File channel="(\d*)".*?filename="(.*?)"', 'tokens');
    ch = cellfun(@(x) str2double(x{1}),chImgNames);
    Ch1ImageNames = chImgNames(ch == 1);
    Ch2ImageNames = chImgNames(ch == 2);
    imgNamesCh1 = cellfun(@(x) x{2},Ch1ImageNames, 'UniformOutput', 0);
    
    imgNamesCh2 = cellfun(@(x) x{2},Ch2ImageNames, 'UniformOutput', 0);
    % Clear xml text form memory, not needed for final img import
    clear xmlfile;

    %%TURN TIFF FILES INTO MATRICES%%%%
    % Read individual tif files  img
    for n = 1:numel(imgNamesCh1)
        I1=imread([pname '\' imgNamesCh1{n}]);
        I2=imread([pname '\' imgNamesCh2{n}]);
        movieCh1{n} = I1;
        movieCh2{n} = I2;
    end

    % [rI, cI] = size(I);
    % stdProjection = zeros(rI, cI);
    % 
    % for i = 1:rI
    %     for k = 1:cI
    %         array = [];
    %         for l = 1:n
    %             array = [array, double(movieM{l}(i, k))];
    %         end
    %         stdProjection(i, k) = std(array);
    %     end
    % end
    % 
    % projImage = mat2gray(stdProjection);
    % imshow(projImage);

    [rows cols] = size(movieCh1{1,1});


    %%%% IMPORT IMAGEJ ROIS%%%%
    [rois] = ReadImageJROI([folder '\RoiSet0' num2str(movieNum) '.zip']);

    cellSignalsCh1 = zeros(length(movieCh1), length(rois));
    cellSignalsCh2 = zeros(length(movieCh2), length(rois));
    coords = zeros(length(rois),4);

    for i = 1:length(rois)
        coords(i,:) = rois{1,i}.vnRectBounds;
        if coords(i, 3) > rows 
            coords(i,3) = rows;
        elseif coords(i,4) > cols
            coords(i,4) = cols;
        end
    end

    coords(coords < 1) = 1;
    coords(coords(:,3) > rows) = rows;
    coords(coords(:,4) > cols) = cols;

    %%%%%EXTRACT MEAN VALUE FOR EACH ROI IN EVERY FRAME%%%%%
    for i=1:length(movieCh1)
        for k= 1:length(rois)
            region1 = movieCh1{1,i}(coords(k,1):coords(k,3), coords(k,2): coords(k,4));
            region2 = movieCh2{1,i}(coords(k,1):coords(k,3), coords(k,2): coords(k,4));
            cellSignalsCh1(i,k) = mean(mean(region1));
            cellSignalsCh2(i,k) = mean(mean(region2));
        end
    end
    cellSignalsCh1 = max(cellSignalsCh1);
    cellSignalsCh2 = max(cellSignalsCh2);
    
    cellSignalsCh1 = cellSignalsCh1/cellSignalsCh1(end);
    cellSignalsCh2 = cellSignalsCh2/cellSignalsCh2(end);
    
    ratioCh1{1,q} = cellSignalsCh1(1:end-1)';
    ratioCh2{1,q} = cellSignalsCh2(1:end-1)';
    
end

figure
set(gcf,'color','white')
for q = 1:size(moviefiles,1)
    colr = colordg(q);
    scatter(ratioCh2{1,q},ratioCh1{1,q},[], colr, 'filled')
    hold on
end

xlabel('f/f0 GCamp6','FontSize',12,'FontWeight','bold')
ylabel('f/f0 TdTomato','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'fontWeight','bold')