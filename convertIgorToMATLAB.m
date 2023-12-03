clear
clc

parent = 'C:\Users\Hector\Box\Postdoc\';

date = '20201209';

folder = [parent date '\Experiment Folder'];

ibw_fns = dir([folder '\w*' '.ibw']);

for i = 2:size(ibw_fns,1)
    D = IBWread([folder '\w' num2str((i-2)) '.ibw']);
        
    start_id = strfind(D.WaveNotes, 'PATTERN');
    
    end_id = strfind(D.WaveNotes, ';H');
    
    pattern = D.WaveNotes((start_id + length('PATTERN:')):(end_id -1));
    
    stim = IBWread([folder '\maf\' pattern '_DA0.ibw']);
    
    D.stim_DA0 = stim.y;
    
    if isempty(strfind(pattern,'membTest')) == 0
        data = filterData(D.y, (D.dx*1000000), 'low_cut',0.1, 'high_cut', 2000, 'order', 1);
        para.deltaV=10;  % deltaV=10mv, the default in the lab
        para.sampInt = D.dx;
        para.sr=1/para.sampInt;
        para.bsl=1:(0.090/para.sampInt);  % the baseline subtraction, taking the first 100ms of data as baseline.
        wsize = 0.45*para.sr;
        wsize = int64(wsize);
        locs = linspace(0,24,25);
        locs = locs*0.45*para.sr + 1;
        locs = int64(locs);
        data = clip(data,locs,wsize);
        data = data{1,1};
        D.clipped = data;
        
        clearvars data locs wsize para
        
    elseif isempty(strfind(pattern,'IVcurve')) == 0
        data = filterData(D.y, (D.dx*1000000), 'low_cut',0.1, 'high_cut', 2000, 'order', 1);
        para.sampInt = D.dx;
        para.sr=1/para.sampInt;
        wsize = 1.2*para.sr;
        wsize = int64(wsize);
        locs = linspace(0,7,8);
        locs = locs*1.2*para.sr + 1;
        locs = int64(locs);
        data = clip(data,locs,wsize);
        data = data{1,1};
        D.clipped = data;
        
        clearvars data locs wsize para

    end
    
    save([folder '\' date '-w' num2str((i-2)) '-' pattern  '.mat'],'D');
    
    clearvars D start_id end_id stim pattern 
end
