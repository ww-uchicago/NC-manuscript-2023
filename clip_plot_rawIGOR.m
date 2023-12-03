clear
clc

parent = 'C:\Users\Hector\Box\Postdoc\';

date = '20201209';

folder = [parent date '\Experiment Folder'];

mat_files = dir([folder '\*' '.mat']);

for i = 1:size(mat_files,1)
    load([folder '\' date 'w' num2str((i-1)) '.mat']);
        
    start_id = strfind(D.WaveNotes, 'PATTERN');
    
    end_id = strfind(D.WaveNotes, ';H');
    
    pattern = D.WaveNotes((start_id + length('PATTERN:')):(end_id -1));
    

    
    clearvars D start_id end_id stim pattern
end
