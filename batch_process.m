function batch_process
% to process multiple .abf files located in one folder. Calls them out in
% an alphabetical order:
list = dir('*.abf');
file_names = {list.name};

for n=1:size(file_names,2)
    filename = file_names(n);
    disp([int2str(n) '. Working on: ' filename{:}])
    APanalysis(filename{:});
end