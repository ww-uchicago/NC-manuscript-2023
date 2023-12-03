function [ pdata] = clip( data, locs, wsize )
% this function would truncate the data according to the locations with
% ceterin window size and return a cell array
if size(data,2)>size(data,1)
    data=data';
end 
pdata=cell(1,size(data,2));
for i=1:size(data,2)
    temp=zeros(wsize,length(locs));
    for j=1:length(locs)
        temp(:,j)=data(locs(j):locs(j)+wsize-1,i);
    end 
    pdata{i}=temp;
end 
