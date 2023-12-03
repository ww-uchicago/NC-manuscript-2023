function [ matrix ] = DataSmooth( matrix,sbins,mode)
if size(matrix,2)>size(matrix,1)
    matrix=matrix';
end
for i=1:size(matrix,2)
    switch mode
        case 'average'
        matrix(:,i)=smooth(matrix(:,i),sbins);
        case 'median'
            matrix(:,i)=medfilt1(matrix(:,i),sbins);
    end
end
end
