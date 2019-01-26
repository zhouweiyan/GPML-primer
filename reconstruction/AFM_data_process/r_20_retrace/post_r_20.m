% correct/delete the outliers
%
% process the experiment data; AFM grid calibration; 
% rastter scanning pattern; scanning rate: 20Hz; sampling rate: 20000Hz
% r20_..._outliers are chosen manually: pre_r_20.m 
% in the figure: scatter3(x, y, z*100,5,z*100,'.')
%
% 'above' represents that the data should equal to the upper limit, i.e., 
% 0.126/100 in this example. 
% 'delete' represents that the data shpuld be deleted
% 
% choose half of the data: trace or retrace
%
% zwy 20190125 night, done

clear
close all
clc
load('r20_outliers_manually.mat')
len = length(z);    

% delete the delete_outliers
%
% 'for i = 1:len' is wrong because len is constant while the length(z) is changing.
% 'for i = 1:len' is still wrong if 'len = len - 1;' is written in the if nest. 
% It is due to the domain of validity. len is not affected by the inside len.
% 'for i = 1:length(z)' is wrong because the terminate of for-loop is unchangable.
r20_delete_outliers = [r20_delete_outliers1;r20_delete_outliers2];
for i = 1:len     
    for j = 1:length(r20_delete_outliers)
        if x(i)==r20_delete_outliers(j,1) && y(i)==r20_delete_outliers(j,2)
            x(i) = []; y(i) = []; z(i) = [];
        end
    end
    if i==length(z)
        break;
    end
end

scatter3(x, y, z*100,5,z*100,'.')

clearvars -except x y z
