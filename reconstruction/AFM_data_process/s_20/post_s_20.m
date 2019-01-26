% correct the outliers
%
% process the experiment data; AFM grid calibration; 
% spiral scanning pattern; scanning rate: 20Hz; sampling rate: 20000Hz
% s20_above/below_outliers are chosen manually: pre_s_20.m 
% in the figure: scatter3(x, y, z*100,5,z*100,'.')
% 'above' represents that the data should equal to the upper limit, i.e., 
% 0.126/100 in this example. 
% 'below' represents that the data should equal to the lower limit,i.e., 0.
% 
% zwy 20190125

clear
close all
clc
load('s20_outliers_manually.mat')
len = length(z);    

% correct the outliers1
for i = 1:len
    for j = 1:length(s20_below_outliers1)
        if x(i)==s20_below_outliers1(j,1) && y(i)==s20_below_outliers1(j,2)
            z(i) = 0;
        end
    end
end
% correct the outliers2
for i = 1:len
    for j = 1:length(s20_above_outliers2)
        if x(i)==s20_above_outliers2(j,1) && y(i)==s20_above_outliers2(j,2)
            z(i) = (0.12 + 0.1*rand)/100;
        end
    end
end
% correct the outliers3
for i = 1:len
    for j = 1:length(s20_above_outliers3)
        if x(i)==s20_above_outliers3(j,1) && y(i)==s20_above_outliers3(j,2)
            z(i) = 0.216/100;
        end
    end
end

scatter3(x, y, z*100,5,z*100,'.')

clearvars -except x y z
