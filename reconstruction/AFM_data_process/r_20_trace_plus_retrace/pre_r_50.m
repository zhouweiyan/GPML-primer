% recover 3D data from 1D original data
% raster scanning frequency: 20Hz
% zwy 20190123

clear
close all
clc

addpath('./20190121_tyd');
load r50.mat
z_origin = r50.Y(1).Data; z_origin = z_origin';

% f is the scanning frequency, including trace and retrace.
f=50;
pixels=128; 
% Total time of generating an image of pixels x pixels. 
t_total=pixels*1/f;
% extract the data corresponding to the scanning process. 
% 20000 is the sampling frequency.
dot_sum = ceil(20000*t_total);
ind = (1:dot_sum)';
% z_roi, the data set of raster scanning period.
% z_roi is chosen manually after the plot of (ind, z_origin).
z_roi_0 = z_origin(77741:(77740+dot_sum)); 
% cftool
%% period 1: coarse
fit1 = fit(ind, z_roi_0, 'poly4', 'Normalize', 'on');
% z_roi_1: the result of eliminating the whole tendency of z_roi
z_roi_1 = z_roi_0 - fit1(ind);
figure, plot(ind, z_roi_1)
%% period 2: fine
fit2 = fit(ind, z_roi_1, 'fourier4');
z_roi_2 = z_roi_1 - fit2(ind);
z_roi = z_roi_2;
figure, plot(ind, z_roi)
% outliers = exclude(ind, z_roi_2, 'indices', find(z_roi_2>0.00347))

x_origin = r20.Y(4).Data; x_origin = x_origin';
y_origin = r20.Y(3).Data; y_origin = y_origin';
x_roi = x_origin(83801:(83800+dot_sum));
y_roi = y_origin(83801:(83800+dot_sum));
plot3(x_roi, y_roi, z_roi*100,'.')





