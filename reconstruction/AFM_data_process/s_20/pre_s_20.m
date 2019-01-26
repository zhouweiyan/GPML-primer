% recover 3D data from 1D original data
% spiral scanning frequency: 20Hz
% zwy 20190124
%
% I. A. Mahmood and S. O. Reza Moheimani, ¡°Fast spiral-scan atomic force microscopy,¡± 
% Nanotechnology, vol. 20, no. 36, p. 365503, Sep. 2009.

clear
close all
% clc
%% surface fitting
load s20_plane_below.txt    
% s20_plane_below.txt is processed with polyworks. Only the points that
% belong to the below plane is saved.
x = s20_plane_below(:,1);
y = s20_plane_below(:,2);
z = s20_plane_below(:,3);
ft = fittype(@(p00,p10,p01,x,y)p00 + p10*x + p01*y, 'independent', {'x', 'y'}, 'dependent', 'z');
f = fit([x, y], z, ft, 'Normalize','on');
% verify the curve fitting results manually
figure; plot(f, [x,y],z)
%% eliminate the whole tendency
clearvars -except f
addpath('E:/OneDrive - hnu.edu.cn/tools/matlabcourse/GPML_matlab/GPML-primer/reconstruction/AFM_data_process/20190121_tyd');
load s20.mat
z_origin = s20.Y(1).Data; z_origin = z_origin';
x_origin = s20.Y(4).Data; x_origin = x_origin';
y_origin = s20.Y(3).Data; y_origin = y_origin';

r_end=10;                   % r_end is the radius, i.e., half of the range
pixels=128; 
P=r_end*2/pixels;           % eq.(3)
fit1=20; w=2*pi*fit1;       % f is the scanning frequency.
t_total=2*pi*r_end/(P*w);   % eq.(5) calculate t_total
% extract the data corresponding to the scanning process
dot_sum = ceil(20000*t_total); 
% _roi is chosen manually after the plot of (ind, z_origin).
x = x_origin(100651:(100650+dot_sum));
y = y_origin(100651:(100650+dot_sum));
z = z_origin(100651:(100650+dot_sum)); 
z = z - f(x, y);
scatter3(x, y, z*100,5,z*100,'.')
clearvars -except x y z











