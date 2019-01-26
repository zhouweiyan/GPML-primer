% generate surface from data points
%
% process the experiment data; AFM grid calibration; 
% raster scanning pattern; scanning rate: 20Hz; sampling rate: 20000Hz
% create the surface with delaunayTriangulation interpolation
% zwy 20190125

%% Period 1: choose the reconstruction boundary without NaN
clear
close all
clc
load('r20_processed_datapoints.mat');

x1_test = linspace(min(x),max(x),128);
x2_test = linspace(min(y),max(y),128);
[X1_test,X2_test] = meshgrid(x1_test,x2_test);
m = griddata(x,y,z,X1_test(:),X2_test(:));
m = reshape(m,length(x2_test),length(x1_test));

% reserve region without NaN manually
m = m(9:120,3:114); 
% sum(sum(isnan(m)))    % 0 is right.
X1_test = X1_test(9:120,3:114); X2_test = X2_test(9:120,3:114);
x1_test = X1_test(1,:); x2_test = X2_test(:,1);
if 0    % visiualization
    figure,surf(X1_test,X2_test,m,'EdgeColor','none')
    colormap(copper)
    axis equal; xlabel('x'); ylabel('y'); view([0 90])
end
side_length = roundn(min(max(x1_test)-min(x1_test),max(x2_test)-min(x2_test)),-2);    
% max(x1_test)-min(x1_test) is small. 
% figure,scatter3(X1_test(:), X2_test(:), m*100,5,m*100,'.')   % verify

%% Period 2: Interpolation with the chosen boundary
clearvars -except side_length x1_test x2_test
% close all
clc
load('r20_processed_datapoints.mat');

x1_test = linspace(min(x1_test),max(x1_test),128);
x2_test = linspace(min(x2_test),min(x2_test)+side_length,128);
[X1_test,X2_test] = meshgrid(x1_test,x2_test);
m = griddata(x,y,z,X1_test(:),X2_test(:));
%% denoise
% waveletSignalDenoiser
% Parameters: level 14; Rule: Soft
load r20_m_denoise.mat
m = reshape(m_denoise,length(x2_test),length(x1_test));
figure,surf(X1_test,X2_test,m,'EdgeColor','none')
colormap(copper)
axis equal; xlabel('x'); ylabel('y'); view([0 90])
%% scale: from signal to reality
%       signal<--> reality
% x/y:  0.04  <--> 3um
% z:    10^-3 <--> 1nm
m = m(:)*10^3;
x = X1_test(:)/0.04*3; y = X2_test(:)/0.04*3;
save r20_processed_delete_NaN.mat x y m








