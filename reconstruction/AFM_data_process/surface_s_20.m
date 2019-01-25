% generate surface from data points
%
% process the experiment data; AFM grid calibration; 
% spiral scanning pattern; scanning rate: 20Hz; sampling rate: 20000Hz
% create the surface with delaunayTriangulation interpolation
% zwy 20190125

clear
% close all
clc
%%
load('s20_processed_datapoints.mat');
x1_test = linspace(min(x),max(x),128);
x2_test = linspace(min(y),max(y),128);
[X1_test,X2_test] = meshgrid(x1_test,x2_test);

% Interpolation
m = griddata(x,y,z,X1_test(:),X2_test(:));

% % Exterpolation
% F = scatteredInterpolant(x,y,z);
% m = F(X1_test(:),X2_test(:));

figure
surf(X1_test,X2_test,reshape(m,length(x2_test),length(x1_test)),'EdgeColor','none')
colormap(copper)
axis equal; xlabel('x'); ylabel('y'); view([0 90])
% figure, scatter3(X1_test(:), X2_test(:), m*100,5,m*100,'.')   % verify

% Why it is circle but not rectangular? Because m includes NaN.
% GRIDDATA is an interpolation algorithm.
% scatter3(X1_test(:), X2_test(:), zeros(length(X1_test(:)),1),5,m*100,'.') %circle
% scatter3(X1_test(:), X2_test(:), zeros(length(X1_test(:)),1),5,'.') % rectangular

%% delete NaN
nan_flags = isnan(m);
m(nan_flags) = []; X1_test_vec = X1_test(:); X2_test_vec = X2_test(:);
X1_test_vec(nan_flags) = []; X2_test_vec(nan_flags) = [];
figure, scatter3(X1_test_vec,X2_test_vec,m,5,m*100,'.');








