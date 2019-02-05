% generate surface from data points
%
% process the experiment data; AFM grid calibration; 
% spiral scanning pattern; scanning rate: 20Hz; sampling rate: 20000Hz
% create the surface with delaunayTriangulation interpolation
% zwy 20190125

clear
close all
clc
%%
load('s20_processed_datapoints.mat');
x1_test = linspace(min(x),max(x),128);
x2_test = linspace(min(y),max(y),128);
[X1_test,X2_test] = meshgrid(x1_test,x2_test);

% Interpolation
m = griddata(x,y,z,X1_test(:),X2_test(:));

% % Exterpolation
% F = scatteredInterpolant(x,y,z); m = F(X1_test(:),X2_test(:));

figure
surf(X1_test,X2_test,reshape(m,length(x2_test),length(x1_test)),'EdgeColor','none')
colormap(copper)
axis equal; xlabel('x'); ylabel('y'); view([0 -90])
% figure, scatter3(X1_test(:), X2_test(:), m*100,5,m*100,'.')   % verify

% Why it is circle but not rectangular? Because m includes NaN.
% GRIDDATA is an interpolation algorithm.
% scatter3(X1_test(:), X2_test(:), zeros(length(X1_test(:)),1),5,m*100,'.') %circle
% scatter3(X1_test(:), X2_test(:), zeros(length(X1_test(:)),1),5,'.') % rectangular

%% replace NaN with right value manually
nan_flags = isnan(m);   % 2D logical matrix, work as a circular mask 
X1_test_vec = X1_test(:); X2_test_vec = X2_test(:);
for i = 1:length(nan_flags)
   flag = nan_flags(i);
   if flag
      if  X1_test_vec(i)<0||(X1_test_vec(i)>0.02&&X1_test_vec(i)<0.04)||...
              (X1_test_vec(i)>0.06&&X1_test_vec(i)<0.08)||...
              (X1_test_vec(i)>0.10&&X1_test_vec(i)<0.12)
          m(i) = 0.216/100*rand;
      else
          m(i) = 0.1/100*(rand-1);
      end
   end
end

figure, %scatter3(X1_test_vec,X2_test_vec,m,5,m*100,'.');
scatter3(X1_test(:), X2_test(:), m,5,m*100,'.') % rectangular

%% denoise
% waveletSignalDenoiser
% Method: Minimax; Parameters: level 12; Rule: Soft
load s20_m_denoise.mat
m = reshape(m_denoise,length(x2_test),length(x1_test));
m(nan_flags) = NaN; % circular roi
figure,surf(X1_test./0.04*3,X2_test./0.04*3,reshape(m.*10^4,length(x2_test./0.04*3),length(x1_test./0.04*3)),'EdgeColor','none')
colormap(copper)
axis equal; axis('off'); xlabel('x'); ylabel('y'); view([0 -90])
caxis([0 21.6]);
tightfig
set_fig_units_cm(10,10)
%% delete NaN, scale: from signal to reality
%       signal<--> reality
% x/y:  0.04  <--> 3um
% z:    10^-4 <--> 1nm
m = m(:)*10^4;
nan_flags = isnan(m);   % 1D NaN detector
m(nan_flags) = []; X1_test_vec(nan_flags) = []; X2_test_vec(nan_flags) = [];
x = X1_test_vec/0.04*3; y = X2_test_vec/0.04*3;
save s20_processed_delete_NaN.mat x y m






