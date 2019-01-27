% reconstruct s20 data points after delete outliers
%
% zwy 20190126, night

clear
% close all
clc

addpath('E:/OneDrive - hnu.edu.cn/tools/matlabcourse/GPML_matlab/GPML-primer/reconstruction/AFM_data_process/s_20');
load s20_processed_datapoints.mat
%          signal <--> reality
% x/y:      0.04  <--> 3um
% z:        10^-3 <--> 1nm
x = x/0.04*3; y = y/0.04*3;
z = z*10e3;
figure,scatter3(x,y,z,5,z,'.')

%% data preparation
% train set
x1_train_spi = x(1:10:end); x2_train_spi = y(1:10:end);
X_train_spi = [x1_train_spi, x2_train_spi];
y_train_spi = z(1:10:end);

% test set
ox = mean([max(x),min(x)]);
oy = mean([max(y),min(y)]);
r_end = roundn(mean([max(x)-min(x),max(y)-min(y)])/2,-2);
pixels = 128;
P=r_end*2/pixels;
w = 20; f = 10000; T = 1/f;
T_spiral = 2*pi*r_end/(P*w);
t = 0:T:T_spiral;
xs = P*w/(2*pi)*t.*cos(w*t)+ox; xs = xs';
ys = P*w/(2*pi)*t.*sin(w*t)+oy; ys = ys';
hold on; plot(xs,ys,'*');xlabel('x1'); ylabel('x2');

X_test = [xs,ys];
meanfunc={@meanNN,[5+((-3*2):3:(3*2)) 5+2+(-3*3:3:(3*1));oy*ones(1,10)]',...
    [21.6*ones(1,5) zeros(1,5)]};
hyp.mean=[];
covfunc={'covNoise'};hyp.cov=log(0.01);
likfunc='likGauss';hyp.lik=log(10);    % log(sn)
hyp=minimize(hyp,@gp,-80,@infGaussLik,meanfunc,covfunc,likfunc,X_train_spi,y_train_spi);
[m,s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train_spi,y_train_spi,X_test);
figure,scatter3(xs,ys,m,5,m,'.');xlabel('x'); ylabel('y');
%% point cloud to surface
% prefix f_ represents face.
f_x1_test = linspace(min(xs),max(xs),128);
f_x2_test = linspace(min(ys),max(ys),128);
[f_X1_test,f_X2_test] = meshgrid(f_x1_test,f_x2_test);
% Interpolation
f_m = griddata(xs,ys,m,f_X1_test(:),f_X2_test(:));
figure
surf(f_X1_test,f_X2_test,reshape(f_m,length(f_x2_test),length(f_x1_test)),'EdgeColor','none')
colormap(copper)
axis equal; xlabel('x'); ylabel('y'); view([0 -90])






