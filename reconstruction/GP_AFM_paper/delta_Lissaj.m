clear
close all
% clc
load('array4_GP.mat')
set(0,'DefaultFigureWindowStyle','docked') 

%% Lissajous scan trajectory
X0=64;Y0=64;
A=128;B=128;    % A,B represent the range of X,Y respectively
% fx=41;fy=43;
fx=5;fy=4;
f=5000;
T=1;  %T = lcm(fx,fy)/(fx*fy)
t=(0:1/f:T)';
x1_train_lisaj=X0+A/2*sin(2*pi*fx*t);
x2_train_lisaj=Y0+B/2*sin(2*pi*fy*t);
figure;
plot(x1_train_lisaj,x2_train_lisaj,'.')
axis equal
axis ([0 128 0 128])
grid on
hold on;plot(x1_train_lisaj,x2_train_lisaj)

%% data preparation
% train set
X_train_lisaj=[x1_train_lisaj(:),x2_train_lisaj(:)];
[y_train_lisaj,y_t_s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_train_lisaj);

% approximate path
path_len=sum(sqrt((x1_train_lisaj(2:end)-x1_train_lisaj(1:(end-1))).^2+(x2_train_lisaj(2:end)-x2_train_lisaj(1:(end-1))).^2));
delta = path_len/(2*128*128)