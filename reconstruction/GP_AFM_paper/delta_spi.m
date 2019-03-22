clear
close all
% clc
load('array4_GP.mat')
opt='GP';
% opt='DT';
set(0,'DefaultFigureWindowStyle','docked')  % keep figures tidy
%% Constant Angular Velocity(CAV)
r_end=64*1.6;
% P=r_end*2/(60-1);
P=r_end*2/(12-1);
%P=spiral radius¡Á2/(number of curves-1)
%eq.(4)
w=120;
f=2000;
T=1/f;
T_spiral=2*pi*r_end/(P*w);	%CAV eq.(6)
t=0:T:T_spiral;
xs=P*w/(2*pi)*t.*cos(w*t)+64;	%eq.(17)
ys=P*w/(2*pi)*t.*sin(w*t)+64;	%eq.(18)
figure;
plot(xs,ys,'.')
axis equal
% axis ([0 128 0 128])
grid on
hold on;plot(xs,ys)

%% data preparation
% train set
x1_train_spi=xs';x2_train_spi=ys';
X_train_spi=[x1_train_spi(:),x2_train_spi(:)];
[y_train_spi,y_t_s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_train_spi);

% approximate path
path_len=sum(sqrt((x1_train_spi(2:end)-x1_train_spi(1:(end-1))).^2+(x2_train_spi(2:end)-x2_train_spi(1:(end-1))).^2));
delta = path_len/(2*128*128)