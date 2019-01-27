% demonstrate the character of an individual covariance function
% zhouweiyan 20180917
clear
% close all
clc
%% initialization
D=2;
seed=0;
rand('state',seed);
randn('state',seed);
%% latent individual covariance
cov={'covSEard'};hyp=log([5;6;5]);
covfunc=cov;
%% visualization
set(0,'DefaultFigureWindowStyle','docked') ;
% 1) query th enumber of parameters
feval(cov{:})
% 2) evaluate the function on x, x and xs to get cross-term
% [K,dK]=feval(cov{:},hyp,x)  % K: n¡Án; kss: ns¡Á1; Ks: n¡Áns
% [kss,dkss]=feval(cov{:},hyp,xs,'diag')
% [Ks,dKs]=feval(cov{:},hyp,x,xs)

% 3) plot a draw from the kernel
x1_test = linspace(0,10,128)';
x2_test = linspace(-5,5,128)';

[X1_test,X2_test]=meshgrid(x1_test,x2_test);
xstar=[X1_test(:) X2_test(:)];
K1=feval(cov{:},hyp,xstar);
K1=K1+(1e-5)*eye(size(K1));
n_samples=1;
samples=mvnrnd(zeros(size(X1_test(:))),K1,n_samples)';
figure
surf(X1_test,X2_test,reshape(samples,length(x2_test),length(x1_test)))
colormap(jet)

K0=feval(cov{:},hyp,xstar,[0 0]);
figure
surf(X1_test,X2_test,reshape(K0,length(x2_test),length(x1_test)),'EdgeColor','none',...
    'LineStyle','none','FaceLighting','phong');
xlabel('x1');ylabel('x2');
colormap(jet)


