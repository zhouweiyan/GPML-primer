% demonstrate the character of an individual covariance function
% zhouweiyan 20180917
clear
% close all
clc
%% initialization
seed=0;
rand('state',seed);
randn('state',seed);
%% latent individual covariance
cov={'covSEard'};hyp=log([4;2;3]);
% covfunc={'covSum',{'covNoise',covfunc}};hyp.cov=[log(0.001);hyp.cov];
% covfunc={'covMask',{mask,{'covSEiso'}}};hyp.cov=log([2;40]);  %log([ell;sf])
%% visualization
set(0,'DefaultFigureWindowStyle','docked') ;
% 1) query th enumber of parameters
feval(cov{:})

% 3) plot a draw from the kernel
% n_xstar=71;
% xrange=linspace(-5,5,n_xstar)';
xrange=(1:4:256)';n_xstar=length(xrange);

[a,b]=meshgrid(xrange);
xstar=[a(:) b(:)];
K1=feval(cov{:},hyp,xstar);
K1=K1+(1e-5)*eye(size(K1));
n_samples=1;
samples=mvnrnd(zeros(size(a(:))),K1,n_samples)';
figure
surf(a,b,reshape(samples,n_xstar,n_xstar))
colormap(jet)

% K0=feval(cov{:},hyp,xstar,[0 0]);
% figure
% surf(a,b,reshape(K0,n_xstar,n_xstar),'EdgeColor','none',...
%     'LineStyle','none','FaceLighting','phong');
% colormap(jet)
