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
% cov={'covSEard'};hyp=log([1;2]);
mask = [0,1];
cgi = {'covMaternard',3};
cma = {'covMask',{mask,cgi}}; hypma = log([1.5;1]);%log([ell;sf]);
covfunc={'covSum',{cma,{'covGE','ard',[]}}};
gamma=1.5;
hyp.cov=[hypma;log([1;1;(gamma/(2-gamma))])];%best now
%% visualization
% set(0,'DefaultFigureWindowStyle','docked') ;
% % 1) query th enumber of parameters
% feval(cov{:})
% 
% % 3) plot a draw from the kernel
% n_xstar=71; 
% xrange=linspace(-5,5,n_xstar)';
% if D~=1
%     [a,b]=meshgrid(xrange);
%     xstar=[a(:) b(:)];
%     K1=feval(cov{:},hyp,xstar);
%     K1=K1+(1e-5)*eye(size(K1));
%     n_samples=1;
%     samples=mvnrnd(zeros(size(a(:))),K1,n_samples)';
%     figure
%     surf(a,b,reshape(samples,n_xstar,n_xstar))
%     colormap(jet)
%     
%     K0=feval(cov{:},hyp,xstar,[0 0]);
%     figure
%     surf(a,b,reshape(K0,n_xstar,n_xstar),'EdgeColor','none',...
%         'LineStyle','none','FaceLighting','phong');
%     colormap(jet)
% else
%     K1=feval(cov{:},hyp,xrange);
%     K1=K1+(1e-5)*eye(size(K1));
%     n_samples=1;
%     samples=mvnrnd(zeros(n_xstar,1),K1,n_samples)';
%     figure
%     plot(xrange,samples,'LineWidth',2);
%     K0=feval(cov{:},hyp,xrange,0);
%     figure
%     plot(xrange,K0,'LineWidth',2);
% end
%% visualization
set(0,'DefaultFigureWindowStyle','docked') ;
% 1) query th enumber of parameters
feval(cov{:})
% 2) evaluate the function on x, x and xs to get cross-term
% [K,dK]=feval(cov{:},hyp,x)  % K: n¡Án; kss: ns¡Á1; Ks: n¡Áns
% [kss,dkss]=feval(cov{:},hyp,xs,'diag')
% [Ks,dKs]=feval(cov{:},hyp,x,xs)

% 3) plot a draw from the kernel
n_xstar=64;
xrange=(1:2:128)';
if D~=1
    [a,b]=meshgrid(xrange);
    xstar=[a(:) b(:)];
    K1=feval(cov{:},hyp,xstar);
    K1=K1+(1e-5)*eye(size(K1));
    n_samples=1;
    samples=mvnrnd(zeros(size(a(:))),K1,n_samples)';
    figure
    surf(a,b,reshape(samples,n_xstar,n_xstar))
    colormap(jet)
    
    K0=feval(cov{:},hyp,xstar,[0 0]);
    figure
    surf(a,b,reshape(K0,n_xstar,n_xstar),'EdgeColor','none',...
        'LineStyle','none','FaceLighting','phong');
    xlabel('x1');ylabel('x2');
    colormap(jet)
else
    K1=feval(cov{:},hyp,xrange);
    K1=K1+(1e-5)*eye(size(K1));
    n_samples=1;
    samples=mvnrnd(zeros(n_xstar,1),K1,n_samples)';
    figure
    plot(xrange,samples);
    K0=feval(cov{:},hyp,xrange,0);
    figure
    plot(xrange,K0);
end