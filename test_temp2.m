clear
% close all
clc
%% initialization
n=5;D=2;

seed=0;
rand('state',seed);
randn('state',seed);
x=randn(n,D);
xs=randn(3,D);

meanfunc=[];
csa={'covSEiso'}; hypsa =log([20;40]);  %log([ell;sf])
csi={'covSEiso'};hypsi=log([50;1]);      %log([ell;sf]);
mask=[1,0];
% cov={'covMask',{mask,{'covProd',{csa,csi}}}};hyp=[hypsa;hypsi];
% cov={'covProd',{'covMask',{mask,csa,csi}}};hyp=[hypsa;hypsi];
cov={'covMask',{mask,{csa,csi}}};hyp=[hypsa;hypsi];%sth. error here

%% visualization
set(0,'DefaultFigureWindowStyle','docked') ;
% 1) query th enumber of parameters
feval(cov{:})
% 2) evaluate the function on x, x and xs to get cross-term
[K,dK]=feval(cov{:},hyp,x)  % K: n¡Án; kss: ns¡Á1; Ks: n¡Áns
[kss,dkss]=feval(cov{:},hyp,xs,'diag')
[Ks,dKs]=feval(cov{:},hyp,x,xs)

% 3) plot a draw from the kernel
n_xstar=71;
xrange=linspace(1,256,n_xstar)';
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
    colormap(jet)
else
    K1=feval(cov{:},hyp,xrange);
    K1=K1+(1e-5)*eye(size(K1));
    n_samples=1;
    samples=mvnrnd(zeros(n_xstar,1),K1,n_samples)';
    figure
    plot(xrange,samples);
end