% Make figures of additive kernels in 3 dimensions
%
% D. K. Duvenaud, ¡°Automatic Model Construction  with Gaussian Processes,¡± p. 157, 2014.
%
% zhouweiyan 20180724

clear
close all
clc
addpath('utils/');

%% generate a grid
range=-4:0.1:4;
[a,b]=meshgrid(range);
xstar=[a(:) b(:)];

%% plot a SE kernel
covfunc='covSEiso';
hyp.cov=log([1;1]);
K=feval(covfunc(:),hyp.cov,xstar,[0 0]);
figure(1);
h=surf(a,b,reshape(K,length(range),length(range)),'EdgeColor','none',...
    'LineStyle','none','FaceLighting','phong');
colormap(jet);
title('covfunc=''covSEiso''');
% plot a draw from that prior
seed=0;
randn('state',seed);
rand('state',seed);
n=length(xstar);
K=feval(covfunc(:),hyp.cov,xstar,xstar);
K=K+eye(n)*0.000001;
y=chol(K,'lower')*gpml_randn(rand,n,1);
figure(2)
h=surf(a,b,reshape(y,length(range),length(range)),'EdgeColor','none',...
    'LineStyle','none','FaceLighting','phong');
colormap(jet);
title('covfunc=''covSEiso''');
%% plot a ADD kernel
covfunc={'covADD',{1,'covSEiso'}};
hyp.cov=log([1;1;1;1;sqrt(2)/2])
K=feval(covfunc{:},hyp.cov,xstar,[0 0]);
figure(3);
h=surf(a,b,reshape(K,length(range),length(range)),'EdgeColor','none',...
    'LineStyle','none','FaceLighting','phong');
colormap(jet)
title('covfunc={''covADD'',{1,''covSEiso''}}')
% plot a draw from that prior
seed=0;
randn('state',seed);
rand('state',seed);
n=length(xstar);
K=feval(covfunc{:},hyp.cov,xstar);
K=K+eye(n)*0.000001;
y=chol(K,'lower')*gpml_randn(rand,n,1);
figure(4)
h=surf(a,b,reshape(y,length(range),length(range)),'EdgeColor','none',...
    'LineStyle','none','FaceLighting','phong');
colormap(jet)
title('covfunc={''covADD'',{1,''covSEiso''}}')
%% plot an additive kernel with 2nd order interactions
covfunc={'covADD',{[1,2],'covSEiso'}};
hyp.cov=log(ones(6,1));
K=feval(covfunc{:},hyp.cov,xstar,[0,0]);
figure(5)
h=surf(a,b,reshape(K,length(range),length(range)),'EdgeColor','none',...
    'LineStyle','none','FaceLighting','phong');
colormap(jet);
title('covfunc={''covADD'',{[1,2],''covSEiso''}}')
% plot a draw from that prior
seed=0;
randn('state',seed);
rand('state',seed);
n=length(xstar);
K=feval(covfunc{:},hyp.cov,xstar);
K=K+eye(n)*0.000001;
y=chol(K,'lower')*gpml_randn(rand,n,1);
figure(6)
h=surf(a,b,reshape(y,length(range),length(range)),'EdgeColor','none',...
    'LineStyle','none','FaceLighting','phong');
colormap(jet)
title('covfunc={''covADD'',{[1,2],''covSEiso''}}')
%% plot an additive prior as the sum of two 1d kernels
covfunc={'covADD',{1,'covSEiso'}};
hyp.cov=log([1;1;sqrt(1/2)]);           % 1d
% hyp.cov = log([1,1,1,1,sqrt(1/2)]);   % 2d
K1=feval(covfunc{:},hyp.cov,a(:),0);
K2=feval(covfunc{:},hyp.cov,b(:),0);
figure(7)
h=surf(a,b,reshape(K1,length(range),length(range)),'EdgeColor','none',...
    'LineStyle','none','FaceLighting','phong');
colormap(jet)
title('additive kernel sum p1');
figure(8)
h=surf(a,b,reshape(K2,length(range),length(range)),'EdgeColor','none',...
    'LineStyle','none','FaceLighting','phong');
colormap(jet)
title('additive kernel sum p2');

figure(9)
h=surf(a,b,reshape(K1+K2,length(range),length(range)),'EdgeColor','none',...
    'LineStyle','none','FaceLighting','phong');
colormap(jet)
title('additive kernel');
% plot a draw from those priors
xstar2=[zeros(length(a(:)),1) b(:)];
xstar1=[a(:) zeros(length(b(:)),1)];
hyp.cov = log([1,1,1,1,sqrt(1/2)]);   % 2d
% seed=0;
seed=5;
randn('state',seed);
rand('state',seed);
n=length(xstar1);
K1=feval(covfunc{:},hyp.cov,xstar1);
K2=feval(covfunc{:},hyp.cov,xstar2);
K1=K1+eye(n)*.000001;
K2=K2+eye(n)*.000001;
y1=chol(K1,'lower')*gpml_randn(rand,n,1);
y2=chol(K2,'lower')*gpml_randn(rand,n,1);
figure(10)
h=surf(a,b,reshape(y1,length(range),length(range)),'EdgeColor','none',...
    'LineStyle','none','FaceLighting','phong');
colormap(jet)
title('additive kernel sum p1')
figure(11)
h=surf(a,b,reshape(y2,length(range),length(range)),'EdgeColor','none',...
    'LineStyle','none','FaceLighting','phong');
colormap(jet)
title('additive kernel sum p2');
figure(12)
h=surf(a,b,reshape(y1+y2,length(range),length(range)),'EdgeColor','none',...
    'LineStyle','none','FaceLighting','phong');
colormap(jet)
title('additive kernel');



