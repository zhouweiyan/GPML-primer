% run in different part, the first part is necessary
% zhouweiyan 20180918
clear 
% close all
clc
write_fig = 0; N = 30;
sd = 3; 
rand('seed',sd)
randn('seed',sd)       % set a seed for reproducability
liktyp = input('Which likelihood?\n  (g)aussian, (l)ogistic: ','s');

a = 0.3; b = 1.2; f = @(x) a*x + b + sin(x);               % underlying function
n = 30; sn = 0.5;          % number of training points, noise standard deviation
x = 2*rand(n,1)-1; x = 1+4*x+sign(x); y = f(x)+sn*randn(n,1);      % sample data

cov = {@covSEiso}; sf = 2; ell = 1.0; hyp.cov = log([ell;sf]);
mean = {@meanSum,{@meanLinear,@meanConst}}; hyp.mean = [a;b];
if isequal(liktyp,'g')
  lik = {@likGauss};    hyp.lik = log(sn); inf = @infGaussLik;
else
  lik = {@likLogistic}; hyp.lik = [];      inf = @infLaplace;  y = sign(y);
end
xs = linspace(-8,10,2e3)'; ys = f(xs);                   % exact function values

%% exact dense prediction
fprintf('Optimise hyperparameters.\n')
hyp = minimize(hyp,@gp,-N,inf,mean,cov,lik,x,y);      % optimise hyperparameters
[ymu,ys2] = gp(hyp,inf,mean,cov,lik,x,y,xs);                   
figure
f = [ymu+2*sqrt(ys2); flip(ymu-2*sqrt(ys2),1)];
fill([xs; flip(xs,1)], f, [7 7 7]/8)
hold on;plot(xs,ymu,'LineWidth',2);plot(x,y,'r+'), plot(xs,ys,'r')
xlim([-8,10]), ylim([-3,6])
%% FITC prediction
nu = 10; xu = linspace(-6,8,nu)'; covf = {@apxSparse,cov,xu};  
[ymuf,ys2f] = gp(hyp,inf,mean,covf,lik,x,y,xs); % FITC
figure
f = [ymuf+2*sqrt(ys2f); flip(ymuf-2*sqrt(ys2f),1)];
fill([xs; flip(xs,1)], f, [7 7 7]/8)
hold on;plot(xs,ymuf,'LineWidth',2);plot(x,y,'r+'), plot(xs,ys,'r')
xlim([-8,10]), ylim([-3,6])
%% grid prediction
ng = 40; xg = linspace(-6,8,ng)'; covg = {@apxGrid,{cov},{xg}};
opt.cg_maxit = 500; opt.cg_tol = 1e-5; opt.pred_var = 100;  % parameters
inf = @(varargin) infGrid(varargin{:},opt);
[ymug,ys2g] = gp(hyp,inf,mean,covg,lik,x,y,xs); % grid
figure
f = [ymug+2*sqrt(ys2g); flip(ymug-2*sqrt(ys2g),1)];
fill([xs; flip(xs,1)], f, [7 7 7]/8)
hold on;plot(xs,ymug,'LineWidth',2);plot(x,y,'r+'), plot(xs,ys,'r')
xlim([-8,10]), ylim([-3,6])
%% fast grid prediction
[nlZ,dnlZ] = gp(hyp,inf,mean,cov,lik,x,y)   % marginal likelihood and derivatives
ng = 40; xg = linspace(-6,8,ng)'; covg = {@apxGrid,{cov},{xg}};
opt.cg_maxit = 500; opt.cg_tol = 1e-5; opt.pred_var = 100;  % parameters
inf = @(varargin) infGrid(varargin{:},opt);
[postg,nlZg,dnlZg] = infGrid(hyp,mean,covg,lik,x,y,opt);  
[fmugf,fs2gf,ymugf,ys2gf] = postg.predict(xs);  % fast-grid
fprintf('Sampling estimators for nlZ and dnlZ\n')
opt.ndcovs = 25;       % ask for (additional) sampling-based (exact) derivatives
opt.ldB2_cheby = true;
opt.ldB2_cheby_hutch = 30;
[posts,nlZs,dnlZs] = infGrid(hyp,mean,covg,lik,x,y,opt);
[dnlZ.cov,dnlZg.cov,dnlZs.covs,dnlZs.cov]
[nlZ,nlZg,nlZs]
figure
f = [ymugf+2*sqrt(ys2gf); flip(ymugf-2*sqrt(ys2gf),1)];
fill([xs; flip(xs,1)], f, [7 7 7]/8)
hold on;plot(xs,ys2gf,'LineWidth',2);plot(x,y,'r+'), plot(xs,ys,'r')
xlim([-8,10]), ylim([-3,6])
%% fast grid prediction using LOVE
opt.pred_var = -10;                                 % ask for 10 Lanczos vectors
postg2 = infGrid(hyp,mean,covg,lik,x,y,opt);   % fast grid prediction using LOVE
[fmugf2,fs2gf2,ymugf2,ys2gf2] = postg2.predict(xs);
%% comparison
figure
subplot(211)
plot(xs,ymu,'k','LineWidth',2), hold on
plot(xs,ymuf,'g-.','LineWidth',2)
plot(xs,ymug,'m:','LineWidth',2)
plot(xs,ymugf,'b-.','LineWidth',2)
plot(xs,ymugf2,'c--','LineWidth',2)
legend('exact','FITC','grid','fast-grid','fast-grid/LOVE')
title('Predictive mean')
plot(x,y,'r+'), plot(xs,ys,'r')
plot(xs,ymu+2*sqrt(ys2),'k'), plot(xs,ymu-2*sqrt(ys2),'k')
xlim([-8,10]), ylim([-3,6])

subplot(212)
plot(xs,sqrt(ys2),'k','LineWidth',2), hold on
plot(xs,sqrt(ys2f),'g-.','LineWidth',2)
plot(xs,sqrt(ys2g),'m:','LineWidth',2)
plot(xs,sqrt(ys2gf),'b-.','LineWidth',2)
plot(xs,sqrt(ys2gf2),'c--','LineWidth',2)
legend('exact','FITC','grid','fast-grid','fast-grid/LOVE')
title('Predictive standard dev')
xlim([-8,10])
if write_fig, print -depsc f10.eps; end