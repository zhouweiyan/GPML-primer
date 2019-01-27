clear
close all
clc

meanfunc0 = {@meanSum, {@meanLinear, @meanConst}};
hyp0.mean = [0.5; 1];
covfunc0 = {@covMaterniso, 3};
ell = 1/4; sf = 1; hyp0.cov = log([ell; sf]);
likfunc0 = @likGauss;
sn = 0.1;
hyp0.lik = log(sn);

n = 20;
x = gpml_randn(0.3, n, 1);
K = feval(covfunc0{:}, hyp0.cov, x);
mu = feval(meanfunc0{:}, hyp0.mean, x);
y = chol(K)'*gpml_randn(0.15, n, 1) + mu + exp(hyp0.lik)*gpml_randn(0.2,n,1);
plot(x, y, '+')
nlml0 = gp(hyp0, @infGaussLik, meanfunc0, covfunc0, likfunc0, x, y);
%%
z = linspace(-1.9, 1.9, 101)';
[m, s2] = gp(hyp0, @infGaussLik, meanfunc0, covfunc0, likfunc0, x, y, z);
f = [m+2*sqrt(s2); flip(m-2*sqrt(s2))];
fill([z; flip(z)], f, [7 7 7]/8);
hold on
plot(z, m);
plot(x, y, '+')
axis([-1.9 1.9 -0.9 3.9])

%%
meanfunc1 = [];
covfunc1 = @covSEiso;
hyp1.cov = [0; 0];
likfunc1 = @likGauss;
hyp1.lik = log(0.1);
hyp1 = minimize(hyp1, @gp, -100, @infGaussLik, meanfunc1, covfunc1, likfunc1, x, y);
exp(hyp1.lik)
nlml1 = gp(hyp1, @infGaussLik, meanfunc1, covfunc1, likfunc1, x, y);
[m, s2] = gp(hyp1, @infGaussLik, meanfunc1, covfunc1, likfunc1, x, y, z);
f = [m+2*sqrt(s2); flip(m-2*sqrt(s2))];
figure
fill([z; flip(z)], f, [7 7 7]/8)
hold on
plot(z, m)
plot(x, y, '+')
axis([-1.9 1.9 -0.9 3.9])
%%
meanfunc2 = {@meanSum, {@meanLinear, @meanConst}};
covfunc2 = @covSEiso;
likfunc2 = @likGauss;
hyp2.mean = [0; 0];
hyp2.cov = [0; 0];
hyp2.lik = log(0.1);
hyp2 = minimize(hyp2, @gp, -100, @infGaussLik, meanfunc2, covfunc2, likfunc2, x, y);
nlml2 = gp(hyp2, @infGaussLik, meanfunc2, covfunc2, likfunc2, x, y);
[m, s2] = gp(hyp2, @infGaussLik, meanfunc2, covfunc2, likfunc2, x, y, z);
f = [m+2*sqrt(s2); flip(m-2*sqrt(s2))];
figure
fill([z; flip(z)], f, [7 7 7]/8)
hold on
plot(z, m)
plot(x, y, '+')
axis([-1.9 1.9 -0.9 3.9])
%%
meanfunc3 = {@meanSum, {@meanLinear, @meanConst}};
covfunc3 = @covSEiso;
likfunc3 = @likGauss;
hyp3.mean = [0; 0];
hyp3.cov = [0; 0];
hyp3.lik = log(0.1);
mu = 1.0;
s2 = 0.01^2;
prior.mean = {{@priorGauss, mu, s2}, {'priorLaplace', mu, s2}};
prior.lik = {{@priorDelta}};
inf = {@infPrior, @infGaussLik, prior};
hyp3 = minimize(hyp3, @gp, -100, inf, meanfunc3, covfunc3, likfunc3, x, y);
nlml3 = gp(hyp3, @infGaussLik, meanfunc3, covfunc3, likfunc3, x, y);
[m, s2] = gp(hyp3, @infGaussLik, meanfunc3, covfunc3, likfunc3, x, y, z);
f = [m+2*sqrt(s2); flip(m-2*sqrt(s2))];
figure
fill([z; flip(z)], f, [7 7 7]/8)
hold on
plot(z, m)
plot(x, y, '+')
axis([-1.9 1.9 -0.9 3.9])
%%
meanfunc4 = {@meanSum, {@meanLinear, @meanConst}};
covfunc4 = @covSEiso;
likfunc4 = @likGauss;
hyp4.mean = [0; 0];
hyp4.cov = [0; 0];
hyp4.lik = log(0.1);

nu = fix(n/2);
u = linspace(-1.9, 1.9, 101)';
inf = @(varargin) infGaussLik(varargin{:}, struct('s', 0.0));
hyp4 = minimize(hyp4, @gp, -100, inf, meanfunc4, covfunc4, likfunc4, x, y);
nlml4 = gp(hyp4, @infGaussLik, meanfunc4, covfunc4, likfunc4, x, y);
[m, s2] = gp(hyp4, @infGaussLik, meanfunc4, covfunc4, likfunc4, x, y, z);
f = [m+2*sqrt(s2); flip(m-2*sqrt(s2))];
figure
fill([z; flip(z)], f, [7 7 7]/8)
hold on
plot(z, m)
plot(x, y, '+')
axis([-1.9 1.9 -0.9 3.9])













