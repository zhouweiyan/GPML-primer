% demonstrate the flow of regression
% http://www.gaussianprocess.org/gpml/code/matlab/doc/
% zwy 20190120

%% clear the workspace
clear  
close all
clc
%% 3c) A Practical Example
x = gpml_randn(0.8, 20, 1);
y =sin(3*x) + 0.1*gpml_randn(0.9, 20, 1);
xs = linspace(-3, 3, 61)';

meanfunc = [];
covfunc = @covSEiso;
likfunc = @likGauss;
hyp = struct('mean', [], 'cov', [0 0], 'lik', -1);
hyp2 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, x, y);
[mu, s2] = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc, x, y, xs);
f = [mu+2*sqrt(s2); flip(mu-2*sqrt(s2))];
fill([xs; flip(xs)], f, [7 7 7]/8)
hold on
plot(xs, mu)
plot(x,y,'+')
%% 4a) Simple Regression
% generate some data from a Gaussian process, i.e., the training set
meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5; 1];
covfunc = {@covMaterniso, 3}; ell = 1/4; sf = 1; hyp.cov = log([ell; sf]);
likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);   

n = 20;
x = gpml_randn(0.3, n, 1);
K = feval(covfunc{:}, hyp.cov, x);
mu = feval(meanfunc{:}, hyp.mean, x);
% exp(hyp.lik)is the noise of y. y is generated according to GPML appendix A.2.
y = chol(K)'*gpml_randn(0.15, n, 1) + mu + exp(hyp.lik)*gpml_randn(0.2, n, 1);  

% plot the generated sample data
figure(100)
plot(x, y, '+', 'MarkerSize', 12)
axis([-1.9 1.9 -0.9 3.9])
grid on
xlabel('input, x')
ylabel('output, y')

% calculate the negative log marginal likelihood, i.e., -log(y|X)
nlml = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, x, y)
% generate the testing set
z = linspace(-1.9, 1.9, 101)';

% predict with the underlying gaussian process
[m, s2] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, x, y, z);
% plot the predicted mean and variance
figure(1)
f = [m+2*sqrt(s2); flip(m-2*sqrt(s2))];
fill([z; flip(z)], f, [7 7 7]/8);
hold on; plot(z, m, 'LineWidth', 2); plot(x, y, '+', 'MarkerSize', 12)
axis([-1.9 1.9 -0.9 3.9])
grid on
xlabel('input, x')
ylabel('output, y')

%% regression with the different meanfunc and covfunc
meanfunc = [];
% covfunc = @covSEiso; hyp2.cov = [0; 0]; 
covfunc = {@covSum,{@covLINiso,@covConst,@covSEiso}};hyp2.cov = [0; 0; 0; 0]; 
hyp2.lik = log(0.1);

% optimize the hyperparameters
hyp2 = minimize(hyp2, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, x, y);

exp(hyp2.lik)   % predicted noise level
% calculate the negative log marginal likelihood, i.e., -log(y|X)
nlml2 = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc, x, y)
% prediction
[m, s2] = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc, x, y, z);

% plot the predicted mean and variance
figure(2)
f = [m+2*sqrt(s2); flip(m-2*sqrt(s2))];
fill([z; flip(z)], f, [7 7 7]/8)
hold on; plot(z, m, 'LineWidth', 2); plot(x, y, '+', 'MarkerSize', 12)
grid on
xlabel('input, x')
ylabel('output, y')
axis([-1.9 1.9 -0.9 3.9])
%% regression with the same meanfunc and different covfunc 
meanfunc = {@meanSum, {@meanLinear, @meanConst}};
covfunc = @covSEiso;
likfunc = @likGauss;
hyp.cov = [0; 0]; hyp.mean = [0; 0]; hyp.lik = log(0.1);
hyp = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, x, y);
[m, s2] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, x, y, z);
nlml3 = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, x, y)

% plot the predicted mean and variance
figure(3)
f = [m+2*sqrt(s2); flip(m-2*sqrt(s2))];
fill([z; flip(z)], f, [7 7 7]/8)
hold on; plot(z, m, 'LineWidth', 2); plot(x, y, '+', 'MarkerSize', 12)
grid on
xlabel('input, x')
ylabel('output, y')
axis([-1.9 1.9 -0.9 3.9])
%% Guiding marginal likelihood optimisation with a hyperprior
meanfunc = {@meanSum, {@meanLinear, @meanConst}};
covfunc = @covSEiso;
likfunc = @likGauss;
hyp.cov = [0; 0]; hyp.mean = [0; 0]; hyp.lik = log(0.1);

mu=0.1;s2=0.01^2;
prior.mean={{@priorGauss,mu,s2},};
prior.mean={{@priorGauss,mu,s2};{@priorLaplace,mu,s2}};
prior.lik={{@priorDelta}};
inf = {@infPrior,@infGaussLik,prior};

hyp=minimize(hyp,@gp,-100,inf,meanfunc,covfunc,likfunc,x,y);
[m, s2] = gp(hyp, inf, meanfunc, covfunc, likfunc, x, y, z);
nlml4 = gp(hyp, inf, meanfunc, covfunc, likfunc, x, y)

% plot the predicted mean and variance
figure(4)
f = [m+2*sqrt(s2); flip(m-2*sqrt(s2))];
fill([z; flip(z)], f, [7 7 7]/8)
hold on; plot(z, m, 'LineWidth', 2); plot(x, y, '+', 'MarkerSize', 12)
grid on
xlabel('input, x')
ylabel('output, y')
axis([-1.9 1.9 -0.9 3.9])
%% sparse gaussian regression
meanfunc = {@meanSum, {@meanLinear, @meanConst}};
covfunc = @covSEiso;
likfunc = @likGauss;
hyp.cov = [0; 0]; hyp.mean = [0; 0]; hyp.lik = log(0.1);

nu = fix(n/2); u = linspace(-1.3,1.3,nu)';
covfuncF = {@apxSparse, {covfunc}, u};

hyp = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfuncF, likfunc, x, y);
[mF, s2F] = gp(hyp, @infGaussLik, meanfunc, covfuncF, likfunc, x, y, z);
nlml5 = gp(hyp, @infGaussLik, meanfunc, covfuncF, likfunc, x, y)

% plot the predicted mean and variance
figure(5)
f = [mF+2*sqrt(s2F); flip(mF-2*sqrt(s2F))];
fill([z; flip(z)], f, [7 7 7]/8)
hold on; plot(z, mF, 'LineWidth', 2); plot(x, y, '+', 'MarkerSize', 12)
plot(u,1,'ko', 'MarkerSize', 12)
grid on
xlabel('input, x')
ylabel('output, y')
axis([-1.9 1.9 -0.9 3.9])

