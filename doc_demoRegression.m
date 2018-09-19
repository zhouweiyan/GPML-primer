clear 
close all
clc

%%   generate GP data
meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5; 1];
covfunc = {@covMaterniso, 3}; ell = 1/4; sf = 1; hyp.cov = log([ell; sf]);
likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);

n = 20;
x = gpml_randn(0.3, n, 1);
K = feval(covfunc{:}, hyp.cov, x);
mu = feval(meanfunc{:}, hyp.mean, x);
y = chol(K)'*gpml_randn(0.15, n, 1) + mu + exp(hyp.lik)*gpml_randn(0.2, n, 1);
 
figure(1)
set(gca, 'FontSize', 24)
plot(x, y, '+', 'MarkerSize', 12)
axis([-1.9 1.9 -0.9 3.9])
grid on
xlabel('input, x')
ylabel('output, y')
nlml = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, x, y)

z = linspace(-1.9, 1.9, 101)';
[m,s2] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, x, y, z);
figure(2)
set(gca, 'FontSize', 24)
f = [m+2*sqrt(s2); flip(m-2*sqrt(s2),1)];
fill([z; flip(z,1)], f, [7 7 7]/8);
hold on; plot(z, m, 'LineWidth', 2); plot(x, y, '+', 'MarkerSize', 12)
axis([-1.9 1.9 -0.9 3.9])
grid on
xlabel('input, x')
ylabel('output, y')

%%   GP regression 1
meanfunc=[]; 
covfunc={@covSum,{@covSEiso,@covLIN}}; 
likfunc=@likGauss;
ell=0.5; sf=1;
hyp2.cov =log([ell;sf]); hyp2.lik = log(0.1);
hyp2 = minimize(hyp2, @gp, -100, @infGaussLik, [], covfunc, likfunc, x, y);
disp(['exp(hyp2.lik)=',exp(hyp2.lik)])
nlml2 = gp(hyp2, @infGaussLik, [], covfunc, likfunc, x, y)
[m,s2] = gp(hyp2, @infGaussLik, [], covfunc, likfunc, x, y, z);
figure
set(gca, 'FontSize', 24)
f = [m+2*sqrt(s2); flip(m-2*sqrt(s2),1)];
fill([z; flip(z,1)], f, [7 7 7]/8)
hold on; plot(z, m, 'LineWidth', 2); plot(x, y, '+', 'MarkerSize', 12)
grid on
xlabel('input, x')
ylabel('output, y')
axis([-1.9 1.9 -0.9 3.9])

%%   GP regression 2
meanfunc={@meanSum,{@meanLinear,@meanConst}};
covfunc=@covSEiso; likfunc=@likGauss;
hyp.cov = [0; 0]; hyp.mean = [0; 0]; hyp.lik = log(0.1);
hyp = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, x, y);
nlml3 = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, x, y)
[m s2] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, x, y, z);
figure(4)
set(gca, 'FontSize', 24)
f = [m+2*sqrt(s2); flip(m-2*sqrt(s2),1)];
fill([z; flip(z,1)], f, [7 7 7]/8)
hold on; plot(z, m, 'LineWidth', 2); plot(x, y, '+', 'MarkerSize', 12)
grid on
xlabel('input, x')
ylabel('output, y')
axis([-1.9 1.9 -0.9 3.9])

%% Guiding marginal likelihood optimisation with a hyperprior
% mu=0.1;s2=0.01^2;
% prior.mean={{@priorGauss,mu,s2},}
% prior.mean={{@priorGauss,mu,s2};{@priorLaplace,mu,s2}}
% prior.lik={{@priorDelta}};
% inf={@infPrior,@infGaussLik,prior}
% hyp=minimize(hyp,@gp,-100,inf,meanfunc,covfunc,likfunc,x,y);

%% Sparse Approximation Regression based on inducing inputs
% u=linspace(-1.3,1.3,fix(n/2))';
u = x(randperm(n,fix(n/2)),:);
covfuncF = {@apxSparse, {covfunc}, u};
[mF,s2F] = gp(hyp, @infGaussLik, meanfunc, covfuncF, likfunc, x, y, z);

figure(5)
set(gca, 'FontSize', 24)
f = [mF+2*sqrt(s2F); flip(mF-2*sqrt(s2F),1)];
fill([z; flip(z,1)], f, [7 7 7]/8)
hold on; plot(z, mF, 'LineWidth', 2); plot(x, y, '+', 'MarkerSize', 12)
plot(u,1,'ko', 'MarkerSize', 12)
grid on
xlabel('input, x')
ylabel('output, y')
axis([-1.9 1.9 -0.9 3.9])
