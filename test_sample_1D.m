%   sampling method based on GP regression
%   zhouweiyan   20180713
clear
close all
clc

%   ideal function
x_0=(-pi:0.01:pi)';
y_0=sin(2*x_0);
figure(1)
plot(x_0,y_0)

%   random sampling points
N=1;
x0_len=length(x_0);
r=randi([1,x0_len],N,1);
x=x_0(r);
y=y_0(r);
figure(2)
set(gca, 'FontSize', 24)
plot(x, y, '+', 'MarkerSize', 12)

meanfunc=[];
covfunc=@covSEiso;
likfunc=@likGauss;
hyp.mean=[];
hyp.cov=[0;0];
hyp.lik=0.1;
hyp=minimize(hyp,@gp,-100,@infGaussLik,meanfunc,covfunc,likfunc,x,y);
nlml=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,x,y)
exp(hyp.cov)    %   inferred length scale and standard deviation
exp(hyp.lik)    %   inferred noise standard deviation
xs=linspace(-pi,pi,x0_len)';
[m s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,x,y,xs);
f=[m+2*sqrt(s2);flip(m-2*sqrt(s2))];
figure(3)
fill([xs;flip(xs)],f,[7 7 7]/8);
hold on
plot(xs,m)
plot(x,y,'b+')
title('GP Regression for the 1st time')
[M,I]=max(s2);
ite=1;


while (2*sqrt(M)>0.01)||(ite<100)
    ite=ite+1;
    x=[x;x_0(I)];
    y=[y;y_0(I)];
    hyp.mean=[];
    hyp.cov=[0;0];
    hyp.lik=0.1;
    hyp=minimize(hyp,@gp,-100,@infGaussLik,meanfunc,covfunc,likfunc,x,y);
    nlml=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,x,y)
    exp(hyp.cov)    %   inferred length scale and standard deviation
    exp(hyp.lik)    %   inferred noise standard deviation
    [m s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,x,y,xs);
    f=[m+2*sqrt(s2);flip(m-2*sqrt(s2))];
    [M,I]=max(s2);
end
figure
fill([xs;flip(xs)],f,[7 7 7]/8);
hold on
plot(xs,m)
plot(x,y,'b+')
str=['GP Regression for the ',num2str(ite),'th time'];
title(str);
