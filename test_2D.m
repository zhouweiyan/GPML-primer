close all;clear;clc
x1=-2:0.1:2;
x2=2:0.1:4;
len1=length(x1);
len2=length(x2);
[x1,x2]=meshgrid(x1,x2);
y=-x1.^2-(x2-3).^2;
mesh(x1,x2,y)
xlabel('x1');ylabel('x2');zlabel('y');
X=cell(len2,len1);
for i=1:len2
    for j=1:len1
       X{i,j}=[x2(i,j) x1(i,j)];
    end
end

% z=zeros(size(X));
% for i=1:len2
%     for j=1:len1
%         z(i,j)=-(X{i,j}(2))^2-(X{i,j}(1)-3)^2;
%     end
% end

X=reshape(X,[],1);
X=cell2mat(X);
y=reshape(y,[],1);

%%   test set
% xs1=-3+6*rand(1,11);      %   11 random numbers in [-3,3]
xs1=-2+(2+2)*rand(1,11);
lens1=length(xs1);xs1=sort(xs1);
% xs2=1+4*rand(1,11);       %   11 random numbers in [1,5]
xs2=2+2*rand(1,13);
lens2=length(xs2);xs2=sort(xs2);
[xs1,xs2]=meshgrid(xs1,xs2);
Xs=cell(lens2,lens1);
for i=1:lens2
    for j=1:lens1
       Xs{i,j}=[xs2(i,j) xs1(i,j)];
    end
end
Xs=reshape(Xs,[],1);
Xs=cell2mat(Xs);

%%   Gaussian Processing Regression for 2-D data
% meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5;0.5;0.5];
meanfunc=[];hyp.mean=[];
covfunc = {@covSEiso}; ell = 1; sf = 1; hyp.cov = log([ell; sf]);
likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);
hyp=minimize(hyp,@gp,-100,@infGaussLik,[],covfunc,likfunc,X,y);
nlml = gp(hyp, @infGaussLik, [], covfunc, likfunc, X, y)
[m s2]=gp(hyp,@infGaussLik,[],covfunc,likfunc,X,y,Xs);
m=reshape(m,lens2,lens1);
hold on
mesh(xs1,xs2,m)


