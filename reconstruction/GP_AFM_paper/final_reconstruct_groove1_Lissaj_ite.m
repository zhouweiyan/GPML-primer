clear
close all
clc
load('groove1_GP.mat')
set(0,'DefaultFigureWindowStyle','docked') 

%% Lissajous scan trajectory
X0=64;Y0=64;
A=128;B=128;    % A,B represent the range of X,Y respectively
fx=41;fy=43;
f=4000;
T=1;  %T = lsm(fx,fy)/(fx*fy)
t=(0:1/f:T)';
x1_train_lisaj=X0+A/2*sin(2*pi*fx*t);
x2_train_lisaj=Y0+B/2*sin(2*pi*fy*t);
figure
plot(x1_train_lisaj,x2_train_lisaj,'+')
axis equal
axis ([0 128 0 128])
grid on
hold on;plot(x1_train_lisaj,x2_train_lisaj)
%% data preparation
% train set
X_train_lisaj=[x1_train_lisaj(:),x2_train_lisaj(:)];
[y_train_lisaj,y_t_s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_train_lisaj);
% approximate path
path_len=sum(sqrt((x1_train_lisaj(2:end)-x1_train_lisaj(1:(end-1))).^2+(x2_train_lisaj(2:end)-x2_train_lisaj(1:(end-1))).^2))

% test set
x1_test=2:127;x2_test=2:127;
[X1_test,X2_test]=meshgrid(x1_test,x2_test);
X_test=[X1_test(:),X2_test(:)];
[y_test_ideal,y_test_s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_test);
y_test_ideal=reshape(y_test_ideal,length(x2_test),length(x1_test));

%% GP regression
meanfunc=[];hyp.mean=[];
% mask=[1,0];
covfunc={'covSEard'};hyp.cov=log([10;10;20]);
% covfunc={'covMask',{mask,{'covSEiso'}}};hyp.cov=log([2;40]);  %log([ell;sf])
likfunc='likGauss';hyp.lik=log(0.1);    % log(sn)
[hyp,hyp_warehouse,i_warehouse]=minimize_ite(hyp,@gp,-80,@infGaussLik,meanfunc,covfunc,likfunc,X_train_lisaj,y_train_lisaj);
PSNR_warehouse=zeros(length(i_warehouse),1);

for i=1:length(i_warehouse)
    i
    hyp_i=rewrap_out(hyp,hyp_warehouse(:,i));
    [m,s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train_lisaj,y_train_lisaj,X_test);
    m_a=reshape(m,length(x2_test),length(x1_test));
    PSNR_warehouse(i)=max(y_test_ideal(:))*sqrt(size(y_test_ideal,1)*size(y_test_ideal,2))/sqrt(sum((m_a(:)-y_test_ideal(:)).^2))
    PSNR_warehouse(i)=20*log(PSNR_warehouse(i))/log(10);
end
plot(i_warehouse,PSNR_warehouse,'-o')
xlabel('iteration/times');ylabel('PSNR/dB')