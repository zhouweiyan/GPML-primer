clear
close all
clc
load('groove1_GP.mat')
set(0,'DefaultFigureWindowStyle','docked') 

%% Constant Angular Velocity(CAV)
r_end=64*1.6;
P=r_end*2/(60-1);
%P=spiral radius¡Á2/(number of curves-1)
%eq.(4)
w=120;
f=2000;
T=1/f;
T_spiral=2*pi*r_end/(P*w);	%CAV eq.(6)
t=0:T:T_spiral;
xs=P*w/(2*pi)*t.*cos(w*t)+64;	%eq.(17)
ys=P*w/(2*pi)*t.*sin(w*t)+64;	%eq.(18)
figure
plot(xs,ys,'+')
axis equal
axis ([0 128 0 128])
grid on
hold on;plot(xs,ys)

%plot the outside dashed circle
% hold on
% xita=0:(pi/180):(2*pi);
% x=r_end*cos(xita)+64;
% y=r_end*sin(xita)+64;
% plot(x,y,'--')
%% data preparation
% train set
x1_train_spi=xs';x2_train_spi=ys';
X_train_spi=[x1_train_spi(:),x2_train_spi(:)];
[y_train_spi,y_t_s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_train_spi);
% approximate path
path_len=sum(sqrt((x1_train_spi(2:end)-x1_train_spi(1:(end-1))).^2+(x2_train_spi(2:end)-x2_train_spi(1:(end-1))).^2))

% test set
x1_test=1:128;x2_test=1:128;
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
[hyp,hyp_warehouse,i_warehouse]=minimize_ite(hyp,@gp,-80,@infGaussLik,meanfunc,covfunc,likfunc,X_train_spi,y_train_spi);
PSNR_warehouse=zeros(length(i_warehouse),1);
for i=1:length(i_warehouse)
    i
    hyp_i=rewrap_out(hyp,hyp_warehouse(:,i));
    [m,s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train_spi,y_train_spi,X_test);
    m_a=reshape(m,length(x2_test),length(x1_test));
    m_a=m_a(min(x2_train):max(x2_train),min(x1_train):max(x1_train));
    y_test_ideal_a=y_test_ideal(min(x2_train):max(x2_train),min(x1_train):max(x1_train));
    PSNR_warehouse(i)=max(y_test_ideal_a(:))*sqrt(size(y_test_ideal_a,1)*size(y_test_ideal_a,2))/sqrt(sum((m_a(:)-y_test_ideal_a(:)).^2));
    PSNR_warehouse(i)=20*log(PSNR_warehouse(i))/log(10);
end
plot(i_warehouse,PSNR_warehouse,'-o')
xlabel('iteration/times');ylabel('PSNR/dB')