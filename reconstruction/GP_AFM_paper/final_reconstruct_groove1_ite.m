% isotropic/unisotropic undersampling
% combine delaunayTriangulation interpolation and GPR
% zhouweiyan 20181119
clear
% close all
clc
load groove1.mat
figure;
surf(groove);colormap(jet)
% view([0 -90])
axis equal
xlim([1,128]);ylim([1,128]);zlim([-90,40]);
xlabel('x1');ylabel('x2');view(-45,50)

%% data preparation
% train set
x1_train=2:6:128; x2_train=1:128;
[X1_train,X2_train]=meshgrid(x1_train,x2_train);
X_train=[X1_train(:),X2_train(:)];
y_train=groove(x2_train,x1_train);
y_train=y_train(:);
figure
surf(X1_train,X2_train,reshape(y_train,length(x2_train),length(x1_train)));colormap(jet)
axis equal
xlabel('x1');ylabel('x2');zlim([-90,40]);

% test set
x1_test=min(x1_train):max(x1_train);x2_test=min(x2_train):max(x2_train);
[X1_test,X2_test]=meshgrid(x1_test,x2_test);
X_test=[X1_test(:),X2_test(:)];
y_test_ideal=groove(x2_test,x1_test);

%% GP regression
meanfunc=[];hyp.mean=[];
% mask=[1,0];
covfunc={'covSEard'};hyp.cov=log([10;10;20]);
% covfunc={'covMask',{mask,{'covSEiso'}}};hyp.cov=log([2;40]);  %log([ell;sf])
likfunc='likGauss';hyp.lik=log(0.1);    % log(sn)
[hyp,hyp_warehouse,i_warehouse]=minimize_ite(hyp,@gp,-80,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train);
PSNR_warehouse=zeros(length(i_warehouse),1);
s2_warehouse=zeros(length(i_warehouse),length(X_test));
sigma=zeros(length(i_warehouse),1);
for i=1:length(i_warehouse)
    hyp_i=rewrap_out(hyp,hyp_warehouse(:,i));
    [m,s2]=gp(hyp_i,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_test);
    m_a=reshape(m,length(x2_test),length(x1_test));
    PSNR_warehouse(i)=max(y_test_ideal(:))*sqrt(length(x1_test)*length(x2_test))/sqrt(sum((m_a(:)-y_test_ideal(:)).^2));
    PSNR_warehouse(i)=20*log(PSNR_warehouse(i))/log(10);
%     sigma(i)=sqrt(sum(s2(:))/length(m));
end
plot(i_warehouse,PSNR_warehouse,'-o')
xlabel('iteration/times');ylabel('PSNR/dB')