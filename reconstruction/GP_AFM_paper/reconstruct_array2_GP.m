% isotropic/unisotropic undersampling
% array2: 128��128 pixels
% zhouweiyan 20180925
% done
clear
% close all
clc
load array2.mat
figure
surf(array);colormap(jet)
% view([0 -90])
axis equal
xlim([1,128]);ylim([1,128]);zlim([-50,40]);
xlabel('x1');ylabel('x2');view(3)

%% data preparation
% train set
x1_train=2:4:128; x2_train=1:1:128;
[X1_train,X2_train]=meshgrid(x1_train,x2_train);
X_train=[X1_train(:),X2_train(:)];
y_train=array(x2_train,x1_train);
y_train=y_train(:);
figure
surf(X1_train,X2_train,reshape(y_train,length(x2_train),length(x1_train)));colormap(jet)
axis equal
xlim([1,128]);ylim([1,128]);zlim([-50,40]);
xlabel('x1');ylabel('x2');view(3)

% test set
x1_test=1:128;x2_test=1:128;
[X1_test,X2_test]=meshgrid(x1_test,x2_test);
X_test=[X1_test(:),X2_test(:)];
y_test_ideal=array(x2_test,x1_test);
%% GP regression
meanfunc=[];hyp.mean=[];
mask = [0,1];
cgi = {'covSEiso'};
cma = {'covMask',{mask,cgi{:}}}; hypma = log([1.5;15]);%log([ell;sf]);
covfunc1={'covSum',{cma,'covSEiso'}};hyp1.cov=[hypma;log([2;3])];%relative good
covfunc2={'covMaternard',1};hyp2.cov=log([1;1;4]);
covfunc={'covSum',{covfunc1,covfunc2}};hyp.cov=[hyp1.cov;hyp2.cov];

likfunc='likGauss';hyp.lik=log(10);    % log(sn)
hyp=minimize(hyp,@gp,-80,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train);
[m,s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_test);
figure
surf(X1_test,X2_test,reshape(m,length(x2_test),length(x1_test)))
colormap(jet)
axis equal
xlim([1,128]);ylim([1,128]);zlim([-50,40]);
xlabel('x1');ylabel('x2');view(3)
%% Error analysis
range=max(array(:))-min(array(:));
% NRMSD=sqrt(sum(m-y_test_ideal(:)).^2)/range
% max(m-y_test_ideal(:))

m_a=reshape(m,length(x2_test),length(x1_test));
m_a=m_a(min(x2_train):max(x2_train),min(x1_train):max(x1_train));
y_test_ideal_a=y_test_ideal(min(x2_train):max(x2_train),min(x1_train):max(x1_train));
figure
surf(y_test_ideal_a-m_a);xlabel('x1');ylabel('x2');
colormap(jet)
NRMSD=sqrt(sum((m_a(:)-y_test_ideal_a(:)).^2))/range
max(abs(m_a(:)-y_test_ideal_a(:)))
PSNR=max(y_test_ideal_a(:))*sqrt(size(y_test_ideal_a,1)*size(y_test_ideal_a,2))/sqrt(sum((m_a(:)-y_test_ideal_a(:)).^2))

% scanning pattern
% figure
% plot3(X1_train,X2_train,ones(128,32),'*')