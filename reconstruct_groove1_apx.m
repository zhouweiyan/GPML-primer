% isotropic/unisotropic undersampling
% zhouweiyan 20180915
clear
close all
clc
load groove1.txt
groove=reshape(groove1,256,256);
figure
surf(groove);colormap(jet)
% view([0 -90])
axis equal
xlim([1,256]);ylim([1,256]);%zlim([-90,40]);
xlabel('x1');ylabel('x2');

%% data preparation
% train set
x1_train=2:4:256; x2_train=2:4:256;
[X1_train,X2_train]=meshgrid(x1_train,x2_train);
X_train=[X1_train(:),X2_train(:)];
y_train=groove(x2_train,x1_train);
y_train=y_train(:);
figure
surf(X1_train,X2_train,reshape(y_train,length(x2_train),length(x1_train)));colormap(jet)
axis equal
xlim([1,256]);ylim([1,256]);%zlim([-50,40]);
xlabel('x1');ylabel('x2');view(3)
%%
X_train_apx=apxGrid('expand',{x1_train',x2_train'});
y_train_apx=groove(x2_train,x1_train);
y_train_apx=y_train_apx';
y_train_apx=y_train_apx(:);
% test set
x1_test=1:256;x2_test=1:256;
X_test=apxGrid('expand',{x1_test',x2_test'});
y_test_ideal=groove(x2_test,x1_test);
y_test_ideal=y_test_ideal';
% y_test_ideal=y_test_ideal(:);
%% GP regression apx
meanfunc=[];hyp.mean=[];
% covfunc={'covSEard'};hyp.cov=log([4;2;3]);%1.020571e+04
xg=apxGrid('create',[X_train;X_test],true,[100,120]);
% cov={'covSEiso','covSEiso'};
cov={{'covMaternard',3},{'covMaternard',3}};
covg={@apxGrid,cov,xg};hyp.cov=log([4;2;4;2]);
likfunc={'likGauss'};hyp.lik=log(1);    % log(sn)
opt.cg_maxit=200;opt.cg_tol=5e-3;
infg=@(varargin)infGaussLik(varargin{:},opt);
hyp=minimize(hyp,@gp,-80,infg,meanfunc,covg,likfunc,X_train,y_train);
[post,nlZ,dnlZ]=infGrid(hyp,{'meanZero'},covg,likfunc,X_train,y_train,opt);
[fm,fs2,m,ys2]=post.predict(X_test);
[X1_test,X2_test]=meshgrid(x1_test,x2_test);
m=reshape(m,length(x1_test),length(x2_test))';
figure
surf(X1_test,X2_test,m);m=m(:);
xlabel('x1');ylabel('x2')
colormap(jet)
axis equal
xlim([1,256]);ylim([1,256]);%zlim([-50,40]);
xlabel('x1');ylabel('x2');view(3)
%% Error analysis
range=max(groove(:))-min(groove(:));
% NRMSD=sqrt(sum(m-y_test_ideal(:)).^2)/range
% max(m-y_test_ideal(:))

m_a=reshape(m,length(x2_test),length(x1_test));
% m_a=m_a(1:256,4:252);y_test_ideal_a=y_test_ideal(1:256,4:252);    % cut
% the boundary of x1
m_a=m_a(min(x2_train):max(x2_train),min(x1_train):max(x1_train));
y_test_ideal_a=y_test_ideal(min(x2_train):max(x2_train),min(x1_train):max(x1_train));
% surf(m_a)
NRMSD=sqrt(sum((m_a(:)-y_test_ideal_a(:)).^2))/range
max(m_a(:)-y_test_ideal_a(:))
