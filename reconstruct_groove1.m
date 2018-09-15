clear
close all
clc
load groove1.txt
groove=reshape(groove1,256,256);
figure
surf(groove);colormap(jet)
% view([0 -90])
axis equal
xlim([1,256]);ylim([1,256]);zlim([-90,40]);
xlabel('x1');ylabel('x2');

%% data preparation
% train set
% x1_train=2:2:256;x2_train=2:2:256;
% groove1_2=groove1(x1_train,x2_train);
x1_train=4:4:256;x2_train=4:4:256;
[X1_train,X2_train]=meshgrid(x1_train,x2_train);
X_train=[X1_train(:),X2_train(:)];
y_train=groove(x1_train,x2_train);y_train=y_train(:);
figure
surf(X1_train,X2_train,reshape(y_train,length(x2_train),length(x1_train)));colormap(jet)
axis equal
xlabel('x1');ylabel('x2');zlim([-90,40]);

% test set
x1_test=1:256;x2_test=1:256;
[X1_test,X2_test]=meshgrid(x1_test,x2_test);
X_test=[X1_test(:),X2_test(:)];
y_test_ideal=groove(x1_train,x2_train);y_test_ideal=y_test_ideal(:);
%% GP regression
meanfunc=[];hyp.mean=[];
% mask=[1,0];
covfunc={'covSEard'};hyp.cov=log([10;10;20]);
% covfunc={'covMask',{mask,{'covSEiso'}}};hyp.cov=log([2;40]);  %log([ell;sf])
likfunc='likGauss';hyp.lik=log(0.1);    % log(sn)
hyp=minimize(hyp,@gp,-50,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train);
[m,s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_test);
figure
surf(X1_test,X2_test,reshape(m,length(x2_test),length(x1_test)))
colormap(jet)
axis equal
xlim([1,256]);ylim([1,256]);zlim([-90,40]);
xlabel('x1');ylabel('x2');
view(3)
%% Error analysis
range=max(groove(:))-min(groove(:));
nrmsd=sqrt(sum((m-groove1).^2))/range




