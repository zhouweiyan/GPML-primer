load array1.txt
array=reshape(array1,256,256);

figure;surf(array)
xlim([1,256]);ylim([1,256]);colormap(jet)
h=fspecial('average',3);
a3=imfilter(array,h,'corr','replicate');
figure;surf(a3)
xlim([1,256]);ylim([1,256]);colormap(jet)
array=a3;
%% data preparation
% train set
x1_train=1:2:128; x2_train=1:2:128;
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
% covfunc={'covSEard'};hyp.cov=log([1;0.5;3]);%1.020571e+04
% covfunc={'covMaternard',3};hyp.cov=log([1;0.5;3]);%9e3
% cad = {'covADD',{1,{'covMaternard',3},{'covMaternard',3}}};%1.042608e+04
% covfunc=cad; hyp.cov=log([7;1;5;4;4]);
mask = [0,1];
% cgi = {'covSEiso'};
% cma = {'covMask',{mask,cgi{:}}}; hypma = log([1.5;15]);%log([ell;sf]);
% covfunc={'covSum',{cma,'covSEiso'}};hyp.cov=[hypma;log([2;3])];%best now
cgi = {'covMaternard',3};
cma = {'covMask',{mask,cgi}}; hypma = log([1.5;15]);%log([ell;sf]);
covfunc={'covSum',{cma,'covSEiso'}};hyp.cov=[hypma;log([2;3])];%best now
likfunc='likGauss';hyp.lik=log(1);    % log(sn)
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
% m_a=m_a(1:256,4:252);y_test_ideal_a=y_test_ideal(1:256,4:252);    % cut
% the boundary of x1
m_a=m_a(min(x2_train):max(x2_train),min(x1_train):max(x1_train));y_test_ideal_a=y_test_ideal(min(x2_train):max(x2_train),min(x1_train):max(x1_train));
figure
surf(y_test_ideal_a-m_a);xlabel('x1');ylabel('x2');
colormap(jet)
NRMSD=sqrt(sum((m_a(:)-y_test_ideal_a(:)).^2))/range
max(m_a(:)-y_test_ideal_a(:))