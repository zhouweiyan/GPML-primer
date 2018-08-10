% zhouweiyan 20180730
% reconstruction of sinusoidal and parabola surface
clear
close all
clc
%% target figure
x1=-10:0.1:10;
x2=-10:0.1:10;
[X1,X2]=meshgrid(x1,x2);
y_ideal=truefunc(X1,X2);
y_noise=noisefunc(X1,X2);
set(0,'DefaultFigureWindowStyle','docked'); % keep figures tidy
subplot(121);surf(X1,X2,reshape(y_ideal,length(x1),length(x2)));
title('target surface without noise')
subplot(122);surf(X1,X2,reshape(y_noise,length(x1),length(x2)));
title('target surface with noise');
colormap(jet)

%%% generate data
% train data
x1_train_length=floor(length(x1)/4); x2_train_length=floor(length(x2)/4);
% rng   %function里面的设置影响这里，下面的随机每次会保持相同
rng(1,'v5normal')
x1_train=x1(randperm(length(x1),x1_train_length));
x2_train=x2(randperm(length(x2),x2_train_length));

[X1_train,X2_train]=meshgrid(x1_train,x2_train);
X_train=[X1_train(:),X2_train(:)];
y_train=truefunc(X1_train,X2_train);
% y_train=noisefunc(X1_train,X2_train);

x1_train_temp=x1(sort(randperm(length(x1),x1_train_length)));
x2_train_temp=x2(sort(randperm(length(x2),x2_train_length)));
rng('shuffle')  % 生成不同随机数
[X1_train_temp,X2_train_temp]=meshgrid(x1_train_temp,x2_train_temp);
y_train_temp=truefunc(X1_train_temp,X2_train_temp);
% y_train_temp=noisefunc(X1_train_temp,X2_train_temp);
figure(2)
surf(X1_train_temp,X2_train_temp,reshape(y_train_temp,length(x1_train),length(x2_train)));
title('data used to train GPR')

% test data
x1_test=-10:0.3:10;
x2_test=-10:0.3:10;
ns=[length(x1_test);length(x2_test)];
[X1_test,X2_test]=meshgrid(x1_test,x2_test);
X_test=[X1_test(:) X2_test(:)];
f_test=truefunc(X1_test,X2_test);

%% gaussian regression 1
meanfunc=[];
covfunc={'covSEiso'};
likfunc='likGauss';
opt='exact';
hyp1.mean=[];
ell=0.1; sf=10; sn=0.1; 
% test_sinu_paraboloid_GA.m
% ell=trace(1,end);sf=trace(2,end);sn=trace(3,end);
% test_sinu_paraboloid_PSO.m
% ell=zbest(1);sf=zbest(2);sn=zbest(3);
hyp1.cov=log([ell;sf]);
hyp1.lik=log(sn);
switch opt
    case 'exact'
        tic
        hyp1=minimize(hyp1,@gp,-100,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train);
        % load('test_sinu_paraboloid_hyp1.mat');
        disp(['exp(hyp1.lik)=',num2str(exp(hyp1.lik))])
        nlml=gp(hyp1,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train)
        [m s2]=gp(hyp1,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_test);
        figure
        surf(X1_test,X2_test,reshape(m,length(x2_test),length(x1_test)));
        title('GPR');
%         hold on
%         surf(X1_test,X2_test,reshape(m+2*sqrt(s2),length(x2_test),length(x1_test)),'FaceAlpha',0.1);
%         surf(X1_test,X2_test,reshape(m-2*sqrt(s2),length(x2_test),length(x1_test)),'FaceAlpha',0.1);
        toc   
    case 'apx'
        tic
        % corresponding sparse approximate regression based on inducing inputs
        u = [x1(randperm(length(x1),floor(length(x1)/4)))';x2(randperm(length(x2),floor(length(x2)/4)))'];
        covfuncF={@apxSparse,{covfunc},u};
        inf=@(varargin)infGaussLik(varargin{:},struct('s',0));
        hyp1=minimize(hyp1,@gp,-100,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train);
        disp(['exp(hyp1.lik)=',num2str(exp(hyp1.lik))])
        nlml1apx=gp(hyp1,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train);
        [m s2]=gp(hyp1,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_test);
        figure
        surf(X1_test,X2_test,reshape(m,length(X1_test),length(X2_test)));
        title('GPR');
        toc  
end

%% gaussian regression 2
% meanfunc=[];
% covfunc={'covPER','iso',{'covSEiso'}};
% likfunc='likGauss';
% hyp.mean=[];
% p=20; ell=0.6; sf=10; 
% hyp.cov=log([p;ell;sf]);
% sn=0.1; hyp.lik=log(sn);
% tic
% hyp=minimize(hyp,@gp,-500,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train);
% disp(['exp(hyp1.lik)=',num2str(exp(hyp1.lik))])
% nlml=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train);
% [m s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_test);
% figure
% surf(X1_test,X2_test,reshape(m,length(X1_test),length(X2_test)));
% title('GPR');
% toc  

%% error analysis
error=f_test-m;
figure
surf(X1_test,X2_test,reshape(error,length(x2_test),length(x1_test)));
view([-90 90])
title('error');
colorbar
disp(['sum(abs(error(:)))=',num2str(sum(abs(error(:))))])
disp(['sum(abs(s2(:)))=',num2str(sum(abs(s2(:))))])
%% define functions
function y_ideal=truefunc(X1,X2)
    y_ideal=-0.01*X1(:).^2-0.01*X2(:).^2+0.15*cos(2*X1(:))+0.15*cos(2*X2(:));
end
function y_noise=noisefunc(X1,X2)
    seed=1;
    rng(seed,'v4');
    s_noise=0.1;
    y_noise=-0.01*X1(:).^2-0.01*X2(:).^2+0.15*cos(2*X1(:))+0.15*cos(2*X2(:))...
        +s_noise*randn(length(X1)*length(X2),1);
end


