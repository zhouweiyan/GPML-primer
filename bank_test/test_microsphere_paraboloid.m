% zhouweiyan 20180811
% reconstruction of microsphere and paraboloid surface
clear
close all
clc
%% target figure
x1=-10:0.1:10;
x2=-10:0.1:10;
len=length(x1);
[X1,X2]=meshgrid(x1,x2);
load('test_microsphere_paraboloid_data.mat')
y_ideal=reshape(y_ideal,length(x2),length(x1));
y_noise=reshape(y_noise,length(x2),length(x1));
set(0,'DefaultFigureWindowStyle','docked'); % keep figures tidy
subplot(121);surf(X1,X2,y_ideal,'EdgeColor','none',... 
    'LineStyle','none','FaceLighting','phong');
title('target surface without noise')
subplot(122);surf(X1,X2,y_noise,'EdgeColor','none',...
    'LineStyle','none','FaceLighting','phong');
title('target surface with noise');
colormap(jet)

%% generate data
% train data
x1_train_length=50; x2_train_length=50;
rng(1,'v5normal')
x1_train_index=randperm(length(x1),x1_train_length);
x2_train_index=randperm(length(x2),x2_train_length);
rng('shuffle')  % invalid the previous rng
x1_train=x1(x1_train_index);
x2_train=x2(x2_train_index);
[X1_train,X2_train]=meshgrid(x1_train(:),x2_train(:));
X_train=[X1_train(:),X2_train(:)];
y_train=zeros(x2_train_length,x1_train_length);
for i=1:x1_train_length
    for j=1:x2_train_length
        y_train(j,i)=y_ideal(x2_train_index(j),x1_train_index(i));
        %y_train(j,i)=y_noise(x2_train_index(j),x1_train_index(i));
    end
end
y_train=y_train(:);

x1_train_index_temp=sort(x1_train_index);
x2_train_index_temp=sort(x2_train_index);
x1_train_temp=x1(x1_train_index_temp);
x2_train_temp=x2(x2_train_index_temp);
[X1_train_temp,X2_train_temp]=meshgrid(x1_train_temp,x2_train_temp);
y_train_temp=zeros(x2_train_length,x1_train_length);
for i=1:x1_train_length
    for j=1:x2_train_length
        y_train_temp(j,i)=y_ideal(x1_train_index_temp(j),x2_train_index_temp(i));
        %y_train(j,i)=y_noise(x2_train_index(j),x1_train_index(i));
    end
end
figure(2)
surf(X1_train_temp,X2_train_temp,reshape(y_train_temp,length(x1_train),length(x2_train)));
title('data used to train GPR')

% test data
x1_test_index=1:3:length(x1);
x2_test_index=1:3:length(x2);
x1_test=x1(x1_test_index);
x2_test=x2(x2_test_index);
[X1_test,X2_test]=meshgrid(x1_test,x2_test);
X_test=[X1_test(:),X2_test(:)];
y_test=zeros(length(x2_test_index),length(x1_test_index));
for i=1:length(x1_test_index)
    for j=1:length(x2_test_index)
        y_test(j,i)=y_ideal(x2_test_index(j),x1_test_index(i));
    end
end
y_test=y_test(:);
%% gaussian regression 1
meanfunc={'meanPoly',2};
% covfunc={'covPER','iso',{@covSEiso}};
covfunc={@covSEiso};
likfunc='likGauss';
opt='exact';
hyp1.mean=[1;1;1;1];
% hyp_p=log(2); ell = 0.9; sf = 3; hypgi = log([ell;sf]); hyp1.cov=[hyp_p;hypgi];
hyp1.cov=log([1;3]);
sn=0.1; 
hyp1.lik=log(sn);
switch opt
    case 'exact'
        tic
        hyp1=minimize(hyp1,@gp,-100,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train);
%         load('test_sinu_paraboloid_hyp1.mat');
        disp(['exp(hyp1.lik)=',num2str(exp(hyp1.lik))])
        nlml=gp(hyp1,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train)
        [m,s2]=gp(hyp1,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_test);
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
        [m,s2]=gp(hyp1,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_test);
        figure
        surf(X1_test,X2_test,reshape(m,length(X1_test),length(X2_test)));
        title('GPR');
        toc  
end
%% error analysis
error=y_test-m;
rela_error=abs(y_test-m)./y_test;
figure
surf(X1_test,X2_test,reshape(rela_error,length(x2_test),length(x1_test)));
view([-90 90])
figure
surf(X1_test,X2_test,reshape(error,length(x2_test),length(x1_test)));
view([-90 90])
title('error');
colorbar
disp(['sum(abs(error(:)))=',num2str(sum(abs(error(:))))])
disp(['sum(abs(s2(:)))=',num2str(sum(abs(s2(:))))])
