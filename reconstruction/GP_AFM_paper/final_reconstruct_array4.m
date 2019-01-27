% reconstruct array4
% zwy 20190122

clear
close all
clc
opt='GP';
% opt='DT';
load array4.txt
array=reshape(array4,256,256);
set(0,'DefaultFigureWindowStyle','docked'); % keep figures tidy
array=array(1:128,129:256);
figure
surf(array);colormap(jet)
axis equal
xlim([1,128]); ylim([1,128]); zlim([-5,5]);
xlabel('x1'); ylabel('x2'); view([0 -90])
%% data preparation
% train set
x1_train=2:2:128; x2_train=1:2:128;
[X1_train,X2_train]=meshgrid(x1_train,x2_train);
X_train=[X1_train(:),X2_train(:)];
y_train=array(x2_train,x1_train);
y_train=y_train(:);
figure
surf(X1_train,X2_train,reshape(y_train,length(x2_train),length(x1_train)));
colormap(jet)
axis equal
xlim([1,128]); ylim([1,128]); zlim([-5,5]);
xlabel('x1'); ylabel('x2'); view([0 -90])

% test set
x1_test=min(x1_train):max(x1_train);x2_test=min(x2_train):max(x2_train);
[X1_test,X2_test]=meshgrid(x1_test,x2_test);
X_test=[X1_test(:),X2_test(:)];
y_test_ideal=array(x2_test,x1_test);
%% GP regression/DT linear interplotation
switch opt
    case 'GP'
        meanfunc=[];hyp.mean=[];
        mask = [0,1];
        covfunc = {'covSEiso'}; hyp.cov = log([10,2])
        likfunc='likGauss';hyp.lik=log(0.1);    % log(sn)
        hyp=minimize(hyp,@gp,-80,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train);
        [m,s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_test);
    case 'DT'
        P=X_train;
        DT=delaunayTriangulation(P);
        % triplot(DT)
        V=y_train;
        Pq=X_test;
        
        [ti,bc]=pointLocation(DT,Pq);
        triVals=V(DT(ti,:));
        Vq=dot(bc',triVals')';
        m=Vq;
end
figure,surf(X1_test,X2_test,reshape(m,length(x2_test),length(x1_test)),'FaceLighting','phong')
colormap(jet)
axis equal
xlim([1,128]); ylim([1,128]); zlim([-5,5]);
xlabel('x1'); ylabel('x2'); view([0 -90])
%% Error analysis
range=max(array(:))-min(array(:));
m_a=reshape(m,length(x2_test),length(x1_test));

max(abs(m_a(:)-y_test_ideal(:)));
PSNR=max(y_test_ideal(:))*sqrt(length(x1_test)*length(x2_test))/sqrt(sum((m_a(:)-y_test_ideal(:)).^2));
PSNR=20*log(PSNR)/log(10)





