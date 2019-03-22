% reconstruct complex morphology
% subline
% zwy 20190131

clear
% close all
% clc

% individual example
% Obj = imread('Corrosion.bmp');
f = rgb2gray(double(Obj)/255)*255; f = f(373:500,1:128);
% Obj = imread('LithiumTitanate.bmp'); f = rgb2gray(double(Obj)/255)*255; f = f(1:128,1:128);
% Obj = imread('Bi_SN_alloy_216.bmp'); f = rgb2gray(double(Obj)/255)*255; f = f(end-128+1:end,1:128);       
% Obj = imread('F12H20 on MoS2_128.bmp'); f = rgb2gray(double(Obj)/255)*255; f = f(end-128+1:end,1:128);
% Obj = imread('KFM.bmp'); f = rgb2gray(double(Obj)/255)*255; f = f(160-128+1:160,300-128+1:300);
% Obj = imread('AtomicLattice.bmp'); f = rgb2gray(double(Obj)/255)*255; f = f(1:128,1:128);
% Obj = imread('DiatomsScan_216.bmp'); f = rgb2gray(double(Obj)/255)*255; f = f(51:50+128,81:80+128);
% Obj = imread('ThalassiosiraPseudonana_128.bmp'); f = rgb2gray(double(Obj)/255)*255;
Obj = imread('ThalassiosiraPseudonana_128_filter.bmp'); f = rgb2gray(double(Obj)/255)*255; h = ones(3,3)/9; f = imfilter(f, h, 'replicate');
% Obj = imread('ProbeOxidationonSiliconSubstrate_128.bmp'); f = rgb2gray(double(Obj)/255)*255; 
% Obj = imread('AlSurfaceEtched_216.bmp'); f = rgb2gray(double(Obj)/255)*255; f = f(1:128,end-128+1:end);  


opt='GP';
% opt='DT';
set(0,'DefaultFigureWindowStyle','docked') 
figure, surf(f, 'EdgeColor', 'None'); 
colormap(jet); axis equal; 
xlabel('x1'); ylabel('x2'); view([0 -90])
%% data preparation
% train set
x1_train = 1:2:128; x2_train = 1:2:128;
[X1_train,X2_train]=meshgrid(x1_train,x2_train);
X_train=[X1_train(:),X2_train(:)];
y_train=f(x2_train,x1_train); y_train=y_train(:);
figure
surf(X1_train,X2_train,reshape(y_train,length(x2_train),length(x1_train)),...
    'EdgeColor', 'None'); 
colormap(jet); axis equal; 
xlabel('x1'); ylabel('x2'); view([0 -90])
% test set
x1_test=min(x1_train):max(x1_train);x2_test=min(x2_train):max(x2_train);
[X1_test,X2_test]=meshgrid(x1_test,x2_test);
X_test=[X1_test(:),X2_test(:)];
y_test_ideal=f(x2_test,x1_test);
%% GP regression/DT linear interplotation
switch opt
    case 'GP'
        meanfunc = []; hyp.mean = [];
        covfunc = {'covSEard'}; hyp.cov = log([5;7;9]);        
        % ThalassiosiraPseudonana_128
        meanfunc = {'meanConst'}; hyp.mean = log(120);
        covfunc = {'covMaternard',3}; hyp.cov = log([2;3;4]);

        % Bi_SN_alloy_216; ProbeOxidationonSiliconSubstrate_128
%         covfunc = {'covMaternard',3}; hyp.cov = log([2;3;4]);        
        
        likfunc = {'likGauss'}; hyp.lik = log(0.1);
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
figure
surf(X1_test,X2_test,reshape(m,length(x2_test),length(x1_test)),'EdgeColor', 'None');
colormap(jet); axis equal;
xlabel('x1'); ylabel('x2'); view([0 -90])
%% Error analysis
range=max(f(:))-min(f(:));
m_a=reshape(m,length(x2_test),length(x1_test));
max(abs(m_a(:)-y_test_ideal(:)));
PSNR=max(y_test_ideal(:))*sqrt(length(x1_test)*length(x2_test))/sqrt(sum((m_a(:)-y_test_ideal(:)).^2));
PSNR=20*log(PSNR)/log(10)



