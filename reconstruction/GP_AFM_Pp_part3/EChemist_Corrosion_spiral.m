% reconstruct complex morphology
% spiral
% zwy 20190131

clear
close all
clc
load EChemist_Corrosion_GP.mat
% opt='GP';
opt='DT'; 
figure, surf(f, 'EdgeColor', 'None'); 
colormap(jet); axis equal; 
xlabel('x1'); ylabel('x2'); view([0 -90])
%% data preparation
% train set
r_end = 64*1.6;     % cover the square of 128*128 pixels
P = r_end*2/(64-1); % P = spiral radius¡Á2/(number of curves-1) eq.(4)
w = 120;
f = 3000;           % inversely proportional to the number of sampling points
T = 1/f;
T_spiral = 2*pi*r_end/(P*w);	%CAV eq.(6)
t = 0:T:T_spiral;
xs = P*w/(2*pi)*t.*cos(w*t)+64;	%eq.(17)
ys = P*w/(2*pi)*t.*sin(w*t)+64;	%eq.(18)
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
%% GP regression/DT linear interplotation
switch opt
    case 'GP'
        meanfunc = []; hyp.mean = [];
        covfunc = {'covSEard'}; hyp.cov = log([5;7;9]);
        likfunc = {'likGauss'}; hyp.lik = log(0.1);
        hyp=minimize(hyp,@gp,-80,@infGaussLik,meanfunc,covfunc,likfunc,X_train_spi,y_train_spi);
        [m,s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train_spi,y_train_spi,X_test);
    case 'DT'
        P=X_train_spi;
        DT=delaunayTriangulation(P);
        % triplot(DT)
        V=y_train_spi;
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
m_a=reshape(m,length(x2_test),length(x1_test));
m_a=m_a(min(x2_train):max(x2_train),min(x1_train):max(x1_train));
y_test_ideal_a=y_test_ideal(min(x2_train):max(x2_train),min(x1_train):max(x1_train));
PSNR=max(y_test_ideal_a(:))*sqrt(size(y_test_ideal_a,1)*size(y_test_ideal_a,2))/sqrt(sum((m_a(:)-y_test_ideal_a(:)).^2));
PSNR=20*log(PSNR)/log(10)





