% simulate Lissajous scanning pattern
% Reference: T. Tuma, J. Lygeros, V. Kartik, A. Sebastian, and A. Pantazi, ¡°High-speed multiresolution scanning probe microscopy based on Lissajous scan trajectories,¡± 
% Nanotechnology, vol. 23, no. 18, p. 185501, May 2012.
% https://ww2.mathworks.cn/help/matlab/ref/unique.html?searchHighlight=unique&s_tid=doc_srchtitle
% https://ww2.mathworks.cn/help/matlab/ref/uniquetol.html
% zhouweiyan 20180926
% done

clear
close all
% clc
opt='GP';
% opt='DT';
load('groove1_GP.mat')
set(0,'DefaultFigureWindowStyle','docked') 
figure
surf(groove);colormap(jet)
% view([0 -90])
axis equal
xlim([1,128]);ylim([1,128]);zlim([-90,40]);
xlabel('x1');ylabel('x2');view(3)

%% Lissajous scan trajectory
X0=64;Y0=64;
A=128;B=128;    % A,B represent the range of X,Y respectively
fx=41;fy=43;
f=4000;
T=1;  %T = lsm(fx,fy)/(fx*fy)
t=(0:1/f:T)';
x1_train_lisaj=X0+A/2*sin(2*pi*fx*t);
x2_train_lisaj=Y0+B/2*sin(2*pi*fy*t);
figure
plot(x1_train_lisaj,x2_train_lisaj,'+')
axis equal
axis ([0 128 0 128])
grid on
hold on;plot(x1_train_lisaj,x2_train_lisaj)

%% data preparation
% train set
X_train_lisaj=[x1_train_lisaj(:),x2_train_lisaj(:)];
[y_train_lisaj,y_t_s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_train_lisaj);
% visualization
figure
[x_temp,y_temp,z_temp]=griddata(x1_train_lisaj,x2_train_lisaj,y_train_lisaj,linspace(min(x1_train_lisaj),max(x1_train_lisaj))',linspace(min(x2_train_lisaj),max(x2_train_lisaj)),'v4');
surf(x_temp,y_temp,z_temp)
colormap(jet)
hold on;plot(x1_train_lisaj,x2_train_lisaj,'+')
axis equal
axis ([0 128 0 128])
plot(x1_train_lisaj,x2_train_lisaj)
% approximate path
path_len=sum(sqrt((x1_train_lisaj(2:end)-x1_train_lisaj(1:(end-1))).^2+(x2_train_lisaj(2:end)-x2_train_lisaj(1:(end-1))).^2))

% test set
x1_test=2:127;x2_test=2:127;
[X1_test,X2_test]=meshgrid(x1_test,x2_test);
X_test=[X1_test(:),X2_test(:)];
[y_test_ideal,y_test_s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_test);
y_test_ideal=reshape(y_test_ideal,length(x2_test),length(x1_test));
%% GP regression/DT linear interplotation
switch opt
    case 'GP'
        meanfunc=[];hyp.mean=[];
        % mask=[1,0];
        covfunc={'covSEard'};hyp.cov=log([10;10;20]);
        % covfunc={'covMask',{mask,{'covSEiso'}}};hyp.cov=log([2;40]);  %log([ell;sf])
        likfunc='likGauss';hyp.lik=log(0.1);    % log(sn)
        hyp=minimize(hyp,@gp,-80,@infGaussLik,meanfunc,covfunc,likfunc,X_train_lisaj,y_train_lisaj);
        [m,s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train_lisaj,y_train_lisaj,X_test);
        figure
        surf(X1_test,X2_test,reshape(m,length(x2_test),length(x1_test)))
        colormap(jet)
        axis equal
        xlim([1,128]);ylim([1,128]);zlim([-90,40]);
        xlabel('x1');ylabel('x2');
        view(3)
    case 'DT'
%         [X_train_lisaj_dt,ia,ic]=uniquetol(X_train_lisaj,'ByRows',true);
%         y_train_lisaj_dt=y_train_lisaj(ia);
%         P=X_train_lisaj_dt;
%         DT=delaunayTriangulation(P);
%         % triplot(DT)
%         V=y_train_lisaj_dt;
        P=X_train_lisaj;
        DT=delaunayTriangulation(P);
        % triplot(DT)
        V=y_train_lisaj;
        Pq=X_test;
        
        [ti,bc]=pointLocation(DT,Pq);
        triVals=V(DT(ti,:));
        Vq=dot(bc',triVals')';
        m=Vq;
        figure
        surf(X1_test,X2_test,reshape(m,length(x2_test),length(x1_test)))
        colormap(jet)
        axis equal
        xlim([min(x1_train),max(x1_train)]);ylim([min(x2_train),max(x2_train)]);zlim([-90,40]);
        xlabel('x1');ylabel('x2');
        view(3)
end
%% Error analysis
range=max(groove(:))-min(groove(:))
% NRMSD=sqrt(sum(m-y_test_ideal(:)).^2)/range
% max(m-y_test_ideal(:))

m_a=reshape(m,length(x2_test),length(x1_test));

figure
surf(y_test_ideal-m_a);xlabel('x1');ylabel('x2');zlim([-90,40]);
colormap(jet)
NRMSD=sqrt(sum((m_a(:)-y_test_ideal(:)).^2))/range
max(abs(m_a(:)-y_test_ideal(:)))
PSNR=max(y_test_ideal(:))*sqrt(size(y_test_ideal,1)*size(y_test_ideal,2))/sqrt(sum((m_a(:)-y_test_ideal(:)).^2))
PSNR=20*log(PSNR)/log(10)


