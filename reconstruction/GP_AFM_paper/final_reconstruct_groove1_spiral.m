% simulate spiral scanning pattern
% Reference: I. A. Mahmood, S. O. R. Moheimani, and B. Bhikkaji, ¡°A New Scanning Method for Fast Atomic Force Microscopy,¡± 
% IEEE Transactions on Nanotechnology, vol. 10, no. 2, pp. 203¨C216, Mar. 2011.
% combine delaunayTriangulation interpolation and GPR
% zhouweiyan 20180926
% done

clear
close all
clc
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
%% Constant Angular Velocity(CAV)
r_end=64*1.6;
P=r_end*2/(60-1);
%P=spiral radius¡Á2/(number of curves-1)
%eq.(4)
w=120;
f=2000;
T=1/f;
T_spiral=2*pi*r_end/(P*w);	%CAV eq.(6)
t=0:T:T_spiral;
xs=P*w/(2*pi)*t.*cos(w*t)+64;	%eq.(17)
ys=P*w/(2*pi)*t.*sin(w*t)+64;	%eq.(18)
figure
plot(xs,ys,'+')
axis equal
axis ([0 128 0 128])
grid on
hold on;plot(xs,ys)

%plot the outside dashed circle
% hold on
% xita=0:(pi/180):(2*pi);
% x=r_end*cos(xita)+64;
% y=r_end*sin(xita)+64;
% plot(x,y,'--')
%% data preparation
% train set
x1_train_spi=xs';x2_train_spi=ys';
X_train_spi=[x1_train_spi(:),x2_train_spi(:)];
[y_train_spi,y_t_s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_train_spi);
% visualization
figure
[x_temp,y_temp,z_temp]=griddata(x1_train_spi,x2_train_spi,y_train_spi,linspace(min(x1_train_spi),max(x1_train_spi))',linspace(min(x2_train_spi),max(x2_train_spi)),'v4');
surf(x_temp,y_temp,z_temp)
colormap(jet)
hold on;plot(xs,ys,'+')
axis equal
axis ([0 128 0 128])
plot(xs,ys)
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
        meanfunc=[];hyp.mean=[];
        % mask=[1,0];
        covfunc={'covSEard'};hyp.cov=log([10;10;20]);
        % covfunc={'covMask',{mask,{'covSEiso'}}};hyp.cov=log([2;40]);  %log([ell;sf])
        likfunc='likGauss';hyp.lik=log(0.1);    % log(sn)
        hyp=minimize(hyp,@gp,-80,@infGaussLik,meanfunc,covfunc,likfunc,X_train_spi,y_train_spi);
        [m,s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train_spi,y_train_spi,X_test);
        figure
        surf(X1_test,X2_test,reshape(m,length(x2_test),length(x1_test)))
        colormap(jet)
        axis equal
        xlim([1,128]);ylim([1,128]);zlim([-90,40]);
        xlabel('x1');ylabel('x2');
        view(3)
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
m_a=m_a(min(x2_train):max(x2_train),min(x1_train):max(x1_train));
y_test_ideal_a=y_test_ideal(min(x2_train):max(x2_train),min(x1_train):max(x1_train));
figure
surf(y_test_ideal_a-m_a);xlabel('x1');ylabel('x2');zlim([-90,40]);
colormap(jet)
NRMSD=sqrt(sum((m_a(:)-y_test_ideal_a(:)).^2))/range
max(abs(m_a(:)-y_test_ideal_a(:)))
PSNR=max(y_test_ideal_a(:))*sqrt(size(y_test_ideal_a,1)*size(y_test_ideal_a,2))/sqrt(sum((m_a(:)-y_test_ideal_a(:)).^2))

