% simulate Lissajous scanning pattern
% Reference: T. Tuma, J. Lygeros, V. Kartik, A. Sebastian, and A. Pantazi, ��High-speed multiresolution scanning probe microscopy based on Lissajous scan trajectories,�� 
% Nanotechnology, vol. 23, no. 18, p. 185501, May 2012.
% https://ww2.mathworks.cn/help/matlab/ref/unique.html?searchHighlight=unique&s_tid=doc_srchtitle
% https://ww2.mathworks.cn/help/matlab/ref/uniquetol.html
% zhouweiyan 20190127
% 

clear
% close all
clc
load('array4_GP.mat')
opt='GP';
% opt='DT';
% set(0,'DefaultFigureWindowStyle','docked') 
%%
figure;     % truth
surf(array,'EdgeColor','none');
colormap(jet)
axis equal
zlim([-20,20]);
xlabel('x_1/\mum');ylabel('x_2/\mum');zlabel('y/nm');
% view(3)
view(-45,50)
x_label=cell(1,length(0:32:128));
for i=1:length(0:32:128)
    x_label{i}=32*(i-1)/128*2;
end
set(gca,'xtick',0:32:128,'xticklabel',x_label);
set(gca,'ytick',0:32:128,'yticklabel',x_label);
set(gca,'ztick',[-20,0,20],'zticklabel',[-20,0,20])
set(gca,'FontSize',18,'FontName', 'Times New Roman');
set_fig_units_cm(12,10)
set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1.5 0; 0 0 0 1]);

%% Lissajous scan trajectory
X0=64;Y0=64;
A=128;B=128;    % A,B represent the range of X,Y respectively
fx=41;fy=43;
f=5000;
T=1;  %T = lcm(fx,fy)/(fx*fy)
t=(0:1/f:T)';
x1_train_lisaj=X0+A/2*sin(2*pi*fx*t);
x2_train_lisaj=Y0+B/2*sin(2*pi*fy*t);
% figure;
% plot(x1_train_lisaj,x2_train_lisaj,'.')
% axis equal
% axis ([0 128 0 128])
% grid on
% hold on;plot(x1_train_lisaj,x2_train_lisaj)

%% data preparation
% train set
X_train_lisaj=[x1_train_lisaj(:),x2_train_lisaj(:)];
[y_train_lisaj,y_t_s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_train_lisaj);

% figure; % visualization
% [x_temp,y_temp,z_temp]=griddata(x1_train_lisaj,x2_train_lisaj,y_train_lisaj,linspace(min(x1_train_lisaj),max(x1_train_lisaj))',linspace(min(x2_train_lisaj),max(x2_train_lisaj)),'v4');
% surf(x_temp,y_temp,z_temp)
% colormap(jet)
% hold on;plot(x1_train_lisaj,x2_train_lisaj,'.')
% axis equal
% axis ([0 128 0 128])
% plot(x1_train_lisaj,x2_train_lisaj)

% approximate path
path_len=sum(sqrt((x1_train_lisaj(2:end)-x1_train_lisaj(1:(end-1))).^2+(x2_train_lisaj(2:end)-x2_train_lisaj(1:(end-1))).^2));
delta = path_len/(2*128*128)

% test set
x1_test=2:127;x2_test=2:127;
[X1_test,X2_test]=meshgrid(x1_test,x2_test);
X_test=[X1_test(:),X2_test(:)];
[y_test_ideal,y_test_s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_test);
y_test_ideal=reshape(y_test_ideal,length(x2_test),length(x1_test));

[X_train_lisaj,ia,ic]=uniquetol(X_train_lisaj,'ByRows',true);
y_train_lisaj=y_train_lisaj(ia);
%%
figure;
f=plot3(X_train_lisaj(:,1),X_train_lisaj(:,2),y_train_lisaj(:),'.');
grid on
axis equal
xlim([min(x1_train_lisaj),max(x1_train_lisaj)]);ylim([min(x2_train_lisaj),max(x2_train_lisaj)]);
zlim([-20,20]);
xlabel('x_1/\mum');ylabel('x_2/\mum');zlabel('y/nm');
% view(3)
view(-45,50)
x_label=cell(1,length(0:32:128));
for i=1:length(0:64:128)
    x_label{i}=64*(i-1)/128*2;
end
set(gca,'xtick',0:64:128,'xticklabel',x_label);
set(gca,'ytick',0:64:128,'yticklabel',x_label);
set(gca,'ztick',[-20,0,20],'zticklabel',[-20,0,20])
set(gca,'FontSize',18,'FontName', 'Times New Roman');
set_fig_units_cm(12,10)
set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1.5 0; 0 0 0 1]);

%% GP regression/DT linear interplotation
switch opt
    case 'GP'
        meanfunc=[];hyp.mean=[];
        mask = [0,1];
        cgi = {'covSEiso'};
        cma = {'covMask',{mask,cgi{:}}}; hypma = log([1.5;15]);%log([ell;sf]);
        covfunc1={'covSum',{cma,'covSEiso'}};hyp1.cov=[hypma;log([2;3])];%relative good
        covfunc2={'covMaternard',1};hyp2.cov=log([1;1;4]);
        covfunc={'covSum',{covfunc1,covfunc2}};hyp.cov=[hyp1.cov;hyp2.cov];
        
        likfunc='likGauss';hyp.lik=log(10);    % log(sn)
        hyp=minimize(hyp,@gp,-80,@infGaussLik,meanfunc,covfunc,likfunc,X_train_lisaj,y_train_lisaj);
        [m,s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train_lisaj,y_train_lisaj,X_test);
        
    case 'DT'
        P=X_train_lisaj;
        DT=delaunayTriangulation(P);
        % triplot(DT)
        V=y_train_lisaj;
        Pq=X_test;
        
        [ti,bc]=pointLocation(DT,Pq);
        triVals=V(DT(ti,:));
        Vq=dot(bc',triVals')';
        m=Vq;
end
%%
figure;
surf(X1_test,X2_test,reshape(m,length(x2_test),length(x1_test)),'FaceLighting','phong')
colormap(jet)
axis equal
xlim([min(x1_train),max(x1_train)]);ylim([min(x2_train),max(x2_train)]);zlim([-20,20]);
xlabel('x_1/\mum');ylabel('x_2/\mum');zlabel('y/nm');
% view(3)
view(-45,50)
x_label=cell(1,length(0:32:128));
for i=1:length(0:32:128)
    x_label{i}=32*(i-1)/128*2;
end
set(gca,'xtick',0:32:128,'xticklabel',x_label);
set(gca,'ytick',0:32:128,'yticklabel',x_label);
set(gca,'ztick',[-20,0,20],'zticklabel',[-20,0,20])
set(gca,'FontSize',18,'FontName', 'Times New Roman');
set_fig_units_cm(12,10)
set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);

%% Error analysis
range=max(array(:))-min(array(:))
% NRMSD=sqrt(sum(m-y_test_ideal(:)).^2)/range
% max(m-y_test_ideal(:))

m_a=reshape(m,length(x2_test),length(x1_test));

NRMSD=sqrt(sum((m_a(:)-y_test_ideal(:)).^2))/range;
abs(min(y_test_ideal(:)-m_a(:)))
abs(max(y_test_ideal(:)-m_a(:)))
PSNR=max(y_test_ideal(:))*sqrt(size(y_test_ideal,1)*size(y_test_ideal,2))/sqrt(sum((m_a(:)-y_test_ideal(:)).^2));
PSNR=20*log(PSNR)/log(10)

abs(max(y_test_ideal(:)-m_a(:))-min(y_test_ideal(:)-m_a(:)))
RMS = sqrt(1/(size(y_test_ideal,1)*size(y_test_ideal,2)))*sqrt(sum((m_a(:)-y_test_ideal(:)).^2))

%%
figure;
surf(y_test_ideal-m_a,'EdgeColor','none','LineStyle','none','FaceLighting','phong');   % surf(y_test_ideal_a-m_a,'FaceAlpha',0.1); 
zlim([-20,20]);
colormap(jet)
view([-90 90])
axis tight
xlabel('x_1/\mum');ylabel('x_2/\mum');

set(gca, 'xtick', [], 'xticklabel', []);
set(gca, 'ytick', [], 'yticklabel', []);
set(gca,'FontSize',18,'FontName', 'Times New Roman');

caxis([-2,2])
colorbar('eastoutside')

set_fig_units_cm(12,10)
set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1.2 0 1.5 0; 0 -1.2 0 1.5; 0 0 1 0; 0 0 0 1]);
set(gca,'FontSize',24,'FontName', 'Times New Roman');

