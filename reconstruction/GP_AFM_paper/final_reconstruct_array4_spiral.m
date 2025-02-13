% simulate spiral scanning pattern
% Reference: I. A. Mahmood, S. O. R. Moheimani, and B. Bhikkaji, ��A New Scanning Method for Fast Atomic Force Microscopy,�� 
% IEEE Transactions on Nanotechnology, vol. 10, no. 2, pp. 203�C216, Mar. 2011.
% combine delaunayTriangulation interpolation and GPR
% zhouweiyan 20190127, 20190205
% done

clear
close all
% clc
load('array4_GP.mat')
% opt='GP';
opt='DT';
% set(0,'DefaultFigureWindowStyle','docked')  % keep figures tidy
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
set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);

%% Constant Angular Velocity(CAV)
r_end=64*1.6;
P=r_end*2/(60-1);
%P=spiral radius��2/(number of curves-1)
%eq.(4)
w=120;
f=2000;
T=1/f;
T_spiral=2*pi*r_end/(P*w);	%CAV eq.(6)
t=0:T:T_spiral;
xs=P*w/(2*pi)*t.*cos(w*t)+64;	%eq.(17)
ys=P*w/(2*pi)*t.*sin(w*t)+64;	%eq.(18)
% figure;
% plot(xs,ys,'.')
% axis equal
% % axis ([0 128 0 128])
% grid on
% hold on;plot(xs,ys)

%% data preparation
% train set
x1_train_spi=xs';x2_train_spi=ys';
X_train_spi=[x1_train_spi(:),x2_train_spi(:)];
[y_train_spi,y_t_s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_train_spi);

% figure; % visualization
% [x_temp,y_temp,z_temp]=griddata(x1_train_spi,x2_train_spi,y_train_spi,linspace(min(x1_train_spi),max(x1_train_spi))',linspace(min(x2_train_spi),max(x2_train_spi)),'v4');
% surf(x_temp,y_temp,z_temp)
% colormap(jet)
% hold on;plot(xs,ys,'.')
% axis equal
% axis ([0 128 0 128])
% plot(xs,ys)

% approximate path
path_len=sum(sqrt((x1_train_spi(2:end)-x1_train_spi(1:(end-1))).^2+(x2_train_spi(2:end)-x2_train_spi(1:(end-1))).^2));
delta = path_len/(2*128*128)

% test set
x1_test=1:128;x2_test=1:128;
[X1_test,X2_test]=meshgrid(x1_test,x2_test);
X_test=[X1_test(:),X2_test(:)];
[y_test_ideal,y_test_s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_test);
y_test_ideal=reshape(y_test_ideal,length(x2_test),length(x1_test));
%%
figure;
f=plot3(x1_train_spi(:),x2_train_spi(:),y_train_spi(:),'.');
grid on
axis equal
xlim([min(x1_train_spi),max(x1_train_spi)]);ylim([min(x2_train_spi),max(x2_train_spi)]);
zlim([-27,27]);
xlabel('x_1/\mum');ylabel('x_2/\mum');zlabel('y/nm');
% view(3)
view(-45,50)
x_label=cell(1,length(0:32:128));
for i=1:length(0:64:128)
    x_label{i}=64*(i-1)/128*2;
end
set(gca,'xtick',0:64:128,'xticklabel',x_label);
set(gca,'ytick',0:64:128,'yticklabel',x_label);
set(gca,'ztick',[-27,0,27],'zticklabel',[-20,0,20])
set(gca,'FontSize',18,'FontName', 'Times New Roman');
set_fig_units_cm(12,10)
set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);

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
        hyp=minimize(hyp,@gp,-80,@infGaussLik,meanfunc,covfunc,likfunc,X_train_spi,y_train_spi);
        [m,s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train_spi,y_train_spi,X_test);
        
    case 'DT'
        P=X_train_spi;
        DT=delaunayTriangulation(P);
        % triplot(DT)
        V=y_train_spi;
        % approximate path
        path_len=sum(sqrt((x1_train_spi(2:end)-x1_train_spi(1:(end-1))).^2+(x2_train_spi(2:end)-x2_train_spi(1:(end-1))).^2))
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
range=max(array(:))-min(array(:));
m_a=reshape(m,length(x2_test),length(x1_test));
m_a=m_a(min(x2_train):max(x2_train),min(x1_train):max(x1_train));
y_test_ideal_a=y_test_ideal(min(x2_train):max(x2_train),min(x1_train):max(x1_train));
abs(min(y_test_ideal_a(:)-m_a(:)));
abs(max(y_test_ideal_a(:)-m_a(:)));
PSNR=max(y_test_ideal_a(:))*sqrt(size(y_test_ideal_a,1)*size(y_test_ideal_a,2))/sqrt(sum((m_a(:)-y_test_ideal_a(:)).^2));
PSNR=20*log(PSNR)/log(10)

RMS = sqrt(1/(size(y_test_ideal_a,1)*size(y_test_ideal_a,2))*sum((m_a(:)-y_test_ideal_a(:)).^2))

sigma = sqrt(sum(s2(:))/length(s2(:)))
%%
figure;
surf(y_test_ideal_a-m_a,'EdgeColor','none','LineStyle','none','FaceLighting','phong');   % surf(y_test_ideal_a-m_a,'FaceAlpha',0.1); 
zlim([-50,40]);
% axis equal
colormap(jet)
view([-90 90])
axis tight
xlabel('x_1/\mum');ylabel('x_2/\mum');

set(gca, 'xtick', [], 'xticklabel', []);
set(gca, 'ytick', [], 'yticklabel', []);
set(gca,'FontSize',18,'FontName', 'Times New Roman');

% caxis([-10,10])
caxis([-2,2])
colorbar('eastoutside')

set_fig_units_cm(12,10)
set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1.5 0; 0 -1.2 0 1.5; 0 0 1 0; 0 0 0 1]);
set(gca,'FontSize',24,'FontName', 'Times New Roman');

