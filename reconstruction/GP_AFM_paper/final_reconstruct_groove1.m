% isotropic/unisotropic undersampling
% combine delaunayTriangulation interpolation and GPR
% zhouweiyan 20181019
% done
clear
% close all
clc
% opt='GP';
opt='DT';
load groove1.mat
figure
surf(groove);colormap(jet)
% view([0 -90])
axis equal
xlim([1,128]);ylim([1,128]);zlim([-90,40]);
xlabel('x1');ylabel('x2');view(-45,50)

%% data preparation
% train set
x1_train=2:6:128; x2_train=1:128;
[X1_train,X2_train]=meshgrid(x1_train,x2_train);
X_train=[X1_train(:),X2_train(:)];
y_train=groove(x2_train,x1_train);
y_train=y_train(:);
figure
surf(X1_train,X2_train,reshape(y_train,length(x2_train),length(x1_train)));colormap(jet)
axis equal
xlabel('x1');ylabel('x2');zlim([-90,40]);

% test set
x1_test=min(x1_train):max(x1_train);x2_test=min(x2_train):max(x2_train);
[X1_test,X2_test]=meshgrid(x1_test,x2_test);
X_test=[X1_test(:),X2_test(:)];
y_test_ideal=groove(x2_test,x1_test);

%% GP regression/DT linear interplotation
switch opt
    case 'GP'
        meanfunc=[];hyp.mean=[];
        % mask=[1,0];
        covfunc={'covSEard'};hyp.cov=log([10;10;20]);
        % covfunc={'covMask',{mask,{'covSEiso'}}};hyp.cov=log([2;40]);  %log([ell;sf])
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
%%
figure
surf(X1_test,X2_test,reshape(m,length(x2_test),length(x1_test)),'FaceLighting','phong')
colormap(jet)
axis equal
xlim([min(x1_train),max(x1_train)]);ylim([min(x2_train),max(x2_train)]);zlim([-90,40]);
xlabel('x_1/\mum');ylabel('x_2/\mum');zlabel('y/nm');
view(3)
x_label=cell(1,length(0:32:128));
for i=1:length(0:32:128)
    x_label{i}=32*(i-1)/128*2;
end
set(gca,'xtick',0:32:128);
set(gca,'xticklabel',x_label);
set(gca,'ytick',0:32:128);
set(gca,'yticklabel',x_label);
tightfig
set_fig_units_cm(10,10)
%% Error analysis
range=max(groove(:))-min(groove(:));
% NRMSD=sqrt(sum(m-y_test_ideal(:)).^2)/range
% max(m-y_test_ideal(:))

m_a=reshape(m,length(x2_test),length(x1_test));
% m_a=m_a(1:256,4:252);y_test_ideal_a=y_test_ideal(1:256,4:252);    % cut
% the boundary of x1

NRMSD=sqrt(sum((m_a(:)-y_test_ideal(:)).^2))/range
max(abs(m_a(:)-y_test_ideal(:)))
PSNR=max(y_test_ideal(:))*sqrt(length(x1_test)*length(x2_test))/sqrt(sum((m_a(:)-y_test_ideal(:)).^2))
PSNR=20*log(PSNR)/log(10)
%%
figure
surf(y_test_ideal-m_a,'EdgeColor','none','LineStyle','none','FaceLighting','phong');   % surf(y_test_ideal_a-m_a,'FaceAlpha',0.1); 
zlim([-90,40]);
axis equal
colormap(jet)
view([-90 90])
axis tight
xlabel('x_1/\mum');ylabel('x_2/\mum');
x_label=cell(1,length(0:32:128));
y_label=cell(1,length(0:32:128));
for i=1:length(0:32:128)
    x_label{i}=32*(i-1)/128*2;
end
for i=1:length(0:32:128)
    y_label{length(0:32:128)+1-i}=32*(i-1)/128*2;
end
set(gca,'xtick',0:32:128);
set(gca,'xticklabel',x_label);
set(gca,'ytick',0:32:128);
set(gca,'yticklabel',y_label);
% tightfig

% caxis([-10,10])
caxis([-5,5])
colorbar('eastoutside')
% set(gca,'CLim',[-10,10])