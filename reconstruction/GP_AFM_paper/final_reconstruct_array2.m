% isotropic/unisotropic undersampling
% combine delaunayTriangulation interpolation and GPR
% zhouweiyan 20180926
% done

clear
close all
clc
% opt='GP';
opt='DT';
addpath(genpath('E:\OneDrive - hnu.edu.cn\tools\matlabcourse\GPML_matlab\GPML-primer\reconstruction\GP_AFM_paper'))
load array2.mat
%%
figure;     % truth
surf(array,'EdgeColor','none');
colormap(jet)
axis equal
zlim([-40,40]);
xlabel('x_1/\mum');ylabel('x_2/\mum');zlabel('y/nm');
% view(3)
view(-45,50)
x_label=cell(1,length(0:32:128));
for i=1:length(0:32:128)
    x_label{i}=32*(i-1)/128*2;
end
set(gca,'xtick',0:32:128,'xticklabel',x_label);
set(gca,'ytick',0:32:128,'yticklabel',x_label);
set(gca,'ztick',[-40,0,40],'zticklabel',[-40,0,40])
set(gca,'FontSize',18,'FontName', 'Times New Roman');
set_fig_units_cm(12,10)
set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);

%% data preparation
% train set
x1_train=2:4:128; x2_train=1:1:128;
[X1_train,X2_train]=meshgrid(x1_train,x2_train);
X_train=[X1_train(:),X2_train(:)];
y_train=array(x2_train,x1_train);
y_train=y_train(:);
% figure;   % visualization
% surf(X1_train,X2_train,reshape(y_train,length(x2_train),length(x1_train)));
% colormap(jet)
% axis equal
% xlim([1,128]);ylim([1,128]);zlim([-50,40]);
% xlabel('x1');ylabel('x2');view(3)

% test set
x1_test=min(x1_train):max(x1_train);x2_test=min(x2_train):max(x2_train);
[X1_test,X2_test]=meshgrid(x1_test,x2_test);
X_test=[X1_test(:),X2_test(:)];
y_test_ideal=array(x2_test,x1_test);
%% sampling simulation
figure;
f=plot3(X1_train(:),X2_train(:),y_train(:),'.');
grid on
axis equal
xlim([min(x1_train),max(x1_train)]);ylim([min(x2_train),max(x2_train)]);zlim([-50,40]);
xlabel('x_1/\mum');ylabel('x_2/\mum');zlabel('y/nm');
% view(3)
view(-45,50)
x_label=cell(1,length(0:32:128));
for i=1:length(0:32:128)
    x_label{i}=32*(i-1)/128*2;
end
set(gca,'xtick',0:32:128,'xticklabel',x_label);
set(gca,'ytick',0:32:128,'yticklabel',x_label);
set(gca,'ztick',[-40,0,40],'zticklabel',[-40,0,40])
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
%% reconstruction results
figure
surf(X1_test,X2_test,reshape(m,length(x2_test),length(x1_test)),'FaceLighting','phong','EdgeColor','none');
colormap(jet)
axis equal
xlim([min(x1_train),max(x1_train)]);ylim([min(x2_train),max(x2_train)]);zlim([-50,40]);
xlabel('x_1/\mum');ylabel('x_2/\mum');zlabel('y/nm');
% view(3)
view(-45,50)
x_label=cell(1,length(0:32:128));
for i=1:length(0:32:128)
    x_label{i}=32*(i-1)/128*2;
end
set(gca,'xtick',0:32:128,'xticklabel',x_label);
set(gca,'ytick',0:32:128,'yticklabel',x_label);
set(gca,'ztick',[-40,0,40],'zticklabel',[-40,0,40])
set(gca,'FontSize',18,'FontName', 'Times New Roman');
set_fig_units_cm(12,10)
set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);

%% Error analysis
range=max(array(:))-min(array(:));
% NRMSD=sqrt(sum(m-y_test_ideal(:)).^2)/range
% max(m-y_test_ideal(:))

m_a=reshape(m,length(x2_test),length(x1_test));

abs(min(y_test_ideal(:)-m_a(:)))
abs(max(y_test_ideal(:)-m_a(:)))

NRMSD=sqrt(sum((m_a(:)-y_test_ideal(:)).^2))/range;
abs(max(m_a(:)-y_test_ideal(:)))
abs(min(m_a(:)-y_test_ideal(:)))

PSNR=max(y_test_ideal(:))*sqrt(length(x1_test)*length(x2_test))/sqrt(sum((m_a(:)-y_test_ideal(:)).^2));
PSNR=20*log(PSNR)/log(10)

RMS = sqrt(1/(size(y_test_ideal,1)*size(y_test_ideal,2))*sum((m_a(:)-y_test_ideal(:)).^2))

% sigma = sqrt(sum(s2(:))/length(s2(:)))
%%
figure
surf(y_test_ideal-m_a,'EdgeColor','none','LineStyle','none','FaceLighting','phong');   % surf(y_test_ideal_a-m_a,'FaceAlpha',0.1); 
zlim([-50,40]);
colormap(jet)
view([-90 90])
axis tight
xlabel('x_1/\mum');ylabel('x_2/\mum');
% x_label=cell(1,length(0:32:128));
% y_label=cell(1,length(0:32:128));
% for i=1:length(0:32:128)
%     x_label{i}=32*(i-1)/128*2;
% end
% for i=1:length(0:32:128)
%     y_label{length(0:32:128)+1-i}=32*(i-1)/128*2;
% end
% set(gca,'xtick',0:32:128);
% set(gca,'xticklabel',x_label);
% set(gca,'ytick',0:32:128);
% set(gca,'yticklabel',y_label);
set(gca, 'xtick', [], 'xticklabel', []);
set(gca, 'ytick', [], 'yticklabel', []);
set(gca,'FontSize',18,'FontName', 'Times New Roman');

% caxis([-10,10])
caxis([-5,5])
colorbar('eastoutside')
% set(gca,'CLim',[-10,10])

set_fig_units_cm(12,10)
set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1.5 0; 0 -1.2 0 1.5; 0 0 1 0; 0 0 0 1]);
set(gca,'FontSize',24,'FontName', 'Times New Roman');
