% reconstruct complex morphology
% spiral
% zwy 20190131 

clear
close all
clc

addpath(genpath('E:\OneDrive - hnu.edu.cn\tools\matlabcourse\GPML_matlab\GPML-primer\reconstruction\GP_AFM_Pp_part3'))
% load EChemist_Corrosion_GP.mat
% load EMaterial_Bi_SN_alloy_216_GP.mat     %
% load EMaterial_F12H20_on_MoS2_128_GP.mat  %
% load EMaterial_KFM_GP.mat 
% load MaterialSci_AtomicLattice_GP.mat     %
load LifeSci_DiatomsScan_216_GP.mat       %
% load MaterialSci_AlSurfaceEtched_216_GP.mat

% load EChemist_LithiumTitanate_GP.mat
% load LifeSci_ThalassiosiraPseudonana_128_filter_GP.mat 
% load MaterialSci_ProbeOxidationonSiliconSubstrate_128_GP.mat  %

opt='GP';
% opt='DT'; 
set(0,'DefaultFigureWindowStyle','docked')
% figure, surf(f, 'EdgeColor', 'None');       
% colormap(jet); axis equal; 
% xlabel('x1'); ylabel('x2'); view([0 -90])
%% data preparation
% train set
r_end = 64*1.6;     % cover the square of 128*128 pixels
P = r_end*2/(64-1); % P = spiral radius¡Á2/(number of curves-1) eq.(4)
w = 120;
f = 5000;           % inversely proportional to the number of sampling points
T = 1/f;
T_spiral = 2*pi*r_end/(P*w)	%CAV eq.(6)
t = 0:T:T_spiral;
xs = P*w/(2*pi)*t.*cos(w*t)+64;	%eq.(17)
ys = P*w/(2*pi)*t.*sin(w*t)+64;	%eq.(18)
x1_train_spi=xs';x2_train_spi=ys';
X_train_spi=[x1_train_spi(:),x2_train_spi(:)];
[y_train_spi,y_t_s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_train_spi);
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
figure, surf(y_test_ideal, 'EdgeColor', 'None');        
colormap(jet); axis equal; axis off;
xlabel('x1'); ylabel('x2'); view([0 -90])
tightfig
set_fig_units_cm(10,10)
set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1.5; 0 0 1 0; 0 0 0 1]);
%% GP regression/DT linear interplotation
switch opt
    case 'GP'
        % EChemist_Corrosion_GP; EMaterial_KFM_GP; MaterialSci_AlSurfaceEtched_216_GP
%         meanfunc = []; hyp.mean = [];
%         covfunc = {'covSEard'}; hyp.cov = log([5;7;9]);

        % EMaterial_Bi_SN_alloy_216_GP
%         meanfunc = []; hyp.mean = [];
%         covfunc = {'covMaternard',3}; hyp.cov = log([4;5;6]);
        
        % LifeSci_DiatomsScan_216_GP
        meanfunc = []; hyp.mean = [];
%         covfunc = {'covSum',{'covSEard',{'covMaternard',3}}}; hyp.cov = log([4;6;8;4;6;8]);
        covfunc = {'covMaternard',3}; hyp.cov = log([4;4;45]);
        
        % ProbeOxidationonSiliconSubstrate_128_GP
%         meanfunc = []; hyp.mean = [];
%         hyp.cov = log([4;5;6]);  

        % LifeSci_ThalassiosiraPseudonana_128_GP
%         meanfunc = {'meanConst'}; hyp.mean = log(120);
%         covfunc = {'covMaternard',3}; hyp.cov = log([1.1;1.2;4]);
%         covfunc = {'covSEard'}; hyp.cov = log([5;7;9]);

        % EMaterial_F12H20_on_MoS2_128_GP
%         meanfunc = []; hyp.mean = [];
%         covfunc = {'covSum',{{'covPER','ard',{'covSEard'}},'covSEard'}}; 
%         hyp.cov = log([500;20;1;2;3;4;5;1.1;2.1;3.1]);

        % MaterialSci_AtomicLattice_GP
%         meanfunc = []; hyp.mean = [];
%         covfunc = {'covSEard'}; hyp.cov = log([1;2;3]);

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
%%
figure
surf(X1_test,X2_test,reshape(m,length(x2_test),length(x1_test)),'EdgeColor', 'None');
colormap(jet); axis equal; axis off;
xlabel('x1'); ylabel('x2'); view([0 -90])
tightfig
set_fig_units_cm(10,10)
set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1.5; 0 0 1 0; 0 0 0 1]);
%% Error analysis
m_a = reshape(m,length(x2_test),length(x1_test));
m_a = m_a(min(x2_train):max(x2_train),min(x1_train):max(x1_train));
y_test_ideal_a = y_test_ideal(min(x2_train):max(x2_train),min(x1_train):max(x1_train));
PSNR = max(y_test_ideal_a(:))*sqrt(size(y_test_ideal_a,1)*size(y_test_ideal_a,2))/sqrt(sum((m_a(:)-y_test_ideal_a(:)).^2));
PSNR = 20*log(PSNR)/log(10)
%%
figure
surf(y_test_ideal_a-m_a,'EdgeColor','none','LineStyle','none','FaceLighting','phong');   % surf(y_test_ideal_a-m_a,'FaceAlpha',0.1); 
 zlim([-90,40]); axis equal;axis off;
colormap(jet)
view([-90 90])
axis tight

max(abs(m_a(:)-y_test_ideal_a(:)))
mean(abs(m_a(:)-y_test_ideal_a(:)))
caxis([-10,10])
colorbar('eastoutside', 'Ticks', [-10, -5, 0, 5, 10], ...
    'TickLabels',[-10, -5, 0, 5, 10], 'FontSize', 24,...
    'FontName', 'Times New Roman')%, 'AxisLocation','in')
set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 4 0; 0 -1 0 1.7; 0 0 1 0; 0 0 0 1]);

