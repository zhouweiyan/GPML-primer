% isotropic/unisotropic undersampling
% x1 change horizontally, x2 change vertically,e.g.(x2,x1)
% (1,1) (1,2) (1,3)... 
% (2,1) (2,2) (2,3)...
% zhouweiyan 20180915
clear
% close all
% clc
load groove1.txt
groove=reshape(groove1,256,256);
figure
surf(groove);colormap(jet)
% view([0 -90])
axis equal
xlim([1,256]);ylim([1,256]);zlim([-90,40]);
xlabel('x1');ylabel('x2');

%% data preparation
% train set
x1_train=2:4:256; x2_train=2:4:256;
[X1_train,X2_train]=meshgrid(x1_train,x2_train);
X_train=[X1_train(:),X2_train(:)];
y_train=groove(x2_train,x1_train);
y_train=y_train(:);
figure
surf(X1_train,X2_train,reshape(y_train,length(x2_train),length(x1_train)));colormap(jet)
axis equal
xlabel('x1');ylabel('x2');zlim([-90,40]);

% test set
x1_test=2:254;x2_test=2:254;
[X1_test,X2_test]=meshgrid(x1_test,x2_test);
X_test=[X1_test(:),X2_test(:)];
y_test_ideal=groove(x2_test,x1_test);
%% bicubic interplotation
m=interp2(X1_train,X2_train,groove(x2_train,x1_train),X1_test,X2_test,'cubic');
figure
surf(X1_test,X2_test,reshape(m,length(x2_test),length(x1_test)))
colormap(jet)
axis equal
xlim([min(x1_train),max(x1_train)]);ylim([min(x2_train),max(x2_train)]);zlim([-90,40]);
xlabel('x1');ylabel('x2');
view(3)
%% Error analysis
range=max(groove(:))-min(groove(:));
% NRMSD=sqrt(sum(m-y_test_ideal(:)).^2)/range
% max(m-y_test_ideal(:))

m_a=reshape(m,length(x2_test),length(x1_test));
% m_a=m_a(1:256,4:252);y_test_ideal_a=y_test_ideal(1:256,4:252);    % cut
% the boundary of x1

figure
surf(y_test_ideal-m_a);   % surf(y_test_ideal_a-m_a,'FaceAlpha',0.1); 
zlim([-90,40]);
NRMSD=sqrt(sum((m_a(:)-y_test_ideal(:)).^2))/range
max(abs(m_a(:)-y_test_ideal(:)))
