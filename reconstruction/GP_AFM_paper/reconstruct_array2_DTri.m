% delaunayTriangulation interpolation
% zhouweiyan 20180925
clear
close all
clc
load array2.mat
figure
surf(array);colormap(jet)
% view([0 -90])
axis equal
xlim([1,128]);ylim([1,128]);zlim([-50,40]);
xlabel('x1');ylabel('x2');view(3)

%% data preparation
% train set
x1_train=2:6:128; x2_train=1:1:128;
[X1_train,X2_train]=meshgrid(x1_train,x2_train);
X_train=[X1_train(:),X2_train(:)];
y_train=array(x2_train,x1_train);
y_train=y_train(:);
figure
surf(X1_train,X2_train,reshape(y_train,length(x2_train),length(x1_train)));colormap(jet)
axis equal
xlim([1,128]);ylim([1,128]);zlim([-50,40]);
xlabel('x1');ylabel('x2');view(3)

P=X_train;
DT=delaunayTriangulation(P);
% triplot(DT)
V=y_train;

% test set
x1_test=min(x1_train):max(x1_train);x2_test=min(x2_train):max(x2_train);
[X1_test,X2_test]=meshgrid(x1_test,x2_test);
X_test=[X1_test(:),X2_test(:)];
y_test_ideal=array(x2_test,x1_test);

Pq=X_test;
%% linear interplotation
[ti,bc]=pointLocation(DT,Pq);
triVals=V(DT(ti,:));
Vq=dot(bc',triVals')';
m=Vq;
figure
surf(X1_test,X2_test,reshape(m,length(x2_test),length(x1_test)))
colormap(jet)
axis equal
xlim([min(x1_train),max(x1_train)]);ylim([min(x2_train),max(x2_train)]);zlim([-50,40]);
xlabel('x1');ylabel('x2');
view(3)

%% Error analysis
range=max(array(:))-min(array(:));
% NRMSD=sqrt(sum(m-y_test_ideal(:)).^2)/range
% max(m-y_test_ideal(:))

m_a=reshape(m,length(x2_test),length(x1_test));

figure
surf(y_test_ideal-m_a);xlabel('x1');ylabel('x2');colormap(jet)
zlim([-50,40]);
NRMSD=sqrt(sum((m_a(:)-y_test_ideal(:)).^2))/range
max(abs(m_a(:)-y_test_ideal(:)))
