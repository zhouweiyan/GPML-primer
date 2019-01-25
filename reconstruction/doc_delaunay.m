% delaunayTriangulation learning
% https://ww2.mathworks.cn/help/matlab/interpolation.html
% https://ww2.mathworks.cn/help/matlab/math/interpolation-using-a-specific-delaunay-triangulation.html
% https://ww2.mathworks.cn/help/matlab/ref/delaunaytriangulation.html
% zhouweiyan 20180925

%% visualization
clear
close all
clc
P=-2.5+5*gallery('uniformdata',[9,2],0);    % 2D scatter
DT=delaunayTriangulation(P)
triplot(DT)
hold on
IC=incenter(DT)
figure,plot(IC(:,1),IC(:,2),'r^')
V=P(:,1).^2+P(:,2).^2;
Pq=-2+4*gallery('uniformdata',[3 2],0);
figure,plot(Pq(:,1),Pq(:,2),'*')
ground_truth=Pq(:,1).^2+Pq(:,2).^2   % ground truth
%% nearest neighbor interpolation
vi=nearestNeighbor(DT,Pq);
Vq=V(vi)
%% linear interpolation
[ti,bc]=pointLocation(DT,Pq);
triVals=V(DT(ti,:));
Vq=dot(bc',triVals')'