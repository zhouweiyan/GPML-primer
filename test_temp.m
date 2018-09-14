clear
close all
clc
load groove1.txt
groove1=reshape(groove1,256,256);
figure
surf(groove1);colormap(jet)
% view([0 -90])
axis equal
xlim([1,256]);ylim([1,256]);zlim([-90,40]);
xlabel('x1');ylabel('x2');
% groove1_2=groove1(2:2:256,2:2:256);
% figure
% surf(2:2:256,2:2:256,groove1_2);colormap(jet)
groove1_4=groove1(4:4:256,4:4:256);
figure
surf(4:4:256,4:4:256,groove1_4);colormap(jet)
axis equal
zlim([-90,40]);
xlabel('x1');ylabel('x2');

meanfunc=[];
csa={'covSEard'};L=[100;100];sf=10;
covfunc={'covMask',{[0,1],csa}};
