% reconstruct array4
% zwy 20190122
clear
close all
clc
load array4.txt
array=reshape(array4,256,256);
set(0,'DefaultFigureWindowStyle','docked'); % keep figures tidy


figure
surf(array);colormap(jet)
view([0 -90])
axis equal
xlim([1,256]);ylim([1,256]);
xlabel('x1');ylabel('x2');
array=array(1:128,129:256);
figure
surf(array);colormap(jet)
axis equal
xlim([1,128]); ylim([1,128]); zlim([-5,5]);
xlabel('x1'); ylabel('x2'); view([0 -90])

