% zhouweiyan 20180928
% done
clear
% close all
% clc

load groove1.mat
figure
% surf(groove);
% view(3)

surf(groove,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
view([0 -90])
axis off

colormap(jet)
axis equal
xlim([1,128]);ylim([1,128]);zlim([-90,40]);
xlabel('x1');ylabel('x2');
tightfig
set_fig_units_cm(10,10)