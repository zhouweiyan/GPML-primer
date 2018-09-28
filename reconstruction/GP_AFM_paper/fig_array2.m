% zhouweiyan 20180928
% done
clear
% close all
% clc

load array2.mat
figure
surf(array);
view(3)

% surf(array,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
% view([0 -90])
% axis off

colormap(jet)
axis equal
xlim([1,128]);ylim([1,128]);zlim([-50,40]);
xlabel('x1');ylabel('x2');
tightfig
set_fig_units_cm(14,12)