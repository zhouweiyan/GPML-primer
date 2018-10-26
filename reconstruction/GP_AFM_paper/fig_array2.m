% zhouweiyan 20180928
% done
clear
% close all
% clc

load array2.mat
figure
surf(array);
% view(3)
view(-45,50)

% surf(array,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
% view([0 -90])
% axis off

colormap(jet)
axis equal
xlim([1,128]);ylim([1,128]);zlim([-50,40]);
xlabel('x_1/\mum');ylabel('x_2/\mum');zlabel('y/nm');

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