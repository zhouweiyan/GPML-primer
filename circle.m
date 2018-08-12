% function []=circle(x,y,r)
% zhouweiyan 20180811
% plot 2D circles, (x,y) is the center of a circle
% r is radius

rectangle('Position',[x-r,y-r,2*r,2*r],'Curvature',[1,1])
grid on
axis equal