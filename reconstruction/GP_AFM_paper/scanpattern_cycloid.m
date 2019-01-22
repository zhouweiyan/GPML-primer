% simulate cycloid scanning pattern
% Reference: Y. K. Yong, S. O. R. Moheimani, and I. R. Petersen, ¡°High-speed cycloid-scan atomic force microscopy,¡± 
% Nanotechnology, vol. 21, no. 36, p. 365503, Sep. 2010.
% zhouweiyan 20190111

%% initialization
clear
close all
clc
%% cycloid
f=1;
w=2*pi*f;
n=8;
r=5;
P=4*r/(2*n-1);
erfa=P*w/(2*pi);
t=0:0.01:8;
x=erfa*t+r*sin(w*t);
y=r*cos(w*t);
plot(x,y)
axis equal
% ylim([-6,6])
hold on
plot([-5 -5 5 5 -5]+5,[-5 5 5 -5 -5],'k--','LineWidth',2)
