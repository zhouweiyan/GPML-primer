clear
% close all
clc
%% Lissajous scan trajectory
X0=64;Y0=64;
A=128;B=128;    % A,B represent the range of X,Y respectively
% fx=41;fy=43;
fx=41;fy=42;
f=2000;
T=1;  %T = lcm(fx,fy)/(fx*fy)
t=(0:1/f:T)';
x=X0+A/2*sin(2*pi*fx*t);
y=Y0+B/2*sin(2*pi*fy*t);

path_len=sum(sqrt((x(2:end)-x(1:(end-1))).^2+(y(2:end)-y(1:(end-1))).^2));
delta = path_len/(2*128*128)

set(0,'DefaultFigureWindowStyle','docked') 
figure
% plot(x,y,'--')
hold on
axis equal;xlim([X0-A/2,X0+A/2]);ylim([Y0-B/2,Y0+B/2])
for i=1:length(t) % show the dynamic process 
    plot(x(i),y(i),'b.')
%     if i~=length(t)
%         plot([x(i),x(i+1)],[y(i),y(i+1)],'--');
%     end
    pause(0.01)
end
plot([0 0 128 128 0],[0 128 128 0 0],'k--','LineWidth',2)
axis off
grid off
% tightfig
set_fig_units_cm(10,10)