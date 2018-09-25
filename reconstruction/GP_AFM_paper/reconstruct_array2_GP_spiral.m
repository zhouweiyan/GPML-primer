% zhouweiyan 20180925
clear
close all
clc
%% Constant Angular Velocity(CAV)
r_end=64;
P=r_end*2/(32-1);
%P=spiral radius¡Á2/(number of curves?1)
%eq.(4)
w=188.5;
f=2000;
T=1/f;
T_spiral=2*pi*r_end/(P*w);	%CAV eq.(6)
t=0:T:T_spiral;
xs=P*w/(2*pi)*t.*cos(w*t)+64;	%eq.(17)
ys=P*w/(2*pi)*t.*sin(w*t)+64;	%eq.(18)
plot(xs,ys,'+')
axis equal
% axis ([-64 64 -64 64])
grid on
hold on
%plot the outside dashed circle
xita=0:(pi/180):(2*pi);
x=r_end*cos(xita)+64;
y=r_end*sin(xita)+64;
plot(x,y,'--')