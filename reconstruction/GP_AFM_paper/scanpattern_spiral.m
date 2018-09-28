% simulate spiral scanning pattern
% Reference: I. A. Mahmood, S. O. R. Moheimani, and B. Bhikkaji, ¡°A New Scanning Method for Fast Atomic Force Microscopy,¡± 
% IEEE Transactions on Nanotechnology, vol. 10, no. 2, pp. 203¨C216, Mar. 2011.
% zhouweiyan 20180925

clear
close all
clc
opt=1;
r_end=64*1.6;
P=r_end*2/(60-1);%P=spiral radius¡Á2/(number of curves-1) eq.(4)
%% Constant Angular Velocity(CAV)
if opt==1
    w=120;
    f=2000;
    T=1/f;
    T_spiral=2*pi*r_end/(P*w);	%CAV eq.(6)
    t=0:T:T_spiral;
    xs=P*w/(2*pi)*t.*cos(w*t)+64;	%eq.(17)
    ys=P*w/(2*pi)*t.*sin(w*t)+64;	%eq.(18)
    
    %plot the outside dashed circle
    % xita=0:(pi/180):(2*pi);
    % x=r_end*cos(xita);
    % y=r_end*sin(xita);
    % plot(x,y,'--')
    
end
%% Constant Linear Velocity(CLV)
if opt==2
    v=10;
    f=2;
    T=1/f;
    T_spiral= pi*r_end*r_end/(P*v); %eq.(15)
    t=0:T:T_spiral;
    %eq.(17)(18)
    xs_0=sqrt(P*v*t/pi).*cos(sqrt(4*pi*v*t/P))+64;
    ys_0=sqrt(P*v*t/pi).*sin(sqrt(4*pi*v*t/P))+64;
    xs=fliplr(xs_0);
    ys=fliplr(ys_0);
end
plot(xs,ys,'*')
hold on
plot([1,1,128,128,1],[1,128,128,1,1],'--');
axis equal
axis ([0 128 0 128])
axis off
grid off
% hold on;plot(xs,ys)
set(0,'DefaultFigureWindowStyle','docked')
tightfig
set_fig_units_cm(10,10)