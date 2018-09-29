% simulate spiral scanning pattern
% Reference: I. A. Mahmood, S. O. R. Moheimani, and B. Bhikkaji, ¡°A New Scanning Method for Fast Atomic Force Microscopy,¡± 
% IEEE Transactions on Nanotechnology, vol. 10, no. 2, pp. 203¨C216, Mar. 2011.
% zhouweiyan 20180928

clear
close all
clc
opt=1;
r_end=64*1.6;
P=r_end*2/(60-1);%P=spiral radius¡Á2/(number of curves-1) eq.(4)

w=120;
f=2000;
T=1/f;
T_spiral=2*pi*r_end/(P*w);	%CAV eq.(6)
t=0:T:T_spiral;
xs=P*w/(2*pi)*t.*cos(w*t);	%eq.(17)
ys=P*w/(2*pi)*t.*sin(w*t);	%eq.(18)
path_len=sum(sqrt((xs(2:end)-xs(1:(end-1))).^2+(ys(2:end)-ys(1:(end-1))).^2))

plot(xs,ys);
hold on
plot([-64,-64,64,64,-64],[-64,64,64,-64,-64],'k--','LineWidth',1);
axis equal
% axis ([0 128 0 128])
axis([-r_end r_end -r_end r_end])
axis off
grid off
% hold on;plot(xs,ys)
set(0,'DefaultFigureWindowStyle','docked')
tightfig
set_fig_units_cm(10,10)

