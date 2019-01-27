% I. A. Mahmood and S. O. Reza Moheimani, ¡°Fast spiral-scan atomic force microscopy,¡± 
% Nanotechnology, vol. 20, no. 36, p. 365503, Sep. 2009.
% zhouweiyan 20190117

clear
close all
clc
r_end=10;   % r_end is the radius, i.e., half of the range
pixels=128; 
% eq.(3)
P=r_end*2/pixels;
% f is the scanning frequency.
f=20;
w=2*pi*f;
% eq.(5) calculate t_total
t_total=2*pi*r_end/(P*w)
% f_c is the sampling frequency
f_c=20000;
% eq.(2),(6),(7)
t=0:(1/f_c):t_total;
r=P/(2*pi)*w*t;
x_s=r.*cos(w*t);
y_s=r.*sin(w*t);
plot(x_s,y_s,'.')
axis equal
grid on
%% period 2: go back through a straight line
x_end=x_s(end)
y_end=y_s(end)
hold on
plot(x_end,y_end,'*')
% x_ends are near 0 and y_ends flatuate very subtly around 20um. 
% Thus, 1 second is allocated to go back to the origin.
t_b=1;
t=0:1/f_c:t_b;
% x_b, y_b represent the position sequences when the probe goes back.
y_b=y_end*(t_b-t)/t_b;
x_b=x_end*(t_b-t)/t_b;
hold on
plot(x_b,y_b,'g+')
%% change the scale -10~10 to -0.15~0.15 or 0~0.1 and save the trajectory of two periods
% ratio=0.1/20; % 0~0.1
ratio=0.15/10;
x_o=zeros(1,50000);
y_o=zeros(1,50000);
x=[x_o,[x_s,x_b]*ratio];
y=[y_o,[y_s,y_b]*ratio];
% x = x + 0.05; % 0~0.1
% y = y + 0.05; % 0~0.1
figure
plot(x,y,'.')



