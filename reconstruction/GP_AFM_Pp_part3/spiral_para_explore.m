% [1]I. A. Mahmood and S. O. Reza Moheimani, ¡°Fast spiral-scan atomic force microscopy,¡± Nanotechnology, vol. 20, no. 36, p. 365503, Sep. 2009.
% explore the details of the parameters in [1]
% zwy 20190131

clear
close all
clc
r_end = 64*1.6;     % cover the square of 128*128 pixels
P = r_end*2/(103-1); % P = spiral radius¡Á2/(number of curves-1) eq.(4)
w = 120;
f = 3000;           % inversely proportional to the number of sampling points
T = 1/f;
T_spiral = 2*pi*r_end/(P*w);	%CAV eq.(6)
t = 0:T:T_spiral;
xs = P*w/(2*pi)*t.*cos(w*t)+64;	%eq.(17)
ys = P*w/(2*pi)*t.*sin(w*t)+64;	%eq.(18)
path_len=sum(sqrt((xs(2:end)-xs(1:(end-1))).^2+(ys(2:end)-ys(1:(end-1))).^2));
delta = path_len/(2*128*128);

figure
plot(xs,ys,'.')
axis equal
axis ([0 128 0 128])
grid on
hold on; plot(xs,ys)