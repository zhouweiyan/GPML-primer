% raster scanning pattern
% x-axis: triangle wave; y-axis: step by step
% zhouweiyan 20190117

clear 
close all
clc
% f is the scanning frequency. t_T is the time of finishing a trace and retrace line.
f=10;
t_T=1/f;
pixels=128; 
% Total time of generating an image of pixels x pixels. 
t_total=pixels*t_T
% range(um) x range(um) is the size of the scanning area.
range=20;
% f_c is the sampling frequency
f_c=20000;

t=0:1/f_c:t_total;
x_t=@(t)(mod(floor(2*t/t_T),2)==0).*(2*range/t_T*(t-floor(t/t_T)*t_T))+...
    (mod(floor(2*t/t_T),2)~=0).*(-2*range/t_T*(t-floor(t/t_T)*t_T-t_T));
x_r=x_t(t);
figure(1)
plot(t,x_r,'.') 
title('x(t)')
% % verify the plot of x
% hold on
% plot(t_T/2,x_t(t_T/2),'*')

figure(2)
y_t=@(t)range/(pixels-1)*(floor(t/t_T));
y_r=y_t(t);
plot(t,y_r)
title('y(t)')

figure(3)
plot(x_r,y_r,'.')
title('x,y')
%% period 2: go back through a straight line
x_end=x_r(end);
y_end=y_r(end);
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
% ratio = 0.1/20;   % 0~0.1
ratio = 0.15/10;    % -0.15~0.15
x_o=zeros(1,50000);
y_o=zeros(1,50000);
x=[x_o,[x_r,x_b]*ratio];
y=[y_o,[y_r,y_b]*ratio];
figure
plot(x,y,'.')










