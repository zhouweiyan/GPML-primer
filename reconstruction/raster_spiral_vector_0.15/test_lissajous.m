% T. Tuma, J. Lygeros, V. Kartik, A. Sebastian, and A. Pantazi
%¡°High-speed multiresolution scanning probe microscopy based on Lissajous scan trajectories,¡± 
% Nanotechnology, vol. 23, no. 18, p. 185501, May 2012.
% zhouweiyan 20190117

% the constants in eq.(1),(2)
% (X_0,Y_0) is the initial point. A,B are ranges in x-axis and y-axis respectively.
X_0=0; Y_0=0; A=20; B=20;   
% f_x,f_y are scanning frequencies in x-axis and y-axis respectively.
f_x=60;
f_y=61;
% eq.(3) T is the duration of the scan trajectory.
T=lcm(f_x,f_y)/(f_x*f_y)
% f_c is the sampling rate. It should be smaller when f_x,f_y is large.
f_c=20000; 
t=0:(1/f_c):T;
% eq.(1),(2)
x=X_0+A/2*sin(2*pi*f_x*t);
y=Y_0+B/2*sin(2*pi*f_y*t);
figure
plot(x,y)

%% sampling process
figure
for t=0:1/f_c:T
    x=X_0+A/2*sin(2*pi*f_x*t);
    y=Y_0+B/2*sin(2*pi*f_y*t);
    plot(x,y,'.')
    pause(0.0001)
    hold on
end



