f = imread('Corrosion.bmp');
f = double(f)/255;
f = rgb2gray(f);
surf(f);
surf(f, 'EdgeColor', 'None');
colormap(copper);
colormap(gray)