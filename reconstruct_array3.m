load array3.txt
% error: array3 11556¡Á1 double
array3=reshape(array3,256,256)';
surf(array3);colormap(jet)
view([0 -90])
axis equal
xlim([1,256])
ylim([1,256])
