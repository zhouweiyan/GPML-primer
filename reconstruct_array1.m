load array1.txt
array1=reshape(array1,256,256)';
surf(array1);colormap(jet)
view([0 -90])
axis equal
xlim([1,256])
ylim([1,256])
