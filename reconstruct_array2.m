load array2.txt
array2=reshape(array2,256,256)';
surf(array2);colormap(jet)
view([0 -90])
axis equal
xlim([1,256])
ylim([1,256])
