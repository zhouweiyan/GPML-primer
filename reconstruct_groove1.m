load groove1.txt
groove1=reshape(groove1,256,256)';
surf(groove1);colormap(jet)
view([0 -90])
axis equal
xlim([1,256])
ylim([1,256])
