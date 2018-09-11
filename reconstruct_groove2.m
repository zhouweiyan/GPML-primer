load groove2.txt
groove2=reshape(groove2,256,256)';
surf(groove2);colormap(jet)
view([0 -90])
axis equal
xlim([1,256])
ylim([1,256])
