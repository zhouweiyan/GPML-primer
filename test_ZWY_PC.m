figure
surf(groove,'EdgeColor','none');colormap(jet)
axis equal
xlim([1,128]);ylim([1,128]);zlim([-90,40]);
xlabel('x1');ylabel('x2');view(3)
grid on
axis equal
xlim([min(x1_train),max(x1_train)]);ylim([min(x2_train),max(x2_train)]);zlim([-90,40]);
xlabel('x_1/\mum');ylabel('x_2/\mum');zlabel('y/nm');
% view(3)
view(-45,50)
x_label=cell(1,length(0:32:128));
for i=1:length(0:32:128)
    x_label{i}=32*(i-1)/128*2;
end
set(gca,'xtick',0:32:128);
set(gca,'xticklabel',x_label);
set(gca,'ytick',0:32:128);
set(gca,'yticklabel',x_label);
tightfig
set_fig_units_cm(10,10)