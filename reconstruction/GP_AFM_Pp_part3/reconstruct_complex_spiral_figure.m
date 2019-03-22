close all
%%
figure, surf(y_test_ideal, 'EdgeColor', 'None');        
colormap(jet); axis equal; axis off;
xlabel('x1'); ylabel('x2'); view([0 -90])
tightfig
set_fig_units_cm(10,10)
set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1.5; 0 0 1 0; 0 0 0 1]);
%%
figure
surf(X1_test,X2_test,reshape(m,length(x2_test),length(x1_test)),'EdgeColor', 'None');
colormap(jet); axis equal; axis off;
xlabel('x1'); ylabel('x2'); view([0 -90])
tightfig
set_fig_units_cm(10,10)
set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1.5; 0 0 1 0; 0 0 0 1]);
%%
figure
surf(y_test_ideal_a-m_a,'EdgeColor','none','LineStyle','none','FaceLighting','phong');   % surf(y_test_ideal_a-m_a,'FaceAlpha',0.1); 
zlim([-90,40]); axis equal;axis off;
colormap(jet)
view([-90 90])
caxis([-10,10])
colorbar('eastoutside', 'Ticks', [-10, -5, 0, 5, 10], ...
    'TickLabels',[-10, -5, 0, 5, 10], 'FontSize', 24,...
    'FontName', 'Times New Roman')%, 'AxisLocation','in')
% set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 4 0; 0 -1 0 1.7; 0 0 1 0; 0 0 0 1]);
set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 2.7 0; 0 -1 0 1.5; 0 0 1 0; 0 0 0 1]);
