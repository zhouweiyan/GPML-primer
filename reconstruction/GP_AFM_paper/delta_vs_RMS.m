% Plot the tendency of RMS vs. delta
% array2, array4; GP, DT
% RMS is a 20x1 vector. 

clear; close all; %clc
load('E:\OneDrive - hnu.edu.cn\tools\matlabcourse\GPML_matlab\GPML-primer\reconstruction\GP_AFM_paper\delta_PSNR\delta_RMS.mat')

ax1 = subplot(2, 1, 1); hold on
ln2_spiral_GP = plot(delta_spiral, RMS_array2_spiral(1:10), '-o',...
     'LineWidth',1.5,'Color',[0.8500 0.3250 0.0980],'MarkerSize',5);
ln2_spiral_DT = plot(delta_spiral, RMS_array2_spiral(11:20),  '-.o',...
    'LineWidth',1.5,'Color',[0 0.4470 0.7410],'MarkerSize',5);
grid on; xlim([0.05,0.501]); %ylim([0 50])
title('groove','FontWeight','bold');ylabel('RMS(nm)');
legend({'GP','DT'},'Location','northeast'); legend('boxoff')

ax2 = subplot(2, 1, 2); hold on
ln4_spiral_GP = plot(delta_spiral, RMS_array4_spiral(1:10), '-o',...
     'LineWidth',1.5,'Color',[0.8500 0.3250 0.0980],'MarkerSize',5);
ln4_spiral_DT = plot(delta_spiral, RMS_array4_spiral(11:20),  '-.o',...
    'LineWidth',1.5,'Color',[0 0.4470 0.7410],'MarkerSize',5);
grid on; xlim([0.05,0.501]);
title('array','FontWeight','bold');xlabel('Subsampling ratio \delta');ylabel('RMS(nm)')
legend({'GP','DT'},'Location','northeast'); legend('boxoff')

set(ax1,'FontSize',14,'FontName', 'Times New Roman');
set(ax2,'FontSize',14,'FontName', 'Times New Roman');
set_fig_units_cm(12,12)




