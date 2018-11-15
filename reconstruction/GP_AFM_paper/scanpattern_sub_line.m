% zhouweiyan 20180928

% A=zeros(128,128);
% A(2:4:128,1:128)=1;
% imshow(A);
clear
close all
clc
figure
for i=2:1:128
% for i=1:128
    plot([1,128],[i,i],'Color',[0 0.4470 0.7410]);
    hold on
end
plot([1 1 128 128 1],[1 128 128 1 1],'k--','LineWidth',2)
axis equal
axis([1 128 1 128])
axis off
% set(0,'DefaultFigureWindowStyle','docked') 
% tightfig
% set_fig_units_cm(6,6)