% 清空运行环境
close all
clear
clc
% 速度更新参数
c1=1.49445;
c2=1.49445;
maxgen=300; % 迭代次数
sizepop=20; % 种群规模
% 个体和速度最大最小值
popmax=2;popmin=-2;
Vmax=0.5;Vmin=-0.5;

%% 种群初始化
for i=1:sizepop
    % 随机产生一个种群
    pop(i,:)=2*rand