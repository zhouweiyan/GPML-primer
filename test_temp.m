% ������л���
close all
clear
clc
% �ٶȸ��²���
c1=1.49445;
c2=1.49445;
maxgen=300; % ��������
sizepop=20; % ��Ⱥ��ģ
% ������ٶ������Сֵ
popmax=2;popmin=-2;
Vmax=0.5;Vmin=-0.5;

%% ��Ⱥ��ʼ��
for i=1:sizepop
    % �������һ����Ⱥ
    pop(i,:)=2*rand