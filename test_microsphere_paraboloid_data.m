% zhouweiyan 20180811
% ideal data of microsphere_paraboloid
clear
close all
clc
x1=-10:0.1:10;
x2=-10:0.1:10;
len=length(x1);
[X1,X2]=meshgrid(x1,x2);
R=0.8;
y_base=-0.1*X1(:).^2-0.1*X2(:).^2;
% figure
% surf(X1,X2,reshape(y_base,length(x2),length(x1)),'EdgeColor','none',...
%     'LineStyle','none','FaceLighting','phong');
% colormap(jet)

% microshpere
y_add=zeros(len,len);
xi=-9:2:9;
yi=-9:2:9;
for i=1:len
    for j=1:len
        for p=1:length(xi)
            for q=1:length(yi)
                if (X1(i,j)-xi(p))^2+(X2(i,j)-yi(q))^2<=R^2
                    y_add(i,j)=sqrt(R^2-((X1(i,j)-xi(p))^2+(X2(i,j)-yi(q))^2));
                end
            end
        end
    end
end
% figure
% surf(X1,X2,y_add,'EdgeColor','none',...
%     'LineStyle','none','FaceLighting','phong');
% zlim([0 20])
% colormap(jet)

y_ideal=y_base+y_add(:);
seed=1;
rng(seed,'v4');
s_noise=0.25;
y_noise=y_ideal+s_noise*randn(size(y_ideal));
rng('shuffle')
set(0,'DefaultFigureWindowStyle','docked'); % keep figures tidy
figure
subplot(121);surf(X1,X2,reshape(y_ideal,length(x1),length(x2)),'EdgeColor','none',... 
    'LineStyle','none','FaceLighting','phong');
title('target surface without noise')
subplot(122);surf(X1,X2,reshape(y_noise,length(x1),length(x2)),'EdgeColor','none',...
    'LineStyle','none','FaceLighting','phong');
title('target surface with noise');
colormap(jet)