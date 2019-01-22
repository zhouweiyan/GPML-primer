%% PSO optimize hyp
% zhouweiyan 20180806   将GA和PSO用于GP超参数初值的确定，PSO效果不行
%% PSO参数初始化
%粒子群算法中的两个参数
% c1 = 1.49445;
% c2 = 1.49445;
c1=0.8;c2=0.8;

N=3;            % 待优化的变量的个数,分别为ell,sf,sn
maxgen=3;       % 进化次数  
sizepop=40;    % 种群规模 n

popmax=[2 50 2];popmin=[1e-3 1e-3 1e-3];
Vmax=0.2*(popmax-popmin)./2;Vmin=-Vmax;
% Vmax=0.5;Vmin=-0.5;
ws=0.9;
we=0.4;
%% 产生初始粒子和速度
%预分配    modified by zwy
pop=zeros(sizepop,N);    % X: n×D，n个粒子，每个粒子∈R^D，D为待优化参数的维数
V=zeros(sizepop,N);      % V: n×D
fitness=1e5*ones(sizepop,1);% n×1

for i=1:sizepop
    i
    %随机产生一个种群
    pop(i,:)=popmin+(popmax-popmin).*rand(1,N);    % 初始种群 X: n×D
    V(i,:)=Vmin+(Vmax-Vmin).*rand(1,N);    % 初始化速度 V: n×D
%     V(i,:)=0.5*rands(1,N);
    %计算适应度
    fitness(i)=fun(pop(i,:),1,X_train,y_train,X_test,f_test); % 染色体的适应度 n×1
end

%% 寻找初始极值，个体极值和群体极值
[bestfitness,bestindex]=min(fitness);
zbest=pop(bestindex,:);     %全局最佳 1×D
gbest=pop;                  %个体最佳 n×D
fitnessgbest=fitness;       %个体最佳适应度值 n×1
fitnesszbest=bestfitness;   %全局最佳适应度值 1×1
yy=zeros(1,maxgen);
%% 迭代寻优
yy=zeros(1,maxgen);
for i=1:maxgen
    i
    w=ws-(ws-we)*( 2*i/maxgen-(i/maxgen)^2 );
    % 粒子位置和速度更新
    for j=1:sizepop
        j
        %速度更新
        V(j,:) = w*V(j,:) + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:));
        V(j,find(V(j,:)>Vmax))=Vmax(find(V(j,:)>Vmax));
        V(j,find(V(j,:)<Vmin))=Vmin(find(V(j,:)<Vmin));
        
        %种群更新
        pop(j,:)=pop(j,:)+V(j,:);
        pop(j,find(pop(j,:)>popmax))=popmax(find(pop(j,:)>popmax));
        pop(j,find(pop(j,:)<popmin))=popmin(find(pop(j,:)<popmin));
        
        %适应度值
        fitness(j)=fun(pop(i,:),1,X_train,y_train,X_test,f_test); 
    end
    
    for j=1:sizepop
        
        %个体最优更新
        if fitness(j) < fitnessgbest(j)
            gbest(j,:) = pop(j,:);
            fitnessgbest(j) = fitness(j);
        end
        
        %群体最优更新
        if fitness(j) < fitnesszbest
            zbest = pop(j,:);
            fitnesszbest = fitness(j);
        end
    end 
    yy(i)=fitnesszbest;      
end
% %% 结果分析
plot(yy)
title('最优个体适应度','fontsize',12);
% xlabel('进化代数','fontsize',12);ylabel('适应度','fontsize',12);
fitnesszbest

%% define functions
function ObjV=fun(hyp_popu,sizepop,X_train,y_train,X_test,f_test)
% 用来分别求解种群中各个个体的目标值
% 输入
% hyp_popu：所有个体的初始值，NIND组超参数的初值ell,sf,sn
% X1_test,X2_test:测试样本输入
% y_test:测试样本期望输出

% gaussian regression
ObjV=zeros(sizepop,1); nlml=zeros(sizepop,1);
meanfunc=[]; covfunc={'covSEiso'}; likfunc='likGauss';
ell=hyp_popu(:,1); sf=hyp_popu(:,2); hyp0.cov=log([ell sf]);
sn=hyp_popu(:,3); hyp0.lik=log(sn);
for i=1:sizepop
    
    hyp=[];         % 初始化
    hyp.mean=[];
    hyp.cov=[hyp0.cov(i,1);hyp0.cov(i,2)];
    hyp.lik=hyp0.lik(i);
    hyp=minimize(hyp,@gp,-100,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train);
    disp(['exp(hyp1.lik)=',num2str(exp(hyp.lik))])
    nlml(i)=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train);
    try
        m=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_test);
    catch
        ObjV(i)=1e5;
        continue
    end
    ObjV(i)=sum(abs(f_test(:)-m(:)));  % error
    disp(['sum(abs(error))=',num2str(ObjV(i))])
end
end

