%% GA optimize hyp
% 更改GP模型需要更改 N，FieldD，
% zhouwieyan 20180802
%% 定义遗传算法参数
N=3;            %待优化的变量的个数,分别为ell,sf,sn
NIND=40;        %个体数目
MAXGEN=4;       %最大遗传代数
PRECI=10;       %变量的二进制位数
GGAP=0.95;      %代沟
px=0.7;         %交叉概率
pm=0.01;        %变异概率
trace=zeros(N+1,MAXGEN);                        %寻优结果的初始值
FieldD=[repmat(PRECI,1,N);[[1e-3;2] [1e-3;50] [1e-3;2]];repmat([1;0;1;1],1,N)]; %区域描述器
Chrom=crtbp(NIND,PRECI*N);                      %初始种群
%% 优化
gen=0;                                 %代计数器
hyp_popu=bs2rv(Chrom,FieldD);          %计算初始种群population的十进制转换
ObjV=objfunc(meanfunc,covfunc,likfunc,hyp_popu,NIND,X_train,y_train,X_test);           %计算目标函数值
min0=min(ObjV);
while gen<MAXGEN
   fprintf('%d\n',gen)
   FitnV=ranking(ObjV);                         %分配适应度值
   SelCh=select('sus',Chrom,FitnV,GGAP);        %选择
   SelCh=recombin('xovsp',SelCh,px);            %重组
   SelCh=mut(SelCh,pm);                         %变异
   hyp_popu=bs2rv(SelCh,FieldD);                %子代个体的十进制转换
   ObjVSel=objfunc(meanfunc,covfunc,likfunc,hyp_popu,size(hyp_popu,1),X_train,y_train,X_test);             %计算子代的目标函数值
   [Chrom,ObjV]=reins(Chrom,SelCh,1,[1 1],ObjV,ObjVSel);    %重插入子代到父代，得到新种群
   hyp_popu=bs2rv(Chrom,FieldD);
   gen=gen+1;                                   %代计数器增加
   %ObjV=objfunc(hyp_popu,size(hyp_popu,1),X_train,y_train,X_test,y_test,ns);
   %获取每代的最优解及其序号，Y为最优解,I为个体的序号
   [Y,I]=min(ObjV);
   trace(1:N,gen)=hyp_popu(I,:);                %记下每代的最优值位置
   trace(end,gen)=Y;                            %记下每代的最优值
end
%% 画进化图
figure;
plot(1:MAXGEN,trace(end,:));
grid on
xlabel('遗传代数')
ylabel('误差的变化')
title('进化过程')
bestX=trace(1:end-1,end);
bestObjV=trace(end,end);
fprintf(['最优初始值\nX=',num2str(bestX'),'\n最小误差err=',num2str(bestObjV),'\n'])

%% define functions
function ObjV=objfunc(meanfunc,covfunc,likfunc,hyp_popu,NIND,X_train,y_train,X_test)
% 用来分别求解种群中各个个体的目标值
% 输入
% hyp_popu：所有个体的初始值，NIND组超参数的初值ell,sf,sn

% gaussian regression
ObjV=zeros(NIND,1); 
ell=hyp_popu(:,1); sf=hyp_popu(:,2); hyp0.cov=log([ell sf]);
sn=hyp_popu(:,3); hyp0.lik=log(sn);
for i=1:NIND
    i
    hyp=[];         % 初始化
    hyp.mean=[];
    hyp.cov=[hyp0.cov(i,1);hyp0.cov(i,2)];
    hyp.lik=hyp0.lik(i);
    hyp=minimize(hyp,@gp,-100,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train);
    disp(['exp(hyp1.lik)=',num2str(exp(hyp.lik))])
    try
        [~,s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_test);
    catch
        ObjV(i)=1e15;
        continue
    end
    ObjV(i)=mean(sqrt(s2(:)))+max(sqrt(s2(:))); 
    disp(['sum(abs(s2(:)))=',num2str(sum(abs(s2(:))))])
end
end
