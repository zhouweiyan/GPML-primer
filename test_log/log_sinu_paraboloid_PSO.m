%% PSO optimize hyp
% zhouweiyan 20180806   ��GA��PSO����GP��������ֵ��ȷ����PSOЧ������
%% PSO������ʼ��
%����Ⱥ�㷨�е���������
% c1 = 1.49445;
% c2 = 1.49445;
c1=0.8;c2=0.8;

N=3;            % ���Ż��ı����ĸ���,�ֱ�Ϊell,sf,sn
maxgen=3;       % ��������  
sizepop=40;    % ��Ⱥ��ģ n

popmax=[2 50 2];popmin=[1e-3 1e-3 1e-3];
Vmax=0.2*(popmax-popmin)./2;Vmin=-Vmax;
% Vmax=0.5;Vmin=-0.5;
ws=0.9;
we=0.4;
%% ������ʼ���Ӻ��ٶ�
%Ԥ����    modified by zwy
pop=zeros(sizepop,N);    % X: n��D��n�����ӣ�ÿ�����ӡ�R^D��DΪ���Ż�������ά��
V=zeros(sizepop,N);      % V: n��D
fitness=1e5*ones(sizepop,1);% n��1

for i=1:sizepop
    i
    %�������һ����Ⱥ
    pop(i,:)=popmin+(popmax-popmin).*rand(1,N);    % ��ʼ��Ⱥ X: n��D
    V(i,:)=Vmin+(Vmax-Vmin).*rand(1,N);    % ��ʼ���ٶ� V: n��D
%     V(i,:)=0.5*rands(1,N);
    %������Ӧ��
    fitness(i)=fun(pop(i,:),1,X_train,y_train,X_test,f_test); % Ⱦɫ�����Ӧ�� n��1
end

%% Ѱ�ҳ�ʼ��ֵ�����弫ֵ��Ⱥ�弫ֵ
[bestfitness,bestindex]=min(fitness);
zbest=pop(bestindex,:);     %ȫ����� 1��D
gbest=pop;                  %������� n��D
fitnessgbest=fitness;       %���������Ӧ��ֵ n��1
fitnesszbest=bestfitness;   %ȫ�������Ӧ��ֵ 1��1
yy=zeros(1,maxgen);
%% ����Ѱ��
yy=zeros(1,maxgen);
for i=1:maxgen
    i
    w=ws-(ws-we)*( 2*i/maxgen-(i/maxgen)^2 );
    % ����λ�ú��ٶȸ���
    for j=1:sizepop
        j
        %�ٶȸ���
        V(j,:) = w*V(j,:) + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:));
        V(j,find(V(j,:)>Vmax))=Vmax(find(V(j,:)>Vmax));
        V(j,find(V(j,:)<Vmin))=Vmin(find(V(j,:)<Vmin));
        
        %��Ⱥ����
        pop(j,:)=pop(j,:)+V(j,:);
        pop(j,find(pop(j,:)>popmax))=popmax(find(pop(j,:)>popmax));
        pop(j,find(pop(j,:)<popmin))=popmin(find(pop(j,:)<popmin));
        
        %��Ӧ��ֵ
        fitness(j)=fun(pop(i,:),1,X_train,y_train,X_test,f_test); 
    end
    
    for j=1:sizepop
        
        %�������Ÿ���
        if fitness(j) < fitnessgbest(j)
            gbest(j,:) = pop(j,:);
            fitnessgbest(j) = fitness(j);
        end
        
        %Ⱥ�����Ÿ���
        if fitness(j) < fitnesszbest
            zbest = pop(j,:);
            fitnesszbest = fitness(j);
        end
    end 
    yy(i)=fitnesszbest;      
end
% %% �������
plot(yy)
title('���Ÿ�����Ӧ��','fontsize',12);
% xlabel('��������','fontsize',12);ylabel('��Ӧ��','fontsize',12);
fitnesszbest

%% define functions
function ObjV=fun(hyp_popu,sizepop,X_train,y_train,X_test,f_test)
% �����ֱ������Ⱥ�и��������Ŀ��ֵ
% ����
% hyp_popu�����и���ĳ�ʼֵ��NIND�鳬�����ĳ�ֵell,sf,sn
% X1_test,X2_test:������������
% y_test:���������������

% gaussian regression
ObjV=zeros(sizepop,1); nlml=zeros(sizepop,1);
meanfunc=[]; covfunc={'covSEiso'}; likfunc='likGauss';
ell=hyp_popu(:,1); sf=hyp_popu(:,2); hyp0.cov=log([ell sf]);
sn=hyp_popu(:,3); hyp0.lik=log(sn);
for i=1:sizepop
    
    hyp=[];         % ��ʼ��
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

