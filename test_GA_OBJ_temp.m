%% GA optimize hyp
% ����GPģ����Ҫ���� N��FieldD��
% zhouwieyan 20180802
%% �����Ŵ��㷨����
N=3;            %���Ż��ı����ĸ���,�ֱ�Ϊell,sf,sn
NIND=40;        %������Ŀ
MAXGEN=4;       %����Ŵ�����
PRECI=10;       %�����Ķ�����λ��
GGAP=0.95;      %����
px=0.7;         %�������
pm=0.01;        %�������
trace=zeros(N+1,MAXGEN);                        %Ѱ�Ž���ĳ�ʼֵ
FieldD=[repmat(PRECI,1,N);[[1e-3;2] [1e-3;50] [1e-3;2]];repmat([1;0;1;1],1,N)]; %����������
Chrom=crtbp(NIND,PRECI*N);                      %��ʼ��Ⱥ
%% �Ż�
gen=0;                                 %��������
hyp_popu=bs2rv(Chrom,FieldD);          %�����ʼ��Ⱥpopulation��ʮ����ת��
ObjV=objfunc(meanfunc,covfunc,likfunc,hyp_popu,NIND,X_train,y_train,X_test);           %����Ŀ�꺯��ֵ
min0=min(ObjV);
while gen<MAXGEN
   fprintf('%d\n',gen)
   FitnV=ranking(ObjV);                         %������Ӧ��ֵ
   SelCh=select('sus',Chrom,FitnV,GGAP);        %ѡ��
   SelCh=recombin('xovsp',SelCh,px);            %����
   SelCh=mut(SelCh,pm);                         %����
   hyp_popu=bs2rv(SelCh,FieldD);                %�Ӵ������ʮ����ת��
   ObjVSel=objfunc(meanfunc,covfunc,likfunc,hyp_popu,size(hyp_popu,1),X_train,y_train,X_test);             %�����Ӵ���Ŀ�꺯��ֵ
   [Chrom,ObjV]=reins(Chrom,SelCh,1,[1 1],ObjV,ObjVSel);    %�ز����Ӵ����������õ�����Ⱥ
   hyp_popu=bs2rv(Chrom,FieldD);
   gen=gen+1;                                   %������������
   %ObjV=objfunc(hyp_popu,size(hyp_popu,1),X_train,y_train,X_test,y_test,ns);
   %��ȡÿ�������Ž⼰����ţ�YΪ���Ž�,IΪ��������
   [Y,I]=min(ObjV);
   trace(1:N,gen)=hyp_popu(I,:);                %����ÿ��������ֵλ��
   trace(end,gen)=Y;                            %����ÿ��������ֵ
end
%% ������ͼ
figure;
plot(1:MAXGEN,trace(end,:));
grid on
xlabel('�Ŵ�����')
ylabel('���ı仯')
title('��������')
bestX=trace(1:end-1,end);
bestObjV=trace(end,end);
fprintf(['���ų�ʼֵ\nX=',num2str(bestX'),'\n��С���err=',num2str(bestObjV),'\n'])

%% define functions
function ObjV=objfunc(meanfunc,covfunc,likfunc,hyp_popu,NIND,X_train,y_train,X_test)
% �����ֱ������Ⱥ�и��������Ŀ��ֵ
% ����
% hyp_popu�����и���ĳ�ʼֵ��NIND�鳬�����ĳ�ֵell,sf,sn

% gaussian regression
ObjV=zeros(NIND,1); 
ell=hyp_popu(:,1); sf=hyp_popu(:,2); hyp0.cov=log([ell sf]);
sn=hyp_popu(:,3); hyp0.lik=log(sn);
for i=1:NIND
    i
    hyp=[];         % ��ʼ��
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
