function ObjV=test_objfunc(meanfunc,covfunc,likfunc,hyp_popu,NIND,X_train,y_train,X_test,ns)
% test_objfunc(hyp_popu',1,X_train,y_train,X_test,f_test,[length(x1_test) length(x2_test)])
% �����ֱ������Ⱥ�и��������Ŀ��ֵ
% ����
% hyp_popu:һ�д���һ�鳬���������и���ĳ�ʼֵ��NIND�鳬�����ĳ�ֵell,sf,sn
% X1_test,X2_test:������������
% y_test:���������������

% gaussian regression
ObjV=zeros(NIND,1); nlml=zeros(NIND,1);
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
    nlml(i)=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train);
    try
        [m,s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_test);
        surf(reshape(X_test(:,1),ns(2),ns(1)),reshape(X_test(:,2),ns(2),ns(1)),reshape(m,ns(2),ns(1)));
        surf(reshape(X_test(:,1),ns(2),ns(1)),reshape(X_test(:,2),ns(2),ns(1)),reshape(m+2*sqrt(s2),ns(2),ns(1)),'FaceAlpha',0.1);
        surf(reshape(X_test(:,1),ns(2),ns(1)),reshape(X_test(:,2),ns(2),ns(1)),reshape(m-2*sqrt(s2),ns(2),ns(1)),'FaceAlpha',0.1);
    catch
        ObjV(i)=1e15;
        continue
    end
    ObjV(i)=sum(abs(s2(:)));   % error
    disp(['sum(abs(error))=',num2str(ObjV(i))])
end
end
