% large scale regression exploiting grid structure
% http://www.gaussianprocess.org/gpml/code/matlab/doc/
% zhouweiyan 20180918

clear
clc
x1=linspace(-2,0,5);
x2=linspace(0,6,6);
x=apxGrid('expand',{x1',x2'});
y=sin(x(:,2))+x(:,1)+0.1*gpml_randn(1,[size(x,1),1]);
xs1=linspace(-4,4,100);xs2=linspace(-7,7,110);  % test points
xs=apxGrid('expand',{xs1',xs2'});
xg=apxGrid('create',xs,true,[50,70]);% xg: 2¡Á1 cell,50¡Á1,70¡Á1; xg is relative to xs
cov={'covSEiso','covSEiso'};
covg={'apxGrid',cov,xg};
mean={'meanZero'};lik={'likGauss'};
hyp.cov=zeros(4,1);hyp.mean=[];hyp.lik=log(0.1);
opt.cg_maxit=200;opt.cg_tol=5e-3;
infg=@(varargin)infGaussLik(varargin{:},opt);
hyp=minimize(hyp,@gp,-50,infg,[],covg,[],x,y);
[post,nlZ,dnlZ]=infGrid(hyp,{'meanZero'},covg,{'likGauss'},x,y,opt);
[fm,fs2,ym,ys2]=post.predict(xs);
[X1,X2]=meshgrid(xs1,xs2);%surf(X1',X2',reshape(ym,length(xs1),length(xs2)))
ym=reshape(ym,length(xs1),length(xs2))';
surf(X1,X2,ym);ym=ym(:);
xlabel('x1');ylabel('x2')

