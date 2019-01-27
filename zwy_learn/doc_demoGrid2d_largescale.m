% large scale regression exploiting grid structure
% http://www.gaussianprocess.org/gpml/code/matlab/doc/
% zhouweiyan 20180918

clear
clc
x1=linspace(-2,2,120);
x2=linspace(-9,9,130);
x=apxGrid('expand',{x1',x2'});
% x is different from the result of x=meshgrid(x1,x2);
% solution
% method 1:
% x=apxGrid('expand',{x2',x1'});x_t=x(:,2);x(:,2)=x(:,1);x(:,1)=x_t;
y=sin(x(:,2))+x(:,1)+0.1*gpml_randn(1,[size(x,1),1]);
xs1=linspace(-4,4,100);xs2=linspace(-7,7,110);  % test points
xs=apxGrid('expand',{xs1',xs2'});
xg=apxGrid('create',[x;xs],true,[50,70]);
% xg: 2¡Á1 cell,50¡Á1,70¡Á1; xg is relative to the bigger of x and xs
cov={'covSEiso','covSEiso'};
covg={'apxGrid',cov,xg};
mean={'meanZero'};lik={'likGauss'};
hyp.cov=zeros(4,1);hyp.mean=[];hyp.lik=log(0.1);
opt.cg_maxit=200;opt.cg_tol=5e-3;
infg=@(varargin)infGaussLik(varargin{:},opt);
hyp=minimize(hyp,@gp,-50,infg,[],covg,[],x,y);
[post,nlZ,dnlZ]=infGrid(hyp,{'meanZero'},covg,lik,x,y,opt);
[fm,fs2,ym,ys2]=post.predict(xs);
[X1,X2]=meshgrid(xs1,xs2);%surf(X1',X2',reshape(ym,length(xs1),length(xs2)))
% method two
ym=reshape(ym,length(xs1),length(xs2))';
surf(X1,X2,ym);ym=ym(:);
xlabel('x1');ylabel('x2')
colormap(jet)
