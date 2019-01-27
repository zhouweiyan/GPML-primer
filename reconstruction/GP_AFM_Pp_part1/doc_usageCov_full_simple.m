% demonstrate usage of simple covariance functions
% zhouweiyan 20180919
clear
close all
clc
%% initialization
n=5;D=2;
seed=0;
rand('state',seed);
randn('state',seed);
x=randn(n,D);
xs=randn(3,D);

%% simple covariance library
opt=8;
switch opt
    case 1
%         co={'covOne'};
%         cov=co;hyp=[];
        cc={'covConst'};
        sf=2;hypc=log(sf);cov=cc;hyp=hypc;
    case 2
        cc={'covCos'};  % x belongs to R1
        p=1;sf=1;hypcc=log([p;sf]);   cov=cc;hyp=hypcc;
    case 3
        ce={'covEye'};  % k(x^p,x^q) = \delta(p,q)
        hype=[];cov=ce;hyp=hype;
    case 4
        cgbi={'covGaboriso'};
        ell=5;p=3;hypgbi=log([ell;p]);  cov=cgbi;hyp=hypgbi;
        % verify K(1,2)
        % x1 = x(1,:)'; x2 = x(2,:)'; t = x1 - x2; exp(-t'*t/(2*ell^2))*cos(2*pi*sum(t)/p)

%         cgba={'covGaborard'};
%         ell=[5;10];p=[1;6];hypgba=log([ell;p]);  cov=cgba;hyp=hypgba;
    case 5  % stationary
%         cge={'covGE','eye',[]};
%         gamma=1.5;hypge=log(gamma/(2-gamma)); cov=cge;hyp=hypge;
        % verify K(1,2), Ks(1,2)
        % x1 = x(1,:)'; x2 = x(2,:)'; r2 = (x1-x2)'*(x1-x2); r = sqrt(r2); exp(-r^(gamma))
        % x1 = x(1,:)'; xs2 = xs(2,:)'; r2 = (x1-xs2)'*(x1-xs2); r = sqrt(r2); exp(-r^(gamma))
        
%         cgi={'covGE','iso',[]};
%         gamma=1.5;ell=2;hypgi=log([ell;(gamma/(2-gamma))]); cov=cgi;hyp=hypgi;
        % x1 = x(1,:)'; x2 = x(2,:)'; r2 = (x1-x2)'*ell^(-2)*(x1-x2); r = sqrt(r2); exp(-r^(gamma))

        cga={'covGE','ard',[]};
        gamma=1;L=(1:D)';hypga=log([L;(gamma/(2-gamma))]); 
        cov=cga;hyp=hypga;
        % x1 = x(1,:)'; x2 = x(2,:)'; r2 = (x1-x2)'*inv(diag(L))*(x1-x2); r = sqrt(r2); exp(-r^(gamma))
    case 6  % dot
%         cl={'covLIN'};
%         hypl=[];cov=cl;hyp=hypl;

        cli={'covLINiso'};    % K= x*ell^(-2)*x', Ks = x*ell^(-2)*xs'
        ell=0.5; hypli=log(ell);cov=cli;hyp=hypli;  % slope is inversely propotional to ell 
        % x1 = x(1,:)'; x2 = x(2,:)'; x1'*ell^(-2)*x2

%         cla={'covLINard'};  % K= x*inv(diag(L.^2))*x'
%         L=(D:-1:1)';hypla=log(L);cov=cla;hyp=hypla;
        % x1 = x(1,:)'; x2 = x(2,:)'; x1'*inv(diag(L.^2))*x2
        
%         clo={'covLINone'};    % K=(x*x'+1)/ell^2;
%         ell=0.9;hyplo=log(ell);cov=clo;hyp=hyplo;
    case 7
%         cmi={'covMaterniso',3}; 
%         L=0.5;sf=2;hypmi=log([L;sf]);  cov=cmi;hyp=hypmi;
        % f = @(d, t)(d==1)*t+(d==3)*(1+t)+(d==5)*(1+t+t^2/3); 
        % x1=x(1,:)'; x2=x(2,:)'; r=sqrt((x1-x2)'*L^(-2)*(x1-x2)); d=3; sf^2*f(d,sqrt(d)*r)*exp(-sqrt(d)*r)
        cma={'covMaternard',5};
        L=[4;2];sf=3;hypma=log([L;sf]); cov=cma;hyp=hypma;
        % r=sqrt((x1-x2)'*(diag(L))^(-2)*(x1-x2)); d=5;
    case 8
        cn={'covNoise'};% k(x^p,x^q) = sf^2 * \delta(p,q); [varargout{:}] = covScale({'covEye'},varargin{:});
        sf=0.1;hypn=log(sf); cov=cn;hyp=hypn;
    case 9
        cnn={'covNNone'};   % the outline shape remains unchanged, seems useless for me
        ell=0.1;sf=1;hypnn=log([ell;sf]);   cov=cnn;hyp=hypnn;
        % x1=[1,x(1,:)]'; x2=[1,x(2,:)]'; sf^2 * asin(x1'*ell^(-2)*x2 / sqrt((1+x1'*ell^(-2)*x1)*(1+x2'*ell^(-2)*x2)))
    case 10
        cpe={'covPeriodic'};    % x belongs to R1
        ell=0.5;p=2;sf=1;hyppe=log([ell;p;sf]); cov=cpe;hyp=hyppe;
%         cpn={'covPeriodicNoDC'};
%         ell=0.9;p=2;sf=1;hyppn=log([ell;p;sf]); cov=cpn;hyp=hyppn;
    case 11
%         cp={'covPoly',3};   % K=sf^2*(c+x*x').^3
%         c=2;sf=1;hypp=log([c;sf]);cov=cp;hyp=hypp;

%         cpi={'covPoly','iso',[],2};   % K=sf^2*(c+x*ell^(-2)*x').^2
%         ell=2;c=2;sf=1;cov=cpi;hyp=log([ell;c;sf]);

        cpa={'covPoly','ard',[],3};     % K=sf^2*(c+x*inv(diag(L.^2))*x').^3
        c=1;sf=1;L=[2;4];cov=cpa;hyp=log([L;c;sf]); 
    case 12
%         cpi={'covPPiso',2};   % [varargout{:}] = covScale({'covPP','iso',[],v},varargin{:});
%         ell=2;sf=1;hyppi=log([ell;sf]);   cov=cpi;hyp=hyppi;
        % v=2; D=2; j=floor(D/2)+v+1; x1=x(1,:)'; x2=x(2,:)'; r=sqrt((x1-x2)'*(ell^(-2))*(x1-x2));
        
        cpa={'covPPard',3};
        L=[2;1];sf=1;hyppa=log([L;sf]);   cov=cpa;hyp=hyppa;
        %v=3; D=2; j=floor(D/2)+v+1; x1=x(1,:)'; x2=x(2,:)'; r=sqrt((x1-x2)'*(diag(L.^2))^(-1)*(x1-x2));

%         cpi={'covPP','iso',[],2};
%         ell=2;hyppi=log(ell);   cov=cpi;hyp=hyppi;
%         cpa={'covPP','ard',[],2};
%         L=[2;1];hyppa=log(L);   cov=cpa;hyp=hyppa;
       
        % kpp=@(D,v,r)sf^2*((v==0)*(max(0,(1-r))^j)+(v==1)*((max(0,(1-r))^(j+v)))*((j+1)*r+1)+...
        %     (v==2)*(max(0,(1-r))^(j+v))*(((j^2+4*j+3)*(r^2)+(3*j+6)*r+3))/3+...
        %     (v==3)*(max(0,(1-r))^(j+v))*(((j^3+9*j^2+23*j+15)*(r^3)+(6*j^2+36*j+45)*(r^2)+(15*j+45)*r+15))/15);
        % kpp(D,v,r)
    case 13
%         cri={'covRQiso'};
%         ell=0.9;sf=2;al=2;hypri=log([ell;sf;al]);   cov=cri;hyp=hypri;
        % verify K(1,2)
        % x1=x(1,:)';x2=x(2,:)'; sf^2*(1+1/(2*al*ell^2)*(x1-x2)'*(x1-x2))^(-al)
        
        cra={'covRQard'};   
        L=[0.5;1];sf=2;al=1;hypra=log([L;sf;al]);cov=cra;hyp=hypra;
        % x1=x(1,:)';x2=x(2,:)'; sf^2*(1+1/(2*al)*(x1-x2)'*diag(L)^(-2)*(x1-x2))^(-al)
    case 14
%         cgu={'covSEisoU'};  % k(x,z) = exp(-(x-z)'*inv(P)*(x-z)/2); P=ell^2*I;
%         ell=1/sqrt(2);hypgu=log(ell);cov=cgu;hyp=hypgu;

        cgi={'covSEiso'};   
        ell=0.5;sf=1;hypgi=log([ell;sf]);cov=cgi;hyp=hypgi;
        % x1=x(1,:)';x2=x(2,:)'; sf^2*exp(-1/(2*ell^2)*(x1-x2)'*(x1-x2))

%         cga={'covSEard'};
%         L=[0.5;1];sf=2;hypga=log([L;sf]);cov=cga;hyp=hypga;
        % x1=x(1,:)';x2=x(2,:)';sf^2*exp(-1/2*(x1-x2)'*diag(L)^(-2)*(x1-x2))
        
%         cgf={'covSE','fact',D}; % including complex part, useless recently
%         L=randn(2,D);L=L(:);f=ones(D,1);hypf=log([L;f]);cov=cgf;hyp=hypf;
    case 15 
%         cgv={'covSEvlen',{'meanZero'}};
%         sf=1;hyp_len=[];cov=cgv;hyp=log([hyp_len;sf]);  % equals to cgi={'covSEiso'};ell=1;sf=1;
        
%         cgv={'covSEvlen',{'meanConst'}};
%         sf=1;c=1;cov=cgv;hyp=log([c;sf]);   % equals to cgi={'covSEiso'};ell=c;sf=1; lx=@(x)c;

        cgv={'covSEvlen',{'meanLinear'}};% hyp_len=[1;1]easy;else confused
        sf=1;hyp_len=[1;1];cov=cgv;hyp=log([hyp_len;sf]);
        x1=x(1,:)';x2=x(2,:)';%lx=@(x)exp(x'*flip(hyp_len));% lx is wrong.
        %a=2*lx(x1)*lx(x2);b=lx(x1)^2+lx(x2)^2;
        % mean={@meanLinear};
        % lx1=exp(feval(mean{:},[1;1],x1'));lx2=exp(feval(mean{:},[1;1],x2'));
        % a=2*lx1*lx2;b=lx1^2+lx2^2;
        % sf^2*(a/b)*exp(-(x1-x2)'*(x1-x2)/b)
    case 16
        Q=2;w = ones(Q,1)/Q; m = rand(D,Q); v = rand(D,Q);
        csm={'covSM',Q};    % confused
        hypsm=log([w;m(:);v(:)]);cov=csm;hyp=hypsm;
    case 17
      cov={'covZero'};hyp=[];  
end

%% visualization
set(0,'DefaultFigureWindowStyle','docked') ;
% 1) query th enumber of parameters
feval(cov{:})
% 2) evaluate the function on x, x and xs to get cross-term
[K,dK]=feval(cov{:},hyp,x)              % K(X,X)            n¡Án
[kss,dkss]=feval(cov{:},hyp,xs,'diag')  % k** = k(x*,x*)    n*¡Á1
[Ks,dKs]=feval(cov{:},hyp,x,xs)         % K* = K(X,X*)      n¡Án*

% 3) plot a draw from the kernel
n_xstar=71;
xrange=linspace(-5,5,n_xstar)';
if D~=1
    [a,b]=meshgrid(xrange);
    xstar=[a(:) b(:)];
    
    K0=feval(cov{:},hyp,xstar,[0 0]);    
    figure
    surf(a,b,reshape(K0,n_xstar,n_xstar),'EdgeColor','none',...
        'LineStyle','none','FaceLighting','phong');
    xlabel('x1');ylabel('x2');title('K(X,0)')
    colormap(jet)
    
    K1=feval(cov{:},hyp,xstar);
    K1=K1+(1e-5)*eye(size(K1));
    n_samples=1;
    samples=mvnrnd(zeros(size(a(:))),K1,n_samples)';
    figure
    surf(a,b,reshape(samples,n_xstar,n_xstar))          
    title('the sample based on the covfunc')
    colormap(jet)
else
    K0=feval(cov{:},hyp,xrange,0);
    figure
    plot(xrange,K0,'LineWidth',2);
    
    K1=feval(cov{:},hyp,xrange);
    K1=K1+(1e-5)*eye(size(K1));
    n_samples=1;
    samples=mvnrnd(zeros(n_xstar,1),K1,n_samples)';
    figure
    plot(xrange,samples);
end