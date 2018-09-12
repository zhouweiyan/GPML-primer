% demonstrate usage of simple covariance functions
% zhouweiyan 20180911
clear
% close all
clc
%% initialization
n=5;D=2;
seed=0;
rand('state',seed);
randn('state',seed);
x=randn(n,D);
xs=randn(3,D);

%% simple covariance library
opt=12;
switch opt
    case 1
        co={'covOne'};
        cov=co;hyp=[];
%         cc={'covConst'};
%         sf=1;hypc=log(sf);cov=cc;hyp=hypc;
    case 2
        cc={'covCos'};  % x belongs to R1
        p=1;sf=1;hypcc=log([p;sf]);   cov=cc;hyp=hypcc;
    case 3
        ce={'covEye'};  % k(x^p,x^q) = \delta(p,q)
        hype=[];cov=ce;hyp=hype;
    case 4
%         cgbi={'covGaboriso'};
%         ell=5;p=3;hypgbi=log([ell;p]);  cov=cgbi;hyp=hypgbi;
        cgba={'covGaborard'};
        ell=[5;10];p=[1;6];hypgba=log([ell;p]);  cov=cgba;hyp=hypgba;
    case 5
%         cge={'covGE','eye',[]};
%         gamma=1.5;hypge=log(gamma/(2-gamma)); cov=cge;hyp=hypge;
%         cgi={'covGE','iso',[]};
%         gamma=1.5;ell=1;hypgi=log([ell;(gamma/(2-gamma))]); 
%         cov=cgi;hyp=hypgi;
        cga={'covGE','ard',[]};
        gamma=1.5;L=(1:D)';hypga=log([L;(gamma/(2-gamma))]); 
        cov=cga;hyp=hypga;
    case 6
%         cl={'covLIN'};
%         hypl=[];cov=cl;hyp=hypl;
%         cli={'covLINiso'};
%         ell=0.5; hypli=log(ell);cov=cli;hyp=hypli;  % slope is inversely propotional to ell 
        cla={'covLINard'};
        L=(D:-1:1)';hypla=log(L);cov=cla;hyp=hypla;
%         clo={'covLINone'};    % K=(x*x'+1)/ell^2;
%         ell=0.9;hyplo=log(ell);cov=clo;hyp=hyplo;
    case 7
%         cmi={'covMaterniso',5}; 
%         L=0.5;sf=2;hypmi=log([L;sf]);  cov=cmi;hyp=hypmi;
        cma={'covMaternard',5};
        L=[2;1];sf=2;hypma=log([L;sf]); cov=cma;hyp=hypma;
    case 8
        cn={'covNoise'};% k(x^p,x^q) = sf^2 * \delta(p,q)
        sf=1;hypn=log(sf);cov=cn;hyp=hypn;
    case 9
        cnn={'covNNone'};   % the outline shape remains unchanged, seems useless for me
        ell=1;sf=1;hypnn=log([ell;sf]);   cov=cnn;hyp=hypnn;
    case 10
        cpe={'covPeriodic'};    % x belongs to R1
        ell=0.5;p=2;sf=1;hyppe=log([ell;p;sf]); cov=cpe;hyp=hyppe;
%         cpn={'covPeriodicNoDC'};
%         ell=0.9;p=2;sf=1;hyppn=log([ell;p;sf]); cov=cpn;hyp=hyppn;
    case 11
%         cp={'covPoly',3};   % K=sf^2*(c+x*x').^3
%         c=2;sf=1;hypp=log([c;sf]);cov=cp;hyp=hypp;
%         cpi={'covPoly','iso',[],2};
%         ell=1;c=2;sf=1;cov=cpi;hyp=log([ell;c;sf]);
        cpa={'covPoly','ard',[],3};
        c=1;sf=1;L=[2;4];cov=cpa;hyp=log([L;c;sf]); %k(x,z) = sf^2*(c+s)^d , where s=x*inv(P)*z
    case 12
        cpi={'covPPiso',2};   % [varargout{:}] = covScale({'covPP','iso',[],v},varargin{:});
        ell=10;sf=2;hyppi=log([ell;sf]);   cov=cpi;hyp=hyppi;
%         cpa={'covPPard',3};
%         L=[2;1];sf=1;hyppa=log([L;sf]);   cov=cpa;hyp=hyppa;
%         cpi={'covPP','iso',[],2};
%         ell=2;hyppi=log(ell);   cov=cpi;hyp=hyppi;
%         cpa={'covPP','ard',[],2};
%         L=[2;1];hyppa=log(L);   cov=cpa;hyp=hyppa;
    case 13
%         cri={'covRQiso'};
%         ell=0.9;sf=2;al=1;hypri=log([ell;sf;al]);   cov=cri;hyp=hypri;
        cra={'covRQard'};   % k(x,z)=sf^2*(1+1/(2*al)*(x-z)'*diag(L)^(-2)*(x-z));
        L=[0.5;1];sf=2;al=1;hypra=log([L;sf;al]);cov=cra;hyp=hypra;
    case 14
%         cgu={'covSEisoU'};  % k(x,z) = exp(-(x-z)'*inv(P)*(x-z)/2); P=ell^2*I;
%         ell=1/sqrt(2);hypgu=log(ell);cov=cgu;hyp=hypgu;
        cgi={'covSEiso'};   % k(x,z) = sf^2 * exp(-(x-z)'*inv(P)*(x-z)/2)
        ell=1;sf=0.5;hypgi=log([ell;sf]);cov=cgi;hyp=hypgi;
%         cga={'covSEard'};
%         L=[1;0.5];sf=2;hypga=log([L;sf]);cov=cga;hyp=hypga;
%         cgf={'covSE','fact',D}; % including complex part, useless recently
%         L=randn(2,D);L=L(:);f=ones(D,1);hypf=log([L;f]);cov=cgf;hyp=hypf;
    case 15 
%         cgv={'covSEvlen',{'meanZero'}};
%         sf=1;hyp_len=[];cov=cgv;hyp=log([hyp_len;sf]);  % equals to cgi={'covSEiso'};ell=1;sf=1;
%         cgv={'covSEvlen',{'meanConst'}};
%         sf=1;c=0.5;cov=cgv;hyp=log([c;sf]);   % equals to cgi={'covSEiso'};ell=c;sf=1;
        cgv={'covSEvlen',{'meanLinear'}};% hyp_len=[1;1]easy;else confused
        sf=1;hyp_len=[2;0.1];cov=cgv;hyp=log([hyp_len;sf]);
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
[K,dK]=feval(cov{:},hyp,x)  % K: n¡Án; kss: ns¡Á1; Ks: n¡Áns
[kss,dkss]=feval(cov{:},hyp,xs,'diag')
[Ks,dKs]=feval(cov{:},hyp,x,xs)

% 3) plot a draw from the kernel
n_xstar=71;
xrange=linspace(-5,5,n_xstar)';
if D~=1
    [a,b]=meshgrid(xrange);
    xstar=[a(:) b(:)];
    K1=feval(cov{:},hyp,xstar);
    K1=K1+(1e-5)*eye(size(K1));
    n_samples=1;
    samples=mvnrnd(zeros(size(a(:))),K1,n_samples)';
    figure
    surf(a,b,reshape(samples,n_xstar,n_xstar))
    colormap(jet)
    
    K0=feval(cov{:},hyp,xstar,[0 0]);
    figure
    surf(a,b,reshape(K0,n_xstar,n_xstar),'EdgeColor','none',...
        'LineStyle','none','FaceLighting','phong');
    colormap(jet)
else
    K1=feval(cov{:},hyp,xrange);
    K1=K1+(1e-5)*eye(size(K1));
    n_samples=1;
    samples=mvnrnd(zeros(n_xstar,1),K1,n_samples)';
    figure
    plot(xrange,samples);
end