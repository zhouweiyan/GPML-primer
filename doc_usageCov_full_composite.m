% demonstrate usage of composite covariance functions
% zhouweiyan 20180912
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

%% present composite covariance library
opt=1;
switch opt
    case 1  % ADD
        %cad={'covADD',{1,{'covMaternard',5}}};  
        % D. K. Duvenaud, H. Nickisch, and C. E. Rasmussen, ¡°Additive Gaussian processes,¡± in International Conference on Neural Information Processing Systems, 2011, pp. 226¨C234.
        %cov=cad;L=[2;1];sf=2;hypma=log([L;sf]);hyp=[hypma;log([1;1])];
        cad = {'covADD',{1,'covSEiso','covLINiso'}};
        cov=cad; hyp=[0;0;0.5;0];
    case 2  % Mask
        mask = [1,0];
        cgi = {'covSEiso'}; 
        ell = 0.9; sf = 2;hypgi = log([ell;sf]);
        cma = {'covMask',{mask,cgi{:}}}; hypma = hypgi;
        cov = cma; hyp = hypma;
%         cov={'covSum',{cma,'covSEiso'}};hyp=[hypma;log([1;0.1])];
    case 3  % Scale: sf^2*k(x,z)
        cgu={'covSEisoU'};ell=1;
        csc={'covScale',cgu};   
        sf=0.5;hypsc=log([ell;sf]);  cov=csc;hyp=hypsc;
    case 4  % Sum: K1+K2+...
        cn  = {'covNoise'}; sn = 0;  hypn = log(sn);
        cc  = {'covConst'};   sf = 0;  hypc = log(sf);
        cl  = {'covLIN'};              hypl = [];
        csu = {'covSum',{cn,cc,cl}}; hypsu = [hypn; hypc; hypl];      % sum
        cov=csu;hyp=hypsu;
    case 5  % PER
        cpi={'covPER','ard',{'covSEisoU'}};
        hyp_p=log([2;4]);
        cov=cpi;hyp=[hyp_p;log(4)];
    case 6  % Prod: K1.*K2.*...
        cpa={'covPoly','ard',[],3};
        c=1;sf=1;L=[2;4];hyppa=log([L;c;sf]);
        cci = {'covPPiso',2}; % compact support poly degree 2
        ell = 2; sf = 2;hypcc = log([ell;sf]); 
        cpr = {@covProd,{cpa,cci}};   hyppr = [hyppa; hypcc];       % product
        cov = cpr; hyp = hyppr;
    case 7
         cpr = {'covPref',{'covSEard'}}; hyppr = log([1;2]);
         cov = cpr; hyp = hyppr;   
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