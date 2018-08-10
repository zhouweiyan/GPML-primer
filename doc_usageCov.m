% demonstrate usage of covariance functions
%
% See also covFunctions.m.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2018-06-15.
%                                      File automatically generated using noweb.
clear  
% close all
clc
addpath('utils/');
addpath('gpml/');
% n = 5; D = 3; x = randn(n,D); xs = randn(3,D);  % create a data set
n = 5; D = 2; 
% rand_state=rng;                                 % zwy  
seed=0;
rand('state',seed);
randn('state',seed);
x = randn(n,D); xs = randn(3,D);                
% rng(rand_state);                                % zwy
% x = randn(n,D); xs = randn(3,D);                % zwy
opt = 10
%% set up covariance functions
switch opt
    % set up simple covariance functions
    case 1
        cn  = {'covNoise'}; sn = 50;  hypn = log(sn);  % one hyperparameter
        % 0) specify a covariance function
        cov = cn; hyp = hypn;
        % K= sn^(-2)*eye(n,n)
        % kss=sn^(-2)*ones(ns,1)
        % Ks=zeros(n,ns)
    case 2
        cc  = {@covConst};   sf = 2;  hypc = log(sf); % function handles OK
        cov = cc; hyp = hypc;
        % K= sf^2*ones(n,n)
        % kss=sf^2*ones(ns,1)
        % Ks=ones(n,ns)
    case 3
        ce  = {@covEye};              hype = [];                 % identity
        cov = ce; hyp = hype;
        % K = eye(n);   kss = ones(n,1);
        % Ks = double(sq_dist(x',z')<tol*tol);
    case 4
        cl  = {@covLIN};              hypl = []; % linear is parameter-free
        % K=x*x';  kss=sum(xs.*xs,2);   Ks=x*xs';
        cov = cl; hyp = hypl;
    case 5
        cla = {'covLINard'}; L = rand(D,1); hypla = log(L);  % linear (ARD)
        cov = cla; hyp = hypla;
        % K=x*(x*diag(exp(-2*hyp)))'
        % kss=sum(xs.*(xs*diag(exp(-2*hyp))), 2)
        % Ks=x*(xs*diag(exp(-2*hyp)))'
    case 6
        cli = {'covLINiso'}; ell = rand(1);   hypli = log(ell);    % linear iso
        cov = cli; hyp = hypli;
        % K=x*(x*ell^-2)'  %(l:ell)
        % kss=sum(xs.*(xs*ell^-2)), 2)
        % Ks=x*(xs*l^-2)'
    case 7
        clo = {@covLINone}; ell = 0.9; hyplo = log(ell);  % linear with bias
        cov = clo; hyp = hyplo;
        % K=(x*x'+1)/ell^2;
        % kss=(sum(xs.*xs,2)+1)/ell^2
        % Ks=(x*xs'+1)/ell^2
    case 8
        cp  = {@covPoly,3}; c = 2; sf = 2; hypp = log([c;sf]);   % third order poly
        cov = cp; hyp = hypp;
        % K=sf^2*(c+x*x').^3
        % kss=sf^2*(2+sum(xs.*xs,2)).^3
        % Ks=sf^2*(c+x*xs').^3
    case 9
        cga = {@covSEard};   
        L = rand(D,1);sf = 2; hypga = log([L;sf]);       % Gaussian with ARD
        cov = cga; hyp = hypga;
        % k(x,z)=sf^2*exp(-1/2*(x-z)'*Lamda^(-2)*(x-z));
        % hyp=[ln(lamda_1),ln(lamda_2),...,ln(lamda_D),ln(sf)]'
    case 10
        cgi = {'covSEiso'}; 
        ell = 0.2; sf = 50;hypgi = log([ell;sf]);    % isotropic Gaussian
        cov = cgi; hyp = hypgi;
        % k(x,z)=sf^2*exp(-1/(2*ell^2)*(x-z)'*(x-z));
        % e.g. K(1,2)=sf^2*exp(-1/(2*ell^2)*(x(1,:)-x(2,:))*(x(1,:)-x(2,:))')
        % Ks(1,2)=sf^2*exp(-1/(2*ell^2)*(x(1,:)-xs(2,:))*(x(1,:)-xs(2,:))')
        % hyp=[ln(ell),ln(sf)]'
    case 11
        cgu = {'covSEisoU'}; ell = 0.9;hypgu = log(ell);   % isotropic Gauss no scale
        cov = cgu; hyp = hypgu;
        % k(x,z)=exp(-1/(2*ell^2)*(x-z)'*(x-z));
        % e.g. K(1,2)= exp(-1/(2*ell^2)*(x(1,:)-x(2,:))*(x(1,:)-x(2,:))')
        % Ks(1,2)=exp(-1/(2*ell^2)*(x(1,:)-xs(2,:))*(x(1,:)-xs(2,:))')
        % hyp=[ln(ell)]
    case 12
        cra = {'covRQard'}; 
        L = rand(D,1);sf = 2;al = 2; hypra = log([L;sf;al]); % ration. quad.
        cov=cra; hyp=hypra;
        % k(x,z)=sf^2*(1+1/(2*al)*(x-z)'*diag(L)^(-2)*(x-z));
        % e.g. K(2,4)=sf^2*(1+1/(2*al)*(x(2,:)-x(4,:))*diag(L)^(-2)*(x(2,:)-x(4,:))')^(-al)
        % hyp=[ln(lamda_1),ln(lamda_2),...,ln(lamda_D),ln(sf),ln(al)]'
    case 13
        cri = {@covRQiso};
        ell = 0.9; sf = 2;al = 2;hypri = log([ell;sf;al]);   % isotropic
        cov=cri;hyp=hypri;
        % k(x,z)=sf^2*(1+1/(2*al*ell^2)*(x-z)'*(x-z));
        % e.g. K(2,4)=sf^2*(1+1/(2*al*ell^2)*(x(2,:)-x(4,:))*(x(2,:)-x(4,:))')^(-al)
        % hyp=[ln(ell),ln(sf),ln(al)]'
    case 14
        cma = {@covMaternard,5};  % Matern class d=5
        L = rand(D,1); sf = 2; hypma = log([L;sf]);
        cov = cma; hyp = hypma;
        % d=5; f_5(t)=1+t+t^2/3
        % k(x,z)=sf^2*f_d(r_d)exp(-r_d), r_d=sqrt(d*(x-z)'*diag(Lamda)^(-2)*(x-z))
        % e.g. K(1,2)  rd=sqrt(d*(x(1,:)-x(2,:))*diag(L)^(-2)*(x(1,:)-x(2,:))')
        % sf^2*(1+rd+rd^2/3)*exp(-rd)
    case 15
        cmi = {'covMaterniso',3};  % Matern class d=3
        ell = 0.9;sf = 3;hypmi = log([ell;sf]);
        cov=cmi;hyp=hypmi;
        % d=3; f_3(t)=1+t
        % k(x,z)=sf^2*f_d(r_d)exp(-r_d), r_d=sqrt(d/ell^2*(x-z)'*(x-z))
        % e.g. K(1,2): rd=sqrt(d/ell^2*(x(1,:)-x(2,:))*(x(1,:)-x(2,:))')
        % sf^2*(1+rd)*exp(-rd)
    case 16
        cnn = {'covNNone'};     % neural network
        %L = rand(D,1);sf = 2;hypnn = log([L;sf]); 
        % the first paramater is accepted as log(ell)
        % the second is accepted as log(sf), the redundancy are omitted
        ell = 1;sf = 2;hypnn = log([ell;sf]);   
        cov=cnn;hyp=hypnn;
        % k(x,z) = sf2 * asin(x'*P*z / sqrt[(1+x'*P*x)*(1+z'*P*z)])
        % where the x and z vectors on the right hand side have an added extra bias entry with unit value.
        % P=ell^-2*eye(D+1);
        % K(1,2)=sf^2*asin([x(1,:) 1]*P*[x(2,:) 1]'/sqrt((1+[x(1,:) 1]*P*[x(1,:) 1]')*(1+[x(2,:) 1]*P*[x(2,:) 1]')))
    case 17
        cpe = {'covPeriodic'}; 
        p = 2;ell = 0.9;sf = 2; hyppe = log([ell;p;sf]);   % periodic
        cov=cpe; hyp=hyppe;
        % k(x,z) = sf^2 * exp(-2*sin^2( pi*(x-z)/p )/ell^2)
        % e.g. K(1,2)=sf^2 * exp(-2*(sin( pi*(x(1)-x(2))/p))^2/ell^2)
    case 18
        cpn = {'covPeriodicNoDC'}; 
        p = 2;ell = 0.9;sf = 2; hyppe = log([ell;p;sf]); % w/o DC
        cov=cpn; hyp=hyppe;
    case 19
        cpc = {'covCos'}; p = 2;sf = 3; hypcpc = log([p;sf]);         % cosine cov
        cov=cpc; hyp=hypcpc;
        % k(x,z) = sf^2*cos(2*pi*(x-z)/p)
    case 20
        cca = {'covPPard',3}; 
        % L = rand(D,1);
        L=[3;6];
        sf = 2;hypcc = log([L;sf]);% compact support poly degree 3
        cov=cca; hyp=hypcc;
    case 21
        cci = {'covPPiso',2}; % compact support poly degree 2
        ell = 0.9; sf = 2;hypcc = log([ell;sf]); 
        cov=cci; hyp=hypcc;
    case 22
        cgb = {'covGaboriso'}; ell = 1; p = 1.2; hypgb=log([ell;p]); % Gabor
        cov=cgb; hyp=hypgb;
        % t=x(1,:)-x(2,:);
        % exp(-t*t'/(2*ell^2))*cos(2*pi*t*ones(D,1)/p)
    case 23
        Q = 2; w = ones(Q,1)/Q; m = rand(D,Q); v = rand(D,Q);
        csm = {@covSM,Q}; hypsm = log([w;m(:);v(:)]);    % Spectral Mixture
        cov = csm; hyp = hypsm;
    case 24
        % long time for thinking the underlying function
        cvl = {@covSEvlen,{@meanLinear}}; hypvl = [1;0.1;1]; % var lenscal
        cov = cvl; hyp = hypvl;
    case 25
        s = 12; cds = {@covDiscrete,s};      % discrete covariance function
        L = randn(s); L = chol(L'*L); L(1:(s+1):end) = log(diag(L));
        hypds = L(triu(true(s))); xd = randi([1,s],[n,1]); xsd = [1;3;6];
        cov = cds; hyp = hypds; x = xd; xs = xsd;
    case 26     % 未定义函数或变量 'covSEfact'。
        cfa = {@covSEfact,2}; hypfa = randn(D*2,1);       % factor analysis 
        cov = cfa; hyp = hypfa;
    
    % set up composite i.e. meta covariance functions
    case 27
        cgu = {'covSEisoU'}; ell=3;sf = 0.9;  % isotropic Gauss no scale
        csc = {'covScale',cgu}; hypsc = [log(ell); log(sf)];  % scale by 0.9
        cov = csc; hyp = hypsc;
        % k(x,z)=sf^2*exp(-1/(2*ell^2)*(x-z)'*(x-z));
        % e.g. K(1,2)=sf^2*exp(-1/(2*3^2)*(x(1,:)-x(2,:))*(x(1,:)-x(2,:))')
    case 28
        cn  = {'covNoise'}; sn = 0.1;  hypn = log(sn);
        cc  = {@covConst};   sf = 2;  hypc = log(sf);
        cl  = {@covLIN};              hypl = [];
        csu = {'covSum',{cn,cc,cl}}; hypsu = [hypn; hypc; hypl];      % sum
        cov=csu;hyp=hypsu;
        % K1=feval(cn{:},hypn,x);
        % K2=feval(cc{:},hypc,x);
        % K3=feval(cl{:},hypl,x);
        % K1+K2+K3
    case 29
        cc  = {@covConst};   sf = 2;  hypc = log(sf);
        cci = {'covPPiso',2}; % compact support poly degree 2
        ell = 0.9; sf = 2;hypcc = log([ell;sf]); 
        cpr = {@covProd,{cc,cci}};   hyppr = [hypc; hypcc];       % product
        cov = cpr; hyp = hyppr;
        % K1=feval(cci{:},hypcc,x)
        % K2=feval(cc{:},hypc,x)
        % K2.*K1 or K1.*K2
    case 30
        %mask = [0,1,0]; %   binary mask excluding all but the 2nd
        %component for 3D x
        mask = [0,1];
        cgi = {'covSEiso'}; 
        ell = 0.9; sf = 2;hypgi = log([ell;sf]);
        cma = {'covMask',{mask,cgi{:}}}; hypma = hypgi;
        cov = cma; hyp = hypma;
        % x0=zeros(size(x));
        % x0(:,2)=x(:,2)
        % feval(cov{:},hyp,x0)
    case 31
        % isotropic periodic rational quadratic
        cpi = {'covPERiso',{@covRQiso}};
    case 32
        % periodic Matern with ARD
        cpa = {'covPERard',{@covMaternard,3}};
    case 33
        % ? additive based on SEiso using unary and pairwise interactions
        cad = {'covADD',{[1,2],'covSEiso'}};
        cov=cad; hyp=[0;0;0;0;0;0];
        K= feval(cov{:},hyp,x)
        % Please see the paper "Additive Gaussian Processes" by Duvenaud,
        % Nickisch and Rasmussen, NIPS, 2011 for details.
    case 34
        % preference covariance with squared exponential base covariance
        cpr = {'covPref',{'covSEiso'}}; hyppr = [0;0];
        xp = randn(n,2*D); xsp = randn(3,2*D);
        cov = cpr; hyp = hyppr; x = xp; xs = xsp;
%         cov0={'covSEiso'};
%         K11=feval(cov0{:},hyp,x(:,1:D));
%         K22=feval(cov0{:},hyp,x(:,(D+1):(2*D)),x(:,(D+1):(2*D)));
%         K12=feval(cov0{:},hyp,x(:,1:D),x(:,(D+1):(2*D)));
%         K21=feval(cov0{:},hyp,x(:,(D+1):(2*D)),x(:,1:D));
%         K11+K22-K12-K21
    case 35
%         cpi={'covPER','eye',{@covSEisoU}};
%         hyp_p=[]; ell = 0.9; sf = 2;%hypgi = log([ell;sf]);
%         hypgi = log([ell]);
%         cov=cpi; hyp=[hyp_p hypgi];
%         % p=ones(D,1);
%         % x1=x(1,:); x1=x1./p'; x2=x(2,:); x2=x2./p';
%         % ux=[sin(2*pi*x1) cos(2*pi*x1)];uz=[sin(2*pi*x2) cos(2*pi*x2)];
%         % exp(-1/(2*ell^2)*(ux-uz)*(ux-uz)')

%         cpi={'covPER','eye',{@covSEiso}};
%         hyp_p=[]; ell = 0.9; sf = 2; hypgi = log([ell;sf]);
%         cov=cpi; hyp=[hyp_p hypgi];
%         % p=ones(D,1);
%         % x1=x(1,:); x1=x1./p'; x2=x(2,:); x2=x2./p';
%         % ux=[sin(2*pi*x1) cos(2*pi*x1)];uz=[sin(2*pi*x2) cos(2*pi*x2)];
%         % sf^2*exp(-1/(2*ell^2)*(ux-uz)*(ux-uz)')

        cpi={'covPER','iso',{@covSEiso}};
        %cpi={'covPERiso',{@covSEiso}};
        hyp_p=log(2); ell = 0.9; sf = 2; hypgi = log([ell;sf]);
        cov=cpi; hyp=[hyp_p;hypgi];
        % p = exp(hyp_p)*ones(D,1);
        % x1=x(1,:); x1=x1./p'; x2=x(2,:); x2=x2./p';
        % ux=[sin(2*pi*x1) cos(2*pi*x1)];uz=[sin(2*pi*x2) cos(2*pi*x2)];
        % sf^2*exp(-1/(2*ell^2)*(ux-uz)*(ux-uz)')
end

%% 1) query the number of parameters
feval(cov{:})
% 2) evaluate the function on x
[K,dK] = feval(cov{:},hyp,x)    % K: n×n; kss: ns×1; Ks: n×ns

% plot a draw from the kernel
n_samples=1;
n_xstar=71;
xrange=linspace(-5,5,n_xstar)';
[a b]=meshgrid(xrange);
K=feval(cov{:},hyp,[a(:) b(:)]);
K=K+(1e-5)*eye(size(K));
samples=mvnrnd(zeros(size(a(:))),K,n_samples)';
figure
surf(a,b,reshape(samples,n_xstar,n_xstar))
colormap(jet)

% 3) evaluate the function on x and xs to get cross-terms
[kss,dkss] = feval(cov{:},hyp,xs,'diag')
% [kss,dkss] = feval(cov{:},hyp,xs,'diag')    % zwy
[Ks, dKs ]  = feval(cov{:},hyp,x,xs)
