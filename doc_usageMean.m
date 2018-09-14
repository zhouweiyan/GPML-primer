% demonstrate usage of mean functions
%
% See also meanFunctions.m.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2018-06-15.
%                                      File automatically generated using noweb.
% zhouweiyan 20180913
clear
% close all
clc
seed=0;
rand('state',seed);
randn('state',seed);
n = 5; D = 2; x = randn(n,D);            % create a random data set
% n = 5; D = 3; x = randn(n,D);            % create a random data set

%% set up mean functions
opt=12;
switch opt
    % set up simple mean functions
    case 1
        mc = {@meanConst};  hypc = 2;       % also function handles are possible
        mean = mc;  hyp = hypc;
    case 2
        s = 12; hypd = randn(s,1);           % discrete mean with 12 hypers
        md = {'meanDiscrete',s};
        mean = md; hyp = hypd; x = randi([1,s],n,1);
    case 3
        hyp.cov = [0;0]; hypg = [];                    % GP predictive mean
        xt = randn(2*n,D); yt = sign(xt(:,1)-xt(:,2));      % training data
        mg = {@meanGP,hyp,@infEP,@meanZero,@covSEiso,@likErf,xt,yt};
        mean = mg; hyp = hypg;
    case 4
        hype = [0;0; log(0.1)];             % regression GP predictive mean
        xt = randn(2*n,D); yt = xt(:,1).*xt(:,2);           % training data
        me = {@meanGPexact,@meanZero,@covSEiso,xt,yt};
        mean = me; hyp = hype;
    case 5
        ml = {@meanLinear};
        hypl = [2;3];                       % m(x) = 2*x1 + 3*x2; 2*x(:,1)+3*x(:,2)
        % hypl = [2;3;1];
        mean = ml;  hyp = hypl;             % ell
    case 6
        mn = {@meanNN,[0,2;1,2;2,2],[0.9,0.5,0.7]}; hypn = [];  % nearest neighbor
        mean = mn; hyp = hypn;
    case 7
        m1 = {'meanOne'};   hyp1 = [];      % no hyperparameters are needed
        mean = m1;  hyp = hyp1;
    case 8
        mp = {@meanPoly,3}; hypp = [1;1;2;3;1;1];   % m(x) = x1+x2+2*x1^2+3*x2^2;
        mean = mp;  hyp = hypp;                 % x(:,1)+x(:,2)+2*x(:,1).^2+3*x(:,2).^2
    case 9
        d=3;mean={'meanWSPC',d};hyp=[1;2;1;1;1;2;1;1;1;2;1;1];
    case 10
        m0 = {'meanZero'};  hyp0 = [];      % no hyperparameters are needed
        mean = m0;  hyp = hyp0;             % 0) specify mean function

        
    % set up composite mean functions
    case 17 % Scale
        m1 = {'meanOne'};   hyp1 = [];
        msc = {'meanScale',m1};      hypsc = [3; hyp1];      % scale by 3
        mean = msc; hyp = hypsc;
    case 11 % Sum
%         m0 = {'meanZero'};  hyp0 = [];
%         mc = {@meanConst};  hypc = 2;
%         ml = {@meanLinear}; hypl = [2;3]; 
%         msu = {'meanSum',{m0,mc,ml}};  hypsu = [hyp0; hypc; hypl];    % sum
%         mean = msu; hyp = hypsu;
        ml = {@meanLinear};hypl = [2;3];
        mp = {@meanPoly,3}; hypp = [1;1;2;3;1;1];
        mean={'meanSum',{ml,mp}};hyp=[hypl;hypp];
    case 12 % Prod
        mc = {@meanConst};  hypc = 2;
        ml = {@meanLinear}; hypl = [0.2;0.3];
        mpr = {@meanProd,{mc,ml}};     hyppr = [hypc; hypl];      % product
        mean = mpr; hyp = hyppr;
    case 13
        m0 = {'meanZero'};  hyp0 = [];
        mc = {@meanConst};  hypc = 2;
        ml = {@meanLinear}; hypl = [2;3]; 
        msu = {'meanSum',{m0,mc,ml}};  hypsu = [hyp0; hypc; hypl];    % sum
        mpo = {'meanPow',3,msu};       hyppo = hypsu;         % third power
        mean = mpo; hyp = hyppo;
    case 14
        ml = {@meanLinear}; hypl = [2;3];
        mask = [false,true];     % mask excluding all but the 2nd component
        mma = {'meanMask',mask,ml};    hypma = hypl(mask);
        mean = mma; hyp = hypma;
    case 15
        ml = {@meanLinear}; 
        mpf = {@meanPref,ml};          hyppf = 2;  % linear pref with slope
        mean = mpf; hyp = hyppf;
        % 2*(x(:,1)-x(:,2))
    case 16
        ml = {@meanLinear};
        mwp = {@meanWarp,ml,@sin,@cos};
        hypwp = [2;1]; % sin of linear, meanWarp(mean, g, dg, hyp, x)
        mean = mwp; hyp = hypwp;
end

%% 1) query the number of parameters
feval(mean{:})

% plot a draw from the kernel
n_samples=1;
n_xstar=51;
xrange=linspace(-3,3,n_xstar)';
[a,b]=meshgrid(xrange);
m=feval(mean{:},hyp,[a(:) b(:)]);
figure
surf(a,b,reshape(m,n_xstar,n_xstar))
colormap(jet)

% 2) evaluate the function on x
feval(mean{:},hyp,x)

% 3) compute the derivatives w.r.t. to hyperparameter i
% i = 2; feval(mean{:},hyp,x,i)
