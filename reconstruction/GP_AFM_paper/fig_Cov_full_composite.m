% demonstrate usage of composite covariance functions
% zhouweiyan 20180928
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
opt=4;
switch opt
    case 1  % ADD
        %cad={'covADD',{1,{'covMaternard',5}}};  
        % D. K. Duvenaud, H. Nickisch, and C. E. Rasmussen, ��Additive Gaussian processes,�� in International Conference on Neural Information Processing Systems, 2011, pp. 226�C234.
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
%         csu = {'covSum',{'covSEiso','covSEiso'}}; hypsu = [0;0;0.5;0];      % sum
%         cov=csu;hyp=hypsu;
        mask = [0,1];
        cgi = {'covSEiso'};
        cma = {'covMask',{mask,cgi{:}}}; hypma = log([1.5;15]);%log([ell;sf]);
        covfunc1={'covSum',{cma,'covSEiso'}};hyp1.cov=[hypma;log([2;3])];%relative good
        covfunc2={'covMaternard',1};hyp2.cov=log([1;1;4]);
        cov={'covSum',{covfunc1,covfunc2}};hyp=[hyp1.cov;hyp2.cov];
    case 5  % PER
        cpi={'covPER','ard',{'covSEard'}};
        hyp_p=log([500;20;1;1;1;1;1]);
        cov=cpi;hyp=hyp_p;
    case 6  % Prod: K1.*K2.*...
%         cpa={'covPoly','ard',[],3};
%         c=1;sf=1;L=[2;4];hyppa=log([L;c;sf]);
%         cci = {'covPPiso',2}; % compact support poly degree 2
%         ell = 2; sf = 2;hypcc = log([ell;sf]); 
%         cpr = {'covProd',{cpa,cci}};   hyppr = [hyppa; hypcc];       % product
%         cov = cpr; hyp = hyppr;
%         csa={'covSEiso'}; hypsa =log([20;40]);
%         csi={'covSEiso'};hypsi=log([50;1]);  
        csa={'covSEard'}; hypsa =log([20;20;40]);
        csi={'covSEiso'};hypsi=log([5;1]);
        cpr = {'covProd',{csa,csi}};   hyppr = [hypsa; hypsi];       % product
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
[K,dK]=feval(cov{:},hyp,x)  % K: n��n; kss: ns��1; Ks: n��ns
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
    surf(a,b,reshape(samples,n_xstar,n_xstar),'EdgeColor','none',...
        'LineStyle','none','FaceLighting','phong');
    colormap(jet)
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    set(gca,'ZTick',[]);
    tightfig
    set_fig_units_cm(6,4)
    
    K0=feval(cov{:},hyp,xstar,[0 0]);
    figure
    surf(a,b,reshape(K0,n_xstar,n_xstar),'EdgeColor','none',...
        'LineStyle','none','FaceLighting','phong');
    xlabel('x1');ylabel('x2')
    colormap(jet)
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    set(gca,'ZTick',[]);
    tightfig
    set_fig_units_cm(6,4)
else
    K1=feval(cov{:},hyp,xrange);
    K1=K1+(1e-5)*eye(size(K1));
    n_samples=1;
    samples=mvnrnd(zeros(n_xstar,1),K1,n_samples)';
    figure
    plot(xrange,samples);
end