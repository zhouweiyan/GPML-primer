% Generate some plots showing addictive GP regression on synthetic
% datasets with 1st-order interactions and censored data.
% D. K. Duvenaud, ¡°Automatic Model Construction  with Gaussian Processes,¡± p. 157, 2014.
%
% zhouweiyan 20180724

addpath(genpath('utils'));
addpath(genpath('gpml'));
savefig=false;

% fixing the seed of the random generators
seed=6;
randn('state',seed);
rand('state',seed);

% rendering setup
n_1d=100;   % fineness of grid
% regression setup
n=100;  % number of observations
s=2;    % number of relevent variables
noise_std=0.02; % standard deviation of noise
X=rand(n,2)*4-2;
noiseless_Y=truefunc(X);
Y=noiseless_Y+randn(n,1)*noise_std;

xlims=[-2.3,6];
% xlims=[-12 12];

% censor data
censor_threshold=-1.5;
trainset=or(X(:,1)<censor_threshold,X(:,2)<censor_threshold);
X=X(trainset,:);
y=Y(trainset);
noiseless_Y=noiseless_Y(trainset);
[N,D]=size(X);

% generate a grid on which to render predictive surface
range=linspace(xlims(1),xlims(2),n_1d);
[a,b]=meshgrid(range);
xstar=[a(:) b(:)];

% set up model
likfunc='likGauss';
sn=0.1;hyp.lik=log(sn);
inference=@infExact;
meanfunc={'meanConst'};
hyp.mean=0;
% train additive model
R=2;
covfunc={'covADD',{1:2,'covSEiso'}};    % construct an additive kernel
hyp.cov=log([ones(1,2*D),ones(1,R)]');   % set hyperparameters
hyp_add=minimize(hyp,@gp,-20,inference,meanfunc,covfunc,likfunc,X,y);
% make addictive predictions
predictions_add=gp(hyp,inference,meanfunc,covfunc,likfunc,X,y,xstar);

figure(1);clf;
nice_plot_surface(a,b,predictions_add,length(range),xlims);
title('sum of SE kernels')

% train ARD model
covfunc={'covSEard'};
hyp.cov=log([ones(1,D) 1]');
hyp=minimize(hyp,@gp,-20,inference,meanfunc,covfunc,likfunc,X,y);
prediction_ard=gp(hyp,inference,meanfunc,covfunc,likfunc,X,y,xstar);

figure(2);clf;
nice_plot_surface(a,b,prediction_ard,length(range),xlims);
title('product of SE kernels')

figure
y=sin(xstar(:,1).*2-1.1)+sin(xstar(:,2).*2-1.1);
surf(a,b,reshape(y,length(a),length(b)),'EdgeColor','none','LineStyle','none','FaceLighting','phong');
colormap(jet);

function y=truefunc(x)
    y=sin(x(:,1).*2-1.1)+sin(x(:,2).*2-1.1);
end

function nice_plot_surface(a,b,Y,ell,xlims)
    surf(a,b,reshape(Y,ell,ell), 'EdgeColor','none','LineStyle','none','FaceLighting','phong');
    colormap(jet)
    xlim(xlims);ylim(xlims);zlim('auto');
    set(gcf,'color','white');
    set(get(gca,'XLabel'),'Rotation',0,'Interpreter','tex','Fontsize',12,'Fontname','Times New Roman');
    set(get(gca,'YLabel'),'Rotation',0,'Interpreter','tex','Fontsize',12,'Fontname','Times New Roman');
    set(get(gca,'ZLabel'),'Rotation',0,'Interpreter','tex','Fontsize',12,'Fontname','Times New Roman');
%     set(gca,'xTickLabel',[]);
%     set(gca,'yTickLabel',[]);
%     set(gca,'zTickLabel',[]);
    xlabel('x_1');
    ylabel('y_1');
    zlabel('f(x)');
    view([-30,42]);
end
