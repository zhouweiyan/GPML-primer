% A script to make a series of plots demonstracting how structure can be
% reflected in a kernel.
% This version updated to be more amenable to tiny plots.
% D. K. Duvenaud, ¡°Automatic Model Construction  with Gaussian Processes,¡± p. 157, 2014.
% 
% zhouweiyan 20180726

close all
clear
clc
seed=1;
randn('state',seed);
rand('state',seed);
set(0,'DefaultFigureWindowStyle','docked')

%% make up some data
X=[1 2 3]'.*2;
y=[1 2 4]';
y=y-mean(y);
N=length(X);

n_samples=2;
n_xstar=201;
xrange=linspace(-10,10,n_xstar)';
post_xrange=linspace(-3,10,n_xstar)';
x0=1;
numerical_noise=1e-5;
model_noise=2;

% SEiso
se_length_scale=2.5;
se_output_var=2;
se_kernel=@(x,y)se_output_var*exp(-1/(2*se_length_scale^2)*(x-y).^2);

% LINone
lin_output_var=0.5;
lin_kernel=@(x,y)lin_output_var*(x.*y+1);

% various SEiso
longse_length_scale=20;
longse_output_var=10;
longse_kernel=@(x,y)longse_output_var*exp(-1/(2*longse_length_scale^2)*(x-y).^2);
shortse_length_scale = 0.5;
shortse_output_var = 0.5;
shortse_kernel = @(x,y) shortse_output_var*exp(-1/(2*shortse_length_scale^2)*(x-y).^2);
medse_length_scale = 10;
medse_output_var = 5;
medse_kernel = @(x,y) medse_output_var*exp(-1/(2*medse_length_scale^2)*(x-y).^2);

% Periodic, x belongs to R^1
per_length_scale=1;
per_period=4;
per_output_var=1.1;
per_kernel=@(x,y)per_output_var*exp(-2/per_length_scale^2*(sin(pi*(x-y)./per_period).^2));

% Cos, x belongs to R^1
cos_period=4;
cos_output_var=1.1;
cos_kernel=@(x,y)cos_output_var*cos(2*pi*(x-y)./cos_period);

% Noise
wn_output_var=1.1;
wn_kernel=@(x,y)wn_output_var*double(x==y);

% RQiso
rq_length_scale=2.5;
rq_output_var=2;
rq_alpha=1.1;
rq_kernel=@(x,y)rq_output_var*(1+(x-y).^2./(2*rq_alpha*rq_length_scale^2)).^(-rq_alpha);

% One
c_kernel=@(x,y)ones(size(x*y));

se_plus_lin=@(x,y)se_kernel(x,y)+lin_kernel(x,y);
se_plus_per=@(x,y)se_kernel(x,y)+per_kernel(x,y);
se_times_lin=@(x,y)se_kernel(x,y).*lin_kernel(x,y);
se_times_per=@(x,y)se_kernel(x,y).*per_kernel(x,y);
lin_plus_per=@(x,y)lin_kernel(x,y)+per_kernel(x,y);
lin_times_per=@(x,y)lin_kernel(x,y).*per_kernel(x,y);
lin_times_lin=@(x,y)lin_kernel(x,y).*lin_kernel(x,y);
longse_plus_per=@(x,y)longse_kernel(x,y)+per_kernel(x,y);
longse_times_per=@(x,y)longse_kernel(x,y).*per_kernel(x,y);
longse_times_lin=@(x,y)longse_kernel(x,y).*lin_kernel(x,y);
longse_plus_se=@(x,y)longse_kernel(x,y)+se_kernel(x,y);
shortse_plus_medse=@(x,y)shortse_kernel(x,y)+medse_kernel(x,y);

%% calculate kernels
% kernel_names={'se_kernel', 'lin_kernel', 'longse_kernel','shortse_kernel','medse_kernel', ...
%             'per_kernel','cos_kernel','wn_kernel','rq_kernel','c_kernel'};
kernel_names={'se_plus_lin', 'se_plus_per', 'se_times_lin', 'se_times_per', ...
           'lin_plus_per','lin_times_per',  'lin_times_lin', ...
            'longse_plus_per','longse_times_per', 'longse_times_lin', 'longse_plus_se',...
            'shortse_plus_medse'};

% automatically build kernel names from function name
for i=1:numel(kernel_names)
    kernels{i}=eval(kernel_names{i});
end

color_ixs=repmat(10,1,numel(kernels));
if 1
    % plot the kernel
    for k=1:numel(kernels)
        figure(k);clf;
        cur_kernel=kernels{k};
        kvals=bsxfun(cur_kernel,xrange,x0);
        kernel_plot(xrange,kvals,color_ixs(k));
        title(kernel_names{k},'Interpreter','none');
    end
end
% plot draws from kernels
for k=1:numel(kernels)
    figure;clf;
    K=bsxfun(kernels{k},xrange',xrange)+eye(n_xstar).*numerical_noise;
    samples=mvnrnd(zeros(size(xrange)),K,n_samples)';
    samples_plot(xrange,samples,[1:n_samples]);
    title(['draws from ',kernel_names{k}],'Interpreter','none');
end

function kernel_plot(xrange,vals,color_ix)
    % figure settings
    lw=2;
    fontsize=10;
    plot(xrange,vals,'Color',colorbrew(color_ix),'LineWidth',lw);
    hold on
    if all(vals>=0)
        lowlim=0;
    else
        lowlim=min(vals);
    end
    % make plot prettier
    xlim([min(xrange),max(xrange)]);
    ylim([min(lowlim,0),max(vals)*1.05]);
%     set(gca,'XTick',0);
%     set(gca,'YTick',0);
%     set(gca,'XTickLabel','');
    set(get(gca,'XLabel'),'Rotation',0,'Interpreter','latex','Fontsize',fontsize);
    set(get(gca,'YLabel'),'Rotation',90,'Interpreter','latex','Fontsize',fontsize);
    set(gca,'Box','off');         % delete the right and up line
    set(gca,'color','white');
    set(gca,'YGrid','off');
end

function samples_plot(xrange,samples,color_ix)
    % figure setting
    lw=2;
    fontsize=10;
    for i=1:size(samples,2)
        plot(xrange,samples(:,i),'Color',colorbrew(color_ix(i)),'LineWidth',lw);
        hold on
    end
    % make plot prettier
    xlim([min(xrange),max(xrange)]);
%     set(gca,'XTick',0);
%     set(gca,'xTickLabel',[]);
%     set(gca,'YTick',[]);
%     set(gca,'yTickLabel',[]);
    set(get(gca,'XLabel'),'Rotation',0,'Interpreter','latex', 'Fontsize', fontsize);
    set(get(gca,'YLabel'),'Rotation',90,'Interpreter','latex', 'Fontsize', fontsize);
    set(gca,'Box','off');
    set(gcf,'color','white');
    set(gca,'YGrid','off');
end











        








