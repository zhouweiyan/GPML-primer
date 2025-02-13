% reconstruct complex morphology
% spiral
% to enlarge the visual contrast
% zwy 20190223 give up

clear
% close all
% clc

%% Constant Angular Velocity(CAV)
r_end=64*1.6;
num_curves = [62;52;42;32;22];
PSNR = zeros(1,length(num_curves));
for i = 1:1%length(num_curves)
    i
    load('E:\OneDrive - hnu.edu.cn\tools\matlabcourse\GPML_matlab\GPML-primer\reconstruction\GP_AFM_Pp_part3\result\Copy_of_EChemist_Corrosion_GP.mat')
    opt='GP';
%     opt='DT';
    if i == length(num_curves)
       r_end = 64*1.7;
    end
    P=r_end*2/(num_curves(i)-1);
    %P=spiral radius��2/(number of curves-1)
    %eq.(4)
    w=120;
    f=3000;
    T=1/f;
    T_spiral=2*pi*r_end/(P*w);	%CAV eq.(6)
    t=0:T:T_spiral;
    xs=P*w/(2*pi)*t.*cos(w*t)+64;	%eq.(17)
    ys=P*w/(2*pi)*t.*sin(w*t)+64;	%eq.(18)
    %% data preparation
    % train set
    x1_train_spi=xs';x2_train_spi=ys';
    X_train_spi=[x1_train_spi(:),x2_train_spi(:)];
    [y_train_spi,y_t_s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_train_spi);

    % approximate path
    path_len=sum(sqrt((x1_train_spi(2:end)-x1_train_spi(1:(end-1))).^2+(x2_train_spi(2:end)-x2_train_spi(1:(end-1))).^2));
    delta = path_len/(2*128*128);
    
    % test set
    x1_test=1:128;x2_test=1:128;
    [X1_test,X2_test]=meshgrid(x1_test,x2_test);
    X_test=[X1_test(:),X2_test(:)];
    [y_test_ideal,y_test_s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_test);
    y_test_ideal=reshape(y_test_ideal,length(x2_test),length(x1_test));
    
    %% GP regression/DT linear interplotation
    switch opt
        case 'GP'
            meanfunc=[];hyp.mean=[];
            mask = [0,1];
            cgi = {'covSEiso'};
            cma = {'covMask',{mask,cgi{:}}}; hypma = log([1.5;15]);%log([ell;sf]);
            covfunc1={'covSum',{cma,'covSEiso'}};hyp1.cov=[hypma;log([2;3])];%relative good
            covfunc2={'covMaternard',1};hyp2.cov=log([1;1;4]);
            covfunc={'covSum',{covfunc1,covfunc2}};hyp.cov=[hyp1.cov;hyp2.cov];
            
            likfunc='likGauss';hyp.lik=log(10);    % log(sn)
            hyp=minimize(hyp,@gp,-80,@infGaussLik,meanfunc,covfunc,likfunc,X_train_spi,y_train_spi);
            [m,s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train_spi,y_train_spi,X_test);
            
        case 'DT'
            P=X_train_spi;
            DT=delaunayTriangulation(P);
            % triplot(DT)
            V=y_train_spi;
            % approximate path
            path_len=sum(sqrt((x1_train_spi(2:end)-x1_train_spi(1:(end-1))).^2+(x2_train_spi(2:end)-x2_train_spi(1:(end-1))).^2))
            Pq=X_test;
            
            [ti,bc]=pointLocation(DT,Pq);
            triVals=V(DT(ti,:));
            Vq=dot(bc',triVals')';
            m=Vq;
    end
    
    %% Error analysis
    m_a=reshape(m,length(x2_test),length(x1_test));
    m_a=m_a(min(x2_train):max(x2_train),min(x1_train):max(x1_train));
    y_test_ideal_a=y_test_ideal(min(x2_train):max(x2_train),min(x1_train):max(x1_train));
    PSNR(i)=max(y_test_ideal_a(:))*sqrt(size(y_test_ideal_a,1)*size(y_test_ideal_a,2))/sqrt(sum((m_a(:)-y_test_ideal_a(:)).^2));
    PSNR(i)=20*log(PSNR(i))/log(10)
    %%
    figure
    surf(X1_test,X2_test,reshape(m,length(x2_test),length(x1_test)),'EdgeColor', 'None');
    colormap(jet); axis equal; axis off;
    xlabel('x1'); ylabel('x2'); view([0 -90])
end
% fileID = fopen('E:\OneDrive - hnu.edu.cn\tools\matlabcourse\GPML_matlab\GPML-primer\reconstruction\GP_AFM_Pp_part3\result\PSNR_f2000_EChemist_Corrosion.txt','a');
% fprintf(fileID, '%8f\n', PSNR);
% fclose(fileID);