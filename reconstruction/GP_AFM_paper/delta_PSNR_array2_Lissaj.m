% simulate Lissajous scanning pattern
% Reference: T. Tuma, J. Lygeros, V. Kartik, A. Sebastian, and A. Pantazi, ¡°High-speed multiresolution scanning probe microscopy based on Lissajous scan trajectories,¡± 
% Nanotechnology, vol. 23, no. 18, p. 185501, May 2012.
% https://ww2.mathworks.cn/help/matlab/ref/unique.html?searchHighlight=unique&s_tid=doc_srchtitle
% https://ww2.mathworks.cn/help/matlab/ref/uniquetol.html
% zhouweiyan 20190223
% 

clear
close all
% clc

%% Lissajous scan trajectory
X0=64;Y0=64;
A=128;B=128;    % A,B represent the range of X,Y respectively
Fxy = [41,43; 38,39; 33,35; 29,31; 25,26; 21,22; 17,18; 12,13; 8,9; 4,5];
Fx = Fxy(:, 1); Fy = Fxy(:, 2);
PSNR = zeros(1,size(Fxy, 1));
for i = 1:size(Fxy, 1)
    i
    load('Copy_of_array2_GP.mat')
%     opt='GP';
    opt='DT';
    fx=Fx(i); fy=Fy(i);
    f=5000;
    T=1;  %T = lcm(fx,fy)/(fx*fy)
    t=(0:1/f:T)';
    x1_train_lisaj=X0+A/2*sin(2*pi*fx*t);
    x2_train_lisaj=Y0+B/2*sin(2*pi*fy*t);
    
    %% data preparation
    X_train_lisaj=[x1_train_lisaj(:),x2_train_lisaj(:)];
    [y_train_lisaj,y_t_s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_train_lisaj);
    % approximate path
    path_len=sum(sqrt((x1_train_lisaj(2:end)-x1_train_lisaj(1:(end-1))).^2+(x2_train_lisaj(2:end)-x2_train_lisaj(1:(end-1))).^2));
    delta = path_len/(2*128*128);
    % test set
    x1_test=2:127;x2_test=2:127;
    [X1_test,X2_test]=meshgrid(x1_test,x2_test);
    X_test=[X1_test(:),X2_test(:)];
    [y_test_ideal,y_test_s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train,y_train,X_test);
    y_test_ideal=reshape(y_test_ideal,length(x2_test),length(x1_test));
    
    [X_train_lisaj,ia,ic]=uniquetol(X_train_lisaj,'ByRows',true);
    y_train_lisaj=y_train_lisaj(ia);
    
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
            hyp=minimize(hyp,@gp,-80,@infGaussLik,meanfunc,covfunc,likfunc,X_train_lisaj,y_train_lisaj);
            [m,s2]=gp(hyp,@infGaussLik,meanfunc,covfunc,likfunc,X_train_lisaj,y_train_lisaj,X_test);
            
        case 'DT'
            P=X_train_lisaj;
            DT=delaunayTriangulation(P);
            % triplot(DT)
            V=y_train_lisaj;
            Pq=X_test;
            
            [ti,bc]=pointLocation(DT,Pq);
            triVals=V(DT(ti,:));
            Vq=dot(bc',triVals')';
            m=Vq;
    end
    
    %% Error analysis
    m_a=reshape(m,length(x2_test),length(x1_test));
    PSNR(i)=max(y_test_ideal(:))*sqrt(size(y_test_ideal,1)*size(y_test_ideal,2))/sqrt(sum((m_a(:)-y_test_ideal(:)).^2));
    PSNR(i)=20*log(PSNR(i))/log(10)
end
fileID = fopen('E:\OneDrive - hnu.edu.cn\tools\matlabcourse\GPML_matlab\GPML-primer\reconstruction\GP_AFM_paper\delta_PSNR\PSNR_array2_Lissaj.txt','a');
fprintf(fileID, '%8f\n', PSNR);
fclose(fileID);