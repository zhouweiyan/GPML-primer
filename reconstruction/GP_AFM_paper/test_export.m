% https://ww2.mathworks.cn/help/matlab/ref/fprintf.html?searchHighlight=fprintf&s_tid=doc_srchtitle#btf98dm
% zwy 20190222

clear
clc
A1 = [9.9, 9900]
A2 = [8.8, 7.7; 8800, 7700]

formatSpec = 'X is %4.2f meters or %8.3f mm\n'
fprintf(formatSpec, A1, A2)
%%
x = 0:0.1:1;
A = [x; exp(x)];
fileID = fopen('E:\OneDrive - hnu.edu.cn\tools\matlabcourse\GPML_matlab\GPML-primer\reconstruction\GP_AFM_paper\delta_PSNR\exp.txt', 'a')
fprintf(fileID, '%6s %12s\n', 'x', 'exp(x)')
fprintf(fileID, '%6.2f %12.8f\n', A)
fclose(fileID)





