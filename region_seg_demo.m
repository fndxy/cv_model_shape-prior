% Demo of "Region Based Active Contours"
% Example:% seg_demo
% Coded by: Shawn Lankton (www.shawnlankton.com)
clear all;clc;close all


 I=imread('Fig12.1-1.jpg');
 I=I(:,:,1); I=imresize(I,0.5);
 size(I)
m=imread('Fig12.1.jpg');m=m(:,:,1);m=imresize(m,0.5);
size(m);
% m=im2bw(m);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%手动选择初始轮廓
% subplot(2,2,1); imshow(I); title('原图像');
% subplot(2,2,2); imshow(I); title('初始轮廓');
% hold on;
% mask = roipoly;  %设置初始轮廓线 mask 内部为1 外部为0
% contour(mask,[0 0],'r');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%自动设置初始轮廓为矩形
k1=size(I,1)/3; 
k2=size(I,2)/3;
mask = -1*ones(size(I,1),size(I,2));
mask(20:size(I,1)-20,20:size(I,2)-20) = 1;
% subplot(2,2,1); imshow(I); title('原图像');
% subplot(2,2,2); imshow(mask); title('初始轮廓');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%下面开始迭代  
% subplot(2,2,3);

tic;
seg = region_seg2(I, mask, 800,m); %--只做简单的减约束
toc;
% subplot(2,2,4); imshow(seg); title('最终分割结果');

