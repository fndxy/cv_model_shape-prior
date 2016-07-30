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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�ֶ�ѡ���ʼ����
% subplot(2,2,1); imshow(I); title('ԭͼ��');
% subplot(2,2,2); imshow(I); title('��ʼ����');
% hold on;
% mask = roipoly;  %���ó�ʼ������ mask �ڲ�Ϊ1 �ⲿΪ0
% contour(mask,[0 0],'r');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�Զ����ó�ʼ����Ϊ����
k1=size(I,1)/3; 
k2=size(I,2)/3;
mask = -1*ones(size(I,1),size(I,2));
mask(20:size(I,1)-20,20:size(I,2)-20) = 1;
% subplot(2,2,1); imshow(I); title('ԭͼ��');
% subplot(2,2,2); imshow(mask); title('��ʼ����');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���濪ʼ����  
% subplot(2,2,3);

tic;
seg = region_seg2(I, mask, 800,m); %--ֻ���򵥵ļ�Լ��
toc;
% subplot(2,2,4); imshow(seg); title('���շָ���');

