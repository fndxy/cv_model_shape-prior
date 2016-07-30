function seg = region_seg2(I,init_mask,max_its,m,alpha,const)
% 本函数加入了先验形状约束，但是只是两个形状相减没有涉及到仿射变换。
% 
%
%         Inputs: I           原图
%         init_mask   初始轮廓，1为前景，0为背景
%         max_its     迭代次数
%         alpha       控制平滑项的参数，越大越平滑，默认为0.2
  if(~exist('alpha','var')) 
    alpha = .5; 
  end
  if(~exist('const','var')) 
    const = .2; 
  end
  %图像预处理
  [a, b, c] = size(I);
  if(isfloat(I)) % double型  Determine whether input is floating-point array
    if(c==3) 
      I = rgb2gray(uint8(I)); 
    end
  else           % int型
    if(c==3) 
      I = rgb2gray(I); 
    end
    I = double(I);
  end
  
%  g=stopfunction1(I);% 停止函数g
 primask=1-m;
%  prishape = bwdist(primask)-bwdist(1-primask)+im2double(primask)-.5; % SDF 距离函数  
 phi = init_mask;%计算图像中点到初始轮廓的欧几里得距离(Euclidean)
 load priorshape phi0;
 prishape = phi0;
 H_shape = Heavisidef(prishape,0.1);
  %--main loop
  for its = 1:max_its  
      H = Heavisidef(phi,0.1);
      D = Dirac(phi,0.1);
      K = curvature(phi);
      u = sum(sum(I.*H))/sum(sum(H));
      v = sum(sum(I.*(1-H)))/sum(sum(1-H));
      mu = 20;
      lambda = 5;
      eshape = H - H_shape;
      num1 = normalization(I-u);
      num2 = normalization(I-v);
      num3 = normalization(eshape);
      phi = phi + 0.1*D.*(K - num1.^2+num2.^2-num3);
      if mod(its,1)==0
          figure(1)
          imshow(I);
          hold on
          contour(phi,[0,0],'r');
          pause(0.1);
          hold off
      end
  end
      
  %-- 最终结果输出
    showCurveAndPhi(I,phi,its);
    seg = phi<=0;
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%-- AUXILIARY FUNCTIONS ----------------------------------------------
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%-- compute curvature along SDF 这里不能使用gradient 函数求横纵向的梯度在求曲率  因为这里的窄带是一个环形
%%内部和外部的值也是有特定值的。
% function curvature = get_curvature(phi,idx)  %  
%     [dimy, dimx] = size(phi);        
%     [y x] = ind2sub([dimy,dimx],idx);  
%     %-- get subscripts of neighbors
%     ym1 = y-1; xm1 = x-1;             
%     yp1 = y+1; xp1 = x+1;              
%     %-- bounds checking  
%     ym1(ym1<1) = 1; xm1(xm1<1) = 1;              
%     yp1(yp1>dimy)=dimy; xp1(xp1>dimx) = dimx;
%     %-- get indexes for 8 neighbors
%     idup = sub2ind(size(phi),yp1,x);  
%     iddn = sub2ind(size(phi),ym1,x);   
%     idlt = sub2ind(size(phi),y,xm1);
%     idrt = sub2ind(size(phi),y,xp1);
%     idul = sub2ind(size(phi),yp1,xm1);
%     idur = sub2ind(size(phi),yp1,xp1);
%     iddl = sub2ind(size(phi),ym1,xm1);
%     iddr = sub2ind(size(phi),ym1,xp1);    
%     %-- get central derivatives of SDF at x,y
%     phi_x  = -phi(idlt)+phi(idrt);                 %这里是求x方向的一阶差分
%     phi_y  = -phi(iddn)+phi(idup);                 %这里是求y方向的一阶差分
%     phi_xx = phi(idlt)-2*phi(idx)+phi(idrt);       %这里是求x方向的二阶差分
%     phi_yy = phi(iddn)-2*phi(idx)+phi(idup);       %这里是求x方向的二阶差分
%     phi_xy = -0.25*phi(iddl)-0.25*phi(idur)...
%              +0.25*phi(iddr)+0.25*phi(idul);       %这里是中心差分
%     phi_x2 = phi_x.^2;
%     phi_y2 = phi_y.^2;    
%     %-- compute curvature (Kappa) 计算曲率
%     curvature = ((phi_x2.*phi_yy + phi_y2.*phi_xx - 2*phi_x.*phi_y.*phi_xy)./...
%               (phi_x2 + phi_y2 +eps).^(3/2)).*(phi_x2 + phi_y2).^(1/2);      



          
function D = sussman(D, dt)
  % forward/backward differences
  a = D - shiftR(D); % backward
  b = shiftL(D) - D; % forward
  c = D - shiftD(D); % backward
  d = shiftU(D) - D; % forward
  
  a_p = a;  a_n = a; % a+ and a-
  b_p = b;  b_n = b;
  c_p = c;  c_n = c;
  d_p = d;  d_n = d;
  
  a_p(a < 0) = 0;
  a_n(a > 0) = 0;
  b_p(b < 0) = 0;
  b_n(b > 0) = 0;
  c_p(c < 0) = 0;
  c_n(c > 0) = 0;
  d_p(d < 0) = 0;
  d_n(d > 0) = 0;
  
  dD = zeros(size(D));
  D_neg_ind = find(D < 0);
  D_pos_ind = find(D > 0);
  dD(D_pos_ind) = sqrt(max(a_p(D_pos_ind).^2, b_n(D_pos_ind).^2) ...
                       + max(c_p(D_pos_ind).^2, d_n(D_pos_ind).^2)) - 1;
  dD(D_neg_ind) = sqrt(max(a_n(D_neg_ind).^2, b_p(D_neg_ind).^2) ...
                       + max(c_n(D_neg_ind).^2, d_p(D_neg_ind).^2)) - 1;
  
  D = D - dt .* sussman_sign(D) .* dD;
  
%-- whole matrix derivatives

function shift = shiftD(M)
  shift = shiftR(M')';

function shift = shiftL(M)
  shift = [ M(:,2:size(M,2)) M(:,size(M,2)) ];

function shift = shiftR(M)
  shift = [ M(:,1) M(:,1:size(M,2)-1) ];

function shift = shiftU(M)
  shift = shiftL(M')';
  
function S = sussman_sign(D)
  S = D ./ sqrt(D.^2 + 1);    
  
  
%-- 显示图像及轮廓线
function showCurveAndPhi(I, phi, i)
  imshow(I,'initialmagnification',200,'displayrange',[0 255]); hold on;
  contour(phi, [0 0], 'g'); hold off; title([num2str(i) '次迭代']); drawnow;

  
function H = Heavisidef(phi,sigma)
[m,n] = size(phi);
H = zeros(m,n);

for i =1:m
    for j = 1:n
        H(i,j) = 0.5*(1+2/pi*atan(phi(i,j)/sigma));
    end
end
function D = Dirac(phi,sigma)
[m,n] = size(phi);
D = zeros(m,n);

for i =1:m
    for j = 1:n
        D(i,j) = 1/pi*(sigma/(sigma^2+phi(i,j)^2));
    end
end 
function D = normalization(A)
t1 = max(A(:));
t2 = min(A(:));
D = (A - t2)/(t1-t2);
  




