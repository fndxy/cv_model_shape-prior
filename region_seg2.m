function seg = region_seg2(I,init_mask,max_its,m,alpha,const)
% ������������������״Լ��������ֻ��������״���û���漰������任��
% 
%
%         Inputs: I           ԭͼ
%         init_mask   ��ʼ������1Ϊǰ����0Ϊ����
%         max_its     ��������
%         alpha       ����ƽ����Ĳ�����Խ��Խƽ����Ĭ��Ϊ0.2
  if(~exist('alpha','var')) 
    alpha = .5; 
  end
  if(~exist('const','var')) 
    const = .2; 
  end
  %ͼ��Ԥ����
  [a, b, c] = size(I);
  if(isfloat(I)) % double��  Determine whether input is floating-point array
    if(c==3) 
      I = rgb2gray(uint8(I)); 
    end
  else           % int��
    if(c==3) 
      I = rgb2gray(I); 
    end
    I = double(I);
  end
  
%  g=stopfunction1(I);% ֹͣ����g
 
%  prishape = bwdist(primask)-bwdist(1-primask)+im2double(primask)-.5; % SDF ���뺯��  
 phi = init_mask;%����ͼ���е㵽��ʼ������ŷ����þ���(Euclidean)
 load priorshape phi0;
 prishape = phi0;

 
 a = 0;
 b = 0;
 r = 1;
 theta = 5/180*pi;
 [f_fx,f_fy] = gradient(prishape);
 [m,n] = size(phi0);
 
 tmp = 2*(I-min(I(:)))/(max(I(:))-min(I(:))) - 1;
 I = tmp;
  %--main loop
  for its = 1:max_its  
      prishape = transformation(phi0,a,b,r,theta);  %���α任
      
      H = Heavisidef(phi,0.1);
      H_shape = Heavisidef(prishape,0.1);
      
      D = Dirac(phi,0.1);
      D_shape = Dirac(prishape,0.1);
      
      K = curvature(phi);
      u = sum(sum(I.*H))/sum(sum(H));
      v = sum(sum(I.*(1-H)))/sum(sum(1-H));
      
      eshape = H - H_shape;
      if(its==1)
          eshape
      end
      num1 = normalization(I-u);
      num2 = normalization(I-v);
      num3 = normalization(eshape);
     
 
      
      a_sum = zeros(m,n);
      b_sum = zeros(m,n);
      theta_sum = zeros(m,n);
      r_sum = zeros(m,n);
 
      labmda = 10;
      mu = 10;
      etaa = 20;
       phi = phi + 0.05*D.*( K - num1.^2+num2.^2-1*eshape);
      
           if mod(its,1)==0
          figure(1)
          imshow(I);
          hold on
          contour(phi,[0,0],'r');
         contour(prishape,[0,0],'g');
          pause(0.01);
         
          hold off
           end
      
           
                for tm = 1:1
      tmp_E_shape = H_shape - H;
%       tmp_E_shape = normalization(tmp_E_shape);
      for i = 1:m
          for j = 1:n
         
              x = r*cos(theta)*i+r*sin(theta)*j-m*r/2*cos(theta)-n*r/2*sin(theta)+a+m/2;  %x*\y*
              y = -1*r*sin(theta)*i+r*cos(theta)*j+m*r/2*sin(theta)-n*r*cos(theta)/2+b+n/2;
              x = round(x);
              y = round(y);
              if x>0 && x<=m && y>0 && y<=n
                  D_s = D_shape(i,j);
                  
                  E_shap = tmp_E_shape(i,j);
%                   a_sum(x,y) = E_shap*D_s*(f_fx(x,y));
%                   b_sum(x,y) = E_shap*D_s*(f_fy(x,y));
%                     tmp_rx = cos(theta)*i+sin(theta)*j-m/2*cos(theta)-n/2*sin(theta); %dx*/dr
                    tmp_rx = (x-a-m/2)/r;
%                     tmp_ry = -sin(theta)*i+cos(theta)*j+m/2*sin(theta)-n/2*cos(theta); %dy*/dr
                    tmp_ry = (y-b-n/2)/r;
%                     r_sum(x,y) = E_shap*D_s*(f_fx(x,y)*tmp_rx+f_fy(x,y)*tmp_ry);
               

                  tmp_thetax = y-b-n/2; %dx*/dtheta
                  tmp_thetay = -(x-a-m/2);
                  theta_sum(x,y) = E_shap*D_s*(f_fx(x,y)*tmp_thetax+f_fy(x,y)*tmp_thetay);
              end
          end
      end
 
%       a = a-0.01*sum(sum(a_sum));
%       b = b-0.01*sum(sum(b_sum));
%         r = r+0.002*sum(sum(r_sum));
      theta = theta+0.002*sum(sum(theta_sum));
      end
           
           
           
  end
      
  %-- ���ս�����
%     showCurveAndPhi(I,phi,its);
    seg = phi<=0;
    
    figure(2)
    imshow(I)
    hold on
    [c,h] = contour(phi,[0,0],'r');
    hold off
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%-- AUXILIARY FUNCTIONS ----------------------------------------------
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%-- compute curvature along SDF ���ﲻ��ʹ��gradient �������������ݶ���������  ��Ϊ�����խ����һ������
%%�ڲ����ⲿ��ֵҲ�����ض�ֵ�ġ�
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
%     phi_x  = -phi(idlt)+phi(idrt);                 %��������x�����һ�ײ��
%     phi_y  = -phi(iddn)+phi(idup);                 %��������y�����һ�ײ��
%     phi_xx = phi(idlt)-2*phi(idx)+phi(idrt);       %��������x����Ķ��ײ��
%     phi_yy = phi(iddn)-2*phi(idx)+phi(idup);       %��������x����Ķ��ײ��
%     phi_xy = -0.25*phi(iddl)-0.25*phi(idur)...
%              +0.25*phi(iddr)+0.25*phi(idul);       %���������Ĳ��
%     phi_x2 = phi_x.^2;
%     phi_y2 = phi_y.^2;    
%     %-- compute curvature (Kappa) ��������
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
  
  
%-- ��ʾͼ��������
function showCurveAndPhi(I, phi, i)
  imshow(I,'initialmagnification',200,'displayrange',[0 255]); hold on;
  contour(phi, [0 0], 'g'); hold off; title([num2str(i) '�ε���']); drawnow;

  
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
  




