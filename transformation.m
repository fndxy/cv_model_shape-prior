function[A] = transformation(B,a,b,r,alpha)
[m,n] = size(B);
theata = alpha;

T2 = [1,0,0;0,1,0;-m/2,-n/2,1];  %x、y轴平移值原点  
T3 = [1,0,0;0,1,0;m/2,n/2,1];    %x、y轴反平移  
A = -1*ones(m,n);
T1 = [r,0,0;0,r,0;0,0,1];   %对应的比例系数矩阵 
T = T2*T1*T3;
T1 = [cos(theata),-sin(theata),0;sin(theata),cos(theata),0;0,0,1];
T = T*(T2*T1*T3);
T1 = [1 0 0;0 1 0;a b 1];
T = T*T1;
Trans = T^-1;
for i = 1:m
    for j = 1:n
        p = round([i,j,1]*Trans);
        x = p(1);
        y = p(2);
        if x>0 && x<=m && y>0 && y<=n
            A(i,j) = B(x,y);
        end
    end
end
A = (A);
% function A = transformation(B,theta)
% [m,n] = size(B);
% w1 = n;
% h1 = ;
% w2 = w1*cos(theta)+h1*sin(theta);
% h2 = w1*sin(theta)+h1*cos(theta);
% A = zeros(m,n);
% for i = 1:m
%     for j = 1:n
%        x = i*cos(theta) - j*sin(theta) +(-0.5*w1*cos(theta)+0.5*h1*sin(theta)+0.5*w2);
%        y = i*sin(theta) + j*sin(theta) +(-0.5*w1*cos(theta)-0.5*h1*sin(theta)+0.5*h2);
%        x = round(x);
%        y = round(y);
%        if x>0 && x<=m && y>0 && y<=n
%            A(i,j) = B(x,y);
%        end
%     end
% end
% A = uint8(A);
%     
