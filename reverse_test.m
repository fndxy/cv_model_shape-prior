function[A,p] = reverse_test(B,alpha,a,b,r)
[m,n] = size(B);
theata = alpha/180*pi;
A = -1*ones(m,n);
% T2 = [1,0,0;0,1,0;-m/2,-n/2,1];  %x、y轴平移值原点  
% T3 = [1,0,0;0,1,0;m/2,n/2,1];    %x、y轴反平移  
% T1 = [1/r,0,0;0,1/r,0;0,0,1];   %对应的比例系数矩阵 
% T = T2*T1*T3;
% T1 = [cos(theata),sin(theata),0;-sin(theata),cos(theata),0;0,0,1];
% T = T*T2*T1*T3;
% T1 = [1 0 0;0 1 0;-a -b 1];
% T = T*T1;
% Trans = T^-1;
% for i = 1:m
%     for j = 1:n
%         p = round([i,j,1]*Trans);
%         x = p(1);
%         y = p(2);
%         if x>0 && x<=m && y>0 && y<=n
%             A(i,j) = B(x,y);
%         end
%     end
% end
% A = uint8(A);
p = zeros(m*n,4);
index=1;
for i = 1:m
    for j = 1:n
        x = ((i+a)*cos(theata)+(j+b)*sin(theata))*r;
        y = (-(i+a)*sin(theata)+(j+b)*cos(theata))*r;
        
        x =round(x);
        y = round(y);
        if x>0 && x<=m && y>0 && y<=n
            A(i,j) = B(x,y);
            p(index,:) = [x,y,i,j];
            index = index+1;
        end
    end
end
A = (A);
p = p(1:index-1,:);

        