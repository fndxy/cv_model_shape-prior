function R = reverse(B,a,b,r,alpha)
[m,n] = size(B);
theta =alpha;
R = -1*ones(m,n);

for i = 1:m
    for j = 1:n
        x = r*cos(theta)*i+r*sin(theta)*j-m*r/2*cos(theta)-n*r/2*sin(theta)+a+m/2;  %x*\y*
        y = -1*r*sin(theta)*i+r*cos(theta)*j+m*r/2*sin(theta)-n*r*cos(theta)/2+b+n/2;
        x = round(x);
        y = round(y);
        if x>0 && x<=m && y>0 && y<=n
            R(i,j) = B(x,y);
        end
    end
end
R = (R);
